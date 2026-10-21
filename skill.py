"""Score every archived StormSense Ensemble run against the water levels observed after it,
with Derek's skill formula, and pick each station's primary guidance per lead band (phase 2).

    python skill.py <out_dir> <work_dir> [run_tarball ...]

Inputs: the run tarballs on this repo's ensemble-YYYY-MM releases (the last RESCORE_DAYS of runs
are downloaded and rescored each time, plus any tarballs named on the command line), and the
cumulative scores of the earlier runs (skill_scores.csv.gz and skill_runs.csv on the newest
release that holds them).

Truth: the observed_ft column of the archived ens72 CSVs, stitched across runs (each run carries
the 24 hours before its last observation; the newest run wins where they overlap). These are the
observations the ensemble lines its members up with; CBBT, MNPT, DAHL and WADC are NOAA CO-OPS.

Each run x station x member x lead band (hours after the station's last observation at run time)
is scored once observations cover MIN_COVERAGE of the band:

    S = 0.25 [R2 + e^(-|Bias|/sigma) + e^(-MAE/sigma) + e^(-RMSE/sigma)],  sigma = population SD of
    the observations in the band (the same formula and code as the Tidewatch ASA, asa/metrics.py)

Adaptive pick: at each update, a member's recent skill at a station and band is the mean of its
scores over the last WINDOW_DAYS, weighted by exp(-age / EFOLD_DAYS). The leader must beat the
current pick by MARGIN for CONFIRM updates running before the pick switches. A station and band
stay on P-ETSS until READY_DAYS of verified runs exist and the member has MIN_RUNS scores.

Outputs in out_dir: scores.csv, summary.csv, summary.json, pick.json, picks_history.csv,
index.html; in work_dir: skill_scores.csv.gz and skill_runs.csv for the release.
"""
import glob
import html as H
import io
import json
import os
import sys
import tarfile
import time
import urllib.request
from datetime import datetime, timedelta, timezone

import numpy as np
import pandas as pd

REPO = os.environ.get("GITHUB_REPOSITORY", "DerekLoftis/StormSense-ensemble")
UA = {"User-Agent": "stormsense-ensemble (VIMS)"}
# (label, lower, upper) in hours after the station's last observation at run time
BANDS = (("0-12", 0, 12), ("12-36", 12, 36), ("36-72", 36, 72), ("0-36", 0, 36), ("0-72", 0, 72))
PICK_BANDS = ("0-12", "12-36", "36-72")
STEP_MIN = 6
MIN_COVERAGE = 0.85
RESCORE_DAYS = 5        # runs newer than this are rescored each time (72 h + slack for late runs)
WINDOW_DAYS = 10        # recent-skill window
EFOLD_DAYS = 3.0        # recency weighting
READY_DAYS = 7          # verified history needed before the pick may leave P-ETSS
MIN_RUNS = 12           # scores a member needs in the window to be picked
MARGIN = 0.01           # skill a challenger must lead by ...
CONFIRM = 2             # ... for this many updates running
REPLAY_DAYS = 30        # the pick is replayed from P-ETSS through this many days of updates
PUBLISH_DAYS = 31       # runs in the published scores.csv (the release keeps them all)
DEFAULT = "petss"
NOT_MEMBERS = {"observed", "astro_tide"}
# Members that can drive the primary forecast. Tidewatch is the other brand and the baseline and
# the tide alone are references every forecast should beat; they are scored but never picked.
# stormsense_raw is StormSense before its site corrector, kept to check that the corrector helps.
REFERENCE = {"tidewatch", "baseline", "tide", "stormsense_raw"}
# Tidewatch publishes another gauge's observations under these ids (handoff aliasing check,
# Oct 2026; USBB is held in ~18-min steps): scored for the record, never used to pick guidance.
SUSPECT_OBS = {"USSC", "USAR", "USCB", "USNR", "USSB", "GS01", "USBB"}
NAMES = {"petss": "NWS P-ETSS", "stofs2d": "NOAA STOFS-2D-Global", "stofs3d": "NOAA STOFS-3D-Atlantic",
         "cbofs": "NOAA CBOFS", "nwm_short": "NWM coastal (short range)", "nwm_medium": "NWM coastal (medium range)",
         "etss": "NWS ETSS", "psurge": "NHC P-Surge (50%)", "nws": "NWS official forecast",
         "tidewatch": "VIMS Tidewatch", "baseline": "Tide + persistence", "tide": "Astronomical tide only",
         "stormsense": "StormSense", "stormsense_raw": "StormSense (uncorrected)"}
SCORE_COLS = ["cycle", "station", "member", "member_run", "band", "t0", "lead_end", "N", "R2", "Bias", "MAE",
              "RMSE", "SDobs", "Skill", "PeakError", "ConeHit"]


# --- Derek's statistics (same as handoff/skill/tidewatch-asa/asa/metrics.py) -------------------
def statistics(x, y):
    """x observed, y forecast: R2 (squared Pearson r), Bias, MAE, RMSE, SDobs (population), PeakError."""
    err = y - x
    r = np.corrcoef(x, y)[0, 1] if x.std() > 0 and y.std() > 0 else np.nan
    return {"R2": r * r, "Bias": err.mean(), "MAE": np.abs(err).mean(), "RMSE": np.sqrt((err ** 2).mean()),
            "SDobs": x.std(ddof=0), "PeakError": y.max() - x.max()}


def skill_score(s):
    """S = 0.25 [R2 + e^(-|Bias|/sd) + e^(-MAE/sd) + e^(-RMSE/sd)], sd the population SD of obs."""
    sd = s["SDobs"]
    if not (sd > 0) or np.isnan(s["R2"]):
        return np.nan
    return 0.25 * (float(np.clip(s["R2"], 0, 1)) + np.exp(-abs(s["Bias"]) / sd)
                   + np.exp(-s["MAE"] / sd) + np.exp(-s["RMSE"] / sd))


# --- archive access ------------------------------------------------------------------------------
def api(path):
    head = dict(UA)
    tok = os.environ.get("GH_TOKEN") or os.environ.get("GITHUB_TOKEN")
    if tok:
        head["Authorization"] = f"Bearer {tok}"
    for i in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(path, headers=head), timeout=120) as r:
                return r.read()
        except Exception:
            if i == 2:
                raise
            time.sleep(5 * (i + 1))


def releases():
    """{tag: {asset name: download url}} for the ensemble-YYYY-MM releases, newest month first."""
    out, page = {}, 1
    while True:
        js = json.loads(api(f"https://api.github.com/repos/{REPO}/releases?per_page=100&page={page}"))
        for r in js:
            if r["tag_name"].startswith("ensemble-"):
                out[r["tag_name"]] = {a["name"]: a["browser_download_url"] for a in r["assets"]}
        if len(js) < 100:
            break
        page += 1
    return dict(sorted(out.items(), reverse=True))


def cycle_of(name):
    return datetime.strptime(os.path.basename(name)[9:23], "%Y-%m-%dT%HZ").replace(tzinfo=timezone.utc)


def gather(work, local, now):
    """Download the cached scores and the recent run tarballs. Returns (scores, runs, tarballs)."""
    tars = {cycle_of(p): p for p in local}
    scores, runs = pd.DataFrame(columns=SCORE_COLS), pd.DataFrame(columns=["cycle", "updated", "station", "t0"])
    try:
        rel = releases()
    except Exception as e:
        print("releases unavailable:", e)
        rel = {}
    cached = next((a for a in rel.values() if "skill_scores.csv.gz" in a and "skill_runs.csv" in a), None)
    if cached:
        scores = pd.read_csv(io.BytesIO(api(cached["skill_scores.csv.gz"])), compression="gzip")
        runs = pd.read_csv(io.BytesIO(api(cached["skill_runs.csv"])))
        print(f"cached scores: {len(scores)} rows, {runs['cycle'].nunique()} runs")
    since = now - timedelta(days=RESCORE_DAYS) if cached else datetime(2000, 1, 1, tzinfo=timezone.utc)
    os.makedirs(os.path.join(work, "tars"), exist_ok=True)
    for assets in rel.values():
        for name, url in assets.items():
            if not (name.startswith("ensemble_") and name.endswith(".tar.gz")):
                continue
            c = cycle_of(name)
            if c >= since and c not in tars:
                p = os.path.join(work, "tars", name)
                if not os.path.exists(p):
                    open(p, "wb").write(api(url))
                tars[c] = p
    print(f"{len(tars)} run tarballs to score")
    return scores, runs, dict(sorted(tars.items()))


def read_run(path):
    """(updated, {station: (t0, DataFrame of the MLLW ens72 CSV indexed by time)})."""
    out = {}
    with tarfile.open(path) as t:
        latest = json.load(t.extractfile("latest.json"))
        for m in t.getmembers():
            n = os.path.basename(m.name)
            if not (m.name.startswith("ens72/") and n.endswith("_MLLW.csv")):
                continue
            st = n[:-len("_ens72_MLLW.csv")]
            info = latest["stations"].get(st)
            if not info or not info.get("t0"):
                continue
            df = pd.read_csv(t.extractfile(m))
            df.index = pd.to_datetime(df.pop("date_time_z"), utc=True)
            out[st] = (pd.Timestamp(info["t0"]), df, info.get("members", {}))
    updated = pd.Timestamp(datetime.strptime(latest["updated"], "%Y-%m-%d %H:%M UTC").replace(tzinfo=timezone.utc))
    return updated, out


# --- scoring ---------------------------------------------------------------------------------------
def score_runs(tars):
    """Score every run in tars against the observations stitched from all of them."""
    loaded = {c: read_run(p) for c, p in tars.items()}
    obs = {}
    for c in sorted(loaded):  # newer runs overwrite older ones
        for st, (_, df, _) in loaded[c][1].items():
            o = df["observed_ft"].dropna()
            obs.setdefault(st, {}).update(o.to_dict())
    obs = {st: pd.Series(v).sort_index() for st, v in obs.items()}
    rows, runs = [], []
    for c, (updated, stations) in loaded.items():
        for st, (t0, df, meta) in stations.items():
            runs.append({"cycle": f"{c:%Y-%m-%dT%HZ}", "updated": updated.isoformat(), "station": st, "t0": t0.isoformat()})
            truth = obs.get(st)
            if truth is None or truth.empty:
                continue
            cols = {k[:-3]: k for k in df.columns if k.endswith("_ft") and not k.endswith(("_p10_ft", "_p90_ft"))}
            cols = {("tide" if k == "astro_tide" else k): v for k, v in cols.items() if k != "observed"}
            for band, lo, hi in BANDS:
                if t0 + pd.Timedelta(hours=hi) > truth.index[-1]:
                    continue  # not verified yet
                sel = (df.index > t0 + pd.Timedelta(hours=lo)) & (df.index <= t0 + pd.Timedelta(hours=hi))
                win = df.loc[sel]
                x_all = truth.reindex(win.index)
                need = MIN_COVERAGE * (hi - lo) * 60 / STEP_MIN
                for m, col in cols.items():
                    ok = x_all.notna() & win[col].notna()
                    if ok.sum() < need:
                        continue
                    x, y = x_all[ok].to_numpy(float), win.loc[ok, col].to_numpy(float)
                    s = statistics(x, y)
                    cone = np.nan
                    if f"{m}_p10_ft" in win and f"{m}_p90_ft" in win:
                        lo_c, hi_c = win.loc[ok, f"{m}_p10_ft"].to_numpy(float), win.loc[ok, f"{m}_p90_ft"].to_numpy(float)
                        good = ~np.isnan(lo_c) & ~np.isnan(hi_c)
                        if good.any():
                            cone = float(np.mean((x[good] >= lo_c[good]) & (x[good] <= hi_c[good])))
                    rows.append({"cycle": f"{c:%Y-%m-%dT%HZ}", "station": st, "member": m,
                                 "member_run": (meta.get(m) or {}).get("run", ""), "band": band,
                                 "t0": t0.isoformat(), "lead_end": (t0 + pd.Timedelta(hours=hi)).isoformat(),
                                 "N": int(ok.sum()), **{k: s[k] for k in ("R2", "Bias", "MAE", "RMSE", "SDobs")},
                                 "Skill": skill_score(s), "PeakError": s["PeakError"], "ConeHit": cone})
    return pd.DataFrame(rows, columns=SCORE_COLS), pd.DataFrame(runs)


def merge(old, new, rescored, order=("cycle", "station", "member", "band")):
    """Old rows for runs not rescored now, plus the new rows."""
    if old is not None and not old.empty:
        new = pd.concat([old[~old["cycle"].isin(rescored)], new], ignore_index=True)
    return new.sort_values(list(order)).reset_index(drop=True)


# --- adaptive pick -------------------------------------------------------------------------------
def recent(scores, at, bands=PICK_BANDS):
    """Recent skill per station x band x member as known at time `at` (verified by then):
    exponentially weighted means over the last WINDOW_DAYS. No look-ahead."""
    s = scores[scores["band"].isin(bands)].copy()
    s["lead_end"] = pd.to_datetime(s["lead_end"], utc=True)
    s = s[(s["lead_end"] <= at) & (s["lead_end"] > at - pd.Timedelta(days=WINDOW_DAYS)) & s["Skill"].notna()]
    if s.empty:
        return pd.DataFrame(columns=["station", "band", "member", "n", "skill_ew", "rmse_ew", "mae_ew", "bias_ew"])
    s["w"] = np.exp(-(at - s["lead_end"]).dt.total_seconds() / 86400 / EFOLD_DAYS)
    g = s.assign(ws=s["w"] * s["Skill"], wr=s["w"] * s["RMSE"], wm=s["w"] * s["MAE"], wb=s["w"] * s["Bias"])
    a = g.groupby(["station", "band", "member"]).agg(n=("w", "size"), w=("w", "sum"), ws=("ws", "sum"),
                                                     wr=("wr", "sum"), wm=("wm", "sum"), wb=("wb", "sum")).reset_index()
    for k, c in (("skill_ew", "ws"), ("rmse_ew", "wr"), ("mae_ew", "wm"), ("bias_ew", "wb")):
        a[k] = a[c] / a["w"]
    return a[["station", "band", "member", "n", "skill_ew", "rmse_ew", "mae_ew", "bias_ew"]]


def replay(scores, runs):
    """Run the pick with hysteresis through every update in order, so the current pick needs no
    saved state. Returns (current picks, history rows)."""
    sc = scores.copy()
    sc["lead_end_t"] = pd.to_datetime(sc["lead_end"], utc=True)
    first = sc[sc["band"].isin(PICK_BANDS)].groupby(["station", "band"])["lead_end_t"].min()
    ups = runs.assign(updated=pd.to_datetime(runs["updated"], utc=True)).groupby("cycle")["updated"].max().sort_index()
    ups = ups[ups > ups.max() - pd.Timedelta(days=REPLAY_DAYS)] if len(ups) else ups
    sc = sc[sc["lead_end_t"] > ups.min() - pd.Timedelta(days=WINDOW_DAYS)] if len(ups) else sc
    state, hist = {}, []
    for cyc, at in ups.items():
        rec = recent(sc, at)
        stations = sorted(runs.loc[runs["cycle"] == cyc, "station"])
        for st in stations:
            for band in PICK_BANDS:
                key = (st, band)
                cur = state.setdefault(key, {"member": DEFAULT, "since": cyc, "challenger": None, "count": 0})
                r = rec[(rec["station"] == st) & (rec["band"] == band)].set_index("member")
                since = first.get(key)
                ready = since is not None and at - since >= pd.Timedelta(days=READY_DAYS) and st not in SUSPECT_OBS
                elig = r[(r["n"] >= MIN_RUNS) & ~r.index.isin(REFERENCE)]
                leader = elig["skill_ew"].idxmax() if not elig.empty else None
                cur_skill = elig["skill_ew"].get(cur["member"], np.nan)
                if not ready:
                    cur.update(member=DEFAULT, challenger=None, count=0)
                elif leader is not None and leader != cur["member"] and (
                        np.isnan(cur_skill) or elig.loc[leader, "skill_ew"] >= cur_skill + MARGIN):
                    cur["count"] = cur["count"] + 1 if cur["challenger"] == leader else 1
                    cur["challenger"] = leader
                    if cur["count"] >= CONFIRM:
                        cur.update(member=leader, since=cyc, challenger=None, count=0)
                else:
                    cur.update(challenger=None, count=0)
                hist.append({"cycle": cyc, "station": st, "band": band, "pick": cur["member"], "ready": ready,
                             "leader": leader, "leader_skill": round(float(elig.loc[leader, "skill_ew"]), 3) if leader else None,
                             "pick_skill": None if np.isnan(elig["skill_ew"].get(cur["member"], np.nan))
                             else round(float(elig["skill_ew"][cur["member"]]), 3),
                             "challenger": cur["challenger"], "count": cur["count"]})
    return state, pd.DataFrame(hist)


# --- outputs -----------------------------------------------------------------------------------------
def r3(v):
    return None if v is None or (isinstance(v, float) and np.isnan(v)) else round(float(v), 3)


def write(out, scores, runs, now):
    os.makedirs(out, exist_ok=True)
    cut = f"{now - pd.Timedelta(days=PUBLISH_DAYS):%Y-%m-%dT%HZ}"
    scores[scores["cycle"] >= cut].to_csv(os.path.join(out, "scores.csv"), index=False, float_format="%.4f")
    state, hist = replay(scores, runs)
    hist.to_csv(os.path.join(out, "picks_history.csv"), index=False)
    rec = recent(scores, now, bands=[b for b, _, _ in BANDS])
    rec.to_csv(os.path.join(out, "summary.csv"), index=False, float_format="%.4f")
    last = hist[hist["cycle"] == hist["cycle"].max()] if not hist.empty else hist
    picks = {}
    for _, h in last.iterrows():
        cur = state[(h["station"], h["band"])]
        picks.setdefault(h["station"], {})[h["band"]] = {
            "member": h["pick"], "since": cur["since"], "ready": bool(h["ready"]), "skill": h["pick_skill"],
            "leader": h["leader"], "leader_skill": h["leader_skill"],
            "challenger": h["challenger"], "challenger_updates": int(h["count"])}
    rule = (f"Per station and lead band (hours after the last observation), the member with the best recent skill "
            f"S = 0.25[R2 + e^(-|Bias|/sd) + e^(-MAE/sd) + e^(-RMSE/sd)] (mean over {WINDOW_DAYS} days, weighted "
            f"exp(-age/{EFOLD_DAYS:g} d)); it must lead the current pick by {MARGIN} for {CONFIRM} updates running. "
            f"P-ETSS until {READY_DAYS} days of verified runs exist and the member has {MIN_RUNS} scores; "
            f"Tidewatch, the baseline and the tide alone are never picked; stations with another gauge's "
            f"observations ({', '.join(sorted(SUSPECT_OBS))}) stay on P-ETSS.")
    doc = {"updated": now.strftime("%Y-%m-%d %H:%M UTC"), "default": DEFAULT, "bands_h": list(PICK_BANDS),
           "rule": rule, "names": NAMES, "stations": picks,
           "ready_stations": sorted(st for st, b in picks.items() if any(v["ready"] for v in b.values()))}
    json.dump(doc, open(os.path.join(out, "pick.json"), "w"), indent=1)
    summ = {}
    for _, r in rec.iterrows():
        summ.setdefault(r["station"], {}).setdefault(r["band"], {})[r["member"]] = {
            "n": int(r["n"]), "skill": r3(r["skill_ew"]), "rmse_ft": r3(r["rmse_ew"]),
            "mae_ft": r3(r["mae_ew"]), "bias_ft": r3(r["bias_ew"])}
    json.dump({"updated": doc["updated"], "window_days": WINDOW_DAYS, "efold_days": EFOLD_DAYS, "stations": summ},
              open(os.path.join(out, "summary.json"), "w"), indent=1)
    write_page(out, scores, rec, picks, doc)
    return doc


def write_page(out, scores, rec, picks, doc):
    """Skill tables: all stations by member and band, then per station."""
    order = [m for m in NAMES if m in set(rec["member"])] + sorted(set(rec["member"]) - set(NAMES))
    bands = [b for b, _, _ in BANDS]
    # overall: mean over stations (stations with suspect observations left out)
    ok = rec[~rec["station"].isin(SUSPECT_OBS)]
    head = "".join(f"<th>{b} h</th>" for b in bands)
    rows = []
    for m in order:
        cells = []
        for b in bands:
            v = ok[(ok["member"] == m) & (ok["band"] == b)]
            cells.append(f"<td>{v['skill_ew'].mean():.3f} <span class='n'>({len(v)} st, {v['rmse_ew'].mean():.2f} ft)</span></td>"
                         if len(v) else "<td></td>")
        rows.append(f"<tr><th>{H.escape(NAMES.get(m, m))}</th>{''.join(cells)}</tr>")
    overall = f"<table><thead><tr><th>Member</th>{head}</tr></thead><tbody>{''.join(rows)}</tbody></table>"
    blocks = []
    for st in sorted(set(rec["station"])):
        r = rec[rec["station"] == st]
        trs = []
        for m in order:
            cells = []
            for b in bands:
                v = r[(r["member"] == m) & (r["band"] == b)]
                if v.empty:
                    cells.append("<td></td>")
                    continue
                v = v.iloc[0]
                best = v["skill_ew"] == r[(r["band"] == b) & ~r["member"].isin(REFERENCE)]["skill_ew"].max()
                cells.append(f"<td{' class=best' if best else ''}>{v['skill_ew']:.3f} <span class='n'>"
                             f"{v['rmse_ew']:.2f} ft, n={int(v['n'])}</span></td>")
            if any(c != "<td></td>" for c in cells):
                trs.append(f"<tr><th>{H.escape(NAMES.get(m, m))}</th>{''.join(cells)}</tr>")
        p = picks.get(st, {})
        pick = "; ".join(f"{b} h: {H.escape(NAMES.get(v['member'], v['member']))}{'' if v['ready'] else ' (not ready)'}"
                         for b, v in p.items())
        flag = " <span class='n'>(observations are another gauge's: not used for the pick)</span>" if st in SUSPECT_OBS else ""
        blocks.append(f"<h3>{st}{flag}</h3><p class='n'>Primary guidance: {pick or 'P-ETSS'}</p>"
                      f"<div class='wrap'><table><thead><tr><th>Member</th>{head}</tr></thead><tbody>{''.join(trs)}</tbody></table></div>")
    ncyc = scores["cycle"].nunique() if not scores.empty else 0
    page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>StormSense Ensemble Skill</title>
<style>
body{{font-family:Arial,'Liberation Sans',sans-serif;margin:0 auto;padding:16px;max-width:1000px;color:#1b2533;background:#fff}}
header img{{width:min(520px,100%);height:auto;display:block}}
h1{{font-size:1.5rem;margin:12px 0 4px;color:#14285a}} h3{{margin:18px 0 2px;color:#14285a}}
p{{margin:6px 0;line-height:1.45}} .disc,.n{{color:#555;font-size:.88rem}}
.wrap{{overflow-x:auto;margin-top:8px}} table{{border-collapse:collapse;font-size:.9rem}}
th,td{{border:1px solid #d5dbe3;padding:4px 8px;text-align:left;white-space:nowrap}}
thead th{{background:#eef6fa;text-align:center}} td.best{{background:#e0fbfc;font-weight:bold}} a{{color:#0d6e8a}}
</style></head>
<body><header><img src="../stormsense_logo.png" alt="StormSense"></header>
<h1>StormSense Ensemble Skill</h1>
<p>Recent skill of each forecast against the observations that came in after it (updated {doc['updated']};
{ncyc} archived runs scored). Lead bands are hours after the station's last observation at run time.</p>
<p class="disc">Skill S = 0.25[R&sup2; + e<sup>&minus;|Bias|/&sigma;</sup> + e<sup>&minus;MAE/&sigma;</sup> + e<sup>&minus;RMSE/&sigma;</sup>],
&sigma; the population SD of the observations in the band: 1 is perfect, about 0.45 is no skill. Each cell is the mean over the
last {WINDOW_DAYS} days weighted toward the newest runs (e-folding {EFOLD_DAYS:g} days), with the weighted RMSE and the number of
runs scored. Highlighted: the best member that can drive the primary forecast.</p>
<p class="disc">{H.escape(doc['rule'])}</p>
<p class="disc">Downloads: <a href="scores.csv">every score (last {PUBLISH_DAYS} days)</a>, <a href="summary.csv">recent skill (csv)</a>,
<a href="summary.json">(json)</a>, <a href="pick.json">guidance pick</a>, <a href="picks_history.csv">pick history</a>.</p>
<p class="disc">Disclaimer: StormSense&trade; Predictive Plots are computer-generated operational storm tide forecasts;
the user assumes any and all risk.</p>
<h3>All stations</h3><p class="n">Mean recent skill over stations (stations, mean RMSE).</p><div class="wrap">{overall}</div>
{''.join(blocks)}
</body></html>
"""
    open(os.path.join(out, "index.html"), "w").write(page)


def main(out="site/skill", work="work/skill", *local):
    now = pd.Timestamp(datetime.now(timezone.utc))
    old_scores, old_runs, tars = gather(work, [p for a in local for p in glob.glob(a)], now)
    new_scores, new_runs = score_runs(tars)
    rescored = set(new_runs["cycle"]) if not new_runs.empty else set()
    scores = merge(old_scores, new_scores, rescored)
    runs = merge(old_runs, new_runs, rescored, ("cycle", "station"))
    os.makedirs(work, exist_ok=True)
    scores.to_csv(os.path.join(work, "skill_scores.csv.gz"), index=False, float_format="%.4f", compression="gzip")
    runs.to_csv(os.path.join(work, "skill_runs.csv"), index=False)
    doc = write(out, scores, runs, now)
    moved = {f"{st} {b}": v["member"] for st, bb in doc["stations"].items() for b, v in bb.items() if v["member"] != DEFAULT}
    print(f"skill: {len(scores)} scores from {scores['cycle'].nunique()} runs; "
          f"ready at {len(doc['ready_stations'])} stations; picks off P-ETSS: {moved or 'none'}")


if __name__ == "__main__":
    main(*sys.argv[1:])
