"""Build the StormSense Ensemble site: every available forecast for each station on
one 72-hour chart, with CSV and JSON, plus an archive tarball of the run.

    python ensemble.py [site_dir] [work_dir]

Members (each where it has a point for the station): NWS P-ETSS (mean and 10-90% cone),
STOFS-2D-Global, STOFS-3D-Atlantic, CBOFS, NWM coastal short and medium range, NWS ETSS
(surge bulletin on the P-ETSS tide), NHC P-Surge (50% level, only while NHC has tropical
advisories out), the official NWS forecast, VIMS Tidewatch, and a tide-plus-persistence baseline.

All members are put on the station's MLLW (Tidewatch's basis; NAVD88 charts subtract
Tidewatch's MLLW-NAVD88 offset). The model guidance (P-ETSS, STOFS, CBOFS, NWM, ETSS, P-Surge) is then
lined up with the gauge: the median difference between the observations and the run over
the overlap in the last 24 hours is added to the whole run, which removes datum offsets
and standing biases. The official NWS and Tidewatch forecasts are already tied to the
gauge by their makers and are shown as issued.
"""
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import time
import urllib.request
from datetime import datetime, timedelta, timezone

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import members as M  # noqa: E402
import pages as P  # noqa: E402
import petss  # noqa: E402
import noaa_override  # noqa: E402
import etss  # noqa: E402
import psurge  # noqa: E402
import stormsense_member as SM  # noqa: E402

TIDEWATCH_URL = os.environ.get("TIDEWATCH_URL", "https://web4.vims.edu/physical/tidewatch/")
NAME = re.compile(r"^([A-Z0-9-]+)_\d{4}_(f36|3d)(?:_(MLLW|NAVD))?(_hc\.html|\.csv)$")
PICKLE = re.compile(r"^[A-Z0-9-]+_\d{4}_f36\.pypickle$")
DROP = re.compile(r"^(NN\d\d|VB\d\d|USNW|USPF|USSV)_")  # retired or non-tidal (Derek, 2026-10-15)
HOURS = 72
STEP = timedelta(minutes=6)
SCHEDULE = "00:40, 06:40, 12:40 and 18:40 UTC"  # keep in step with .github/workflows/ensemble.yml
CONE = "rgba(20, 38, 96, 0.5)"
# key: (label, color, dash, aligned to the gauge)
MEMBERS = {
    "petss": ("P-ETSS (mean)", "#14285a", "Solid", True),
    "stofs2d": ("STOFS-2D-Global", "#ff7f0e", "Solid", True),
    "stofs3d": ("STOFS-3D-Atlantic", "#bcbd22", "Solid", True),
    "cbofs": ("CBOFS", "#8c564b", "Solid", True),
    "nwm_short": ("NWM coastal (short range)", "#9467bd", "Solid", True),
    "nwm_medium": ("NWM coastal (medium range)", "#e377c2", "Solid", True),
    "etss": ("NWS ETSS", "#17becf", "Solid", True),
    "psurge": ("NHC P-Surge (50%)", "#d62728", "Solid", True),
    "nws": ("NWS official", "#000000", "Solid", False),
    "tidewatch": ("VIMS Tidewatch", "#2ca02c", "Solid", False),
    "baseline": ("Tide + persistence", "#7f7f7f", "ShortDash", False),
    "stormsense": ("StormSense model", "#e6007e", "Solid", False),
}
OWN = ["stormsense", "stormsense_p10", "stormsense_p90"]  # appended after the existing CSV columns
UA = {"User-Agent": "stormsense-ensemble (VIMS)"}


def get(url, tries=3):
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60) as r:
                return r.read()
        except Exception:
            if i == tries - 1:
                raise
            time.sleep(5 * (i + 1))


def fetch_tidewatch(raw):
    """Tidewatch f36 and 3d files (charts as templates, CSVs for observations and offsets)."""
    listing = get(TIDEWATCH_URL).decode("utf-8", "replace")
    names = sorted(n for n in set(re.findall(r'href="([^"?/]+)"', listing))
                   if (NAME.match(n) or PICKLE.match(n)) and not DROP.match(n))
    if not names:
        sys.exit("no Tidewatch files found in the listing")
    for n in names:
        prod = "f36" if PICKLE.match(n) else NAME.match(n).group(2)
        os.makedirs(os.path.join(raw, prod), exist_ok=True)
        try:
            data = get(TIDEWATCH_URL + n)
        except Exception as e:
            print(f"  failed {n}: {e}")
            continue
        if data:
            open(os.path.join(raw, prod, n), "wb").write(data)


def series_csv(path, col):
    """(times, values) from a Tidewatch CSV column."""
    try:
        head, rows = P.read_csv(path)
    except OSError:
        return None
    i = head.index(col)
    pairs = [(P.parse_t(r[0]), float(r[i])) for r in rows if len(r) > i and r[i] not in ("", "nan")]
    if len(pairs) < 2:
        return None
    t, v = zip(*pairs)
    return np.array(t), np.array(v)


def sec(times):
    return np.array([t.timestamp() for t in times])


def on_grid(series, grid):
    """Linear interpolation onto grid, NaN outside the series' own span."""
    t, v = series
    s, g = sec(t), sec(grid)
    out = np.interp(g, s, v)
    out[(g < s[0]) | (g > s[-1])] = np.nan
    return out


def alignment(obs, series, t_end):
    """Median of observation minus member over their overlap in the last 24 h (None if under 3 h)."""
    ot, ov = obs
    sel = (ot > t_end - timedelta(hours=24)) & (ot <= t_end)
    t = ot[sel]
    if not len(t):
        return None
    m = on_grid(series, t)
    ok = ~np.isnan(m)
    if ok.sum() < 30:  # 30 six-minute observations = 3 hours
        return None
    return float(np.median(ov[sel][ok] - m[ok]))


def fetch_members(now, stations, jobs=None):
    """{member key: (run time, {station: (times, values ft, native datum)})} for the newest runs."""
    lid_st = {s["nws_lid"]: s["station"] for s in stations.values() if s["nws_lid"]}
    latlon = {s["station"]: (float(s["lat"]), float(s["lon"])) for s in stations.values() if s["lat"]}
    out = {}
    jobs = jobs or {"stofs2d": lambda: M.stofs2d(now, set(lid_st)), "stofs3d": lambda: M.stofs3d(now, set(lid_st)),
            "nwm_short": lambda: M.nwm(now, set(lid_st), "short"),
            "nwm_medium": lambda: M.nwm(now, set(lid_st), "medium"),
            "cbofs": lambda: M.cbofs(now, latlon), "nws": lambda: M.nws(set(lid_st)),
            "psurge": lambda: psurge.fetch(now, stations)}  # only while NHC has advisories out
    for key, job in jobs.items():
        try:
            run, ser, datum, prev = job()
        except Exception as e:
            print(f"{key}: unavailable ({e})")
            continue
        by_st = {}
        for k, (t, v) in ser.items():
            st = lid_st.get(k, k)
            if st in stations:
                d = stations[st]["nws_datum"] if datum == "per-gauge" else datum
                d = "NAVD88" if d.startswith("NAVD") else d
                if k in prev:  # the run before, ahead of this run's start: overlaps the latest observations
                    pt, pv = prev[k]
                    early = pt < t[0]
                    by_st[st] = (t, v, d, (np.concatenate([pt[early], t]), np.concatenate([pv[early], v])))
                else:
                    by_st[st] = (t, v, d, (t, v))
        out[key] = (run, by_st)
        print(f"{key}: run {run:%Y-%m-%d %HZ} at {len(by_st)} stations" if run else f"{key}: {len(by_st)} stations")
    return out


def build_station(st, obs, off, raw_members, pet, pet_cycle, tw_px30):
    """Members on the station's MLLW basis, for t0-24 h to t0+72 h. Returns (t0, grid, cols, notes)."""
    t0 = obs[0][-1]
    grid = np.array([t0 - timedelta(hours=24) + k * STEP for k in range((HOURS + 24) * 10 + 1)])
    cols = {"observed": on_grid(obs, grid)}
    cols["observed"][grid > t0] = np.nan
    notes = {}
    tide = None
    if pet is not None:
        try:
            f = petss.forecast(pet, t0, pet_cycle, hours=HOURS)
            hourly = (np.array([t for t in pet["time"]]), pet["TWL"])
            ok = ~np.isnan(hourly[1])
            adj = alignment(obs, (hourly[0][ok], hourly[1][ok]), t0) or 0.0
            for k, name in (("SSforecast", "petss"), ("low", "petss_p10"), ("high", "petss_p90"),
                            ("astro_tide", "astro_tide")):
                cols[name] = on_grid((f["time"], f[k] + adj), grid)
            tide_all = pet["TIDE"]
            okt = ~np.isnan(tide_all)
            tide = on_grid((pet["time"][okt], tide_all[okt] + adj), grid)
            cols["astro_tide"] = tide.copy()
            notes["petss"] = {"run": pet_cycle.strftime("%Y-%m-%dT%HZ"), "aligned_ft": round(adj, 3)}
        except ValueError as e:
            print(f"  {st}: P-ETSS {e}")
    for key, (run, by_st) in raw_members.items():
        if st not in by_st:
            continue
        t, v, datum, (st_t, st_v) = by_st[st]
        if datum == "NAVD88":
            if off is None:
                continue
            v, st_v = v + off, st_v + off
        aligned = MEMBERS[key][3]
        adj = 0.0
        if aligned:  # the newest run, preceded by the run before, against the last 24 h of observations
            adj = alignment(obs, (st_t, st_v), t0)
            if adj is None:
                if datum == "MSL":  # no fixed tie to MLLW without an overlap
                    continue
                adj = 0.0
        series = on_grid((t, v + adj), grid)
        series[grid < t0 - timedelta(hours=24)] = np.nan
        if np.isnan(series[grid >= t0]).all():
            continue
        cols[key] = series
        notes[key] = {"run": run.strftime("%Y-%m-%dT%H:%MZ") if run else None,
                      "aligned_ft": round(adj, 3) if aligned else None, "native_datum": datum}
    if tw_px30 is not None:
        cols["tidewatch"] = on_grid(tw_px30, grid)
        notes["tidewatch"] = {"run": tw_px30[0][0].strftime("%Y-%m-%dT%H:%MZ"), "aligned_ft": None}
    if tide is not None:  # tide plus the last 3 hours' residual, decaying over a day
        recent = (grid > t0 - timedelta(hours=3)) & (grid <= t0)
        r0 = np.nanmedian(cols["observed"][recent] - tide[recent])
        if not np.isnan(r0):
            dt = (sec(grid) - t0.timestamp()) / 3600.0
            base = tide + r0 * np.exp(-np.clip(dt, 0, None) / 24.0)
            base[grid < t0] = np.nan
            cols["baseline"] = base
            notes["baseline"] = {"residual_ft": round(float(r0), 3), "e_folding_h": 24}
    return t0, grid, cols, notes


def chart_series(grid, cols, shift, t0):
    keep = grid >= t0 - timedelta(hours=24)
    g = grid[keep]
    s = [{"data": P.pts(g, cols["observed"][keep] + shift), "type": "line", "color": "red", "name": "Observed",
          "lineWidth": 3, "zIndex": 5}]
    if "astro_tide" in cols:
        s.append({"data": P.pts(g, cols["astro_tide"][keep] + shift), "type": "line", "color": "blue",
                  "name": "Astronomic Tide", "lineWidth": 1, "dashStyle": "Dot"})
    if "petss_p10" in cols:
        lo, hi = cols["petss_p10"][keep] + shift, cols["petss_p90"][keep] + shift
        s.append({"data": [[P.ms(t), round(float(a), 3), round(float(b), 3)] for t, a, b in zip(g, lo, hi)
                           if not (np.isnan(a) or np.isnan(b))], "type": "arearange", "color": CONE,
                  "fillColor": CONE, "lineWidth": 0, "zIndex": 0, "name": "P-ETSS 10%–90%"})
    for key, (label, color, dash, _) in MEMBERS.items():
        if key in cols:
            s.append({"data": P.pts(g, cols[key][keep] + shift), "type": "line", "color": color, "name": label,
                      "dashStyle": dash, "lineWidth": 2})
    return s


LEVEL_COLS = (["observed", "astro_tide"] + [k for k in MEMBERS if k not in OWN] + ["petss_p10", "petss_p90"]
              + OWN)
HEAD = ["date_time_z"] + [f"{k}_ft" for k in LEVEL_COLS]


def write_outputs(site, st, datum, grid, cols, shift, notes, updated):
    rows = []
    for i, t in enumerate(grid):
        rows.append([P.iso(t)] + [P.fmt(cols[k][i] + shift) if k in cols else "" for k in LEVEL_COLS])
    base = os.path.join(site, "ens72", f"{st}_ens72_{datum}")
    P.write_csv(base + ".csv", HEAD, rows)
    doc = {"station": st, "product": "ens72", "description": "StormSense Ensemble: 72-hour storm tide forecasts "
           "from every available source", "datum": "NAVD88" if datum == "NAVD" else "MLLW", "units": "feet",
           "time_zone": "UTC (Z)", "interval_minutes": 6, "updated": updated, "members": notes, "fields": HEAD,
           "data": [{h: (r[j] if j == 0 else (float(r[j]) if r[j] else None)) for j, h in enumerate(HEAD)} for r in rows]}
    json.dump(doc, open(base + ".json", "w"), separators=(",", ":"))


def main(site="site", work="work"):
    shutil.rmtree(work, ignore_errors=True)
    raw = os.path.join(work, "raw")
    os.makedirs(raw)
    now = datetime.now(timezone.utc)
    updated = now.strftime("%Y-%m-%d %H:%M UTC")
    fetch_tidewatch(raw)
    f36, d3 = os.path.join(raw, "f36"), os.path.join(raw, "3d")
    # CBBT, MNPT, DAHL, WADC: observations straight from NOAA CO-OPS (Tidewatch's are stale or aliased)
    noaa_notes = noaa_override.apply({"f36": f36, "3d": d3}, f36)
    stations = {r["station"]: r for r in csv.DictReader(open(os.path.join(HERE, "stations.csv")))
                if os.path.exists(os.path.join(f36, f"{r['station']}_2018_f36_MLLW.csv"))}
    print(f"{len(stations)} stations with Tidewatch files")
    raw_members = fetch_members(now, stations)
    pet, pet_cycle = {}, None
    try:
        pet_cycle, a, b = petss.newest_cycle(now)
        pet = petss.load(pet_cycle, a, b, {s["nws_lid"] for s in stations.values() if s["nws_lid"]})
        print(f"petss: run {pet_cycle:%Y-%m-%d %HZ} at {len(pet)} gauges")
    except Exception as e:
        print("petss: unavailable", e)
    # VIMS gauges P-ETSS doesn't forecast (BRDH): a nearby station's P-ETSS surge on their own tide
    proxy, _ = noaa_override.petss_proxies(pet)
    # ETSS's bulletin is surge only: it rides on the P-ETSS (NOAA harmonic) tide
    raw_members.update(fetch_members(now, stations, {
        "etss": lambda: etss.fetch(now, stations, etss.tides_from_petss(pet, stations))}))
    pages_dir = os.path.join(work, "pages")
    os.makedirs(pages_dir)
    os.makedirs(os.path.join(site, "ens72"), exist_ok=True)
    summary, raw_dump = {}, {}
    skill = SM.load_skill()
    for st, info in sorted(stations.items()):
        obs = series_csv(os.path.join(d3, f"{st}_2018_3d_MLLW.csv"), "feet")
        if obs is None or now - obs[0][-1] > timedelta(hours=12):
            print(f"  {st}: no recent observations; skipped")
            continue
        off = P.offset(f36, st)
        tw = series_csv(os.path.join(f36, f"{st}_2018_f36_MLLW.csv"), "px30")
        t0, grid, cols, notes = build_station(st, obs, off, raw_members, pet.get(info["nws_lid"]) or proxy.get(st), pet_cycle, tw)
        SM.add(st, grid, t0, cols, notes, skill.get(st))
        members_here = [MEMBERS[k][0] for k in MEMBERS if k in cols]
        guide = f"StormSense Ensemble: {len(members_here)} forecasts; model guidance lined up with the gauge's last 24 hours"
        for datum, template in (("MLLW", f"{st}_2018_f36_MLLW_hc.html"), ("NAVD", f"{st}_2018_f36_hc.html")):
            if datum == "NAVD" and off is None:
                continue
            shift = 0.0 if datum == "MLLW" else -off
            write_outputs(site, st, datum, grid, cols, shift, notes, updated)
            tpath = os.path.join(f36, template)
            if not os.path.exists(tpath) or os.path.getsize(tpath) == 0:
                continue
            html, _ = P.read_page(tpath)
            if "var data = " not in html:
                continue
            P.write_page(os.path.join(pages_dir, template), html, chart_series(grid, cols, shift, t0),
                         title=f"72 Hour StormSense Ensemble for {st}", guidance=guide)
        summary[st] = {"t0": P.iso(t0), "members": notes}
        raw_dump[st] = {k: {"t": [P.iso(t) for t in by[st][0]], "v": [round(float(x), 3) for x in by[st][1]],
                            "datum": by[st][2]} for k, (run, by) in raw_members.items() if st in by}
    style(raw, work, pages_dir)
    place(work, pages_dir, site)
    finish(site, work, summary, raw_dump, raw_members, pet_cycle, updated, now, noaa_notes)


def style(raw, work, pages_dir):
    thresholds = os.path.join(work, "thresholds.json")
    r = subprocess.run([sys.executable, os.path.join(HERE, "make_thresholds.py"), thresholds], capture_output=True, text=True)
    if r.returncode or not os.path.exists(thresholds):
        shutil.copy(os.path.join(HERE, "thresholds.json"), thresholds)
    out = os.path.join(work, "styled")
    r = subprocess.run([sys.executable, os.path.join(HERE, "styler", "style_pages.py"), "--product", "f36",
                        "--pages", pages_dir, "--csv", os.path.join(raw, "f36"), "--pickles", os.path.join(raw, "f36"),
                        "--out", out, "--jsdir=../jslib", "--brand=stormsense", "--thresholds=" + thresholds,
                        "--navdlines=mllw", "--nom30"], capture_output=True, text=True)
    print((r.stdout + r.stderr).strip()[-1500:])
    if r.returncode == 0:
        for n in os.listdir(out):
            shutil.copy(os.path.join(out, n), os.path.join(pages_dir, n))


def place(work, pages_dir, site):
    for n in os.listdir(pages_dir):
        m = re.match(r"^([A-Z0-9-]+)_2018_f36_(MLLW_)?hc\.html$", n)
        if m:
            shutil.copy(os.path.join(pages_dir, n),
                        os.path.join(site, "ens72", f"{m.group(1)}_ens72_{'MLLW' if m.group(2) else 'NAVD'}.html"))
    os.makedirs(os.path.join(site, "jslib"), exist_ok=True)
    for f in os.listdir(os.path.join(HERE, "styler", "jslib")):
        shutil.copy(os.path.join(HERE, "styler", "jslib", f), os.path.join(site, "jslib", f))
    shutil.copy(os.path.join(HERE, "styler", "stormsense_logo.png"), os.path.join(site, "jslib"))
    shutil.copy(os.path.join(HERE, "stormsense_logo.png"), site)


def finish(site, work, summary, raw_dump, raw_members, pet_cycle, updated, now, noaa_notes=None):
    runs = {k: (r.strftime("%Y-%m-%dT%H:%MZ") if r else None) for k, (r, _) in raw_members.items()}
    runs["petss"] = pet_cycle.strftime("%Y-%m-%dT%HZ") if pet_cycle else None
    latest = {"updated": updated, "runs": runs, "stations": summary, "noaa_observations": noaa_notes or {}}
    json.dump(latest, open(os.path.join(site, "latest.json"), "w"), indent=1)
    write_index(site, summary, runs, updated)
    open(os.path.join(site, ".nojekyll"), "w").close()
    # archive: this run's CSV/JSON plus every member's native series, one tarball
    cyc = now.replace(minute=0, second=0, microsecond=0)
    cyc -= timedelta(hours=cyc.hour % 6)
    os.makedirs(os.path.join(work, "archive"), exist_ok=True)
    json.dump({"updated": updated, "runs": runs, "series": raw_dump},
              open(os.path.join(work, "members_native.json"), "w"), separators=(",", ":"))
    tar = os.path.join(work, "archive", f"ensemble_{cyc:%Y-%m-%dT%HZ}.tar.gz")
    with tarfile.open(tar, "w:gz") as t:
        t.add(os.path.join(site, "latest.json"), "latest.json")
        t.add(os.path.join(work, "members_native.json"), "members_native.json")
        for n in sorted(os.listdir(os.path.join(site, "ens72"))):
            if n.endswith((".csv", ".json")):
                t.add(os.path.join(site, "ens72", n), f"ens72/{n}")
    print(f"{len(summary)} stations; archive {os.path.basename(tar)} ({os.path.getsize(tar) // 1024} KB)")


def write_index(site, summary, runs, updated):
    import html as H
    rows = []
    for st in sorted(summary):
        links = []
        for d in ("MLLW", "NAVD"):
            parts = [f'<a href="ens72/{st}_ens72_{d}.{e}">{lab}</a>' for e, lab in (("html", "chart"), ("csv", "csv"), ("json", "json"))
                     if os.path.exists(os.path.join(site, "ens72", f"{st}_ens72_{d}.{e}"))]
            links.append(f"<td>{' '.join(parts)}</td>")
        n = len(summary[st]["members"])
        rows.append(f"<tr><th>{st}</th>{''.join(links)}<td>{n}</td></tr>")
    run_list = "".join(f"<li>{H.escape(MEMBERS[k][0])}: {v or 'unavailable'}</li>" for k, v in runs.items() if k in MEMBERS)
    page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>StormSense Ensemble</title>
<style>
body{{font-family:Arial,'Liberation Sans',sans-serif;margin:0 auto;padding:16px;max-width:1000px;color:#1b2533;background:#fff}}
header img{{width:min(520px,100%);height:auto;display:block}}
h1{{font-size:1.5rem;margin:12px 0 4px;color:#14285a}}
p,li{{margin:6px 0;line-height:1.45}}
.disc{{color:#555;font-size:.92rem}}
.wrap{{overflow-x:auto;margin-top:14px}}
table{{border-collapse:collapse;font-size:.92rem}}
th,td{{border:1px solid #d5dbe3;padding:4px 8px;text-align:left;white-space:nowrap}}
thead th{{background:#eef6fa;text-align:center}}
a{{color:#0d6e8a}}
</style></head>
<body><header><img src="stormsense_logo.png" alt="StormSense"></header>
<h1>StormSense Ensemble Forecasts</h1>
<p>Every available 72-hour storm tide forecast for each sensor on one chart, updated every 6 hours (updated {updated};
about {SCHEDULE}). Links keep the same address and are overwritten on each update; each run is also archived.</p>
<p class="disc">Disclaimer: StormSense&trade; Predictive Plots are computer-generated operational storm tide forecasts;
the user assumes any and all risk.</p>
<p class="disc">Model guidance (P-ETSS, STOFS-2D, STOFS-3D, CBOFS, NWM, ETSS, P-Surge) is lined up with each gauge: the median
difference from the observations over the last 24 hours is added to the whole run, removing datum offsets and standing
bias. The official NWS and VIMS Tidewatch forecasts are shown as issued. Tide + persistence is the astronomical tide plus
the last 3 hours' residual, fading over a day: the baseline a forecast should beat. Click a legend entry to hide or show it.</p>
<p class="disc"><a href="skill/">Skill of each forecast</a> against the observations that came in after it, and which
guidance drives the primary StormSense forecast at each station.</p>
<p class="disc">Newest runs used (UTC):</p><ul class="disc">{run_list}</ul>
<div class="wrap"><table><thead><tr><th>Station</th><th>MLLW</th><th>NAVD88</th><th>Forecasts</th></tr></thead>
<tbody>
{chr(10).join(rows)}
</tbody></table></div></body></html>
"""
    open(os.path.join(site, "index.html"), "w").write(page)


if __name__ == "__main__":
    main(*sys.argv[1:])
