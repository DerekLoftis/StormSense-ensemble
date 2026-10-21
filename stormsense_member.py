"""The StormSense model as an ensemble member (`stormsense`).

The model lives in DerekLoftis/StormSense_forecastmodel, which ensemble.yml checks out at
forecastmodel/ (SSMODEL_DIR overrides). For each station it takes the members this run has
already lined up with the gauge, StormSense's own harmonic tide (forecastmodel/harmonics,
falling back to this run's P-ETSS tide where no fit exists yet) and the skill scoring's recent
RMSE per member and lead band (skill/summary.csv from the last published run; equal weights
without it). That blend is the predictor; the site-specific corrector then removes the bias the
blend has shown at this gauge over the last week, learned by re-running it on the archived runs
(see forecastmodel/ssmodel/corrector.py). `stormsense_raw` keeps the uncorrected blend so the
skill scoring can check the corrector. A failure only drops this member; the run carries on.
"""
import csv
import io
import os
import sys
import urllib.request
from datetime import datetime, timezone

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MODEL_DIR = os.environ.get("SSMODEL_DIR", os.path.join(HERE, "forecastmodel"))
SKILL_URL = os.environ.get("SKILL_SUMMARY_URL", "https://derekloftis.github.io/StormSense-ensemble/skill/summary.csv")
UA = {"User-Agent": "stormsense-ensemble (VIMS)"}


def _model():
    if MODEL_DIR not in sys.path:
        sys.path.insert(0, MODEL_DIR)
    from ssmodel import archive, corrector, surge, tide
    return archive, corrector, surge, tide


def load_skill(url=SKILL_URL):
    """{station: {member: {band: recent RMSE ft}}} from the skill summary; {} if unavailable."""
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=30) as r:
            text = r.read().decode()
    except Exception as e:
        print(f"stormsense: no skill summary ({e}); equal weights")
        return {}
    out = {}
    for row in csv.DictReader(io.StringIO(text)):
        try:
            rmse = float(row["rmse_ew"])
        except (KeyError, ValueError):
            continue
        out.setdefault(row["station"], {}).setdefault(row["member"], {})[row["band"]] = rmse
    return out


def load_past(work, now=None):
    """(runs, truth): the last week of archived ensemble runs and their stitched observations,
    for the corrector (None if unavailable)."""
    try:
        archive = _model()[0]
        paths = archive.download_recent(os.path.join(work, "past_runs"), now or datetime.now(timezone.utc))
        runs = archive.read_runs(paths)
        print(f"stormsense: {len(runs)} archived runs for the corrector")
        return runs, archive.truth_series(runs)
    except Exception as e:
        print(f"stormsense: no archived runs ({e}); no correction")
        return None


def _tide(st, grid, cols):
    tide = _model()[3]
    folder = os.path.join(MODEL_DIR, "harmonics")
    if os.path.exists(os.path.join(folder, f"{st}.json")):
        return tide.predict(tide.load(st, folder), list(grid)), "StormSense harmonic fit"
    if "astro_tide" in cols:
        return cols["astro_tide"], "P-ETSS tide (no StormSense fit yet)"
    return None, None


def add(st, grid, t0, cols, notes, skill=None, past=None):
    """Add stormsense (corrected), its 10-90% band and stormsense_raw to cols (ft MLLW), and a note."""
    try:
        _, corrector, surge, _ = _model()
        own, source = _tide(st, grid, cols)
        if own is None:
            return
        out, note = surge.forecast(grid, t0, cols, own, skill or None)
        if np.isnan(out["stormsense"]).all():
            return
        raw = out["stormsense"].copy()
        info = {"applied": False, "runs": 0}
        if past:
            def predict(g, t, c):
                td, _ = _tide(st, g, c)
                return (None if td is None else (surge.forecast(g, t, c, td, skill or None)[0]["stormsense"], td))
            runs, truth = past
            here = dict(truth.get(st, {}))  # this run's observations win
            here.update((t, v) for t, v in zip(grid, cols["observed"]) if not np.isnan(v))
            coefs, info = corrector.correction_for(st, runs, {st: here}, t0, predict)
            corrector.apply(grid, t0, out, own, coefs)
        cols.update(out)
        cols["stormsense_raw"] = raw
        notes["stormsense"] = dict(note, tide=source, corrector=info, run=t0.strftime("%Y-%m-%dT%H:%MZ"),
                                   aligned_ft=None)
    except Exception as e:
        print(f"  {st}: stormsense unavailable ({e})")
