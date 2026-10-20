"""The StormSense model as an ensemble member (`stormsense`).

The model lives in DerekLoftis/StormSense_forecastmodel, which ensemble.yml checks out at
forecastmodel/ (SSMODEL_DIR overrides). For each station it takes the members this run has
already lined up with the gauge, StormSense's own harmonic tide (forecastmodel/harmonics,
falling back to this run's P-ETSS tide where no fit exists yet) and the skill scoring's recent
RMSE per member and lead band (skill/summary.csv from the last published run; equal weights
without it). A failure only drops this member; the rest of the run carries on.
"""
import csv
import io
import os
import sys
import urllib.request

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MODEL_DIR = os.environ.get("SSMODEL_DIR", os.path.join(HERE, "forecastmodel"))
SKILL_URL = os.environ.get("SKILL_SUMMARY_URL", "https://derekloftis.github.io/StormSense-ensemble/skill/summary.csv")
UA = {"User-Agent": "stormsense-ensemble (VIMS)"}


def _model():
    if MODEL_DIR not in sys.path:
        sys.path.insert(0, MODEL_DIR)
    from ssmodel import surge, tide
    return surge, tide


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


def add(st, grid, t0, cols, notes, skill=None):
    """Add stormsense, stormsense_p10 and stormsense_p90 to cols (ft MLLW) and a note."""
    try:
        surge, tide = _model()
        if os.path.exists(os.path.join(MODEL_DIR, "harmonics", f"{st}.json")):
            own = tide.predict(tide.load(st, os.path.join(MODEL_DIR, "harmonics")), list(grid))
            source = "StormSense harmonic fit"
        elif "astro_tide" in cols:
            own, source = cols["astro_tide"], "P-ETSS tide (no StormSense fit yet)"
        else:
            return
        out, note = surge.forecast(grid, t0, cols, own, skill or None)
        if np.isnan(out["stormsense"]).all():
            return
        cols.update(out)
        notes["stormsense"] = dict(note, tide=source, run=t0.strftime("%Y-%m-%dT%H:%MZ"), aligned_ft=None)
    except Exception as e:
        print(f"  {st}: stormsense unavailable ({e})")
