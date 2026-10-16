"""Fetch the newest run of each ensemble member as point time series.

Each fetcher returns (run time, {LID or station: (times, values ft)}, datum, the same
for the run before) where
datum says what the values are referenced to: "MLLW", "NAVD88", "MSL" or "per-gauge"
(NWS: each gauge's own datum, from stations.csv). ensemble.py puts them all on the
station's MLLW basis and lines them up with the gauge.

    STOFS-2D-Global   SHEF, ft MLLW, 30 min, 4 runs a day, 180 h    (NOAA AWS bucket)
    STOFS-3D-Atlantic SHEF, ft MLLW, 12Z only, 96 h                 (NOAA AWS bucket)
    NWM coastal       SHEF, ft NAVD88, short range 18 h / medium range 10 days (NOAA AWS bucket)
    CBOFS             station netCDF, m above MSL, 48 h              (NOAA AWS bucket)
    NWS official      NWPS stage forecast JSON per gauge             (api.water.noaa.gov)
P-ETSS comes from petss.py and Tidewatch from its own CSVs.
"""
import io
import json
import re
import time
import urllib.request
from datetime import datetime, timedelta, timezone

import numpy as np

UA = {"User-Agent": "stormsense-ensemble (VIMS)"}
S3 = "https://{bucket}.s3.amazonaws.com/{key}"
NWPS = "https://api.water.noaa.gov/nwps/v1/gauges/{lid}/stageflow"


def get(url, tries=2, timeout=120):
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout) as r:
                return r.read()
        except Exception:
            if i == tries - 1:
                raise
            time.sleep(3)


def cycles(now, hours=(0, 6, 12, 18), back=36):
    """Candidate run times, newest first, over the last `back` hours."""
    c = now.replace(minute=0, second=0, microsecond=0)
    out = []
    for _ in range(back + 1):
        if c.hour in hours:
            out.append(c)
        c -= timedelta(hours=1)
    return out


def newest(bucket, key, now, hours=(0, 6, 12, 18), runs=2):
    """[(run time, bytes)] for the newest `runs` runs of an S3 key template, newest first.
    The run before the newest one is kept to line the newest up with the gauge: it
    overlaps the latest observations, which the newest run often starts after."""
    out = []
    for c in cycles(now, hours, back=60):
        k = key.format(ymd=c.strftime("%Y%m%d"), hh=c.strftime("%H"), y=c.strftime("%Y"), m=c.strftime("%m"),
                       d=c.strftime("%d"))
        try:
            out.append((c, get(S3.format(bucket=bucket, key=k), tries=1)))
        except Exception:
            continue
        if len(out) == runs:
            break
    if not out:
        raise RuntimeError(f"no run of {bucket}/{key} in the last 60 hours")
    return out


def _shef_runs(found, lids):
    """(newest run, {LID: series}, {LID: series of the run before})."""
    ser = [parse_shef(b.decode("latin-1"), lids) for _, b in found]
    return found[0][0], ser[0], (ser[1] if len(ser) > 1 else {})


def parse_shef(text, lids, pe_prefix="HM"):
    """{LID: (times, values)} from a SHEF .E bulletin, first water-level block per gauge."""
    out, cur = {}, None
    for line in text.splitlines():
        if line.startswith(".E "):
            cur = None
            parts = line[3:].strip().split("/")
            head = parts[0].split()
            lid, date = head[0], head[1]
            if lid not in lids or lid in out:
                continue
            dh = di = pe = None
            vals = []
            for tok in head[3:] + [t.strip() for t in parts[1:]]:
                if not tok:
                    continue
                if tok.startswith("DH"):
                    dh = tok[2:]
                elif tok.startswith("DI"):
                    di = tok[2:]
                elif tok.startswith("DC"):
                    pass
                elif re.fullmatch(r"[A-Z]{2}[A-Z0-9]{3,5}", tok):
                    pe = tok
                else:
                    vals.append(tok)
            if not pe or not pe.startswith(pe_prefix) or not dh or not di:
                continue
            t0 = datetime(int(date[:4]), int(date[4:6]), int(date[6:8]), int(dh[:2]), int(dh[2:4] or 0),
                          tzinfo=timezone.utc)
            step = timedelta(minutes=int(di[1:])) if di[0] == "N" else timedelta(hours=int(di[1:]))
            cur = {"t0": t0, "step": step, "vals": vals}
            out[lid] = cur
        elif re.match(r"^\.E\d", line) and cur is not None:
            cur["vals"] += [t.strip() for t in line.split(None, 1)[1].split("/") if t.strip()]
    series = {}
    for lid, b in out.items():
        v = np.array([float(x) if re.fullmatch(r"-?\d+(\.\d+)?", x) else np.nan for x in b["vals"]])
        v[v < -900] = np.nan
        t = np.array([b["t0"] + i * b["step"] for i in range(len(v))])
        ok = ~np.isnan(v)
        if ok.sum() >= 2:
            series[lid] = (t[ok], v[ok])
    return series


def stofs2d(now, lids):
    run, ser, prev = _shef_runs(newest("noaa-gestofs-pds", "stofs_2d_glo.{ymd}/stofs_2d_glo.t{hh}z.points.cwl.shef",
                                       now), lids)
    return run, ser, "MLLW", prev


def stofs3d(now, lids):
    run, ser, prev = _shef_runs(newest("noaa-nos-stofs3d-pds",
                                       "STOFS-3D-Atl/stofs_3d_atl.{ymd}/stofs_3d_atl.t{hh}z.points.cwl.shef",
                                       now, hours=(12,)), lids)
    return run, ser, "MLLW", prev


def nwm(now, lids, rng="short"):
    key = {"short": "nwm.{ymd}/short_range_coastal_atlgulf/nwm.t{hh}z.short_range_coastal.total_water.atlgulf.shef",
           "medium": "nwm.{ymd}/medium_range_coastal_atlgulf_mem1/nwm.t{hh}z.medium_range_coastal.total_water.atlgulf.shef"}[rng]
    run, ser, prev = _shef_runs(newest("noaa-nwm-pds", key, now), lids)
    return run, ser, "NAVD88", prev


def cbofs(now, latlon, max_km=3.0):
    """CBOFS station zeta at the station nearest each gauge (within max_km)."""
    import netCDF4
    found = newest("noaa-nos-ofs-pds", "cbofs/netcdf/{y}/{m}/{d}/cbofs.t{hh}z.{ymd}.stations.forecast.nc", now)
    res = []
    for _, data in found:
        with netCDF4.Dataset("cbofs.nc", memory=data) as d:
            t = netCDF4.num2date(d["ocean_time"][:], d["ocean_time"].units, only_use_cftime_datetimes=False)
            t = np.array([datetime(x.year, x.month, x.day, x.hour, x.minute, tzinfo=timezone.utc) for x in t])
            zeta = d["zeta"][:].filled(np.nan)
            lon, lat = d["lon_rho"][:].filled(np.nan), d["lat_rho"][:].filled(np.nan)
        out = {}
        for st, (la, lo) in latlon.items():
            dist = np.hypot((lat - la) * 111.0, (lon - lo) * 111.0 * np.cos(np.radians(la)))
            j = int(np.nanargmin(dist))
            if dist[j] <= max_km:
                v = zeta[:, j] / 0.3048
                ok = ~np.isnan(v)
                out[st] = (t[ok], v[ok])
        res.append(out)
    return found[0][0], res[0], "MSL", (res[1] if len(res) > 1 else {})


def nws(lids):
    """Official NWS forecasts, by gauge; run time is the newest issuance."""
    out, issued = {}, []
    for lid in sorted(lids):
        try:
            fc = json.loads(get(NWPS.format(lid=lid), tries=1, timeout=30)).get("forecast") or {}
        except Exception:
            continue
        pts = [p for p in fc.get("data", []) if p.get("primary") not in (None, -999)]
        if len(pts) < 2:
            continue
        t = np.array([datetime.strptime(p["validTime"], "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc) for p in pts])
        out[lid] = (t, np.array([float(p["primary"]) for p in pts]))
        if fc.get("issuedTime"):
            issued.append(datetime.strptime(fc["issuedTime"][:19], "%Y-%m-%dT%H:%M:%S").replace(tzinfo=timezone.utc))
    return (max(issued) if issued else None), out, "per-gauge", {}
