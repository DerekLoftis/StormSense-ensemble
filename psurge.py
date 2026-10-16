"""NHC P-Surge (Probabilistic Hurricane Storm Surge, NWS MDL) station guidance.

P-Surge runs only while NHC issues advisories for a tropical cyclone that threatens the US
coast (normally 00, 06, 12 and 18Z, posted about 3 h later). Each run goes to NOMADS:

    https://nomads.ncep.noaa.gov/pub/data/nccf/com/psurge/prod/psurge.YYYYMMDD/
        psurge.tHHz.<storm>.e10|e50|e90.conus.shef        storm surge + tide at gauges (SHEF)
        psurge.tYYYYMMDDHHz.<storm>_eNN_inc_dat.h102.conus_625m.grib2   gridded, above NAVD88
        psurge_YYYYMMDDHH_<storm>.meta, YYYYMMDDHH.go, <storm>_YYYYMMDDHH_wlevel.dat
    (<storm> is e.g. al092026; NOMADS keeps about a week of date folders.)

The SHEF bulletins give hourly storm surge + tide in feet above MLLW, 0-102 h, at about 100
gauges keyed by NWS id (PE code HMIFZ), for the 10%, 50% and 90% exceedance levels. fetch()
reads them. P-Surge has no ensemble mean, so the default level is e50 (the median); pass
level="e10" for the 10%-exceedance height NHC uses for surge watches and warnings.

fetch_grid() reads the 625 m GRIB2 grid instead (hourly maximum over each hour, metres above
NAVD88, eNN_inc_dat) at each station's nearest wet cell within 3 km: it covers gauges that
have no SHEF point. It needs pygrib (pip install pygrib; the wheel bundles ecCodes) and
downloads about 27 MB per run.

Both return (run time, {station: (times, values ft)}, datum, {same for the run before}) and
(None, {}, datum, {}) when no run was posted in the last 24 hours, which is the normal case.
"""
import re
from datetime import datetime, timedelta, timezone

import numpy as np

import members as M

BASE = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/psurge/prod/"
MAX_AGE = timedelta(hours=24)


def runs(now, max_age=MAX_AGE):
    """[(cycle, storm id)] posted in the last max_age, newest first, from the date folders.
    A folder that does not exist (404) means no run that day."""
    found = set()
    for d in (now, now - timedelta(days=1)):
        try:
            html = M.get(f"{BASE}psurge.{d:%Y%m%d}/", tries=1, timeout=60).decode("latin-1")
        except Exception:
            continue
        for hh, storm in re.findall(r'href="psurge\.t(\d\d)z\.(\w+)\.e50\.conus\.shef"', html):
            c = datetime(d.year, d.month, d.day, int(hh), tzinfo=timezone.utc)
            if now - max_age <= c <= now:
                found.add((c, storm))
    return sorted(found, reverse=True)


def _by_station(series, stations):
    lid_st = {s["nws_lid"]: st for st, s in stations.items() if s.get("nws_lid")}
    return {lid_st[k]: v for k, v in series.items() if k in lid_st}


def _merge(a, b):
    """Two storms in the same cycle: keep, per station and hour, the higher water level."""
    out = dict(a)
    for st, (t, v) in b.items():
        if st not in out:
            out[st] = (t, v)
            continue
        ta, va = out[st]
        tt = np.array(sorted(set(ta) | set(t)))
        x = np.full(len(tt), np.nan)
        for ts, vs in ((ta, va), (t, v)):
            idx = np.searchsorted(tt, ts)
            x[idx] = np.fmax(x[idx], vs)
        out[st] = (tt, x)
    return out


def read_shef(text, stations):
    """{station: (times, ft MLLW)} from one P-Surge SHEF bulletin; missing (dry) hours dropped."""
    lids = {s["nws_lid"] for s in stations.values() if s.get("nws_lid")}
    return _by_station(M.parse_shef(text, lids, pe_prefix="HM"), stations)


def fetch(now, stations, level="e50"):
    """Newest P-Surge run at the gauges in its SHEF bulletin, ft above MLLW."""
    found = runs(now)
    if not found:
        return None, {}, "MLLW", {}
    cycles = sorted({c for c, _ in found}, reverse=True)
    res = []
    for c in cycles[:2]:
        ser = {}
        for cc, storm in found:
            if cc != c:
                continue
            url = f"{BASE}psurge.{c:%Y%m%d}/psurge.t{c:%H}z.{storm}.{level}.conus.shef"
            try:
                ser = _merge(ser, read_shef(M.get(url, tries=2).decode("latin-1"), stations))
            except Exception as e:
                print(f"psurge: {url} unavailable ({e})")
        res.append((c, ser))
    res = [r for r in res if r[1]]
    if not res:
        return None, {}, "MLLW", {}
    return res[0][0], res[0][1], "MLLW", (res[1][1] if len(res) > 1 else {})


def read_grib(data, stations, max_km=3.0):
    """{station: (times, ft NAVD88)} from an eNN_inc_dat GRIB2 file: the hourly maximum at the
    station's nearest grid cell within max_km that is wet (unmasked) nearly every hour.
    Times are the end of each hour."""
    import pygrib
    import tempfile
    with tempfile.NamedTemporaryFile(suffix=".grib2") as f:
        f.write(data)
        f.flush()
        g = pygrib.open(f.name)
        msgs = list(g)
        lat, lon = msgs[0].latlons()
        pts = {}
        for st, s in stations.items():
            try:
                la, lo = float(s["lat"]), float(s["lon"])
            except (TypeError, ValueError, KeyError):
                continue
            # a window of +-0.05 deg (5 km) around the station keeps this cheap
            win = (np.abs(lat - la) < 0.05) & (np.abs(lon - lo) < 0.06)
            if not win.any():
                continue
            ii, jj = np.nonzero(win)
            dist = np.hypot((lat[ii, jj] - la) * 111.0, (lon[ii, jj] - lo) * 111.0 * np.cos(np.radians(la)))
            pts[st] = (ii, jj, dist)
        times, win_vals = [], {st: [] for st in pts}
        for m in msgs:
            v = np.ma.filled(np.ma.masked_greater(m.values.astype(float), 9000), np.nan)
            times.append(m.analDate.replace(tzinfo=timezone.utc) + timedelta(hours=int(m["endStep"])))
            for st, (ii, jj, _) in pts.items():
                win_vals[st].append(v[ii, jj])
        g.close()
    t = np.array(times)
    out = {}
    for st, rows in win_vals.items():
        a = np.array(rows)  # hours x window cells
        dist = pts[st][2]
        frac = np.where(dist <= max_km, np.mean(~np.isnan(a), axis=0), 0.0)  # share of hours wet
        if frac.max() == 0:
            continue
        # the nearest cell that is wet nearly all the time (water, not intermittently flooded
        # land), else the cell wet most often
        good = frac >= 0.9 * frac.max()
        x = a[:, np.argmin(np.where(good, dist, np.inf))] / 0.3048
        ok = ~np.isnan(x)
        if ok.sum() >= 2:
            out[st] = (t[ok], x[ok])
    return out


def fetch_grid(now, stations, level="e50", max_km=3.0):
    """Newest P-Surge run from the 625 m grid (hourly maxima), ft above NAVD88."""
    found = runs(now)
    if not found:
        return None, {}, "NAVD88", {}
    c = found[0][0]
    ser = {}
    for cc, storm in found:
        if cc != c:
            continue
        url = (f"{BASE}psurge.{c:%Y%m%d}/psurge.t{c:%Y%m%d%H}z.{storm}_{level}_inc_dat"
               f".h102.conus_625m.grib2")
        try:
            ser = _merge(ser, read_grib(M.get(url, tries=2, timeout=300), stations, max_km))
        except Exception as e:
            print(f"psurge grid: {url} unavailable ({e})")
    return (c if ser else None), ser, "NAVD88", {}
