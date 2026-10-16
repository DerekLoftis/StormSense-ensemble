"""NWS ETSS (Extra-Tropical Storm Surge, MDL; GFS-driven, deterministic) station guidance.

NOMADS refuses ETSS from GitHub Actions (HTTP 403 on .../com/etss/prod/ and on every file
under it), so this reads the East Coast text bulletin the NWS sends out each cycle
(FQUS23 KWNO, AWIPS id MRPSSE), which tgftp keeps the newest copy of:

    https://tgftp.nws.noaa.gov/data/raw/fq/fqus23.kwno.mrp.sse.txt

Four cycles a day (00, 06, 12, 18Z), posted about 4.5-5 h after the cycle. For each of about
75 gauges (named, not numbered) it gives the storm surge only, in tenths of a foot, for hour 0
(at the end of the station's name line) and hours 1-96 (four lines of 24). "NOT VALID FOR
TROPICAL STORMS". The surge has no datum: it is the water level change from wind and pressure.
Tidewatch adds it to its own tide and anomaly; here a tide series (ft MLLW) can be passed in,
for example P-ETSS's harmonic tide (tides_from_petss), to get a water level on MLLW that
ensemble.py can line up with the gauge like the other members.
"""
import re
from datetime import datetime, timedelta, timezone

import numpy as np

import members as M

URL = "https://tgftp.nws.noaa.gov/data/raw/fq/fqus23.kwno.mrp.sse.txt"
MAX_AGE = timedelta(hours=24)

# bulletin name -> NWS id in stations.csv. The ETSS points are NOS gauges; where the bulletin
# point is the same water a few km from one of ours, it is mapped to ours (noted).
NAMES = {
    "HAMPTON ROADS/ VA": "SWPV2",       # 8638610 Sewells Point
    "KIPTOPEKE BEACH/ VA": "KPTV2",     # 8632200
    "LEWISETTA/ VA": "LWTV2",           # 8635750
    "CHESAPEAKE BBT/ VA": "CHBV2",      # on the bridge-tunnel, as ours (8638901)
    "BISHOPS HEAD/ MD": "BISM2",        # 8571421
    "SOLOMONS/ MD": "SLIM2",            # 8577330
    "WASHINGTON/ DC": "WASD2",          # 8594900
    "WACHAPREAGUE/ VA": "WAHV2",        # 8631044
    "GLOUCESTER POINT/ VA": "YKTV2",    # 8637624, across the York River from USCG Yorktown 8637689
    "OLD PT. COMFORT/ VA": "FTMV2",     # Old Point Comfort = the USGS Fort Monroe gauge (USFM)
    "POTOMAC RV ALEX/ VA": "AXTV2",     # Potomac River at Alexandria
}


def read_bulletin(text, now):
    """(cycle, {bulletin name: (times, surge ft)}) from an MRPSSE bulletin."""
    lines = text.splitlines()
    head = next(l for l in lines if re.match(r"^FQUS\d\d KWNO \d{6}", l))
    day, hh, mm = (int(x) for x in re.findall(r"(\d\d)(\d\d)(\d\d)$", head.strip())[0])
    issued = None
    for back in range(0, 3):  # the header has only the day of the month
        y, m = now.year, now.month - back
        while m < 1:
            y, m = y - 1, m + 12
        try:
            t = datetime(y, m, day, hh, mm, tzinfo=timezone.utc)
        except ValueError:
            continue
        if t <= now + timedelta(hours=1):
            issued = t
            break
    if issued is None:
        raise ValueError("cannot date bulletin header " + head)
    label = next(l for l in lines if re.match(r"^\s*\d\dZ\s+\d\dZ", l))
    first = int(label.split()[0][:2])  # hour of the first column = cycle + 1 h
    cycle = (issued - timedelta(hours=(issued.hour - (first - 1)) % 24)).replace(minute=0)
    out, name, vals = {}, None, []

    def done():
        if name and len(vals) >= 2:
            v = np.array(vals, dtype=float) / 10.0
            t = np.array([cycle + timedelta(hours=i) for i in range(len(v))])
            ok = ~np.isnan(v)
            out[name] = (t[ok], v[ok])

    for line in lines:
        m = re.match(r"^ ([A-Z][A-Z0-9 .'/&-]*?/ ?[A-Z]{2})\s+(-?\d+)\s*$", line)
        if m:
            done()
            name, vals = re.sub(r"\s+", " ", m.group(1).strip()), [float(m.group(2))]
            continue
        if name and re.fullmatch(r"[\s\d-]+", line) and line.strip():
            body = line.rstrip()
            for i in range(0, len(body), 3):  # 3 columns per value; "-9-10-11" runs together
                tok = body[i:i + 3].strip()
                vals.append(float(tok) if re.fullmatch(r"-?\d+", tok) else np.nan)
    done()
    return cycle, out


def tides_from_petss(pet, stations):
    """{station: (times, tide ft MLLW)} from petss.load()'s {NWS id: station dict}."""
    lid_st = {s["nws_lid"]: st for st, s in stations.items() if s.get("nws_lid")}
    out = {}
    for lid, d in (pet or {}).items():
        ok = ~np.isnan(d["TIDE"])
        if lid in lid_st and ok.sum() >= 2:
            out[lid_st[lid]] = (d["time"][ok], d["TIDE"][ok])
    return out


def fetch(now, stations, tide=None):
    """Newest ETSS run. With tide ({station: (times, ft MLLW)}), returns surge + tide in ft MLLW
    for the stations that have a tide over the forecast; without it, the surge alone with
    datum "surge"."""
    cycle, by_name = read_bulletin(M.get(URL, tries=2, timeout=60).decode("latin-1"), now)
    if now - cycle > MAX_AGE:
        raise RuntimeError(f"newest ETSS bulletin is from {cycle:%Y-%m-%d %HZ}")
    lid_st = {s["nws_lid"]: st for st, s in stations.items() if s.get("nws_lid")}
    surge = {lid_st[NAMES[n]]: v for n, v in by_name.items() if n in NAMES and NAMES[n] in lid_st}
    if tide is None:
        return cycle, surge, "surge", {}
    sec = lambda a: np.array([x.timestamp() for x in a])
    out = {}
    for st, (t, s) in surge.items():
        if st not in tide:
            continue
        tt, tv = tide[st]
        inside = (sec(t) >= sec(tt[:1])[0]) & (sec(t) <= sec(tt[-1:])[0])
        if inside.sum() >= 2:
            out[st] = (t[inside], s[inside] + np.interp(sec(t[inside]), sec(tt), tv))
    return cycle, out, "MLLW", {}
