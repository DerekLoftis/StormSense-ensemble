# New members: P-Surge and ETSS

Both modules use the members.py return shape:
`(run time or None, {station: (times, values ft)}, datum, {station: series of the run before})`.
Both return station IDs as keys, so `lid_st.get(k, k)` in `ensemble.fetch_members` passes them through.

## psurge.py: NHC P-Surge

- **Source:** NOMADS `https://nomads.ncep.noaa.gov/pub/data/nccf/com/psurge/prod/psurge.YYYYMMDD/`.
  The fetcher reads the date-folder listing (a missing day returns 404) to find
  `psurge.tHHz.<storm>.e10|e50|e90.conus.shef`, for example `psurge.t00z.al092026.e50.conus.shef`.
  The `prod/` root listing times out from Actions, but the date folders and the files list fine.
- **Runs:** only while NHC has advisories for a US-threatening storm, at 00/06/12/18Z, posted
  about 3 h after the cycle.
- **`fetch(now, stations, level="e50")`:** hourly storm surge + tide, ft above **MLLW**, 0-102 h,
  keyed by NWS id. The PE code is HMIFZD1/D5/D9 for e10/e50/e90, and `members.parse_shef` reads it.
  - Coverage: about 100 gauges per run. 17 of ours: BISH CBBT DAHL KIPT LEWI MNPT NSWV SLIM SWPT USHR USJT USLI USPT WACH WADC WMPT YRCG.
  - Hours when a gauge's cell is dry are `M` and get dropped.
  - Default level: P-Surge publishes no mean, so the default is e50, the median.
    `level="e10"` gives NHC's 10%-exceedance height (warning-level guidance). Calling it with e10 and e90 too would give a band like P-ETSS's.
  - No run in the last 24 h returns `(None, {}, "MLLW", {})` without printing anything. This is the normal case.
  - Two storms in the same cycle are combined by taking the higher value each hour.
  - The run before is returned as prev.
- **`fetch_grid(now, stations, level="e50")` (optional):** fills in the gauges with no SHEF point.
  - Reads `psurge.tYYYYMMDDHHz.<storm>_eNN_inc_dat.h102.conus_625m.grib2`, about 27 MB: 102 hourly maxima in metres above **NAVD88**, on a 625 m Lambert grid.
  - Takes the nearest cell within 3 km that is wet in at least 90% as many hours as the wettest nearby cell.
  - Datum "NAVD88"; covered 38 of our stations in the test.
  - Needs `pygrib` (add `pygrib` to requirements.txt; its manylinux wheel bundles ecCodes).
  - Takes about 5 min per run to decode.
  - Use it as a separate member key (for example `psurge_grid`), because its datum differs from the SHEF output.

## etss.py: ETSS (deterministic, GFS)

- **NOMADS is blocked from Actions.** Every path under `nomads.ncep.noaa.gov/pub/data/nccf/com/etss/` returns HTTP 403, including the `prod/` listing, `etss/v2.3/` and a file URL.
  `ftp.ncep.noaa.gov` and `www.ftp.ncep.noaa.gov` time out (curl 000).
- **Source used:** the NWS text bulletin FQUS23 KWNO / MRPSSE, "GFS BASED STORM SURGE", East Coast:
  `https://tgftp.nws.noaa.gov/data/raw/fq/fqus23.kwno.mrp.sse.txt`.
  - It holds only the newest copy, so prev is always `{}`.
  - 4 cycles a day, posted about 4h50m after the cycle (the 00Z run was posted at 04:49).
- **Format:** 74 named points. Surge only, in tenths of a foot, for hour 0 plus hours 1-96, hourly. The text says "NOT VALID FOR TROPICAL STORMS".
  - `etss.NAMES` maps bulletin names to our NWS ids. This gives 11 stations: SWPT (Hampton Roads), KIPT, LEWI, CBBT, BISH, SLIM, WADC, WACH, YRCG (Gloucester Pt), USFM (Old Pt Comfort), AXTV (Potomac Rv Alex).
- **`fetch(now, stations, tide=None)`:**
  - With `tide={station: (times, ft MLLW)}`, it returns surge + tide in ft **MLLW**. `etss.tides_from_petss(pet, stations)` builds that from `petss.load()`.
  - Without a tide, it returns the surge alone with datum `"surge"`, which ensemble.py doesn't handle yet.
  - It raises an error if the newest bulletin is more than 24 h old.
  - Integration note: `fetch_members` currently runs before `petss.load`. Either move the P-ETSS load earlier, or call `etss.fetch` after it.
- **Possible other source (untested):** MDL's ETSS 2.0 page serves station JSON (`ts, pred, ss, anom, obs, twl`) from
  `https://slosh.nws.noaa.gov/etsurge2.0/fixed/php/getData.php`. It only accepts POST `st=<NOS id>`; a GET returns an empty body.
  The station list is `https://slosh.nws.noaa.gov/etsurge2.0/datafiles/master.csv` and the cycle is in `PDY.txt`.
  It is an experimental site with 8x5 support.
