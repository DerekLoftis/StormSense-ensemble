# StormSense Ensemble

Every 6 hours (00:40, 06:40, 12:40, 18:40 UTC) GitHub Actions fetches each surge guidance
source for the Tidewatch stations, lines the model guidance up with each gauge and publishes one
72-hour chart per station and datum, with CSV and JSON, to GitHub Pages:
https://derekloftis.github.io/StormSense-ensemble/

## Members

| Key | Source | Native datum | Aligned to gauge |
|---|---|---|---|
| petss | NWS P-ETSS mean, with the 10-90% cone | MLLW (via petss.py) | yes |
| stofs2d | STOFS-2D-Global (AWS) | MLLW | yes |
| stofs3d | STOFS-3D-Atlantic, 12Z run (AWS) | MLLW | yes |
| cbofs | CBOFS station output within 3 km (AWS) | MSL | yes (dropped without overlap) |
| nwm_short, nwm_medium | NWM coastal total water (AWS) | NAVD88 | yes |
| etss | NWS ETSS East Coast surge bulletin (tgftp FQUS23 KWNO), added to the P-ETSS tide | MLLW | yes |
| psurge | NHC P-Surge 50% storm tide (NOMADS SHEF), only while NHC has tropical advisories out | MLLW | yes |
| nws | NWS official stage forecast (NWPS) | each gauge's datum | no |
| tidewatch | Tidewatch f36 px30 | MLLW | no |
| baseline | astro tide plus the last 3 h residual, decaying over 24 h | MLLW | n/a |
| stormsense | StormSense's own forecast ([StormSense_forecastmodel](https://github.com/DerekLoftis/StormSense_forecastmodel)): its own harmonic tide plus the gauge's current residual handed over to a skill-weighted blend of the members above, with a 10-90% band (`stormsense_p10_ft`, `stormsense_p90_ft`, after the other columns) | MLLW | anchored to the gauge |

Alignment: the median of (observed minus member) over the last 24 h of overlap (at least 30
six-minute observations), using the newest run with the previous run stitched in front of it.
NAVD88 charts use Tidewatch's MLLW-NAVD88 offset for each station.

## Outputs

- `ens72/<ST>_ens72_<MLLW|NAVD>.html|csv|json` (columns `date_time_z, observed_ft, astro_tide_ft,
  <member>_ft..., petss_p10_ft, petss_p90_ft`)
- `latest.json`: each member's run time, and per station which members are present and their
  alignment shift in feet
- Archive: one `ensemble_<cycle>.tar.gz` per run (latest.json, every member's native series and the
  CSV/JSON) on that month's release, `ensemble-YYYY-MM`

## Files

CBBT, MNPT, DAHL and WADC observations come straight from NOAA CO-OPS (`noaa_override.py`), since Tidewatch's are stale or from another gauge.

`ensemble.py` builds everything (`python ensemble.py site work`), `members.py` fetches the
models, `petss.py` reads P-ETSS, `stations.csv` lists the stations and which members cover them,
`styler/` applies the StormSense chart style and NWS flood thresholds (`thresholds.json`).

## Skill and the guidance pick (phase 2)

`skill.py` runs after each archive upload (`python skill.py site/skill work/skill <new tarball>`). It scores every
archived run against the observations that came in after it, per station, member and lead band (0-12, 12-36,
36-72 h after the station's last observation, plus 0-36 and 0-72), once the band has passed, with Derek's skill
S = 0.25[R² + e^(-|Bias|/σ) + e^(-MAE/σ) + e^(-RMSE/σ)], σ the population SD of the observations in the band.
Truth is `observed_ft` stitched across the archived runs. The cumulative scores (`skill_scores.csv.gz`,
`skill_runs.csv`) are kept on the month's release; the last 5 days of runs are rescored each time.

Published under `skill/`: `index.html` (recent skill tables), `scores.csv` (every score, last 31 days),
`summary.csv`/`summary.json` (10-day recent skill, weighted exp(-age/3 d)), `pick.json` and `picks_history.csv`.
The pick per station and band is the member with the best recent skill; it must lead the current pick by 0.01 for
two updates running before it switches, and stays on P-ETSS until 7 days of verified runs exist and the member has
12 scores. Tidewatch, the baseline and the tide alone are scored but never picked, nor are stations whose Tidewatch
observations are another gauge's (USSC, USAR, USCB, USNR, USSB, GS01, USBB). StormSense Latest reads `pick.json`
when its `guidance.json` mode is `adaptive`.

This is phase 1 of the design in the project's `stormsense-latest/ensemble_design.md`; skill-based
selection of the primary StormSense guidance (phase 2) is above.
