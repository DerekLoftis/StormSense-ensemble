# StormSense Ensemble

https://derekloftis.github.io/StormSense-ensemble/

Every available 72-hour storm tide forecast for each VIMS StormSense station on one chart, updated every
6 hours: NWS P-ETSS (with its 10–90% range), NWS ETSS, NOAA STOFS-2D-Global and STOFS-3D-Atlantic, NOAA
CBOFS, NWM coastal (short and medium range), the NWS official forecast, NHC P-Surge (during tropical
advisories), VIMS Tidewatch, tide + persistence, the StormSense model, and the published StormSense Latest
forecast.

- `ens72/<ST>_ens72_<MLLW|NAVD>.html|csv|json`: the chart with three skill tables under it, and one
  column per forecast in the csv and json
- `skill/`: skill of every forecast against the observations that followed (SS-Skill, R², bias, MAE,
  RMSE, per 6 hours, peak-height and timing errors, hits, misses and false alarms at each flood stage),
  and which forecast drives StormSense Latest at each station (`pick.json`)

The published files live on the [`site` branch](../../tree/site); `main` only holds this README and the
workflow that deploys the site to GitHub Pages. Run archives from 2026-10-16 to 27 stay on this
repository's `ensemble-2026-10` release.

Disclaimer: StormSense™ Predictive Plots are computer-generated operational storm tide forecasts; the
user assumes any and all risk.
