# LAIC pipeline audit report

## Verdict

The corrected, deterministic pipeline ran on the supplied archive. Under the explicitly defined exploratory rule (15-day trailing baseline, 1.5×IQR threshold, no missing-value interpolation), it generated **0 candidate alert days across the analysis period** and **0 candidate alert days in the 5-day pre-event window (2025-01-02 to 2025-01-06)**. This is not a validated precursor detector and should not be presented as an earthquake prediction result.

## Key source checks

- Event origin time and location use USGS reviewed metadata: 2025-01-07 01:05:16 UTC; 28.604°N, 87.378°E; Mw 7.1.
- The 36-day gas table is read only from `04_Gas_S5P/multi_gas_s5p.csv`; `s5p_co.csv` is not consulted, avoiding `glob()` order dependence. No gas values are interpolated.
- CH4 has 23/36 valid days (63.89%); SO2 contains 10 negative daily values. The CSV lacks quality flags, pixel counts, retrieval uncertainties and source-product metadata, so it cannot support pixel-validity, SNR or spatial percentile claims.
- The GNSS parser uses the NGL `TROP/SOLUTION` epoch field and the documented PWV/WVAPOR value at record field 10 (1-based); it does not guess a plausible number. It retained 33/36 daily PWV values. Missing/insufficient coverage dates: 2025-01-10, 2025-01-11, 2025-01-12.
- The OMNI parser reads the final three data columns as Kp*10, Dst and F10.7 based on their value patterns. The archive has no `.fmt` sidecar, so verify the original OMNIWeb parameter-selection metadata before journal submission.
- VTEC values are bilinearly interpolated from CODE global ionosphere maps at coordinates, with the IONEX geocentric latitude convention taken into account. These are model-grid estimates from the same GIM, not independent epicentral, Dhanbad CORS or control-station observations. No Dhanbad RINEX processing is in the archive.

## Jan 4 CO check

For Jan 4, 2025, full-precision CO is **0.018714612** (CSV units) and the trailing-IQR upper bound under W=15, k=1.5 is **0.018917023**. The value is therefore **not flagged** by this exact rule.

## Space-weather gate

During Jan 2–6, the daily maximum Kp-like index and daily minimum Dst are checked directly from the supplied OMNI file. The quiet gate is not satisfied on dates with Kp ≥ 3 or Dst ≤ −30 nT. This can suppress an alert by design; it must be documented rather than described as proof that every non-seismic cause has been eliminated.

## Statistical interpretation

The parameter-sensitivity grid is exploratory and shows how the output changes for W ∈ {10, 15, 20} days and k ∈ {1.0, 1.25, 1.5, 1.75, 2.0}. The circular-shift file enumerates 1,295 alternative independent circular alignments of the CO and PWV flag streams while holding VTEC and the quiet gate fixed. Here, an absence of alert in the selected pre-event window makes a p-like comparison uninformative; the shift frequency is not a real-world false-alarm rate. A single retrospectively selected earthquake cannot support claims such as 100% precision, F1=1.00, 0% false alarms or validated lead time.

## Required before a journal-ready revision

1. Recover the Google Earth Engine scripts, collection IDs, product versions, QA filters, region geometry and gas units. Re-export pixel-level statistics if valid-pixel counts, uncertainty/SNR or spatial percentiles are to be reported.
2. Verify the OMNI parameter selection from the missing `.fmt`/query metadata.
3. Obtain and process actual RINEX/CORS data before referring to Dhanbad or any control coordinate as a station validation.
4. Pre-register the rule and evaluate multiple earthquakes plus matched quiet periods using the same procedure.
5. Revise the manuscript to describe a retrospective case study and candidate signals, not validated precursors.

## Outputs

`daily_analysis.csv`, `satellite_gas_quality.csv`, `gnss_quality_summary.csv`, `site_distance_audit.csv`, `gim_point_correlations.csv`, `parameter_sensitivity.csv`, `circular_shift_null_distribution.csv`, `circular_shift_null_summary.csv`, `diagnostic_timeseries.png`, and `summary.json`.

## External documentation

- USGS reviewed event page: https://earthquake.usgs.gov/earthquakes/eventpage/us6000pi9w/executive
- NASA OMNIWeb documentation: https://omniweb.gsfc.nasa.gov/omniweb/html/ow_data.html
- International GNSS Service GIM product description (grid resolution and cadence): https://igs.org/products/
- Sentinel-5P methane product documentation (quality flags/uncertainties): https://sentinels.copernicus.eu/data-products/-/asset_publisher/fp37fc19FN8F/content/tropomi-level-2-methane
