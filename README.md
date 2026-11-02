# A False-Alarm-Calibrated Benchmark of GNSS Water-Vapour Anomalies Before Large Earthquakes (2015-2025)

**Author:** Parth Jichkar (Department of Electronics and Communication Engineering & COSMOS Club, MIT World Peace University, Pune)

## Overview
Retrospective studies of Lithosphere-Atmosphere-Ionosphere Coupling (LAIC) frequently report pre-seismic anomalies in Global Navigation Satellite System (GNSS) Precipitable Water Vapor (PWV) within 1-30 days before major earthquakes using short-window heuristic thresholds ($\mu \pm 2\sigma$). However, such methods are rarely audited across multi-year non-seismic control periods to quantify their empirical False Alarm Rate (FAR) or tested against numerical weather reanalysis to separate synoptic meteorology from crustal processes.

This repository provides a self-contained, leak-free, false-alarm-calibrated statistical benchmark (`step13_final_pipeline.py`) evaluating 11-year (2015-2025) continuous 5-minute GNSS tropospheric solutions from the Nevada Geodetic Laboratory (NGL) and ERA5 total-column water vapour (`TCWV`) reanalysis across five $M_w \ge 7.1$ earthquakes ($N = 1,589$ catalog-certified 10-day quiet windows).

## Key Methodological & Scientific Contributions
1. **Quantifying the 34.0% False-Alarm Trap of Standard $2\sigma$ Literature Rules:** Replicating the standard 30-day moving $\mu + 2\sigma$ precursor rule across certified non-seismic quiet periods reveals a mean quiet-period False Alarm Rate of **34.0% ($\pm 5.0\%$)** (`Literature_vs_Rigorous_Benchmark.csv`).
2. **Leak-Free Guarded Baseline & All-Window Placebo Calibration:** Using a $\pm 45$-day temporal guard and leap-safe $\pm 15$-day day-of-year median/MAD baseline, held-out calibration locks the quiet-window false-alarm rate at **5.48%** (`PASS` across all 53 year boundaries and $1,589$ quiet windows).
3. **Physical Detectability Curve in Millimeters:** Synthetic anomaly injection across $1,589$ quiet windows establishes the exact 80% detection power threshold for 3-day pre-seismic anomalies: **$3.61\sigma \approx 22.2\text{ mm}$** in raw GNSS PWV and **$2.52\sigma \approx 3.8\text{ mm}$** in floored GNSS-minus-ERA5 residuals (`Step13_Injection_raw.csv`, `Step13_Injection_residual.csv`).
4. **Empirical Findings (2015-2025):**
   - **Primary Family (`Raw PWV`, `held_out_all`, `1-sided(+)`):** **0 of 3 core earthquakes** (`2024 Noto Mw 7.5`, `2024 Hyuga-nada Mw 7.1`, `2025 Tibet-Nepal Mw 7.1`; 95% Clopper-Pearson CI $0\text{-}71\%$, smallest Holm-adjusted $p = 0.231$) and **0 of 5 total earthquakes** (adding ultra-near-field `P595` at $10\text{ km}$ from `2019 Ridgecrest Mw 7.1` and `MERS` at $257\text{ km}$ from `2023 Turkey Mw 7.8`; 95% CI $0\text{-}52\%$, smallest Holm-adjusted $p = 0.385$) exceeded the locked threshold (`Step13_Event_Level_Primary.csv`).
   - **Synoptic Weather Reproduction:** Deseasonalized GNSS PWV anomalies correlate with ERA5 total-column water vapour reanalysis at **$r_Z = 0.940\text{-}0.989$** (RMS $1.09\text{-}3.24\text{ mm}$), demonstrating that pre-event spikes at `TSK2`, `MIZU`, and `IITK` are synoptic weather fronts (`Step13_Residual_Audit.csv`).

## Quick Start & Reproducibility
```bash
# 1. Run fast offline unit tests (29 checks)
python step13_final_pipeline.py --selftest

# 2. Run end-to-end offline smoke test on synthetic data
python step13_final_pipeline.py --smoketest

# 3. Execute the full 2015-2025 benchmark (fetches & verifies NGL/USGS/ERA5 data)
python step13_final_pipeline.py

# 4. Byte-verify local NGL archives against fresh server downloads
python step13_final_pipeline.py --verify-provenance
```

## Repository Contents
- `step13_final_pipeline.py` - Self-contained benchmark pipeline (zero external module dependencies).
- `config.yaml` & `DEVIATIONS_AND_SELECTION_RULE.md` - Locked configuration parameters and transparent deviations log.
- `Step13_*.csv` - Complete statistical output tables (Event-Level Primary, 84-row Multiverse Variant Table, Residual Audit, Yearly Bias, Floor Sensitivity, Event Sensitivity Grid, Injection Curves, Boundary Checks, Raw 10-Day Window Dumps).
- `03_Troposphere/PROVENANCE.csv` & `Step13_run_manifest.json` - SHA-256 hashes and URL provenance for all station-year archives.
