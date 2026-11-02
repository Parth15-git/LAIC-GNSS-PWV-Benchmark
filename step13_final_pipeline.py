#!/usr/bin/env python3
"""
step13_final_pipeline.py  -  FINAL, SELF-CONTAINED PIPELINE (replaces step10, step11, step12)
===========================================================================================
Needs only: python 3.9+, numpy, pandas, matplotlib (scipy optional). It does NOT import step3/step6/step9.

    python step13_final_pipeline.py --selftest            # fast unit tests, offline
    python step13_final_pipeline.py --smoketest           # whole pipeline on fake data, offline (about 1-2 min)
    python step13_final_pipeline.py --verify-provenance   # really verify every local NGL zip against a fresh download
    python step13_final_pipeline.py --list-candidates     # USGS list of Mw>=7.1 events (for the event-selection rule)
    python step13_final_pipeline.py                       # full run

WHAT IS NEW COMPARED WITH step12 (each item is a mistake that step12 made or hid)
 1. Provenance: step12 labelled every local file 'verified_local' without verifying anything, so 'ORIGIN_UNKNOWN = 0'
    was meaningless. Now local files are 'local_unverified' until --verify-provenance proves they are byte-identical
    to a fresh NGL download; old fake labels are relabelled automatically.
 2. Download failures were hidden. Now every failure is printed with its reason (HTTP 404 = year not available,
    anything else = 'uncertain, network').
 3. Checks that step11 had and step12 deleted are back: placebo, injection curve, boundary check (it is CALLED now),
    run manifest, config.yaml, deviations file, self-test, candidate listing, catalog epicentres.
 4. The floored residual is now tested: placebo FAR, floor sensitivity (1.0/1.5/2.0/3.0 mm), minimum detectable
    anomaly in mm for raw and residual series.
 5. Event-level statistics: stations of one earthquake are collapsed (any-of, Bonferroni inside the event), Holm is
    applied across events, and the sensitivity grid uses a Poisson-binomial expectation that accounts for two
    stations at one event.
 6. ERA5 cache is no longer used silently when it does not cover the needed years (loud warning, or fetch).
 7. Yearly GNSS-minus-ERA5 bias table: detects a product step (for example IGS14 vs IGS20) between years.
 8. Catalogs (USGS ComCat M>=5.5 within 800 km) and mainshock epicentres are fetched and cached by this script.
 9. The 'negative controls' (Taiwan, Guam) are labelled honestly: different climate, NOT controls.
10. Docstring claims are replaced by a printed FINAL CHECKLIST computed from the actual run.

WHAT THIS SCRIPT CANNOT DO FOR YOU (printed again at the end of every run)
 - The 5 events are a convenience sample. The selection rule is NOT applied automatically.
 - VTEC and trace gas are not part of this pipeline: the study is PWV-only.
 - IITK has only the years NGL provides (2022-2025 at the time of writing): that event stays weak.
"""
import os
import sys
import io
import gc
import json
import math
import gzip
import glob
import shutil
import hashlib
import platform
import argparse
import zipfile
import tempfile
import urllib.parse
import urllib.request
import urllib.error
import datetime as _dt
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

DIR_TROP = "03_Troposphere"
DIR_CAT = "00_Catalogs"
DIR_CACHE = "05_Local_ERA5_Cache"
ORIGINS_FILE = os.path.join(DIR_TROP, "VERIFIED_ORIGINS.json")

# ======================================================================================
# SINGLE SOURCE OF TRUTH (config.yaml is generated from this dict)
# ======================================================================================
CONFIG = {
    "observable_primary": "GNSS-derived PWV (NGL WVAPOR, mm; NGL uses VMF1 weather-model inputs)",
    "observable_secondary": "GNSS horizontal gradient magnitude G = sqrt(TGN^2 + TGE^2) (invariant to the N/E column swap)",
    "time_standard": "UTC days; a day needs >= 200 of 288 five-minute epochs",
    "day_zero": "UTC calendar date of mainshock origin time; pre-event window = UTC days [-W, -1]",
    "doy_half_window_days": 15,
    "guard_days": 45,
    "min_ref_values": 30,
    "window_len": 10,
    "min_valid_frac": 0.8,
    "target_far": 0.05,
    "primary_series": "Raw PWV",
    "primary_calibration": "held_out_all",
    "primary_direction": "one_sided_positive",
    "held_out_exclusion_days": 180,
    "season_half_width_days": 45,
    "season_exclusion_days": 90,
    "min_cal_windows": 30,
    "quiet_buffer_days": 10,
    "quiet_m_near": 5.5,
    "quiet_r_near_km": 300.0,
    "quiet_m_far": 6.5,
    "quiet_r_far_km": 800.0,
    "event_exclusion_days": 30,
    "catalog_min_mag": 5.5,
    "catalog_radius_km": 800,
    "n_bootstrap": 1000,
    "placebo_gate_pct": [3.0, 7.0],
    "boundary_gate_pct": 20.0,
    "residual_floor_mm": 1.50,
    "residual_floor_sensitivity_mm": [1.0, 1.5, 2.0, 3.0],
    "window_lengths_tested": [7, 10, 14],
    "target_fars_tested": [0.02, 0.05, 0.10],
    "min_corr_gnss_vs_era5_tcwv": 0.85,
    "yearly_bias_flag_mm": 1.0,
    "detect_level_pct": 80.0,
    "trop_base_urls": [
        "https://geodesy.unr.edu/gps_timeseries/IGS20/trop",
        "https://geodesy.unr.edu/gps_timeseries/trop",
    ],
    "random_seed": 20261031,
    "event_selection_status": "NOT APPLIED (convenience sample of 5 events)",
    "event_selection_proposed_rule": "USGS ComCat, Mw >= 7.1, depth <= 70 km, 2019-01-01..2025-03-31, >=1 NGL station "
                                     "within 750 km with >=8/10 valid pre-event days and >=30 quiet windows; include ALL that pass",
}
DIRECTIONS = ("one_sided_positive", "two_sided_abs")
S_RAW, S_LOG, S_GRAD, S_RES = "Raw PWV", "Log(PWV)", "GNSS Gradient (G)", "Col-Weather-Res (Floored)"

STATIONS_DEF = [
    dict(ev="2024 Noto (Mw 7.5)", stn="TSK2", lat=36.06, lon=140.09, date="2024-01-01", years=(2015, 2025), role="Core", region="Japan", mw=7.5, epi_fb=(37.50, 137.27)),
    dict(ev="2024 Noto (Mw 7.5)", stn="MIZU", lat=39.13, lon=141.13, date="2024-01-01", years=(2015, 2025), role="Core", region="Japan", mw=7.5, epi_fb=(37.50, 137.27)),
    dict(ev="2024 Hyuga-nada (Mw 7.1)", stn="AIRA", lat=31.82, lon=130.60, date="2024-08-08", years=(2015, 2025), role="Core", region="Japan", mw=7.1, epi_fb=(31.74, 131.53)),
    dict(ev="2025 Tibet-Nepal (Mw 7.1)", stn="IITK", lat=26.52, lon=80.23, date="2025-01-07", years=(2015, 2025), role="Core", region="Himalaya", mw=7.1, epi_fb=(28.64, 87.36)),
    dict(ev="2019 Ridgecrest (Mw 7.1)", stn="P595", lat=35.69, lon=-117.55, date="2019-07-06", years=(2015, 2025), role="Extended", region="California", mw=7.1, epi_fb=(35.77, -117.60)),
    dict(ev="2023 Turkey (Mw 7.8)", stn="MERS", lat=36.56, lon=34.25, date="2023-02-06", years=(2015, 2025), role="Extended", region="Anatolia", mw=7.8, epi_fb=(37.17, 37.03)),
]
EXTRA_STATIONS = [  # context only, NOT controls
    dict(stn="USUD", lat=36.13, lon=138.36, years=(2015, 2025), desc="Usuda, Japan", region="Japan"),
    dict(stn="TWTF", lat=24.95, lon=121.16, years=(2015, 2025), desc="Taoyuan, Taiwan", region="Japan"),
    dict(stn="GUAM", lat=13.59, lon=144.87, years=(2015, 2025), desc="Dededo, Guam", region="Japan"),
]


# ======================================================================================
# SMALL UTILITIES
# ======================================================================================
def _haversine_km(lat1, lon1, lat2, lon2):
    R = 6371.0088
    p1 = np.radians(float(lat1))
    p2 = np.radians(np.asarray(lat2, dtype=float))
    dphi = p2 - p1
    dl = np.radians(np.asarray(lon2, dtype=float) - float(lon1))
    a = np.sin(dphi / 2.0) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2.0) ** 2
    return 2.0 * R * np.arcsin(np.sqrt(a))


def norm_doy(index):
    """Leap-safe day of year on a 365-day circle (Feb 29 -> Feb 28; later days of a leap year shifted by one)."""
    idx = pd.DatetimeIndex(index)
    doy = idx.dayofyear.values.astype(int)
    leap = np.asarray(idx.is_leap_year)
    month, day = idx.month.values, idx.day.values
    doy = np.where(leap & (month > 2), doy - 1, doy)
    doy = np.where(leap & (month == 2) & (day == 29), 59, doy)
    return doy


def circ_dist(a, b, period=365):
    d = np.abs(np.asarray(a) - np.asarray(b))
    return np.minimum(d, period - d)


def wilson_ci(k, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    den = 1.0 + z * z / n
    c = (p + z * z / (2.0 * n)) / den
    h = z * math.sqrt(p * (1.0 - p) / n + z * z / (4.0 * n * n)) / den
    return (max(0.0, c - h), min(1.0, c + h))


def clopper_pearson(k, n, alpha=0.05):
    if n == 0:
        return (float("nan"), float("nan"))
    try:
        from scipy.stats import beta
        lo = 0.0 if k == 0 else float(beta.ppf(alpha / 2.0, k, n - k + 1))
        hi = 1.0 if k == n else float(beta.ppf(1.0 - alpha / 2.0, k + 1, n - k))
        return (lo, hi)
    except ImportError:  # exact fallback for the k == 0 and k == n cases, Wilson otherwise
        if k == 0:
            return (0.0, 1.0 - (alpha / 2.0) ** (1.0 / n))
        if k == n:
            return ((alpha / 2.0) ** (1.0 / n), 1.0)
        return wilson_ci(k, n)


def binom_tail(k, n, p):
    return float(sum(math.comb(n, i) * (p ** i) * ((1.0 - p) ** (n - i)) for i in range(k, n + 1)))


def poisson_binomial_tail(k, ps):
    """P(X >= k) where X is a sum of independent Bernoulli(p_i)."""
    dist = np.zeros(len(ps) + 1)
    dist[0] = 1.0
    for p in ps:
        new = np.zeros_like(dist)
        new[0] = dist[0] * (1.0 - p)
        new[1:] = dist[1:] * (1.0 - p) + dist[:-1] * p
        dist = new
    return float(dist[k:].sum())


def holm_adjust(pvals):
    p = np.asarray(pvals, dtype=float)
    m = len(p)
    if m == 0:
        return p
    order = np.argsort(p)
    adj = np.empty(m)
    running = 0.0
    for rank, idx in enumerate(order):
        running = max(running, (m - rank) * p[idx])
        adj[idx] = min(1.0, running)
    return adj


def k_star_from_scores(scores, target_far):
    """Smallest threshold k with mean(scores > k) <= target_far (exact order statistic)."""
    q = np.sort(np.asarray(scores, dtype=float))[::-1]
    m = int(math.floor(target_far * len(q) + 1e-9))
    return float(q[m]) if m < len(q) else float(q[-1])


def sha256_file(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def now_utc():
    return _dt.datetime.now(_dt.timezone.utc).isoformat()


# ======================================================================================
# NETWORK (single entry point, so tests can replace it)
# ======================================================================================
def http_get(url, timeout=120):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (LAIC research script)"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def _classify_error(e):
    if isinstance(e, urllib.error.HTTPError):
        return f"HTTP {e.code}", (e.code in (404, 410))
    if isinstance(e, urllib.error.URLError):
        return f"network: {e.reason}", False
    return f"{type(e).__name__}: {e}", False


# ======================================================================================
# NGL ZIP PARSER (header driven; same layout as your real files)
# ======================================================================================
def _valid_zip(path):
    try:
        if os.path.getsize(path) <= 10_000 or not zipfile.is_zipfile(path):
            return False
        with zipfile.ZipFile(path) as z:
            return z.testzip() is None
    except Exception:
        return False


def parse_station_trop_zips(stn, years, trop_dir=None):
    """
    Returns daily means (needs >= 200 of 288 epochs): ZWD_mean, PWV_mean, PWV_sig_med, G_mean, count.
    Header: *SITE ___EPOCH____ TROTOT _SIG TRWET TGETOT _SIG TGNTOT _SIG WVAPOR _SIG MTEMP
    Columns are located by NAME. The gradient magnitude does not depend on the N/E column swap.
    """
    trop_dir = DIR_TROP if trop_dir is None else trop_dir
    cols_out = ["ZWD_mean", "PWV_mean", "PWV_sig_med", "G_mean", "count"]
    recs = []
    for yr in years:
        zpath = os.path.join(trop_dir, f"{stn}.{yr}.trop.zip")
        if not os.path.exists(zpath) or not _valid_zip(zpath):
            continue
        with zipfile.ZipFile(zpath, "r") as z:
            for member in sorted(m for m in z.namelist() if not m.endswith("/")):
                raw = z.read(member)
                text = gzip.decompress(raw).decode("ascii", errors="ignore") if raw[:2] == b"\x1f\x8b" else raw.decode("ascii", errors="ignore")
                in_sol, header = False, []
                ix = {}
                for line in text.splitlines():
                    s = line.strip()
                    if not s:
                        continue
                    if s.startswith("+TROP/SOLUTION"):
                        in_sol = True
                        continue
                    if s.startswith("-TROP/SOLUTION"):
                        in_sol = False
                        continue
                    if in_sol and s.startswith("*SITE"):
                        header = s.lstrip("*").split()
                        ix = {}
                        for pos, tok in enumerate(header):
                            tu = tok.upper()
                            if "EPOCH" in tu:
                                ix["epoch"] = pos
                            elif tu in ("TGNTOT", "TGETOT", "TRWET", "WVAPOR"):
                                ix[tu] = pos
                        if "WVAPOR" in ix:
                            ix["WVSIG"] = ix["WVAPOR"] + 1
                        continue
                    need = {"epoch", "TGNTOT", "TGETOT", "TRWET", "WVAPOR", "WVSIG"}
                    if not (in_sol and header and need.issubset(ix) and not s.startswith(("*", "#", "%"))):
                        continue
                    parts = s.split()
                    if len(parts) < len(header):
                        continue
                    try:
                        ep = parts[ix["epoch"]].split(":")
                        yy, doy, sod = int(ep[0]), int(ep[1]), int(float(ep[2]))
                        year_full = 2000 + yy if yy < 80 else 1900 + yy
                        dtm = pd.Timestamp(year=year_full, month=1, day=1) + pd.Timedelta(days=doy - 1, seconds=sod)
                        zwd, pwv = float(parts[ix["TRWET"]]), float(parts[ix["WVAPOR"]])
                        wvsig = float(parts[ix["WVSIG"]])
                        g1, g2 = float(parts[ix["TGNTOT"]]), float(parts[ix["TGETOT"]])
                    except (ValueError, IndexError):
                        continue
                    if 0.0 < pwv <= 120.0 and -50.0 < zwd < 800.0 and abs(g1) < 50.0 and abs(g2) < 50.0:
                        recs.append((dtm, zwd, pwv, wvsig, math.sqrt(g1 * g1 + g2 * g2)))
    if not recs:
        return pd.DataFrame(columns=cols_out)
    df5 = pd.DataFrame(recs, columns=["datetime", "ZWD_mm", "PWV_mm", "WVSIG_mm", "G_mm"]).drop_duplicates("datetime")
    df5 = df5.set_index("datetime").sort_index()
    daily = df5.resample("1D").agg(ZWD_mean=("ZWD_mm", "mean"), PWV_mean=("PWV_mm", "mean"),
                                   PWV_sig_med=("WVSIG_mm", "median"), G_mean=("G_mm", "mean"), count=("PWV_mm", "count"))
    return daily[daily["count"] >= 200].copy()


# ======================================================================================
# DOWNLOADS WITH HONEST PROVENANCE
# ======================================================================================
def load_origins():
    if not os.path.exists(ORIGINS_FILE):
        return {}
    try:
        with open(ORIGINS_FILE, "r", encoding="utf-8") as f:
            d = json.load(f)
    except Exception:
        return {}
    for k, v in d.items():  # step12 wrote a made-up 'verified_local' label: never trust it
        if "verified_local" in str(v.get("origin", "")):
            v["origin"] = "local_unverified (relabelled: the old 'verified_local' label was not a real check)"
    return d


def save_origins(d):
    os.makedirs(DIR_TROP, exist_ok=True)
    with open(ORIGINS_FILE, "w", encoding="utf-8") as f:
        json.dump(d, f, indent=2)


def ensure_trop_zip(stn, year, force=False, prov=None, origins=None, trop_dir=None, verbose=True):
    """Returns the path of a valid zip or None. Every failure is reported with its reason."""
    trop_dir = DIR_TROP if trop_dir is None else trop_dir
    os.makedirs(trop_dir, exist_ok=True)
    fname = f"{stn}.{year}.trop.zip"
    path = os.path.join(trop_dir, fname)
    origins = {} if origins is None else origins
    if os.path.exists(path) and not force and _valid_zip(path):
        sha = sha256_file(path)
        rec = origins.get(fname)
        if rec and rec.get("sha256") == sha:
            origin = rec["origin"]
        else:
            origin = "local_unverified"
            origins[fname] = {"origin": origin, "sha256": sha}
        if prov is not None:
            prov.append({"file": fname, "origin": origin, "bytes": os.path.getsize(path), "sha256": sha, "utc": now_utc()})
        return path
    failures, only_404 = [], True
    for base in CONFIG["trop_base_urls"]:
        url = f"{base}/{stn}/{fname}"
        tmp = path + ".part"
        try:
            data = http_get(url, timeout=180)
            with open(tmp, "wb") as f:
                f.write(data)
            if not _valid_zip(tmp):
                failures.append((url, "not a valid zip"))
                only_404 = False
                os.remove(tmp)
                continue
            os.replace(tmp, path)
            sha = sha256_file(path)
            origins[fname] = {"origin": url, "sha256": sha}
            if verbose:
                print(f"    downloaded {fname} from {url}")
            if prov is not None:
                prov.append({"file": fname, "origin": url, "bytes": os.path.getsize(path), "sha256": sha, "utc": now_utc()})
            return path
        except Exception as e:
            why, is_404 = _classify_error(e)
            only_404 &= is_404
            failures.append((url, why))
            if os.path.exists(tmp):
                os.remove(tmp)
    if verbose:
        kind = "NOT AVAILABLE on NGL" if only_404 else "UNCERTAIN (network/other error, not proof the year is missing)"
        print(f"    {fname}: {kind}: " + "; ".join(f"{u.split('/gps_timeseries/')[1].split('/')[0]}: {w}" for u, w in failures))
    if prov is not None:
        prov.append({"file": fname, "origin": "NOT_AVAILABLE" if only_404 else "UNCERTAIN_FETCH_ERROR", "bytes": 0, "sha256": "", "utc": now_utc()})
    return None


def verify_provenance(trop_dir=None):
    """Compare every local zip (years before the current year) with a fresh download. Only identical files are certified."""
    trop_dir = DIR_TROP if trop_dir is None else trop_dir
    origins = load_origins()
    cur = _dt.datetime.now(_dt.timezone.utc).year
    rows = []
    for path in sorted(glob.glob(os.path.join(trop_dir, "*.trop.zip"))):
        fname = os.path.basename(path)
        parts = fname.split(".")
        stn, year = parts[0], int(parts[1])
        if year >= cur:
            rows.append((fname, "SKIPPED (current year is updated daily by NGL)", ""))
            continue
        sha_local = sha256_file(path)
        result, used = "UNAVAILABLE (could not fetch a reference copy)", ""
        for base in CONFIG["trop_base_urls"]:
            url = f"{base}/{stn}/{fname}"
            try:
                data = http_get(url, timeout=180)
                used = url
                if hashlib.sha256(data).hexdigest() == sha_local:
                    result = "IDENTICAL"
                    origins[fname] = {"origin": f"verified_identical_to {url}", "sha256": sha_local}
                else:
                    result = "DIFFERENT (reprocessed by NGL, or a different product)"
                break
            except Exception:
                continue
        rows.append((fname, result, used))
        print(f"  {fname:<24} {result}")
    save_origins(origins)
    n_id = sum(1 for r in rows if r[1] == "IDENTICAL")
    n_diff = sum(1 for r in rows if r[1].startswith("DIFFERENT"))
    print(f"\nVerified identical: {n_id} | different: {n_diff} | unavailable/skipped: {len(rows) - n_id - n_diff}")
    if n_diff:
        print("DIFFERENT files mean the local copy is NOT the same product as the current NGL file: re-run with the")
        print("redownload option for those files before trusting any baseline built from them.")
    return rows


# ======================================================================================
# CATALOGS (USGS ComCat, fetched and cached here)
# ======================================================================================
_USGS = "https://earthquake.usgs.gov/fdsnws/event/1/query?"


def usgs_query(params, cache_name):
    os.makedirs(DIR_CAT, exist_ok=True)
    path = os.path.join(DIR_CAT, cache_name)
    cols = ["time", "latitude", "longitude", "depth", "mag", "place"]
    if os.path.exists(path) and os.path.getsize(path) > 0:
        df = pd.read_csv(path)
    else:
        url = _USGS + urllib.parse.urlencode({**params, "format": "csv", "orderby": "time-asc"})
        text = http_get(url, timeout=180).decode("utf-8", errors="ignore")
        if not text.strip():
            text = ",".join(cols) + "\n"
        with open(path, "w", encoding="utf-8") as f:
            f.write(text)
        df = pd.read_csv(io.StringIO(text))
    for c in cols:
        if c not in df.columns:
            df[c] = np.nan
    df["time_dt"] = pd.to_datetime(df["time"], utc=True, errors="coerce").dt.tz_localize(None)
    return df[["time_dt", "latitude", "longitude", "depth", "mag", "place"]].dropna(subset=["time_dt"]).reset_index(drop=True)


def station_catalog(c):
    y0, y1 = c["years"]
    p = {"latitude": c["lat"], "longitude": c["lon"], "maxradiuskm": CONFIG["catalog_radius_km"],
         "minmagnitude": CONFIG["catalog_min_mag"], "starttime": f"{y0}-01-01", "endtime": f"{y1 + 1}-01-01"}
    return usgs_query(p, f"usgs_{c['stn']}_{y0}_{y1}_m{CONFIG['catalog_min_mag']}_r{CONFIG['catalog_radius_km']}.csv")


def mainshock(c):
    d = pd.Timestamp(c["date"])
    p = {"starttime": (d - pd.Timedelta(days=1)).strftime("%Y-%m-%d"), "endtime": (d + pd.Timedelta(days=2)).strftime("%Y-%m-%d"),
         "minmagnitude": 6.0}
    df = usgs_query(p, f"usgs_main_{c['date']}.csv")
    if df.empty:
        raise RuntimeError("no mainshock found")
    r = df.loc[df["mag"].astype(float).idxmax()]
    return float(r["latitude"]), float(r["longitude"]), float(r["mag"]), r["time_dt"]


def list_candidate_events():
    p = {"starttime": "2019-01-01", "endtime": "2025-03-31", "minmagnitude": 7.1, "maxdepth": 70}
    df = usgs_query(p, "candidate_events_mw71_2019_2025.csv")
    print(df[["time_dt", "latitude", "longitude", "depth", "mag", "place"]].to_string(index=False))
    print(f"\n{len(df)} candidate events saved in {DIR_CAT}/candidate_events_mw71_2019_2025.csv")
    print("Apply the rule BEFORE looking at any result: for every event, check whether an NGL station within 750 km has")
    print(">= 8/10 valid pre-event days and >= 30 quiet windows. Include ALL that pass; list the failures with reasons.")


# ======================================================================================
# ERA5 TOTAL COLUMN WATER VAPOUR (cache-aware, coverage-checked)
# ======================================================================================
def _read_tcwv_csv(p):
    df = pd.read_csv(p, parse_dates=["date"])
    return df.set_index("date")["tcwv_mm"].astype(float)


def load_era5_tcwv(stn, lat, lon, start, end):
    """Returns (series or None, note). Uses cached files only if they COVER [start, end]; otherwise fetches (3-year chunks)."""
    os.makedirs(DIR_CACHE, exist_ok=True)
    s0, s1 = pd.Timestamp(start), pd.Timestamp(end)
    pieces = []
    for p in sorted(glob.glob(os.path.join(DIR_CACHE, f"tcwv*{stn}*.csv"))):
        try:
            if os.path.getsize(p) > 2000:
                pieces.append(_read_tcwv_csv(p))
        except Exception:
            continue
    have = pd.concat(pieces).sort_index() if pieces else pd.Series(dtype=float)
    have = have[~have.index.duplicated(keep="last")]
    tol = pd.Timedelta(days=45)
    if len(have) and have.index.min() <= s0 + tol and have.index.max() >= s1 - tol:
        return have[(have.index >= s0) & (have.index <= s1)], f"cache ({len(have)} days)"
    chunks, err = [], None
    y = s0.year
    while y <= s1.year:
        y_end = min(y + 2, s1.year)
        url = ("https://archive-api.open-meteo.com/v1/archive?latitude=%s&longitude=%s&start_date=%s&end_date=%s"
               "&hourly=total_column_integrated_water_vapour&timezone=GMT" % (lat, lon, f"{y}-01-01", f"{y_end}-12-31"))
        try:
            js = json.loads(http_get(url, timeout=180).decode("utf-8"))
            h = js["hourly"]
            v = pd.Series(h["total_column_integrated_water_vapour"], index=pd.to_datetime(h["time"]), dtype=float)
            agg = v.resample("1D").agg(["mean", "count"])
            chunks.append(agg.loc[agg["count"] >= 20, "mean"])
        except Exception as e:
            err = e
        y = y_end + 1
    if chunks:
        got = pd.concat(chunks)
        if len(have):
            got = pd.concat([have, got])
        got = got[~got.index.duplicated(keep="last")].sort_index()
        got = got[(got > 0) & (got < 120)]
        got.to_frame("tcwv_mm").rename_axis("date").to_csv(os.path.join(DIR_CACHE, f"tcwv_verified_{stn}_{s0.year}_{s1.year}.csv"))
        return got[(got.index >= s0) & (got.index <= s1)], "fetched from Open-Meteo (ERA5)" + (" [PARTIAL: some chunks failed]" if err else "")
    if len(have):
        print(f"    [WARNING] cached ERA5 for {stn} does NOT cover {start}..{end} and the fetch failed ({err}). Using the partial cache.")
        return have, "PARTIAL CACHE"
    print(f"    [WARNING] ERA5 total-column water vapour not available for {stn}: {err}")
    return None, "unavailable"


# ======================================================================================
# BASELINE, FLOOR, WINDOWS, CALIBRATION
# ======================================================================================
def guarded_baseline(series, half_window=None, guard_days=None, min_ref=None):
    """
    Reference = same (leap-safe) day of year +/- half_window, from days MORE than guard_days away from t (any year).
    Median and 1.4826*MAD. Fewer than min_ref reference days -> NaN (no silent fallback).
    Columns: raw, med, mad (already scaled spread), Z, n_ref.
    """
    half_window = CONFIG["doy_half_window_days"] if half_window is None else half_window
    guard_days = CONFIG["guard_days"] if guard_days is None else guard_days
    min_ref = CONFIG["min_ref_values"] if min_ref is None else min_ref
    s = series.dropna().astype(float)
    cols = ["raw", "med", "mad", "Z", "n_ref"]
    if s.empty:
        return pd.DataFrame(columns=cols)
    idx = pd.DatetimeIndex(s.index).normalize()
    vals = s.values
    n = len(vals)
    doy = norm_doy(idx)
    day_ints = (idx - pd.Timestamp("1970-01-01")).days.values.astype(np.int64)
    med_a, spr_a, z_a = np.full(n, np.nan), np.full(n, np.nan), np.full(n, np.nan)
    n_ref = np.zeros(n, dtype=int)
    for i in range(n):
        mask = (circ_dist(doy, doy[i]) <= half_window) & (np.abs(day_ints - day_ints[i]) > guard_days)
        ref = vals[mask]
        n_ref[i] = len(ref)
        if len(ref) < min_ref:
            continue
        med = float(np.median(ref))
        spread = float(1.4826 * np.median(np.abs(ref - med)))
        if spread < 1e-6:
            spread = float(np.std(ref))
        if spread < 1e-6:
            continue
        med_a[i], spr_a[i], z_a[i] = med, spread, (vals[i] - med) / spread
    return pd.DataFrame({"raw": vals, "med": med_a, "mad": spr_a, "Z": z_a, "n_ref": n_ref}, index=idx)


def apply_floor(base_df, floor_mm):
    """Re-scale an UNFLOORED baseline with spread = max(spread, floor). Needed because a residual of ~1 mm has a
    spread smaller than the instrument noise, which would turn sub-millimetre differences into large z-scores."""
    out = base_df.copy()
    spr = out["mad"].values
    new = np.where(np.isnan(spr), np.nan, np.maximum(spr, floor_mm))
    out["mad"] = new
    out["Z"] = (out["raw"].values - out["med"].values) / new
    return out


def boundary_continuity(z_df, gate_pct=None):
    gate_pct = CONFIG["boundary_gate_pct"] if gate_pct is None else gate_pct
    out = []
    if z_df is None or z_df.empty:
        return out
    for y in sorted(set(z_df.index.year))[:-1]:
        seg = z_df.loc[f"{y}-12-20":f"{y + 1}-01-10"].dropna(subset=["mad"])
        if len(seg) < 10:
            continue
        ch = seg["mad"].pct_change().abs().dropna() * 100.0
        same = (seg[["n_ref", "med", "mad"]].diff().abs().sum(axis=1).iloc[1:] == 0).values
        run = best = 0
        for s_ in same:
            run = run + 1 if s_ else 0
            best = max(best, run)
        out.append({"year_boundary": f"{y}/{y + 1}", "max_spread_change_pct": float(ch.max()), "where": str(ch.idxmax().date()),
                    "longest_identical_run": int(best), "pass": bool(ch.max() < gate_pct and best < 3)})
    return out


def _catalog_times(cat_df):
    t = pd.to_datetime(cat_df["time_dt"])
    if getattr(t.dt, "tz", None) is not None:
        t = t.dt.tz_convert("UTC").dt.tz_localize(None)
    return t


def build_quiet_windows(stn_lat, stn_lon, z_df, cat_df, years, exclude_dates, win_len=None, min_valid_frac=None):
    """Non-overlapping windows. Quiet = not within +/-event_exclusion_days of any event of the region AND no catalog
    event M>=m_near within r_near km or M>=m_far within r_far km between (start - buffer) and (end + buffer)."""
    C = CONFIG
    win_len = C["window_len"] if win_len is None else win_len
    min_valid_frac = C["min_valid_frac"] if min_valid_frac is None else min_valid_frac
    if z_df is None or z_df.empty or "Z" not in z_df.columns:
        return []
    min_valid = int(math.ceil(win_len * min_valid_frac))
    big_times = np.array([], dtype="datetime64[ns]")
    if cat_df is not None and len(cat_df):
        d = _haversine_km(stn_lat, stn_lon, cat_df["latitude"].values, cat_df["longitude"].values)
        m = cat_df["mag"].values.astype(float)
        big = ((m >= C["quiet_m_near"]) & (d <= C["quiet_r_near_km"])) | ((m >= C["quiet_m_far"]) & (d <= C["quiet_r_far_km"]))
        big_times = np.sort(_catalog_times(cat_df)[big].values.astype("datetime64[ns]"))
    ex = [pd.Timestamp(e).normalize() for e in exclude_dates]
    buf, exd = pd.Timedelta(days=C["quiet_buffer_days"]), pd.Timedelta(days=C["event_exclusion_days"])
    curr = pd.Timestamp(f"{min(years)}-01-01")
    end_dt = pd.Timestamp(f"{max(years)}-12-31")
    out = []
    while curr + pd.Timedelta(days=win_len - 1) <= end_dt:
        w0, w1 = curr, curr + pd.Timedelta(days=win_len - 1)
        curr = curr + pd.Timedelta(days=win_len)
        if any((w0 <= e + exd) and (w1 >= e - exd) for e in ex):
            continue
        if len(big_times) and np.any((big_times >= np.datetime64(w0 - buf)) & (big_times <= np.datetime64(w1 + buf))):
            continue
        zs = z_df.loc[w0:w1, "Z"].dropna().values
        if len(zs) < min_valid:
            continue
        srt, srt_abs = np.sort(zs), np.sort(np.abs(zs))
        out.append({"start": w0, "end": w1, "mid_ts": w0 + pd.Timedelta(days=win_len // 2), "year": int(w0.year),
                    "doy": int(norm_doy([w0])[0]), "n_valid": int(len(zs)), "z_values": zs,
                    "z2_pos": float(srt[-2]), "z2_abs": float(srt_abs[-2]), "zmax_pos": float(srt[-1])})
    return out


def common_windows(qw_dict):
    if not qw_dict:
        return qw_dict
    keep = set.intersection(*[set(w["start"] for w in lst) for lst in qw_dict.values()])
    return {name: [w for w in lst if w["start"] in keep] for name, lst in qw_dict.items()}


def select_cal_pool(qw_all, eq_ts, cal_mode):
    C = CONFIG
    if cal_mode == "in_sample_all":
        return list(qw_all)
    if cal_mode == "held_out_all":
        return [w for w in qw_all if abs((w["mid_ts"] - eq_ts).days) > C["held_out_exclusion_days"]]
    if cal_mode == "held_out_season45d":
        ev_doy = int(norm_doy([eq_ts])[0])
        return [w for w in qw_all if abs((w["mid_ts"] - eq_ts).days) > C["season_exclusion_days"]
                and circ_dist(w["doy"], ev_doy) <= C["season_half_width_days"]]
    raise ValueError(f"unknown cal_mode {cal_mode}")


def block_bootstrap_k(pool, key, target_far, n_boot, rng):
    if n_boot <= 0:
        return (float("nan"), float("nan"))
    by_year = {}
    for w in pool:
        by_year.setdefault(w["year"], []).append(w[key])
    ys = list(by_year)
    boots = []
    for _ in range(n_boot):
        if len(ys) >= 3:
            pick = rng.integers(0, len(ys), len(ys))
            sample = np.concatenate([by_year[ys[i]] for i in pick])
        else:
            allv = np.concatenate([np.asarray(v) for v in by_year.values()])
            sample = rng.choice(allv, size=len(allv), replace=True)
        boots.append(k_star_from_scores(sample, target_far))
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return (float(lo), float(hi))


def evaluate_variant(z_df, qw_all, eq_dt, cal_mode="held_out_all", direction="one_sided_positive",
                     target_far=None, win_len=None, n_boot=None, rng=None):
    C = CONFIG
    target_far = C["target_far"] if target_far is None else target_far
    win_len = C["window_len"] if win_len is None else win_len
    n_boot = C["n_bootstrap"] if n_boot is None else n_boot
    rng = np.random.default_rng(C["random_seed"]) if rng is None else rng
    if not qw_all or z_df is None or z_df.empty:
        return None
    eq_ts = pd.Timestamp(eq_dt).normalize()
    key = "z2_pos" if direction == "one_sided_positive" else "z2_abs"
    pool = select_cal_pool(qw_all, eq_ts, cal_mode)
    n_cal = len(pool)
    if n_cal < 5:
        return None
    scores = np.array([w[key] for w in pool])
    k_star = k_star_from_scores(scores, target_far)
    far_ach = float(np.mean(scores > k_star))
    ci = block_bootstrap_k(pool, key, target_far, n_boot, rng)
    pre = z_df.loc[eq_ts - pd.Timedelta(days=win_len): eq_ts - pd.Timedelta(days=1), "Z"].dropna().values
    if len(pre) < int(math.ceil(win_len * C["min_valid_frac"])):
        return None
    vals = pre if direction == "one_sided_positive" else np.abs(pre)
    ev_stat = float(np.sort(vals)[-2])
    days_over = int(np.sum(vals > k_star))
    p = float((np.sum(scores >= ev_stat) + 1) / (n_cal + 1))
    reliable = n_cal >= C["min_cal_windows"]
    hit = days_over >= 2
    return {"N_cal": n_cal, "k_star": k_star, "FAR_achieved": far_ach, "k_CI95": ci, "ev_stat": ev_stat,
            "days_over": days_over, "raw_p": p, "hit": bool(hit), "reliable": bool(reliable),
            "verdict": ("HIT" if hit else "MISS") if reliable else "UNRELIABLE(N<%d)" % C["min_cal_windows"]}


# ======================================================================================
# PLACEBO / INJECTION
# ======================================================================================
def heldout_thresholds(all_qw, direction, target_far=None):
    C = CONFIG
    target_far = C["target_far"] if target_far is None else target_far
    key = "z2_pos" if direction == "one_sided_positive" else "z2_abs"
    excl_ns = C["held_out_exclusion_days"] * 86400.0 * 1e9
    out = []
    for stn, wins in all_qw.items():
        if not wins:
            continue
        scores = np.array([w[key] for w in wins])
        mids = np.array([w["mid_ts"].value for w in wins], dtype=float)
        for i, w in enumerate(wins):
            keep = np.abs(mids - mids[i]) > excl_ns
            if keep.sum() < C["min_cal_windows"]:
                continue
            out.append((stn, w, k_star_from_scores(scores[keep], target_far)))
    return out


def far_with_ci(thr, direction, rng, n_boot=2000):
    key = "z2_pos" if direction == "one_sided_positive" else "z2_abs"
    hits = np.array([w[key] > k for _, w, k in thr], dtype=bool)
    n = len(hits)
    if n == 0:
        return None
    clusters = {}
    for (stn, w, _), h in zip(thr, hits):
        c = clusters.setdefault((stn, w["year"]), [0, 0])
        c[0] += int(h)
        c[1] += 1
    arr = np.array(list(clusters.values()), dtype=float)
    boots = []
    for _ in range(n_boot):
        s = arr[rng.integers(0, len(arr), len(arr))]
        boots.append(s[:, 0].sum() / s[:, 1].sum())
    lo, hi = np.percentile(boots, [2.5, 97.5])
    wl, wh = wilson_ci(int(hits.sum()), n)
    return {"n": n, "hits": int(hits.sum()), "far_pct": 100.0 * hits.mean(), "wilson": (100 * wl, 100 * wh), "cluster": (100 * lo, 100 * hi)}


def gate_passed(res):
    lo_g, hi_g = CONFIG["placebo_gate_pct"]
    t = 100.0 * CONFIG["target_far"]
    return bool(res["cluster"][0] <= t <= res["cluster"][1] and lo_g <= res["far_pct"] <= hi_g)


def injection_curve(thr, direction, amplitudes, durations, rng, reps=5):
    rows = []
    for dur in durations:
        for amp in amplitudes:
            rec = tot = 0
            for _, w, k in thr:
                z = w["z_values"]
                for _r in range(reps):
                    z2 = z.copy()
                    days = rng.choice(len(z2), size=min(dur, len(z2)), replace=False)
                    z2[days] += amp
                    v = z2 if direction == "one_sided_positive" else np.abs(z2)
                    rec += int(np.sort(v)[-2] > k)
                    tot += 1
            rows.append({"duration_days": dur, "amplitude_sigma": amp, "recovery_pct": 100.0 * rec / max(tot, 1)})
    return pd.DataFrame(rows)


def min_detectable_sigma(inj_df, duration, level):
    sub = inj_df[inj_df["duration_days"] == duration].sort_values("amplitude_sigma")
    a, r = sub["amplitude_sigma"].values, sub["recovery_pct"].values
    if len(a) == 0 or r.max() < level:
        return float("nan")
    i = int(np.argmax(r >= level))
    if i == 0:
        return float(a[0])
    return float(a[i - 1] + (level - r[i - 1]) * (a[i] - a[i - 1]) / max(r[i] - r[i - 1], 1e-9))


# ======================================================================================
# GOVERNANCE FILES (generated, honest)
# ======================================================================================
def _yaml(obj, ind=0):
    pad, lines = "  " * ind, []
    if isinstance(obj, dict):
        for k, v in obj.items():
            if isinstance(v, (dict, list)) and not (isinstance(v, list) and all(not isinstance(x, (dict, list)) for x in v)):
                lines.append(f"{pad}{k}:")
                lines.extend(_yaml(v, ind + 1))
            else:
                lines.append(f"{pad}{k}: {json.dumps(v)}")
    elif isinstance(obj, list):
        for x in obj:
            lines.append(f"{pad}- {json.dumps(x)}")
    return lines


def write_governance(prov_summary, years_info):
    with open("config.yaml", "w", encoding="utf-8") as f:
        f.write("# Generated by step13_final_pipeline.py from its CONFIG dict. Do not edit by hand.\n" + "\n".join(_yaml(CONFIG)) + "\n")
    md = f"""# Event selection, deviations and data status (generated by the script; honest status)

## 1. Event-selection rule: {CONFIG['event_selection_status']}
Proposed rule: {CONFIG['event_selection_proposed_rule']}

The five events (Noto 2024, Hyuga-nada 2024, Tibet-Nepal 2025, Ridgecrest 2019, Turkey 2023) are a CONVENIENCE SAMPLE.
The paper must say so until the rule has been applied to the full USGS list (`--list-candidates`) and every event
that passes has been added, with the failures listed and their reasons.

## 2. Scope of the study
- Observable: GNSS PWV (plus log PWV, gradient magnitude, and a GNSS-minus-ERA5 residual used as a diagnostic).
- VTEC and Sentinel-5P trace gas are NOT part of this pipeline. Title and claims must say PWV only.

## 3. Data status of this run
- NGL archives: {prov_summary}
- Years used per station: {years_info}
- ERA5 total-column water vapour is model-derived. NGL PWV also uses VMF1 weather-model inputs, so GNSS minus ERA5 is a partial,
  not a full, weather removal.

## 4. Deviations from the Version 2 plan
| Item | Version 2 plan | This pipeline |
| --- | --- | --- |
| Events | 3 events | same 3 (4 stations) + 2 added cases (P595, MERS): exploratory |
| Observables | PWV, VTEC, gas | PWV family only |
| Baseline | leave-one-year-out | +/-45-day guard, leap-safe day of year |
| Weather | dewpoint regression | GNSS PWV minus ERA5 TCWV with a {CONFIG['residual_floor_mm']} mm floor (diagnostic) |
| Primary variant | none | {CONFIG['primary_series']}, {CONFIG['primary_calibration']}, {CONFIG['primary_direction']} |
| Placebo | 200 random windows | all quiet windows, held out, with Wilson and year-cluster intervals |

## 5. Known limitations
- 3 independent core earthquakes (Noto has two stations); weak or local effects cannot be excluded.
- The residual floor ({CONFIG['residual_floor_mm']} mm) is a chosen constant; a sensitivity table is printed.
- Distant stations (Taiwan, Guam) lie in different climate regimes: they are context, not controls.
"""
    with open("DEVIATIONS_AND_SELECTION_RULE.md", "w", encoding="utf-8") as f:
        f.write(md)


def write_manifest(extra):
    man = {"utc": now_utc(), "python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
           "config": CONFIG, **extra}
    try:
        man["script_sha256"] = sha256_file(os.path.abspath(__file__))
    except Exception:
        pass
    with open("Step13_run_manifest.json", "w", encoding="utf-8") as f:
        json.dump(man, f, indent=2, default=str)


# ======================================================================================
# MAIN PIPELINE
# ======================================================================================
def _hdr(t):
    print("\n" + "=" * 105 + f"\n{t}\n" + "=" * 105)


def main(force_redownload=False, years_override=None, n_boot=None, verbose_download=True):
    C = CONFIG
    if n_boot is not None:
        C["n_bootstrap"] = n_boot
    rng = np.random.default_rng(C["random_seed"])
    summary = {}
    _hdr("STEP 13: FINAL PIPELINE (guarded baseline, ERA5 residual with floor, event-level statistics, all gates)")

    STATIONS = [dict(c) for c in STATIONS_DEF]
    if years_override:
        for c in STATIONS:
            c["years"] = years_override
    origins = load_origins()
    prov, cats = [], {}

    # ---- epicentres, catalogs, zone membership ----
    for c in STATIONS:
        c["dt"] = pd.Timestamp(c["date"])
        try:
            la, lo, mg, t_utc = mainshock(c)
            c["epi_lat"], c["epi_lon"], c["epi_mag"], c["epi_src"] = la, lo, mg, f"USGS {t_utc:%Y-%m-%d %H:%M} UTC"
            if t_utc.normalize() != c["dt"]:
                raise RuntimeError(f"{c['ev']}: USGS origin time {t_utc} is on a different UTC date than {c['date']}")
        except RuntimeError:
            raise
        except Exception as e:
            c["epi_lat"], c["epi_lon"] = c["epi_fb"]
            c["epi_mag"], c["epi_src"] = c["mw"], f"APPROXIMATE FALLBACK (USGS unavailable: {e}) - verify!"
        c["dist_km"] = float(_haversine_km(c["epi_lat"], c["epi_lon"], c["lat"], c["lon"]))
        c["radius_km"] = 10 ** (0.43 * c["mw"])
        try:
            cats[c["stn"]] = station_catalog(c)
        except Exception as e:
            raise RuntimeError(f"USGS catalog for {c['stn']} could not be fetched or read from {DIR_CAT}: {e}")
    print("\nEPICENTRE / DISTANCE CHECK (computed)")
    for c in STATIONS:
        zone = "inside" if c["dist_km"] <= c["radius_km"] else "OUTSIDE"
        print(f"  {c['ev']:<28} {c['stn']}: {c['dist_km']:7.0f} km from epicentre [{c['epi_src']}], Dobrovolsky R = {c['radius_km']:,.0f} km -> {zone} zone")
    region_dates = {}
    for c in STATIONS:
        region_dates.setdefault(c["region"], set()).add(c["dt"])

    ZS, RESBASE, QW, TC, WVSIG, DAILY = {}, {}, {}, {}, {}, {}
    avail_map, diag_rows, bias_rows, boundary_rows = {}, [], [], []
    variant_rows, raw_dump = [], []

    for c in STATIONS:
        stn = c["stn"]
        y0, y1 = c["years"]
        _hdr(f"{stn}  ({c['ev']})  requested years {y0}-{y1}")
        avail = []
        for y in range(y0, y1 + 1):
            p = ensure_trop_zip(stn, y, force=force_redownload, prov=prov, origins=origins, verbose=verbose_download)
            if p is not None:
                avail.append(y)
        save_origins(origins)
        avail_map[stn] = avail
        print(f"  years available for {stn}: {len(avail)} {avail}")
        if not avail:
            print(f"  [ERROR] no data for {stn}; skipped")
            continue
        daily = parse_station_trop_zips(stn, avail)
        if daily.empty:
            print(f"  [ERROR] no valid daily records for {stn}; skipped")
            continue
        DAILY[stn] = daily
        pwv = daily["PWV_mean"].astype(float)
        WVSIG[stn] = float(daily["PWV_sig_med"].median())
        floor_eff = max(2.0 * WVSIG[stn], C["residual_floor_mm"])
        z = {S_RAW: guarded_baseline(pwv), S_LOG: guarded_baseline(np.log(pwv.clip(lower=0.5))),
             S_GRAD: guarded_baseline(daily["G_mean"].astype(float))}
        for b in boundary_continuity(z[S_RAW]):
            boundary_rows.append({"station": stn, **b})
        bad = [b for b in boundary_continuity(z[S_RAW]) if not b["pass"]]
        print(f"  baseline boundary check: {len(boundary_continuity(z[S_RAW])) - len(bad)} pass, {len(bad)} fail"
              + ("".join(f"\n    FAIL {b['year_boundary']}: spread change {b['max_spread_change_pct']:.1f}% at {b['where']}, identical run {b['longest_identical_run']}" for b in bad)))

        tc, note = load_era5_tcwv(stn, c["lat"], c["lon"], f"{min(avail)}-01-01", f"{max(avail)}-12-31")
        print(f"  ERA5 TCWV: {note}")
        if tc is not None:
            common = pwv.dropna().index.intersection(tc.dropna().index)
            if len(common) >= 300:
                diff = pwv.loc[common] - tc.loc[common]
                bias, rms = float(diff.median()), float(np.sqrt(np.mean(diff ** 2)))
                spread_mm = float(1.4826 * np.median(np.abs(diff - bias)))
                corr_raw = float(pwv.loc[common].corr(tc.loc[common]))
                z_era = guarded_baseline(tc.loc[common])
                cz = z[S_RAW]["Z"].dropna().index.intersection(z_era["Z"].dropna().index)
                corr_anom = float(z[S_RAW].loc[cz, "Z"].corr(z_era.loc[cz, "Z"]))
                if corr_raw < C["min_corr_gnss_vs_era5_tcwv"]:
                    print(f"  [WARNING] raw GNSS-ERA5 correlation {corr_raw:.3f} < {C['min_corr_gnss_vs_era5_tcwv']}: residual series SKIPPED")
                else:
                    base = guarded_baseline(diff)
                    RESBASE[stn] = base
                    z[S_RES] = apply_floor(base, floor_eff)
                    TC[stn] = (tc, z_era)
                    pre = lambda zz: zz.loc[c["dt"] - pd.Timedelta(days=10): c["dt"] - pd.Timedelta(days=1), "Z"].dropna().values
                    pu, pf = pre(base), pre(z[S_RES])
                    pdiff = diff.loc[c["dt"] - pd.Timedelta(days=10): c["dt"] - pd.Timedelta(days=1)]
                    diag_rows.append({"Station": stn, "N_Days": len(common), "Raw_Corr_r": corr_raw, "Anomaly_Corr_rZ": corr_anom,
                                      "Median_Bias_mm": bias, "RMS_Diff_mm": rms, "Residual_Spread_mm": spread_mm,
                                      "Effective_Floor_mm": floor_eff, "Max_PreEvent_AbsDiff_mm": float(pdiff.abs().max()) if len(pdiff) else np.nan,
                                      "Unfloored_Z2": float(np.sort(pu)[-2]) if len(pu) >= 2 else np.nan,
                                      "Floored_Z2": float(np.sort(pf)[-2]) if len(pf) >= 2 else np.nan})
                    print(f"  residual audit: N={len(common)}d raw r={corr_raw:.4f} anomaly r={corr_anom:.4f} RMS={rms:.2f} mm "
                          f"spread={spread_mm:.2f} mm floor={floor_eff:.2f} mm")
                    yb = diff.groupby(diff.index.year).median()
                    bias_rows.append({"Station": stn, **{str(k): float(v) for k, v in yb.items()},
                                      "Range_mm": float(yb.max() - yb.min()),
                                      "FLAG": bool((yb.max() - yb.min()) > C["yearly_bias_flag_mm"])})
        ZS[stn] = z
        qw = {name: build_quiet_windows(c["lat"], c["lon"], zdf, cats[stn], avail, region_dates[c["region"]]) for name, zdf in z.items()}
        QW[stn] = common_windows(qw)
        print(f"  quiet windows (common across series): {len(QW[stn][S_RAW])}")
        for name, zdf in z.items():
            for cal in ["held_out_all", "held_out_season45d", "in_sample_all"]:
                for dr in DIRECTIONS:
                    if name in (S_LOG, S_GRAD) and (cal != "held_out_all" or dr != C["primary_direction"]):
                        continue
                    r = evaluate_variant(zdf, QW[stn][name], c["dt"], cal, dr, rng=rng)
                    if r is None:
                        continue
                    primary = (name == C["primary_series"] and cal == C["primary_calibration"] and dr == C["primary_direction"])
                    variant_rows.append({"Event": c["ev"], "Station": stn, "Role": c["role"], "Series": name, "Calibration": cal,
                                         "Direction": "1-sided(+)" if dr == "one_sided_positive" else "2-sided(|Z|)",
                                         "Primary": "PRIMARY" if primary else "exploratory", "N_cal": r["N_cal"], "k_star": r["k_star"],
                                         "FAR_ach": r["FAR_achieved"], "k_CI95": f"[{r['k_CI95'][0]:.2f}, {r['k_CI95'][1]:.2f}]",
                                         "Stat_Z2": r["ev_stat"], "Days>k": r["days_over"], "Emp_p": r["raw_p"], "Verdict": r["verdict"]})
        zr = z[S_RAW].loc[c["dt"] - pd.Timedelta(days=10): c["dt"] - pd.Timedelta(days=1)].copy()
        zr.insert(0, "event", c["ev"])
        zr.insert(1, "station", stn)
        if stn in TC:
            tcs, zer = TC[stn]
            zr["ERA5_TCWV_mm"] = tcs.reindex(zr.index)
            zr["ERA5_Z"] = zer["Z"].reindex(zr.index)
            zr["GNSS_minus_ERA5_mm"] = zr["raw"] - zr["ERA5_TCWV_mm"]
            zr["Z_Res_Floored"] = z[S_RES]["Z"].reindex(zr.index)
        raw_dump.append(zr.reset_index().rename(columns={"index": "date_utc"}))
        gc.collect()

    if not ZS:
        raise RuntimeError("No station could be processed. Check network access and the 03_Troposphere folder.")

    # ---------------- GATES ----------------
    _hdr("GATE 1: BASELINE CONTINUITY AT EVERY YEAR BOUNDARY (Raw PWV)")
    bdf = pd.DataFrame(boundary_rows)
    n_b_fail = int((~bdf["pass"]).sum()) if len(bdf) else 0
    print(f"  boundaries checked: {len(bdf)} | failures: {n_b_fail} | worst spread change: {bdf['max_spread_change_pct'].max():.1f}% (gate {C['boundary_gate_pct']}%)")
    bdf.to_csv("Step13_Boundary_Checks.csv", index=False)
    summary["boundary_failures"] = n_b_fail

    _hdr("GATE 2: HELD-OUT PLACEBO FALSE-ALARM RATE (all quiet windows; PASS = target inside the year-cluster interval AND inside 3-7%)")
    gate_res = {}
    for name in (S_RAW, S_RES, S_GRAD, S_LOG):
        for dr in DIRECTIONS:
            all_qw = {s: QW[s][name] for s in QW if name in QW[s]}
            if not all_qw:
                continue
            res = far_with_ci(heldout_thresholds(all_qw, dr), dr, rng)
            if res is None:
                continue
            ok = gate_passed(res)
            gate_res[(name, dr)] = ok
            print(f"  {name:<26} {dr:<19} FAR {res['far_pct']:5.2f}% (hits {res['hits']}/{res['n']}; Wilson {res['wilson'][0]:.1f}-{res['wilson'][1]:.1f}%, "
                  f"year-cluster {res['cluster'][0]:.1f}-{res['cluster'][1]:.1f}%) -> {'PASS' if ok else 'FAIL'}")
    summary["placebo"] = {f"{k[0]}|{k[1]}": v for k, v in gate_res.items()}

    _hdr("GATE 3: INJECTION SENSITIVITY (recovery % of +A sigma on D random days) AND MINIMUM DETECTABLE ANOMALY IN mm")
    amps = [0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6]
    durs = [1, 2, 3, 5, 7]
    mdet = {}
    for name in (S_RAW, S_RES):
        all_qw = {s: QW[s][name] for s in QW if name in QW[s]}
        if not all_qw:
            continue
        thr = heldout_thresholds(all_qw, C["primary_direction"])
        inj = injection_curve(thr, C["primary_direction"], amps, durs, rng)
        inj.to_csv(f"Step13_Injection_{'raw' if name == S_RAW else 'residual'}.csv", index=False)
        print(f"\n  {name}: recovery % (rows = amplitude in sigma, columns = duration in days)")
        print(inj.pivot(index="amplitude_sigma", columns="duration_days", values="recovery_pct").to_string(float_format=lambda v: f"{v:5.1f}"))
        spreads = [float(ZS[s][name]["mad"].median()) for s in ZS if name in ZS[s]]
        sp_mm = float(np.median(spreads))
        a3 = min_detectable_sigma(inj, 3, C["detect_level_pct"])
        mdet[name] = (a3, a3 * sp_mm)
        print(f"  -> {C['detect_level_pct']:.0f}% recovery for a 3-day anomaly needs about {a3:.2f} sigma = about {a3 * sp_mm:.1f} mm "
              f"(typical baseline spread {sp_mm:.2f} mm)")
    summary["min_detectable_3day"] = {k: {"sigma": v[0], "mm": v[1]} for k, v in mdet.items()}

    # ---------------- RESULTS ----------------
    _hdr("RESIDUAL AUDIT (GNSS PWV vs ERA5 total-column water vapour)")
    df_diag = pd.DataFrame(diag_rows)
    if len(df_diag):
        print(df_diag.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
        df_diag.to_csv("Step13_Residual_Audit.csv", index=False)
    _hdr("YEARLY GNSS-minus-ERA5 BIAS (a jump between years would mean a product step; FLAG = range above %.1f mm)" % C["yearly_bias_flag_mm"])
    if bias_rows:
        dfb = pd.DataFrame(bias_rows)
        print(dfb.to_string(index=False, float_format=lambda v: f"{v:6.2f}"))
        dfb.to_csv("Step13_Yearly_Bias.csv", index=False)
    summary["yearly_bias_flags"] = [b["Station"] for b in bias_rows if b["FLAG"]]

    df = pd.DataFrame(variant_rows)
    _hdr("MASTER VARIANT TABLE (PRIMARY = the pre-declared test; everything else is exploratory)")
    print(df.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    df.to_csv("Step13_Variant_Table.csv", index=False)

    prim = df[df["Primary"] == "PRIMARY"].copy()
    _hdr("PRIMARY FAMILY: STATION LEVEL AND EVENT LEVEL (stations of one earthquake are collapsed)")
    ev_rows = []
    for ev, sub in prim.groupby("Event", sort=False):
        n = len(sub)
        ev_rows.append({"Event": ev, "Role": sub["Role"].iloc[0], "Stations": ",".join(sub["Station"]), "n_stations": n,
                        "min_p_station": sub["Emp_p"].min(), "p_event": min(1.0, sub["Emp_p"].min() * n),
                        "hit_event": bool((sub["Verdict"] == "HIT").any())})
    dfe = pd.DataFrame(ev_rows)
    fam_out = {}
    for label, mask in (("Core events", dfe["Role"] == "Core"), ("All events", dfe["Role"].notna())):
        sub = dfe[mask].copy()
        if sub.empty:
            continue
        sub["Holm_p"] = holm_adjust(sub["p_event"].values)
        hits = int(sub["hit_event"].sum())
        lo, hi = clopper_pearson(hits, len(sub))
        print(f"\n  {label}: hits {hits}/{len(sub)}, 95% CI for the event hit rate {100 * lo:.0f}-{100 * hi:.0f}%, "
              f"min event p {sub['p_event'].min():.3f}, min Holm p {sub['Holm_p'].min():.3f}")
        print(sub[["Event", "Stations", "min_p_station", "p_event", "Holm_p", "hit_event"]].to_string(index=False, float_format=lambda v: f"{v:.3f}"))
        fam_out[label] = {"hits": hits, "n": len(sub), "ci": [lo, hi], "min_holm_p": float(sub["Holm_p"].min())}
    dfe.to_csv("Step13_Event_Level_Primary.csv", index=False)
    summary["families"] = fam_out
    print("\n  Station-level primary rows:")
    print(prim[["Event", "Station", "N_cal", "k_star", "k_CI95", "Stat_Z2", "Days>k", "Emp_p", "Verdict"]].to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    reliable = df[~df["Verdict"].str.startswith("UNRELIABLE")]
    nh = int((reliable["Verdict"] == "HIT").sum())
    print(f"\n  MULTIVERSE (descriptive, variants are correlated): {nh} HIT among {len(reliable)} reliable variants = "
          f"{100 * nh / max(len(reliable), 1):.1f}% (nominal {100 * C['target_far']:.0f}%); unreliable variants: {len(df) - len(reliable)}")
    summary["multiverse"] = {"hits": nh, "reliable": int(len(reliable))}

    _hdr("RESIDUAL FLOOR SENSITIVITY (is the residual conclusion an artefact of the floor?)")
    frows = []
    for fl in C["residual_floor_sensitivity_mm"]:
        hit_s, minp, all_q = [], 1.0, {}
        for c in STATIONS:
            stn = c["stn"]
            if stn not in RESBASE:
                continue
            zf = apply_floor(RESBASE[stn], max(2.0 * WVSIG[stn], fl))
            q = build_quiet_windows(c["lat"], c["lon"], zf, cats[stn], avail_map[stn], region_dates[c["region"]])
            all_q[stn] = q
            r = evaluate_variant(zf, q, c["dt"], "held_out_all", C["primary_direction"], n_boot=0, rng=rng)
            if r:
                minp = min(minp, r["raw_p"])
                if r["hit"] and r["reliable"]:
                    hit_s.append(stn)
        res = far_with_ci(heldout_thresholds(all_q, C["primary_direction"]), C["primary_direction"], rng, n_boot=500) if all_q else None
        frows.append({"Floor_mm": fl, "Stations_hit": ",".join(hit_s) or "-", "Min_p": minp,
                      "Placebo_FAR_pct": res["far_pct"] if res else np.nan,
                      "Placebo_gate": ("PASS" if res and gate_passed(res) else "FAIL") if res else "n/a"})
    dff = pd.DataFrame(frows)
    print(dff.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    dff.to_csv("Step13_Floor_Sensitivity.csv", index=False)

    _hdr("SENSITIVITY GRID BY EVENT (an event hits if ANY of its stations hits; expectation uses 1-(1-FAR)^stations per event)")
    ev_names = []
    for c in STATIONS:
        if c["stn"] in ZS and c["ev"] not in ev_names:
            ev_names.append(c["ev"])
    n_per_ev = {e: sum(1 for c in STATIONS if c["ev"] == e and c["stn"] in ZS) for e in ev_names}
    srows = []
    for wl in C["window_lengths_tested"]:
        qwl = {c["stn"]: build_quiet_windows(c["lat"], c["lon"], ZS[c["stn"]][S_RAW], cats[c["stn"]], avail_map[c["stn"]],
                                              region_dates[c["region"]], win_len=wl) for c in STATIONS if c["stn"] in ZS}
        for tf in C["target_fars_tested"]:
            hit_stn, hit_ev = [], set()
            for c in STATIONS:
                if c["stn"] not in ZS:
                    continue
                r = evaluate_variant(ZS[c["stn"]][S_RAW], qwl[c["stn"]], c["dt"], "held_out_all", C["primary_direction"], tf, wl, n_boot=0, rng=rng)
                if r and r["reliable"] and r["hit"]:
                    hit_stn.append(c["stn"])
                    hit_ev.add(c["ev"])
            ps = [1.0 - (1.0 - tf) ** n_per_ev[e] for e in ev_names]
            srows.append({"Window_days": wl, "Target_FAR_pct": 100 * tf, "Events": len(ev_names), "Event_hits": len(hit_ev),
                          "Expected_hits": round(sum(ps), 2), "p_PoissonBinomial": round(poisson_binomial_tail(len(hit_ev), ps), 3),
                          "Stations_hit": ",".join(hit_stn) or "-"})
    dfs = pd.DataFrame(srows)
    print(dfs.to_string(index=False))
    dfs.to_csv("Step13_Sensitivity_Grid.csv", index=False)
    print("  The 10-day / 5% row is the pre-declared one. Other rows are sensitivity checks, not extra chances to find a hit.")

    # ---------------- WEATHER EXPLANATION TABLE ----------------
    _hdr("WEATHER EXPLANATION: do ERA5 anomalies reproduce the GNSS anomalies in the 10 days before each event?")
    for c in STATIONS:
        stn = c["stn"]
        if stn not in TC:
            continue
        tcs, zer = TC[stn]
        z_raw = ZS[stn][S_RAW]
        w = z_raw.loc[c["dt"] - pd.Timedelta(days=10): c["dt"] - pd.Timedelta(days=1)]
        tab = pd.DataFrame({"GNSS_PWV_mm": w["raw"], "GNSS_Z": w["Z"], "ERA5_TCWV_mm": tcs.reindex(w.index), "ERA5_Z": zer["Z"].reindex(w.index)})
        tab["GNSS-ERA5_mm"] = tab["GNSS_PWV_mm"] - tab["ERA5_TCWV_mm"]
        tab["Z_resid_floored"] = ZS[stn][S_RES]["Z"].reindex(w.index)
        tab.index = tab.index.strftime("%Y-%m-%d")
        print(f"\n  {c['ev']} | {stn}")
        print(tab.to_string(float_format=lambda v: f"{v:7.2f}"))

    # ---------------- REGIONAL COINCIDENCE (context, not controls) ----------------
    _hdr("REGIONAL COINCIDENCE AND DISTANCE CONTEXT")
    print("  Stations in other climate regimes (Taiwan, Guam) are NOT negative controls: a mid-latitude front over Japan")
    print("  does not reach them. Weather attribution rests on the ERA5 table above, not on these stations.")
    extra_z, extra_info = {}, []
    for ex in EXTRA_STATIONS:
        e_av = []
        for y in range(ex["years"][0], ex["years"][1] + 1):
            if ensure_trop_zip(ex["stn"], y, prov=prov, origins=origins, verbose=False) is not None:
                e_av.append(y)
        save_origins(origins)
        d = parse_station_trop_zips(ex["stn"], e_av) if e_av else pd.DataFrame()
        if len(d):
            extra_z[ex["stn"]] = guarded_baseline(d["PWV_mean"].astype(float))
            extra_info.append((ex, e_av))
    seen = set()
    for c in STATIONS:
        if c["ev"] in seen or c["region"] != "Japan":
            continue
        seen.add(c["ev"])
        R = c["radius_km"]
        members = {s["stn"]: ZS[s["stn"]][S_RAW] for s in STATIONS if s["region"] == "Japan" and s["stn"] in ZS}
        pos = {s["stn"]: (s["lat"], s["lon"]) for s in STATIONS if s["region"] == "Japan"}
        for ex, e_av in extra_info:
            members[ex["stn"]] = extra_z[ex["stn"]]
            pos[ex["stn"]] = (ex["lat"], ex["lon"])
        print(f"\n  {c['ev']} (R = {R:,.0f} km)")
        cols = {}
        for s, zdf in members.items():
            dist = float(_haversine_km(c["epi_lat"], c["epi_lon"], pos[s][0], pos[s][1]))
            cols[s] = zdf.loc[c["dt"] - pd.Timedelta(days=10): c["dt"] - pd.Timedelta(days=1), "Z"]
            print(f"    {s}: {dist:7.0f} km  {'INSIDE zone' if dist <= R else 'OUTSIDE zone'}")
        tab = pd.concat(cols, axis=1, sort=True)
        tab.index = tab.index.strftime("%Y-%m-%d")
        print(tab.to_string(float_format=lambda v: f"{v:6.2f}"))
        inside = [s for s in members if float(_haversine_km(c["epi_lat"], c["epi_lon"], pos[s][0], pos[s][1])) <= R and s in ZS]
        if len(inside) >= 2:
            allz = pd.concat({s: members[s]["Z"] for s in inside}, axis=1, sort=True)
            reg = allz.mean(axis=1).where(allz.notna().sum(axis=1) >= 2)
            pre = reg.loc[c["dt"] - pd.Timedelta(days=10): c["dt"] - pd.Timedelta(days=1)].dropna()
            exd = pd.Timedelta(days=C["event_exclusion_days"])
            evd = list(region_dates["Japan"])
            null, cur = [], reg.dropna().index.min()
            last = reg.dropna().index.max()
            while cur + pd.Timedelta(days=9) <= last:
                w1 = cur + pd.Timedelta(days=9)
                if not any((cur <= d + exd) and (w1 >= d - exd) for d in evd):
                    seg = reg.loc[cur:w1].dropna()
                    if len(seg) >= 8:
                        null.append(float(np.sort(seg.values)[-2]))
                cur = cur + pd.Timedelta(days=10)
            if len(pre) >= 8 and len(null) >= 20:
                st = float(np.sort(pre.values)[-2])
                print(f"    regional-mean Z(2) over inside-zone stations {inside} = {st:.2f}; empirical p = "
                      f"{(np.sum(np.array(null) >= st) + 1) / (len(null) + 1):.3f} over {len(null)} quiet windows")

    # ---------------- FIGURE ----------------
    fig, axes = plt.subplots(len(ZS), 1, figsize=(9, 2.2 * len(ZS)))
    for ax, c in zip(np.atleast_1d(axes), [s for s in STATIONS if s["stn"] in ZS]):
        zr = ZS[c["stn"]][S_RAW]
        y = c["dt"].year if (c["dt"].month > 1 or c["dt"].day > 20) else c["dt"].year - 1
        seg = zr.loc[f"{y}-12-10":f"{y + 1}-01-20"]
        if len(seg):
            ax.plot(seg.index, seg["mad"], "s-", ms=3, label="baseline spread (mm)")
            ax.plot(seg.index, seg["med"], "--", label="baseline median (mm)")
            ax.axvline(pd.Timestamp(f"{y + 1}-01-01"), color="gray", ls=":")
        ax.set_title(f"{c['stn']}: baseline across {y}/{y + 1}", fontsize=9)
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        ax.legend(fontsize=7)
    plt.tight_layout()
    plt.savefig("Figure_YearBoundary_All_Stations.png", dpi=200)
    plt.close()

    # ---------------- OUTPUT FILES, PROVENANCE, CHECKLIST ----------------
    if raw_dump:
        pd.concat(raw_dump, sort=False).to_csv("Step13_Raw_Window_Dumps.csv", index=False)
    dfp = pd.DataFrame(prov)
    os.makedirs(DIR_TROP, exist_ok=True)
    dfp.to_csv(os.path.join(DIR_TROP, "PROVENANCE.csv"), index=False)
    n_dl = int(dfp["origin"].str.startswith("http").sum()) if len(dfp) else 0
    n_ver = int(dfp["origin"].str.startswith("verified_identical").sum()) if len(dfp) else 0
    n_unv = int(dfp["origin"].str.startswith("local_unverified").sum()) if len(dfp) else 0
    n_unc = int((dfp["origin"] == "UNCERTAIN_FETCH_ERROR").sum()) if len(dfp) else 0
    prov_summary = f"{n_dl} downloaded this run, {n_ver} verified identical to NGL, {n_unv} local files UNVERIFIED, {n_unc} fetches uncertain"
    years_info = "; ".join(f"{s}: {len(v)} yr" for s, v in avail_map.items())
    write_governance(prov_summary, years_info)
    write_manifest({"provenance": prov, "provenance_summary": prov_summary, "years": avail_map, "summary": summary,
                    "stations": [{k: (str(v) if k == "dt" else v) for k, v in c.items()} for c in STATIONS]})

    _hdr("FINAL CHECKLIST (computed from this run; nothing here is hand-written)")
    fam_core = fam_out.get("Core events")
    short = [s for s, v in avail_map.items() if len(v) < 8]
    print(f"  [{'OK ' if n_unv == 0 and n_unc == 0 else 'TODO'}] provenance: {prov_summary}"
          + ("" if n_unv == 0 else "  -> run: python step13_final_pipeline.py --verify-provenance"))
    print(f"  [{'OK ' if n_b_fail == 0 else 'TODO'}] baseline boundary gate: {n_b_fail} failing boundaries")
    for (name, dr), ok in gate_res.items():
        if dr == C["primary_direction"]:
            print(f"  [{'OK ' if ok else 'TODO'}] placebo gate {name}: {'PASS' if ok else 'FAIL -> treat this series as diagnostic only, make no claim from it'}")
    print(f"  [{'OK ' if not short else 'WARN'}] history length: stations with < 8 years: {short or 'none'}")
    print(f"  [{'OK ' if not summary['yearly_bias_flags'] else 'WARN'}] yearly GNSS-ERA5 bias flags: {summary['yearly_bias_flags'] or 'none'}")
    print(f"  [TODO] event-selection rule: {C['event_selection_status']}")
    print("  [INFO] scope: PWV only (no VTEC, no gas). 3 independent core earthquakes.")
    if fam_core:
        print(f"\n  RESULT: {fam_core['hits']} of {fam_core['n']} core events exceeded the locked threshold (primary variant); "
              f"smallest Holm-adjusted event p = {fam_core['min_holm_p']:.3f}; 95% CI for the event hit rate "
              f"{100 * fam_core['ci'][0]:.0f}-{100 * fam_core['ci'][1]:.0f}%.")
    for name, (s3, mm3) in mdet.items():
        print(f"  DETECTABILITY ({name}): a 3-day anomaly is caught {C['detect_level_pct']:.0f}% of the time only above ~{s3:.1f} sigma (~{mm3:.1f} mm).")
    print("  A null result here means: no anomaly above these sizes before these events. It does not prove that weaker or")
    print("  local effects do not exist.")
    return summary


# ======================================================================================
# SELF-TEST (offline)
# ======================================================================================
def _synthetic_daily(seed, start="2015-01-01", end="2025-12-31", storm=None, base=15.0):
    idx = pd.date_range(start, end, freq="D")
    rng = np.random.default_rng(seed)
    doy = norm_doy(idx)
    seasonal = base + 10 * np.sin(2 * np.pi * (doy - 100) / 365.0)
    sigma = 1.5 + 1.0 * (1 + np.sin(2 * np.pi * (doy - 100) / 365.0))
    e = np.zeros(len(idx))
    for i in range(1, len(idx)):
        e[i] = 0.6 * e[i - 1] + rng.normal(0, 1)
    s = pd.Series(np.clip(seasonal + sigma * e * 0.8, 0.5, None), index=idx)
    for d, amp in (storm or {}).items():
        s.loc[d] += amp
    return s


def _fake_trop_zip_bytes(stn, year, days=3, pad=12000, with_gz=True):
    """Builds a small NGL-like archive (header and column order as in the real files) for parser tests."""
    lines = ["%=TRO 2.00 NGL", "+TROP/SOLUTION",
             "*SITE ___EPOCH____ TROTOT _SIG TRWET TGETOT _SIG TGNTOT _SIG WVAPOR _SIG MTEMP"]
    yy = year % 100
    for d in range(1, days + 1):
        for k in range(288):
            sod = k * 300
            lines.append(f"{stn} {yy:02d}:{d:03d}:{sod:05d} 2329.7 4.2 33.3 0.19 0.54 -1.08 0.74 {5.0 + 0.01 * d:.2f} 0.63 265.5")
    lines.append("-TROP/SOLUTION")
    text = ("\n".join(lines) + "\n").encode("ascii")
    member = f"{stn}.{year}.001.trop" + (".gz" if with_gz else "")
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr(member, gzip.compress(text) if with_gz else text)
        z.writestr("padding.bin", os.urandom(pad))
    return buf.getvalue()


def selftest():
    ok_all = True

    def check(name, cond, info=""):
        nonlocal ok_all
        ok_all &= bool(cond)
        print(f"  [{'PASS' if cond else 'FAIL'}] {name} {info}")

    print("SELF-TEST (offline)")
    # leap-year day of year
    d = pd.DatetimeIndex(["2024-12-31", "2025-01-01", "2024-02-29", "2024-03-01", "2023-03-01"])
    nd = norm_doy(d)
    check("Dec 31 2024 vs Jan 1 2025 are 1 day apart on the circle", circ_dist(nd[0], nd[1]) == 1)
    check("Mar 1 has the same normalised day in leap and non-leap years", nd[3] == nd[4])
    # statistics helpers
    check("Holm adjust known values", np.allclose(holm_adjust([0.01, 0.04, 0.03]), [0.03, 0.06, 0.06]))
    check("binomial tail P(X>=2|n=6,p=.05) = 0.0328", abs(binom_tail(2, 6, 0.05) - 0.0328) < 5e-4)
    check("Poisson-binomial equals binomial for equal p", abs(poisson_binomial_tail(2, [0.05] * 6) - binom_tail(2, 6, 0.05)) < 1e-9)
    check("Poisson-binomial with unequal p", abs(poisson_binomial_tail(1, [0.1, 0.2]) - (1 - 0.9 * 0.8)) < 1e-9)
    lo, hi = clopper_pearson(0, 4)
    check("Clopper-Pearson 0/4 upper bound about 0.60", abs(hi - 0.602) < 0.01, f"(hi={hi:.3f})")
    lo, hi = clopper_pearson(0, 3)
    check("Clopper-Pearson 0/3 upper bound about 0.708", abs(hi - 0.708) < 0.01, f"(hi={hi:.3f})")
    sc = np.random.default_rng(0).normal(size=137)
    check("k* gives FAR <= 5%", np.mean(sc > k_star_from_scores(sc, 0.05)) <= 0.05)
    # baseline and floor
    s = _synthetic_daily(1, storm={"2024-12-28": 25.0, "2024-12-29": 12.0})
    z = guarded_baseline(s)
    b24 = [b for b in boundary_continuity(z) if b["year_boundary"] == "2024/2025"][0]
    check("year-boundary continuity (storm on Dec 28 does not leak into January)", b24["pass"], f"(max change {b24['max_spread_change_pct']:.1f}%)")
    check("storm day is flagged", z.loc["2024-12-28", "Z"] > 4)
    check("no fallback: every Z has n_ref >= min_ref", (z.dropna(subset=["Z"])["n_ref"] >= CONFIG["min_ref_values"]).all())
    zf = apply_floor(z, 5.0)
    check("apply_floor never shrinks the spread", (zf["mad"].dropna() >= 5.0 - 1e-12).all() and (zf["mad"].dropna() >= z["mad"].dropna() - 1e-12).all())
    check("apply_floor leaves large spreads unchanged", np.allclose(apply_floor(z, 0.0)["Z"].dropna(), z["Z"].dropna()))
    # quiet windows
    cat = pd.DataFrame({"time_dt": [pd.Timestamp("2021-06-15")], "latitude": [10.0], "longitude": [10.0], "mag": [6.0]})
    yrs = list(range(2015, 2026))
    qw = build_quiet_windows(10.0, 10.5, z, cat, yrs, [pd.Timestamp("2024-08-08")])
    bad_cat = [w for w in qw if w["start"] - pd.Timedelta(days=10) <= pd.Timestamp("2021-06-15") <= w["end"] + pd.Timedelta(days=10)]
    bad_ev = [w for w in qw if abs((w["mid_ts"] - pd.Timestamp("2024-08-08")).days) <= 25]
    check("catalog event excludes windows +/-10 days", len(bad_cat) == 0)
    check("other target events exclude windows +/-30 days", len(bad_ev) == 0)
    check("most windows remain", len(qw) > 330, f"(N={len(qw)})")
    # held-out placebo on null data
    allqw = {f"S{i}": build_quiet_windows(0, 0, guarded_baseline(_synthetic_daily(10 + i)), None, yrs, []) for i in range(3)}
    rng = np.random.default_rng(5)
    thr = heldout_thresholds(allqw, "one_sided_positive")
    res = far_with_ci(thr, "one_sided_positive", rng, n_boot=300)
    check("held-out FAR close to 5% on null data", 3.5 <= res["far_pct"] <= 7.0, f"(FAR={res['far_pct']:.2f}%, n={res['n']})")
    inj = injection_curve(thr, "one_sided_positive", [1, 3, 5], [1, 3, 5], rng, reps=2)
    g = lambda a, d_: inj[(inj.amplitude_sigma == a) & (inj.duration_days == d_)]["recovery_pct"].iloc[0]
    check("injection recovery grows with amplitude and duration", g(1, 3) < g(3, 3) < g(5, 5))
    check("min detectable amplitude is finite and decreases with duration", 0 < min_detectable_sigma(inj, 5, 80) <= min_detectable_sigma(inj, 3, 80) + 1e-9
          or np.isnan(min_detectable_sigma(inj, 3, 80)))
    # event evaluation
    zi = guarded_baseline(_synthetic_daily(77))
    qwi = build_quiet_windows(0, 0, zi, None, yrs, [pd.Timestamp("2023-09-15")])
    r0 = evaluate_variant(zi, qwi, "2023-09-15", n_boot=100)
    check("evaluate_variant works on a null event", r0 is not None)
    s2 = _synthetic_daily(77, storm={"2023-09-07": 12, "2023-09-08": 12, "2023-09-09": 12})
    z2 = guarded_baseline(s2)
    r2 = evaluate_variant(z2, build_quiet_windows(0, 0, z2, None, yrs, [pd.Timestamp("2023-09-15")]), "2023-09-15", n_boot=100)
    check("a planted 3-day anomaly is a HIT", r2["hit"], f"(Z2={r2['ev_stat']:.2f}, k*={r2['k_star']:.2f})")
    # parser on a synthetic archive with the real header layout
    tmp = tempfile.mkdtemp()
    try:
        with open(os.path.join(tmp, "ABCD.2021.trop.zip"), "wb") as f:
            f.write(_fake_trop_zip_bytes("ABCD", 2021, days=3))
        dd = parse_station_trop_zips("ABCD", [2021], trop_dir=tmp)
        check("parser: 3 valid days, PWV by column NAME", len(dd) == 3 and abs(dd["PWV_mean"].iloc[0] - 5.01) < 1e-6, f"(rows={len(dd)})")
        check("parser: sigma and gradient magnitude", abs(dd["PWV_sig_med"].iloc[0] - 0.63) < 1e-9 and abs(dd["G_mean"].iloc[0] - math.hypot(0.19, -1.08)) < 1e-6)
        # downloader with a fake network
        calls = []
        orig = globals()["http_get"]

        def fake_get(url, timeout=120):
            calls.append(url)
            if "/IGS20/trop/AAAA/AAAA.2020" in url:
                return _fake_trop_zip_bytes("AAAA", 2020)
            raise urllib.error.HTTPError(url, 404, "nf", None, None)
        globals()["http_get"] = fake_get
        try:
            origins = {}
            p = ensure_trop_zip("AAAA", 2020, origins=origins, trop_dir=tmp, verbose=False)
            check("downloader: file saved and origin recorded as the URL", p is not None and origins["AAAA.2020.trop.zip"]["origin"].startswith("https://"))
            p2 = ensure_trop_zip("AAAA", 2019, origins=origins, trop_dir=tmp, verbose=False)
            check("downloader: a 404 on all URLs returns None (no crash)", p2 is None)
            prov = []
            ensure_trop_zip("AAAA", 2020, origins=origins, prov=prov, trop_dir=tmp, verbose=False)
            check("existing file keeps its recorded origin", prov[0]["origin"].startswith("https://"))
            prov = []
            with open(os.path.join(tmp, "ZZZZ.2020.trop.zip"), "wb") as f:
                f.write(_fake_trop_zip_bytes("ZZZZ", 2020))
            ensure_trop_zip("ZZZZ", 2020, origins={}, prov=prov, trop_dir=tmp, verbose=False)
            check("unknown local file is labelled local_unverified, never 'verified'", prov[0]["origin"] == "local_unverified")
        finally:
            globals()["http_get"] = orig
        # old fake label is relabelled
        old_dir = os.getcwd()
        os.chdir(tmp)
        try:
            os.makedirs(DIR_TROP, exist_ok=True)
            with open(ORIGINS_FILE, "w") as f:
                json.dump({"x.zip": {"origin": "https://geodesy.unr.edu/gps_timeseries/trop_or_IGS20_verified_local", "sha256": "a"}}, f)
            check("old fake 'verified_local' label is relabelled", "unverified" in load_origins()["x.zip"]["origin"])
        finally:
            os.chdir(old_dir)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print("\nSELF-TEST RESULT:", "ALL PASSED" if ok_all else "SOME FAILED")
    return ok_all


def smoketest():
    """Runs main() end to end on fake data, with the network replaced. Writes into a temporary folder."""
    old = os.getcwd()
    tmp = tempfile.mkdtemp()
    os.chdir(tmp)
    orig_get, orig_parse = globals()["http_get"], globals()["parse_station_trop_zips"]
    saved_cfg = dict(CONFIG)
    try:
        def fake_get(url, timeout=120):
            if "archive-api.open-meteo.com" in url:
                q = urllib.parse.parse_qs(urllib.parse.urlparse(url).query)
                s0, s1 = pd.Timestamp(q["start_date"][0]), pd.Timestamp(q["end_date"][0])
                idx = pd.date_range(s0, s1 + pd.Timedelta(hours=23), freq="h")
                days = pd.DatetimeIndex(idx.normalize().unique())
                base = _synthetic_daily(abs(hash(q["latitude"][0])) % 1000, str(s0.date()), str(s1.date()))
                vals = base.reindex(idx.normalize()).values + np.random.default_rng(1).normal(0, 1.0, len(idx))
                return json.dumps({"hourly": {"time": [t.strftime("%Y-%m-%dT%H:%M") for t in idx],
                                              "total_column_integrated_water_vapour": [float(v) for v in vals]}}).encode()
            if "earthquake.usgs.gov" in url:
                q = urllib.parse.parse_qs(urllib.parse.urlparse(url).query)
                if "latitude" in q:  # station catalog: no big events nearby
                    return b"time,latitude,longitude,depth,mag,place\n"
                a = pd.Timestamp(q["starttime"][0]) + pd.Timedelta(days=1)
                ep = {"2024-01-01": (37.49, 137.27, 7.5), "2024-08-08": (31.76, 131.50, 7.1), "2025-01-07": (28.60, 87.38, 7.1),
                      "2019-07-06": (35.77, -117.60, 7.1), "2023-02-06": (37.22, 37.01, 7.8)}[str(a.date())]
                return (f"time,latitude,longitude,depth,mag,place\n{a:%Y-%m-%d}T03:10:00.000Z,{ep[0]},{ep[1]},10,{ep[2]},test\n").encode()
            if "geodesy.unr.edu" in url:
                if "/IGS20/" in url:
                    return _fake_trop_zip_bytes("TEST", 2020, days=1)
                raise urllib.error.HTTPError(url, 404, "nf", None, None)
            raise urllib.error.URLError("blocked")
        globals()["http_get"] = fake_get

        def fake_parse(stn, years, trop_dir=None):
            ys = sorted(years)
            if stn == "IITK":
                ys = [y for y in ys if y >= 2022]
            seed = abs(hash(stn)) % 1000
            s = _synthetic_daily(seed, f"{ys[0]}-01-01", f"{ys[-1]}-12-31", base=15 + (seed % 10))
            g = _synthetic_daily(seed + 7, f"{ys[0]}-01-01", f"{ys[-1]}-12-31") / 10
            return pd.DataFrame({"ZWD_mean": s * 6, "PWV_mean": s, "PWV_sig_med": 0.65, "G_mean": g, "count": 288}, index=s.index)
        globals()["parse_station_trop_zips"] = fake_parse
        os.makedirs(DIR_TROP, exist_ok=True)
        # IITK years before 2022 must look unavailable: the fake network serves everything, so patch the zip existence
        orig_ensure = globals()["ensure_trop_zip"]

        def fake_ensure(stn, year, force=False, prov=None, origins=None, trop_dir=None, verbose=True):
            if stn == "IITK" and year < 2022:
                if prov is not None:
                    prov.append({"file": f"{stn}.{year}.trop.zip", "origin": "NOT_AVAILABLE", "bytes": 0, "sha256": "", "utc": now_utc()})
                return None
            return orig_ensure(stn, year, force, prov, origins, trop_dir, verbose=False)
        globals()["ensure_trop_zip"] = fake_ensure
        summ = main(n_boot=100, verbose_download=False)
        need = ["Step13_Variant_Table.csv", "Step13_Event_Level_Primary.csv", "Step13_Sensitivity_Grid.csv", "Step13_Floor_Sensitivity.csv",
                "Step13_Injection_raw.csv", "Step13_Injection_residual.csv", "Step13_Boundary_Checks.csv", "Step13_Residual_Audit.csv",
                "Step13_Yearly_Bias.csv", "Step13_Raw_Window_Dumps.csv", "Step13_run_manifest.json", "config.yaml",
                "DEVIATIONS_AND_SELECTION_RULE.md", "Figure_YearBoundary_All_Stations.png", os.path.join(DIR_TROP, "PROVENANCE.csv")]
        missing = [f for f in need if not os.path.exists(f)]
        print("\nSMOKE TEST: missing output files:", missing or "none")
        ok = not missing and "families" in summ
        print("SMOKE TEST RESULT:", "PASSED" if ok else "FAILED")
        return ok
    finally:
        globals()["http_get"], globals()["parse_station_trop_zips"] = orig_get, orig_parse
        globals()["ensure_trop_zip"] = orig_ensure
        CONFIG.clear()
        CONFIG.update(saved_cfg)
        os.chdir(old)
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--smoketest", action="store_true")
    ap.add_argument("--verify-provenance", action="store_true")
    ap.add_argument("--list-candidates", action="store_true")
    ap.add_argument("--force-redownload", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(0 if selftest() else 1)
    if a.smoketest:
        sys.exit(0 if smoketest() else 1)
    if a.verify_provenance:
        verify_provenance()
        sys.exit(0)
    if a.list_candidates:
        list_candidate_events()
        sys.exit(0)
    main(force_redownload=a.force_redownload)
