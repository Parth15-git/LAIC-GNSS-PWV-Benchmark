import os
import glob
import gzip
import re
import urllib.request
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# ==============================================================================
# CONFIGURATION: JAN 7, 2025 TIBET-NEPAL EARTHQUAKE (Mw 7.1)
# ==============================================================================
EQ_DATE = datetime(2025, 1, 7, 1, 5)  # 01:05 UTC
EPI_LAT = 28.50                       # Degrees North
EPI_LON = 87.45                       # Degrees East

DIR_OMNI  = "01_SpaceWeather"
DIR_IONEX = "02_IONEX"
DIR_TROP  = "03_Troposphere"
DIR_GAS   = "04_Gas_S5P"

# ==============================================================================
# 0. ENSURE VALID IONEX FILES EXIST (RE-DOWNLOADS CORRUPT/EMPTY FILES)
# ==============================================================================
def ensure_ionex_files(folder):
    os.makedirs(folder, exist_ok=True)
    targets = [(2024, d) for d in range(343, 367)] + [(2025, d) for d in range(1, 13)]
    
    for yr, doy in targets:
        fname_gz = f"COD0OPSFIN_{yr}{doy:03d}0000_01D_01H_GIM.INX.gz"
        fname_inx = fname_gz[:-3]
        fpath_gz = os.path.join(folder, fname_gz)
        fpath_inx = os.path.join(folder, fname_inx)

        # Check if a valid (>5 KB) file already exists
        valid_exists = False
        for p in [fpath_gz, fpath_inx]:
            if os.path.exists(p) and os.path.getsize(p) > 5000:
                valid_exists = True
                break
            elif os.path.exists(p) and os.path.getsize(p) <= 5000:
                os.remove(p)  # Remove broken XML error file

        if not valid_exists:
            url = f"https://zhw-b.s3.cloud.switch.ch/aiub/CODE/{yr}/{fname_gz}"
            try:
                print(f"  Downloading {fname_gz} ...")
                urllib.request.urlretrieve(url, fpath_gz)
            except Exception as e:
                print(f"  [Warning] Could not fetch {fname_gz}: {e}")

# ==============================================================================
# 1. PARSE NASA OMNIWEB SPACE WEATHER (.lst / .txt)
# ==============================================================================
def load_omni(folder):
    files = [
        f for f in glob.glob(os.path.join(folder, "*"))
        if not f.lower().endswith(".fmt") and os.path.isfile(f)
    ]
    if not files:
        raise FileNotFoundError(f"No OMNI data file found in {folder}")
    
    records = []
    with open(files[0], "r", errors="ignore") as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 5 and parts[0].isdigit() and len(parts[0]) == 4:
                yr, doy, hr = int(parts[0]), int(parts[1]), int(parts[2])
                dt = datetime.strptime(f"{yr}-{doy:03d}-{hr:02d}", "%Y-%j-%H")
                vals = [float(x) for x in parts[3:]]
                kp = vals[0] / 10.0 if vals[0] < 90 else np.nan
                dst = vals[1] if abs(vals[1]) < 900 else np.nan
                records.append({"datetime": dt, "Kp": kp, "Dst": dst})
    if not records:
        raise ValueError(f"Could not parse hourly rows from OMNI file: {files[0]}")
    df = pd.DataFrame(records).set_index("datetime").sort_index()
    return df.resample("1D").agg({"Kp": "max", "Dst": "min"})

# ==============================================================================
# 2. PARSE CODE IONEX (80-COLUMN 16I5 FORMAT) & INTERPOLATE VTEC
# ==============================================================================
def read_lines_auto(filepath):
    with open(filepath, "rb") as fb:
        magic = fb.read(2)
    if magic == b"\x1f\x8b":
        with gzip.open(filepath, "rt", errors="ignore") as f:
            return f.readlines()
    else:
        with open(filepath, "r", errors="ignore") as f:
            return f.readlines()

def parse_ionex_file(filepath, target_lat, target_lon):
    lines = read_lines_auto(filepath)
    records = []
    in_tec_map = False
    current_epoch = None
    lat_grid = {}
    current_lat = None
    lon1, lon2, dlon = -180.0, 180.0, 5.0
    exponent = -1

    for line in lines:
        label = line[60:].strip() if len(line) > 60 else ""
        has_alpha_label = any(c.isalpha() for c in label)

        if has_alpha_label and "EXPONENT" in label:
            exponent = int(line[:60].strip().split()[0])
        elif has_alpha_label and "START OF TEC MAP" in label:
            in_tec_map = True
            lat_grid = {}
            current_lat = None
        elif has_alpha_label and "END OF TEC MAP" in label and in_tec_map:
            in_tec_map = False
            if current_epoch and lat_grid:
                lats = sorted(lat_grid.keys())
                closest_lat = min(lats, key=lambda x: abs(x - target_lat))
                row = [np.nan if x == 9999 else x for x in lat_grid[closest_lat]]
                if len(row) >= 30:
                    lons = np.linspace(lon1, lon2, len(row))
                    vtec_val = np.interp(target_lon, lons, row) * (10 ** exponent)
                    records.append({"datetime": current_epoch, "VTEC": vtec_val})
        elif in_tec_map:
            if has_alpha_label and "EPOCH OF CURRENT MAP" in label:
                p = [int(x) for x in line[:60].split()[:6]]
                if p[3] == 24:
                    current_epoch = datetime(p[0], p[1], p[2]) + pd.Timedelta(days=1)
                else:
                    current_epoch = datetime(p[0], p[1], p[2], p[3], p[4], p[5])
            elif has_alpha_label and "LAT/LON1/LON2/DLON/H" in label:
                current_lat = float(line[2:8])
                lon1 = float(line[8:14])
                lon2 = float(line[14:20])
                dlon = float(line[20:26])
                lat_grid[current_lat] = []
            elif current_lat is not None and not has_alpha_label:
                # Full 80-column numeric row (16I5 format)
                nums = [int(x) for x in re.findall(r"-?\d+", line)]
                lat_grid[current_lat].extend(nums)

    return records

def load_all_ionex(folder, lat, lon):
    ensure_ionex_files(folder)
    # Deduplicate file list on case-insensitive Windows filesystem
    inx_files = sorted(list({
        os.path.abspath(f) for f in glob.glob(os.path.join(folder, "*"))
        if ".inx" in f.lower()
    }))
    if not inx_files:
        raise FileNotFoundError(f"No .INX or .INX.gz files found in {folder}.")

    all_recs = []
    for f in inx_files:
        all_recs.extend(parse_ionex_file(f, lat, lon))

    print(f"  Parsed {len(inx_files)} IONEX files ({len(all_recs)} hourly VTEC maps).")
    if not all_recs:
        raise RuntimeError("0 VTEC maps extracted. Check that the .INX files are valid CODE GIM files.")

    df = pd.DataFrame(all_recs).drop_duplicates("datetime").set_index("datetime").sort_index()
    return df.resample("1D").max()

# ==============================================================================
# 3. PARSE GOOGLE EARTH ENGINE CSVs (GAS & PWV)
# ==============================================================================
def load_gee_csv(folder, col_name):
    csv_files = glob.glob(os.path.join(folder, "*.csv"))
    if not csv_files:
        raise FileNotFoundError(f"No CSV file found inside {folder}")
    filepath = csv_files[0]
    df = pd.read_csv(filepath)
    date_col = df.columns[0]
    val_col = df.columns[1]
    df["datetime"] = pd.to_datetime(df[date_col])
    df[col_name] = pd.to_numeric(df[val_col].astype(str).str.replace(",", ""), errors="coerce")
    df = df[["datetime", col_name]].dropna().set_index("datetime").sort_index()
    return df.resample("1D").mean().interpolate(method="time")

# ==============================================================================
# 4. COMPUTE 15-DAY SLIDING MEDIAN & IQR BOUNDS (M +/- 1.5 * IQR)
# ==============================================================================
def compute_sliding_iqr(series, window=15, k=1.5):
    med = series.rolling(window=window, min_periods=5).median().shift(1)
    q75 = series.rolling(window=window, min_periods=5).quantile(0.75).shift(1)
    q25 = series.rolling(window=window, min_periods=5).quantile(0.25).shift(1)
    iqr = q75 - q25
    ub = med + k * iqr
    lb = med - k * iqr
    return med, ub, lb

# ==============================================================================
# MAIN EXECUTION & 4-PANEL IEEE FIGURE GENERATION
# ==============================================================================
if __name__ == "__main__":
    print("Loading datasets...")
    df_omni = load_omni(DIR_OMNI)
    df_vtec = load_all_ionex(DIR_IONEX, EPI_LAT, EPI_LON)
    df_pwv  = load_gee_csv(DIR_TROP, "PWV")
    df_gas  = load_gee_csv(DIR_GAS, "CO")

    # Combine into unified daily timeline
    df = pd.concat([df_omni, df_gas, df_pwv, df_vtec], axis=1)
    df = df.loc["2024-12-08":"2025-01-12"]

    # Compute IQR Bounds
    for param in ["CO", "PWV", "VTEC"]:
        df[f"{param}_med"], df[f"{param}_ub"], df[f"{param}_lb"] = compute_sliding_iqr(df[param], window=15, k=1.5)

    # Plot 4-Panel Figure
    fig, axes = plt.subplots(4, 1, figsize=(11, 10), sharex=True)
    plt.subplots_adjust(hspace=0.18)

    # Panel A: Space Weather (Kp & Dst)
    ax0 = axes[0]
    ax0_twin = ax0.twinx()
    ax0.bar(df.index, df["Kp"], color="steelblue", alpha=0.6, width=0.6, label=r"Max Daily $K_p$")
    ax0.axhline(3.0, color="navy", linestyle="--", linewidth=0.9)
    ax0_twin.plot(df.index, df["Dst"], color="crimson", lw=1.5, marker="o", ms=3, label=r"Min Daily $D_{st}$ (nT)")
    ax0_twin.axhline(-30, color="crimson", linestyle=":", linewidth=0.9)
    ax0.set_ylabel(r"$K_p$ Index", fontweight="bold", color="steelblue")
    ax0_twin.set_ylabel(r"$D_{st}$ (nT)", fontweight="bold", color="crimson")
    ax0.set_title(r"(a) Geomagnetic Conditions ($K_p$ and $D_{st}$ Indices)", loc="left", fontweight="bold", fontsize=10)
    ax0.grid(True, alpha=0.3)

    # Helper function for Panels B, C, D
    def plot_iqr_panel(ax, col, ylabel, title, color):
        ax.plot(df.index, df[f"{col}_med"], color="gray", linestyle="--", lw=1.2, label="15-Day Sliding Median")
        ax.fill_between(df.index, df[f"{col}_lb"], df[f"{col}_ub"], color="lightgray", alpha=0.5, label=r"$M \pm 1.5\times IQR$")
        ax.plot(df.index, df[col], color=color, lw=2.0, marker="o", ms=4, label=f"Observed {col}")
        
        # Highlight Positive Anomalies
        anom = df[df[col] > df[f"{col}_ub"]]
        ax.scatter(anom.index, anom[col], color="red", s=55, zorder=5, label="Precursor Anomaly (> UB)")
        
        ax.axvline(EQ_DATE, color="black", linestyle="-.", lw=1.8, label=r"Mainshock ($M_w$ 7.1)")
        ax.set_ylabel(ylabel, fontweight="bold")
        ax.set_title(title, loc="left", fontweight="bold", fontsize=10)
        ax.legend(loc="upper left", fontsize=8, ncol=4, framealpha=0.9)
        ax.grid(True, alpha=0.3)

    # Panel B: Lithosphere/Atmosphere CO Gas
    plot_iqr_panel(axes[1], "CO", r"CO ($\mathrm{mol/m^2}$)",
                   "(b) Lithospheric/Atmospheric Sentinel-5P CO Column Density", "darkorange")

    # Panel C: Tropospheric PWV
    plot_iqr_panel(axes[2], "PWV", "PWV (mm)",
                   "(c) Tropospheric Precipitable Water Vapor (PWV)", "teal")

    # Panel D: Ionospheric VTEC
    plot_iqr_panel(axes[3], "VTEC", "Daily Max VTEC (TECU)",
                   "(d) Ionospheric Vertical Total Electron Content (VTEC) over Epicenter", "purple")

    axes[3].xaxis.set_major_formatter(mdates.DateFormatter("%b %d\n%Y"))
    axes[3].xaxis.set_major_locator(mdates.DayLocator(interval=3))
    plt.xlabel("Date (UTC)", fontweight="bold")

    out_file = "Figure2_LAIC_Precursors.png"
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    print(f"\nSuccess! Saved publication figure to: {out_file}")
    
    # Print numerical summary of detected anomalies for your paper text
    print("\n--- DETECTED PRECURSOR ANOMALY DATES (> M + 1.5*IQR) ---")
    for param in ["CO", "PWV", "VTEC"]:
        anom_dates = df[df[param] > df[f"{param}_ub"]].index.strftime("%Y-%m-%d").tolist()
        print(f"{param:4s} Anomalies: {anom_dates}")

    plt.show()