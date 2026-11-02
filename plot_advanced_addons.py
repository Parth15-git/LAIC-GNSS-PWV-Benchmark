import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime

from plot_laic_precursors import (
    load_omni, parse_ionex_file, ensure_ionex_files, load_gee_csv,
    compute_sliding_iqr, EQ_DATE, EPI_LAT, EPI_LON, DIR_OMNI, DIR_IONEX, DIR_TROP, DIR_GAS
)

# ==============================================================================
# 1. CONTINUOUS MORLET WAVELET TRANSFORM (PURE NUMPY IMPLEMENTATION)
# ==============================================================================
def morlet_cwt(signal, dt_hours=1.0, periods_hours=np.linspace(2, 36, 60), w0=6.0):
    sig = np.nan_to_num(signal - np.nanmean(signal))
    n = len(sig)
    scales = (periods_hours * w0) / (2 * np.pi * dt_hours)
    power = np.zeros((len(scales), n))

    t = np.arange(-n // 2, n // 2)
    for i, s in enumerate(scales):
        eta = t / s
        wavelet = (np.pi ** -0.25) * np.exp(1j * w0 * eta) * np.exp(-0.5 * (eta ** 2))
        wavelet /= np.sqrt(s)
        conv = np.convolve(sig, np.conj(wavelet[::-1]), mode="same")
        power[i, :] = np.abs(conv) ** 2
    return periods_hours, power

# ==============================================================================
# 2. LOAD FULL HOURLY VTEC (DEC 8 TO JAN 12) SO DEC 12 HAS NO FLAT WARMUP
# ==============================================================================
def load_hourly_vtec(folder, lat, lon):
    ensure_ionex_files(folder)
    inx_files = sorted(list({
        os.path.abspath(f) for f in glob.glob(os.path.join(folder, "*"))
        if ".inx" in f.lower()
    }))
    all_recs = []
    for f in inx_files:
        all_recs.extend(parse_ionex_file(f, lat, lon))
    df = pd.DataFrame(all_recs).drop_duplicates("datetime").set_index("datetime").sort_index()
    return df.resample("1h").interpolate(method="time").loc["2024-12-08":"2025-01-12"]

if __name__ == "__main__":
    print("1. Computing Hourly VTEC Continuous Wavelet Transform (CWT) Scalogram...")
    df_hr_full = load_hourly_vtec(DIR_IONEX, EPI_LAT, EPI_LON)

    # Compute 15-day sliding hourly median using min_periods=1 (eliminates flat line)
    df_hr_full["hour"] = df_hr_full.index.hour
    df_hr_full["bg_med"] = df_hr_full.groupby("hour")["VTEC"].transform(
        lambda x: x.rolling(15, min_periods=1).median().shift(1).bfill()
    )
    df_hr_full["dVTEC"] = (df_hr_full["VTEC"] - df_hr_full["bg_med"]).fillna(0)
    df_hr = df_hr_full.loc["2024-12-12":"2025-01-12"]

    periods, cwt_power = morlet_cwt(df_hr["dVTEC"].values, dt_hours=1.0)
    cwt_norm = (cwt_power - np.mean(cwt_power)) / (np.std(cwt_power) + 1e-9)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 6.6), sharex=True,
                                   gridspec_kw={"height_ratios": [1.05, 1.35]})
    plt.subplots_adjust(hspace=0.18)

    # Panel (a): Hourly Differential VTEC Residual (dVTEC)
    ax1.plot(df_hr.index, df_hr["dVTEC"], color="indigo", lw=1.3,
             label=r"Hourly Differential Residual $\Delta\mathrm{VTEC}$ (TECU)")
    ax1.axhline(0, color="gray", linestyle="--", lw=0.8)
    ax1.axvline(EQ_DATE, color="red", linestyle="-.", lw=1.8, label=r"Mainshock ($M_w$ 7.1)")
    ax1.set_ylim(-25, 36)
    ax1.annotate("Pre-Seismic LAIC Impulse\nJan 4, 2025 (+24.8 TECU)",
                 xy=(mdates.date2num(datetime(2025, 1, 4, 9)), 24.8),
                 xytext=(mdates.date2num(datetime(2024, 12, 25, 12)), 25.0),
                 arrowprops=dict(facecolor="black", shrink=0.06, width=1.2, headwidth=5),
                 fontsize=8.5, fontweight="bold",
                 bbox=dict(boxstyle="round,pad=0.28", fc="lightyellow", ec="indigo", lw=1.2))
    ax1.set_ylabel(r"$\Delta\mathrm{VTEC}$ (TECU)", fontweight="bold")
    ax1.set_title("(a) Hourly Diurnal-Removed Differential Ionospheric Residual over Epicenter (1-Hour Resolution)",
                  loc="left", fontweight="bold", fontsize=10)
    ax1.legend(loc="upper left", fontsize=8.5, framealpha=0.95)
    ax1.grid(True, alpha=0.3)

    # Panel (b): Morlet Wavelet Scalogram
    T, P = np.meshgrid(mdates.date2num(df_hr.index), periods)
    cf = ax2.contourf(T, P, cwt_norm, levels=20, cmap="magma")
    cbar = plt.colorbar(cf, ax=ax2, orientation="horizontal", pad=0.32, aspect=42)
    cbar.set_label("Normalized Wavelet Spectral Power ($Z$-Score)", fontweight="bold")

    ax2.axvline(EQ_DATE, color="cyan", linestyle="-.", lw=1.8, label=r"Mainshock ($M_w$ 7.1)")
    ax2.set_ylabel("Period (Hours)", fontweight="bold")
    ax2.set_xlabel("Date (UTC)", fontweight="bold", labelpad=6)
    ax2.set_title("(b) Continuous Wavelet Transform (CWT) Time-Frequency Scalogram (Pre-Seismic AGW / LAIC Burst on Jan 4, 2025)",
                  loc="left", fontweight="bold", fontsize=10)
    ax2.xaxis.set_major_formatter(mdates.DateFormatter("%b %d\n%Y"))
    ax2.xaxis.set_major_locator(mdates.DayLocator(interval=3))
    ax2.grid(True, linestyle=":", color="white", alpha=0.35)

    out_w = "Figure5_Wavelet_Scalogram.png"
    plt.savefig(out_w, dpi=300, bbox_inches="tight")
    print(f"Saved clean Wavelet Scalogram to: {out_w}")
    plt.close()

    # ==========================================================================
    # 2. MULTI-SPECIES GAS FINGERPRINT & SEQUENTIAL CASCADE (FIGURE 6)
    # ==========================================================================
    print("2. Computing Multi-Species Lithospheric Degassing Cascade...")
    s5p_co_path = os.path.join(DIR_GAS, "s5p_co.csv")
    raw_co = pd.read_csv(s5p_co_path)
    raw_co["datetime"] = pd.to_datetime(raw_co.iloc[:, 0])
    raw_co["CO"] = pd.to_numeric(raw_co.iloc[:, 1].astype(str).str.replace(",", ""), errors="coerce")
    df_co = raw_co[["datetime", "CO"]].dropna().set_index("datetime").sort_index().resample("1D").mean().interpolate(method="time")
    df_co = df_co.loc["2024-12-08":"2025-01-12"]

    multi_path = os.path.join(DIR_GAS, "multi_gas_s5p.csv")
    mg = pd.read_csv(multi_path)
    mg["datetime"] = pd.to_datetime(mg["date"])
    for col in ["CO", "SO2", "NO2", "CH4"]:
        mg[col] = pd.to_numeric(mg[col].astype(str).str.replace(",", ""), errors="coerce")
    mg = mg[["datetime", "CO", "SO2", "NO2", "CH4"]].set_index("datetime").sort_index()
    df_gases = mg.resample("1D").mean().interpolate(method="time").bfill().ffill().reindex(df_co.index)
    # Use the high-resolution chart CO series for exact consistency with Figure 2b
    df_gases["CO"] = df_co["CO"]

    norm_gases = pd.DataFrame(index=df_gases.index)
    for g in ["CO", "CH4", "SO2", "NO2"]:
        med, ub, lb = compute_sliding_iqr(df_gases[g], window=15, k=1.5)
        iqr = (ub - med) / 1.5
        norm_gases[g] = (df_gases[g] - med) / (iqr + 1e-12)

    # Multi-Gas Composite Index (Root-Sum-Square of positive degassing anomalies)
    pos_sq = norm_gases[["CO", "CH4", "SO2", "NO2"]].clip(lower=0) ** 2
    norm_gases["Composite"] = np.sqrt(pos_sq.sum(axis=1))

    df_sub = norm_gases.loc["2024-12-14":"2025-01-12"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 4.9), gridspec_kw={"width_ratios": [1.45, 1.0]})
    plt.subplots_adjust(wspace=0.26)

    styles = {
        "CO":  ("darkorange", "o", "-",  "CO (Shear Micro-Fracturing: Peak Jan 4)"),
        "CH4": ("seagreen",   "s", "--", r"$\mathrm{CH_4}$ (Early Pore-Fluid Seepage: Peak Dec 31)"),
        "SO2": ("crimson",    "^", "-",  r"$\mathrm{SO_2}$ (Hydrothermal Fault Valve: Peak Jan 5)"),
        "NO2": ("royalblue",  "D", ":",  r"$\mathrm{NO_2}$ (Radon/$p$-Hole Ionization: Peak Co-Seismic)"),
    }
    for g, (col, mk, ls, lbl) in styles.items():
        ax1.plot(df_sub.index, df_sub[g], color=col, marker=mk, linestyle=ls, lw=1.8, ms=4, label=lbl)

    ax1.axhline(1.5, color="red", linestyle="--", lw=1.2, label=r"Anomaly Threshold ($\pm 1.5\times IQR$)")
    ax1.axhline(-1.5, color="red", linestyle="--", lw=1.0)
    ax1.axhspan(-1.5, 1.5, color="lightgray", alpha=0.3)
    ax1.axvline(EQ_DATE, color="black", linestyle="-.", lw=1.8, label=r"Mainshock ($M_w$ 7.1)")
    ax1.set_ylabel(r"Normalized Gas Anomaly ($\frac{X - M}{IQR}$)", fontweight="bold")
    ax1.set_xlabel("Date (UTC)", fontweight="bold")
    ax1.set_title("(a) Sequential Multi-Species Trace Gas Cascade (Dec 14, 2024 – Jan 12, 2025)",
                  fontweight="bold", fontsize=10)
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax1.xaxis.set_major_locator(mdates.DayLocator(interval=5))
    ax1.legend(loc="lower left", fontsize=7.5, framealpha=0.92)
    ax1.grid(True, alpha=0.3)

    # Panel (b): Peak Exceedance across the 10-Day Preparation & Rupture Window (Dec 28 - Jan 9)
    win_co   = df_sub.loc["2024-12-28":"2025-01-06", "CO"].max()
    win_ch4  = df_sub.loc["2024-12-28":"2025-01-06", "CH4"].max()
    win_no2  = df_sub.loc["2024-12-28":"2025-01-09", "NO2"].max()
    win_so2  = df_sub.loc["2024-12-28":"2025-01-06", "SO2"].max()
    win_comp = df_sub.loc["2025-01-03":"2025-01-06", "Composite"].max()

    gases_lbl = [
        r"$\mathrm{SO_2}$ (Day -2)" + "\nElectrochem.",
        r"$\mathrm{CH_4}$ (Day -7)" + "\nNDIR Optical",
        r"$\mathrm{NO_2}/{}^{222}\mathrm{Rn}$" + "\nPIN Photodiode",
        r"$\mathrm{CO}$ (Day -3)" + "\nMOS / Electrochem.",
        "4-Gas RSS\nComposite (Day -3)"
    ]
    vals = [win_so2, win_ch4, win_no2, win_co, win_comp]
    colors = ["crimson", "seagreen", "royalblue", "darkorange", "purple"]

    bars = ax2.barh(gases_lbl, vals, color=colors, edgecolor="black", height=0.55, alpha=0.88)
    ax2.axvline(1.5, color="red", linestyle="--", lw=1.5, label=r"Threshold ($1.5\times IQR$)")
    for bar, v in zip(bars, vals):
        ax2.text(v + 0.05, bar.get_y() + bar.get_height()/2, f"+{v:.2f} IQR",
                 va="center", fontweight="bold", fontsize=8.8)

    ax2.set_xlim(0, max(vals) * 1.28)
    ax2.set_xlabel("Peak Standardized Exceedance (IQR Units)", fontweight="bold")
    ax2.set_title("(b) Multi-Stage Gas Cascade -> Sensor Priority", fontweight="bold", fontsize=10)
    ax2.legend(loc="lower right", fontsize=8.5)
    ax2.grid(True, axis="x", alpha=0.3)

    out_g = "Figure6_MultiGas_Correlation.png"
    plt.savefig(out_g, dpi=300, bbox_inches="tight")
    print(f"Saved updated Multi-Gas Cascade plot to: {out_g}")
    print("\n--- UPDATED MULTI-STAGE GAS CASCADE PEAKS ---")
    for name, v in zip(["SO2 (Day -2)", "CH4 (Day -7)", "NO2/Rn", "CO (Day -3)", "4-Gas Composite"], vals):
        print(f"  {name:16s}: +{v:.2f} x IQR")
    plt.show()