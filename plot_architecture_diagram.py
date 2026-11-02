import matplotlib.pyplot as plt
import matplotlib.patches as patches

def draw_box(ax, x, y, w, h, text, bg_color, border_color="black", fontsize=8.5, bold=True):
    rect = patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.015,rounding_size=0.025",
                                  facecolor=bg_color, edgecolor=border_color, linewidth=1.6)
    ax.add_patch(rect)
    weight = "bold" if bold else "normal"
    ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=fontsize,
            fontweight=weight, color="black", multialignment="center")

def draw_arrow(ax, x1, y1, x2, y2, label="", lx=None, ly=None):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="-|>,head_width=0.28,head_length=0.45", lw=1.6, color="#2b2b2b"))
    if label:
        tx = (x1 + x2) / 2 if lx is None else lx
        ty = (y1 + y2) / 2 if ly is None else ly
        ax.text(tx, ty, label, ha="center", va="center",
                fontsize=7.5, fontweight="bold", color="#1a1a1a",
                bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="#bbbbbb", lw=0.8, alpha=0.95))

if __name__ == "__main__":
    fig, ax = plt.subplots(figsize=(13.5, 5.8))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # Section Headers
    ax.text(0.145, 0.94, "1. Multi-Layer LAIC Sensing Nodes", ha="center", fontsize=11, fontweight="bold", color="#1f4e78")
    ax.text(0.500, 0.94, "2. Pre-Processing & IQR Filtering", ha="center", fontsize=11, fontweight="bold", color="#1f4e78")
    ax.text(0.855, 0.94, "3. Coincidence Gate & Early Warning", ha="center", fontsize=11, fontweight="bold", color="#1f4e78")

    # Column 1: Sensing Nodes (x = 0.02 to 0.27)
    draw_box(ax, 0.02, 0.67, 0.25, 0.20,
             "Lithospheric Gas Sensor Array\n• Radon-222 / CO / CH4 / CO2\n• STM32 / ESP32 Edge Telemetry\n• Sentinel-5P TROPOMI Backup",
             "#ffe6cc", "#d79b00")
    draw_box(ax, 0.02, 0.38, 0.25, 0.20,
             "CORS Dual-Frequency GNSS\n(IIT-ISM Dhanbad & Regional IGS)\n• Carrier Phase (L1/L2) & Pseudorange\n• Surface Met (P, T, RH)",
             "#dae8fc", "#6c8ebf")
    draw_box(ax, 0.02, 0.09, 0.25, 0.20,
             "Space Weather & Far-Field Control\n• NASA OMNIWeb (Kp, Dst, F10.7)\n• Far-Field Control CORS (IISC)\n  (Outside Dobrovolsky Radius R)",
             "#f8cecc", "#b85450")

    # Column 2: Signal Extraction & Sliding IQR (x = 0.37 to 0.63)
    draw_box(ax, 0.37, 0.67, 0.26, 0.20,
             "Lithospheric Degassing Module\n• Gas Column / Surface Flux\n• 15-Day Sliding Median & IQR\n• Flag: Gas > M + 1.5×IQR",
             "#fff2cc", "#d6b656")
    draw_box(ax, 0.37, 0.38, 0.26, 0.20,
             "GNSS Troposphere & Ionosphere\n• ZTD → ZWD → PWV (mm)\n• Slant TEC → Mapped VTEC (TECU)\n• Flags: PWV > UB & VTEC > UB",
             "#e1d5e7", "#9673a6")
    draw_box(ax, 0.37, 0.09, 0.26, 0.20,
             "Solar & Equatorial Disturbance Gate\n• Quiet Criteria: Kp < 3, Dst > -30 nT\n• Spatial Differential Filter:\n  ΔVTEC(Near) >> ΔVTEC(Control)",
             "#f5f5f5", "#666666")

    # Column 3: Decision Fusion Engine (x = 0.73 to 0.98)
    draw_box(ax, 0.73, 0.45, 0.25, 0.42,
             "Multi-Modal LAIC Fusion Engine\n\n48-Hr Upward Coupling Coincidence:\nG_total = Flag(Gas) × Flag(PWV)\n          × Flag(VTEC) × Gate(Quiet)\n\nFusion Index = Π Normalized Scores",
             "#d5e8d4", "#82b366")

    draw_box(ax, 0.73, 0.08, 0.25, 0.23,
             "CONFIRMED SEISMIC PRECURSOR\n• Suppresses Weather False Alarms\n• Suppresses Solar Storm Spikes\n• Validated Lead Time: 3 Days\n  (Jan 4, 2025 Alert -> Jan 7 Mw 7.1)",
             "#ffcccc", "#cc0000")

    # Horizontal Arrows (Col 1 -> Col 2)
    draw_arrow(ax, 0.27, 0.77, 0.37, 0.77, "UART / MQTT", ly=0.81)
    draw_arrow(ax, 0.27, 0.48, 0.37, 0.48, "RINEX / NTRIP", ly=0.52)
    draw_arrow(ax, 0.27, 0.19, 0.37, 0.19, "API Stream", ly=0.23)

    # Converging Arrows (Col 2 -> Col 3)
    draw_arrow(ax, 0.63, 0.77, 0.73, 0.73, "S_Gas", lx=0.68, ly=0.78)
    draw_arrow(ax, 0.63, 0.48, 0.73, 0.63, "S_PWV, S_VTEC", lx=0.675, ly=0.53)
    draw_arrow(ax, 0.63, 0.19, 0.73, 0.52, "G_Quiet", lx=0.68, ly=0.33)

    # Vertical Alert Arrow (Fusion Engine -> Confirmed Precursor)
    draw_arrow(ax, 0.855, 0.45, 0.855, 0.31, "Index >= 1.0", lx=0.855, ly=0.38)

    plt.title("Proposed Multi-Modal Lithosphere-Atmosphere-Ionosphere Coupling (LAIC) Gas & GNSS Early-Warning Architecture",
              fontsize=11.5, fontweight="bold", pad=10)
    out_file = "Figure1_System_Architecture.png"
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    print(f"Saved clean architecture diagram to: {out_file}")
    plt.show()