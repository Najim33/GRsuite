"""Rebuild the validation figures shipped in ``docs/assets``.

Reads ``docs/data/validation_stats.json`` and writes static SVGs. The
palette is colour-vision-safe and the background is explicitly light, so
the figures stay readable on GitHub in both light and dark themes.

Usage::

    python docs/make_figures.py
"""

import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data", "validation_stats.json")
OUT = os.path.join(HERE, "assets")

INK = "#12202A"
INK2 = "#45585F"
MUTED = "#7C8C94"
RULE = "#DDE4E7"
PAPER = "#FFFFFF"
ACCENT = "#0A7EA4"
ALERT = "#C2521F"
MOSS = "#0F8A5F"
CATS = ["#0A7EA4", "#C2521F", "#0F8A5F", "#8B4BB0", "#9A7A10"]
SEQ = ["#D3E7EF", "#9CC9DA", "#5BA6C2", "#1E82A6", "#08536D"]

#: The stats file is produced in French; the published figures are English.
RAIL_LABELS = {
    "Exigence de validation": "Validation requirement",
    "Simulation directe (502 grandeurs)": "Direct simulation (502 quantities)",
    "Débits simulés, validation": "Simulated discharge, validation",
    "Critères NSE / KGE": "NSE and KGE criteria",
    "Paramètres calés": "Calibrated parameters",
    "Epsilon machine (double précision)": "Machine epsilon (float64)",
}
PERF_LABELS = {
    "CemaNeige-GR4J / NSE": "CemaNeige-GR4J / NSE",
}

plt.rcParams.update({
    "figure.facecolor": PAPER,
    "axes.facecolor": PAPER,
    "savefig.facecolor": PAPER,
    "text.color": INK,
    "axes.labelcolor": INK2,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "axes.edgecolor": RULE,
    "font.family": "DejaVu Sans",
    "font.size": 9.5,
    "axes.titlesize": 11,
    "figure.dpi": 110,
    # Sans sel fixe, matplotlib tire des identifiants aleatoires dans le SVG :
    # chaque regeneration produirait un diff, meme sans changement de contenu.
    "svg.hashsalt": "grsuite",
})


def tidy(ax, grid="y"):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.grid(axis=grid, color=RULE, lw=0.8, zorder=0)
    ax.set_axisbelow(True)


def save(fig, name):
    """Write the SVG shipped in the docs, plus a PNG for quick previews.

    ``Date: None`` keeps the SVG byte-identical from one run to the next, so
    regenerating the figures only shows up in ``git diff`` when the numbers
    behind them have actually changed.
    """
    for ext in ("svg", "png"):
        path = os.path.join(OUT, "%s.%s" % (name, ext))
        metadata = {"Date": None} if ext == "svg" else None
        fig.savefig(path, format=ext, bbox_inches="tight", transparent=False,
                    dpi=150, metadata=metadata)
    plt.close(fig)
    print("  ", os.path.relpath(os.path.join(OUT, name + ".svg"), HERE))


# ---------------------------------------------------------------------------


def fig_deviation_scale(d):
    """Where the measured deviations sit, against the 5 % requirement."""
    items = d["rail"]
    fig, ax = plt.subplots(figsize=(9.4, 3.1))

    labels, values, kinds = [], [], []
    for it in items:
        labels.append(RAIL_LABELS.get(it["label"], it["label"]))
        values.append(it["value"])
        kinds.append(it["kind"])
    order = np.argsort(values)
    labels = [labels[i] for i in order]
    values = [values[i] for i in order]
    kinds = [kinds[i] for i in order]
    y = np.arange(len(values))

    requirement = 0.05
    ax.axvspan(1e-17, requirement, color=MOSS, alpha=0.05, zorder=0)
    ax.axvline(requirement, color=ALERT, lw=2, zorder=3)

    for i, (v, k) in enumerate(zip(values, kinds)):
        colour = ALERT if k == "req" else (MUTED if k == "eps" else ACCENT)
        ax.plot([v, requirement], [i, i], color=colour, lw=0.9, ls=":",
                alpha=0.5, zorder=2)
        ax.plot(v, i, "o", ms=7, color=colour, mec=PAPER, mew=1.5, zorder=4)

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9, color=INK)
    ax.set_xscale("log")
    ax.set_xlim(6e-17, 0.35)
    ax.set_ylim(-0.7, len(values) - 0.3)
    ax.set_xlabel("relative deviation, GRsuite vs airGR")
    ax.set_title("Measured agreement against the 5 % requirement", loc="left",
                 color=INK)
    tidy(ax, grid="x")

    for i, v in enumerate(values):
        exponent = int(np.floor(np.log10(v)))
        mant = v / 10 ** exponent
        ax.annotate(r"$%.1f\times10^{%i}$" % (mant, exponent),
                    (v, i), textcoords="offset points", xytext=(9, -11),
                    fontsize=8, color=INK2)
    ax.annotate("5 % requirement", (requirement, -0.5),
                textcoords="offset points", xytext=(-10, 0), ha="right",
                va="center", fontsize=8.5, color=ALERT, weight="bold")
    save(fig, "deviation_scale")


def fig_deviation_histogram(d):
    """Distribution of the deviations over the 500 calibrations."""
    edges = d["hist_param"]["edges"]
    series = [("Calibrated parameters", d["hist_param"]["counts"], CATS[0]),
              ("Discharge, calibration", d["hist_qsim_cal"]["counts"], CATS[2]),
              ("Discharge, validation", d["hist_qsim_val"]["counts"], CATS[4])]

    fig, ax = plt.subplots(figsize=(9.4, 3.4))
    x = np.arange(len(edges) - 1)
    width = 0.26
    for i, (label, counts, colour) in enumerate(series):
        ax.bar(x + (i - 1) * width, counts, width * 0.9, label=label,
               color=colour, zorder=3)

    ax.set_xticks(x)
    ax.set_xticklabels([r"$10^{%i}$" % e for e in edges[:-1]])
    ax.set_xlabel("relative deviation, GRsuite vs airGR")
    ax.set_ylabel("calibrations")
    ax.set_title("500 calibrations: how far apart the two implementations land",
                 loc="left", color=INK)
    ax.legend(frameon=False, loc="upper right", fontsize=9)
    tidy(ax)
    save(fig, "deviation_histogram")


def fig_performance(d):
    """NSE in validation, per model configuration."""
    perf = d["perf"]
    fig, ax = plt.subplots(figsize=(9.4, 3.6))

    for i, p in enumerate(perf):
        colour = CATS[i % len(CATS)]
        pts = np.array(p["nse_val_points"])
        jitter = (np.arange(pts.size) % 11 - 5) / 42.0
        ax.plot(pts, i + jitter, "o", ms=3.4, color=colour, alpha=0.32,
                mec="none", zorder=2)
        b = p["nse_val"]
        ax.plot([b["min"], b["max"]], [i, i], color=colour, lw=1.4, alpha=0.6,
                zorder=3)
        ax.add_patch(plt.Rectangle((b["q1"], i - 0.19), b["q3"] - b["q1"], 0.38,
                                   fill=True, facecolor=PAPER, alpha=0.9,
                                   edgecolor=colour, lw=1.7, zorder=4))
        ax.plot([b["med"], b["med"]], [i - 0.24, i + 0.24], color=colour,
                lw=2.6, zorder=5, solid_capstyle="round")
        ax.annotate("%.3f" % b["med"], (1.005, i), xycoords=("axes fraction",
                    "data"), va="center", fontsize=8.5, color=INK2)

    ax.set_yticks(range(len(perf)))
    ax.set_yticklabels([p["label"] for p in perf], color=INK)
    ax.invert_yaxis()
    ax.set_xlim(0.28, 1.0)
    ax.set_xlabel("NSE over the independent validation period (2010–2019)")
    ax.set_title("Model performance on 100 CAMELS-FR catchments", loc="left",
                 color=INK)
    tidy(ax, grid="x")
    save(fig, "performance")


def fig_basins(d):
    """Hydro-climatic coverage of the reference catchments."""
    basins = d["basins"]
    z = np.array([b["z"] for b in basins], dtype=float)
    ar = np.array([b["ar"] for b in basins], dtype=float)
    sn = np.array([b["sn"] for b in basins], dtype=float)

    bounds = [0.01, 0.025, 0.05, 0.10]
    labels = ["< 1 %", "1–2.5 %", "2.5–5 %", "5–10 %", "> 10 %"]
    idx = np.digitize(sn, bounds)

    fig, ax = plt.subplots(figsize=(9.4, 3.9))
    for k in range(5):
        m = idx == k
        if not m.any():
            continue
        ax.scatter(ar[m], z[m], s=42, color=SEQ[k], edgecolor="#5A6B73",
                   linewidth=0.5, label=labels[k], zorder=3)
    ax.set_xlabel("aridity index  (potential evapotranspiration / precipitation)")
    ax.set_ylabel("median catchment elevation [m]")
    ax.set_title("100 reference catchments, stratified across French hydrology",
                 loc="left", color=INK)
    legend = ax.legend(frameon=False, title="solid precipitation",
                       fontsize=8.5, title_fontsize=8.5, loc="upper right",
                       ncol=2)
    legend.get_title().set_color(MUTED)
    tidy(ax, grid="both")
    save(fig, "reference_basins")


def main():
    os.makedirs(OUT, exist_ok=True)
    with open(DATA, encoding="utf-8") as fh:
        d = json.load(fh)
    print("Writing figures to docs/assets:")
    fig_deviation_scale(d)
    fig_deviation_histogram(d)
    fig_performance(d)
    fig_basins(d)


if __name__ == "__main__":
    main()
