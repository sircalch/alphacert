"""
Figures and tables for the AlphaCert v2 manuscript, built only from validation/results/.

    python validation/make_figures.py

Style: Elsevier double-column width 190 mm, 8 pt sans-serif text, one fixed colour and marker per
series in every figure (validated categorical palette; identity never relies on colour alone).
"""
import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
FIG = os.path.join(HERE, "figures")
TAB = os.path.join(HERE, "tables")
MM = 1 / 25.4
DOUBLE = 190 * MM
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
NEW, OLD, REF = "#2a78d6", "#e87ba4", "#1baf7a"


def setup():
    matplotlib.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "xtick.labelsize": 7,
        "ytick.labelsize": 7, "legend.fontsize": 6.5, "axes.edgecolor": INK2, "axes.labelcolor": INK,
        "xtick.color": INK2, "ytick.color": INK2, "axes.linewidth": 0.6, "axes.spines.top": False,
        "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5,
        "axes.axisbelow": True, "legend.frameon": False, "lines.linewidth": 1.2, "lines.markersize": 4,
        "savefig.dpi": 600, "pdf.fonttype": 42, "ps.fonttype": 42})


def panel(ax, letter):
    ax.text(-0.14, 1.03, f"({letter})", transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", color=INK)


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


def fig_pdb(d):
    fig, (a, b) = plt.subplots(1, 2, figsize=(DOUBLE, 72 * MM))
    x = d.mp_clashscore + 0.5
    a.scatter(x, d.v100_clashscore + 0.5, s=10, marker="^", color=OLD, lw=0.2, ec="white", label="AlphaCert 1.0.0", zorder=2)
    a.scatter(x, d.ac_clashscore + 0.5, s=10, marker="o", color=NEW, lw=0.2, ec="white", label="AlphaCert 1.1", zorder=3)
    lim = (0.4, 3000)
    a.plot(lim, lim, color=INK2, lw=0.7, ls="--", zorder=1)
    rho = spearmanr(d.ac_clashscore, d.mp_clashscore)[0]
    a.set(xscale="log", yscale="log", xlim=(0.4, 200), ylim=lim, xlabel="MolProbity clashscore + 0.5 (wwPDB)",
          ylabel="AlphaCert clash score + 0.5")
    a.text(0.97, 0.05, f"1.1: Spearman $\\rho$ = {rho:.2f}", transform=a.transAxes, ha="right", fontsize=7, color=NEW)
    a.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2)
    panel(a, "a")

    b.scatter(d.mp_rama_out_pct, d.v100_rama_out_pct, s=10, marker="^", color=OLD, lw=0.2, ec="white", label="AlphaCert 1.0.0 (rectangles, inverted sign)", zorder=2)
    b.scatter(d.mp_rama_out_pct, d.ac_rama_out_pct, s=10, marker="o", color=NEW, lw=0.2, ec="white", label="AlphaCert 1.1 (Top8000 contours)", zorder=3)
    b.plot([0, 10], [0, 10], color=INK2, lw=0.7, ls="--", zorder=1)
    b.set(xscale="symlog", yscale="symlog", xlim=(-0.05, 10), ylim=(-0.05, 100),
          xlabel="MolProbity Ramachandran outliers (%)", ylabel="AlphaCert outliers (%)")
    b.set_xticks([0, 0.5, 1, 2, 5, 10], ["0", "0.5", "1", "2", "5", "10"])
    b.set_yticks([0, 1, 10, 100], ["0", "1", "10", "100"])
    exact = float(((d.ac_rama_out_pct - d.mp_rama_out_pct).abs() < 0.01).mean())
    b.text(0.97, 0.05, f"1.1: identical in {exact:.0%} of entries", transform=b.transAxes, ha="right", fontsize=7, color=NEW)
    b.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=1)
    panel(b, "b")
    fig.tight_layout(w_pad=2.5)
    save(fig, "fig2_pdb")


def table_pdb(d):
    rows = [("MolProbity (wwPDB)", d.mp_clashscore, d.mp_rama_out_pct),
            ("AlphaCert 1.1", d.ac_clashscore, d.ac_rama_out_pct),
            ("AlphaCert 1.0.0", d.v100_clashscore, d.v100_rama_out_pct)]
    lines = [r"\begin{tabular}{lrrrr}", r"\toprule",
             r" & \multicolumn{2}{c}{clash score} & \multicolumn{2}{c}{Ramachandran outliers (\%)} \\",
             r"\cmidrule(lr){2-3}\cmidrule(lr){4-5}",
             r" & median [IQR] & Spearman $\rho$ & median [IQR] & Spearman $\rho$ \\", r"\midrule"]
    for name, c, r in rows:
        rc = spearmanr(c, d.mp_clashscore)[0]
        rr = spearmanr(r, d.mp_rama_out_pct)[0]
        lines.append(f"{name} & {c.median():.1f} [{c.quantile(0.25):.1f}, {c.quantile(0.75):.1f}] & {rc:.2f} & "
                     f"{r.median():.2f} [{r.quantile(0.25):.2f}, {r.quantile(0.75):.2f}] & {rr:.3f} " + r"\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    with open(os.path.join(TAB, "table_pdb.tex"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


def main():
    for d in (FIG, TAB):
        os.makedirs(d, exist_ok=True)
    setup()
    pdb = pd.read_csv(os.path.join(RES, "pdb_bench.csv"))
    pdb = pdb[pdb.ok]
    fig_pdb(pdb)
    table_pdb(pdb)
    print("PDB figures written")


if __name__ == "__main__":
    main()
