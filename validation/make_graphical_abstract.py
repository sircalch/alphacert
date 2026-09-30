"""
Graphical abstract for the AlphaCert v2 manuscript (Elsevier: 531 x 1328 px minimum, readable at
5 x 13 cm). Built only from validation/results/.
"""
import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
INK, INK2 = "#0b0b0b", "#52514e"
NEW, OLD = "#2a78d6", "#e87ba4"


def main():
    pdb = pd.read_csv(os.path.join(RES, "pdb_bench.csv"))
    pdb = pdb[pdb.ok]
    afdb = pd.read_csv(os.path.join(RES, "afdb_bench.csv"))
    afdb = afdb[afdb.ok]
    matplotlib.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                                "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                                "axes.edgecolor": INK2, "xtick.color": INK2, "ytick.color": INK2})
    fig = plt.figure(figsize=(13 / 2.54, 5 / 2.54))
    a = fig.add_axes([0.09, 0.2, 0.33, 0.55])
    b = fig.add_axes([0.55, 0.2, 0.33, 0.55])

    a.scatter(pdb.mp_rama_out_pct, pdb.v100_rama_out_pct, s=5, marker="^", color=OLD, lw=0, label="v1.0.0")
    a.scatter(pdb.mp_rama_out_pct, pdb.ac_rama_out_pct, s=5, marker="o", color=NEW, lw=0, label="v1.1")
    a.plot([0, 10], [0, 10], color=INK2, lw=0.6, ls="--")
    a.set(xscale="symlog", yscale="symlog", xlim=(-0.05, 10), ylim=(-0.05, 100))
    a.set_xticks([0, 1, 10], ["0", "1", "10"], fontsize=7)
    a.set_yticks([0, 1, 10, 100], ["0", "1", "10", "100"], fontsize=7)
    a.minorticks_off()
    a.set_xlabel("MolProbity outliers (%)", fontsize=7)
    a.set_ylabel("AlphaCert outliers (%)", fontsize=7)
    a.set_title("Ramachandran, 150 X-ray entries", fontsize=8, color=INK)
    a.legend(fontsize=6, frameon=False, loc="lower right", borderaxespad=0.1, handletextpad=0.2, markerscale=1.8)

    t, c = afdb.ted_domains.clip(upper=4), afdb["ac_domains_r0.5"]
    groups = [afdb[t == k] for k in range(1, 5)]
    fr = np.array([[(g["ac_domains_r0.5"] < g.ted_domains).mean(), (g["ac_domains_r0.5"] == g.ted_domains).mean(),
                    (g["ac_domains_r0.5"] > g.ted_domains).mean()] for g in groups])
    x = np.arange(4)
    bottom = np.zeros(4)
    for j, (lab, col) in enumerate((("fewer", "#8a8984"), ("equal", NEW), ("more", "#eda100"))):
        b.bar(x, fr[:, j], 0.7, bottom=bottom, color=col, label=lab, edgecolor="white", lw=0.8)
        bottom += fr[:, j]
    b.set_xticks(x, ["1", "2", "3", "4+"], fontsize=7)
    b.set_yticks([0, 0.5, 1], ["0", "50%", "100%"], fontsize=7)
    b.set_ylabel("AFDB models", fontsize=7)
    b.set_xlabel("TED structural domains", fontsize=7)
    b.set_title("PAE rigid bodies vs TED, 200 models", fontsize=8, color=INK)
    h, l = b.get_legend_handles_labels()
    b.legend(h[::-1], l[::-1], title="rigid bodies", title_fontsize=6, fontsize=6, frameon=False, loc="center left",
             bbox_to_anchor=(1.0, 0.5), borderaxespad=0.2)
    b.set_xlim(-0.5, 3.5)
    fig.text(0.5, 0.975, "AlphaCert: checking the checks on predicted protein structures",
             ha="center", va="top", fontsize=9, fontweight="bold", color=INK)
    out = os.path.join(HERE, "figures", "graphical_abstract")
    fig.savefig(out + ".png", dpi=400)
    fig.savefig(out + ".pdf")
    from PIL import Image
    w, h = Image.open(out + ".png").size
    print("graphical abstract", w, "x", h, "px")


if __name__ == "__main__":
    main()
