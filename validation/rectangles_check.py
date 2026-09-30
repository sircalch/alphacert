"""
How much of the 1.0.0 Ramachandran failure was the sign error and how much the rectangular regions:
applies the 1.0.0 rectangles (legacy_v100/geometry_v100.py) to correctly signed phi/psi (gemmi) of
the 150 X-ray structures of benchmark_pdb.py, and compares the outlier percentage with MolProbity.

    python validation/rectangles_check.py CACHE_DIR

Adds v100rect_correct_sign_out_pct to results/pdb_bench.csv and prints the summary.
"""
import math
import os
import sys

import gemmi
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from legacy_v100.geometry_v100 import _classify_ramachandran  # noqa: E402


def outlier_pct(pdb_path):
    st = gemmi.read_structure(pdb_path)
    st.setup_entities()
    n = out = 0
    for ch in st[0]:
        pol = ch.get_polymer()
        for i, res in enumerate(pol):
            prev = pol[i - 1] if i > 0 else None
            nxt = pol[i + 1] if i + 1 < len(pol) else None
            if prev is None or nxt is None:
                continue
            phi, psi = gemmi.calculate_phi_psi(prev, res, nxt)
            if math.isnan(phi) or math.isnan(psi):
                continue
            n += 1
            out += _classify_ramachandran(math.degrees(phi), math.degrees(psi)) == "OUTLIER"
    return 100.0 * out / n if n else np.nan


def main(cache):
    p = os.path.join(HERE, "results", "pdb_bench.csv")
    d = pd.read_csv(p)
    d["v100rect_correct_sign_out_pct"] = [outlier_pct(os.path.join(cache, pid, "protein.pdb")) if ok else np.nan
                                          for pid, ok in zip(d.pdb, d.ok)]
    d.to_csv(p, index=False)
    ok = d[d.ok]
    r = ok.v100rect_correct_sign_out_pct
    print(f"n = {len(ok)}; median outliers: MolProbity {ok.mp_rama_out_pct.median():.2f}%, "
          f"1.0.0 rectangles with correct sign {r.median():.2f}%, 1.0.0 as released {ok.v100_rama_out_pct.median():.1f}%, "
          f"Top8000 (1.1) {ok.ac_rama_out_pct.median():.2f}%")
    print(f"mean: MolProbity {ok.mp_rama_out_pct.mean():.2f}%, rectangles {r.mean():.2f}%; "
          f"ratio of means {r.mean() / ok.mp_rama_out_pct.mean():.1f}; Spearman {spearmanr(r, ok.mp_rama_out_pct)[0]:.2f}")


if __name__ == "__main__":
    main(sys.argv[1])
