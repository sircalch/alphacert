"""Adds the AlphaCert 1.0.0 clash score and Ramachandran outlier fraction to the PDB benchmark table
(same protein-only inputs), for comparison with 1.1 and MolProbity.

Usage: python legacy_pdb.py BENCH_CSV CACHE_DIR OUT_CSV
"""
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)
from alphacert.parsers.structure_io import parse_protein_structure   # noqa: E402
from legacy_v100 import geometry_v100                                # noqa: E402


def main(bench, cache, out):
    d = pd.read_csv(bench)
    cl, ro = [], []
    for pid in d.pdb:
        s = parse_protein_structure(os.path.join(cache, pid, "protein.pdb"))
        g = geometry_v100.evaluate_stereochemistry(s["coordinates"], s["atom_elements"], s["backbone_atoms"])
        cl.append(g.clashscore)
        ro.append(100 * g.frac_rama_outliers)
    d["v100_clashscore"] = cl
    d["v100_rama_out_pct"] = ro
    d.to_csv(out, index=False)
    print(d[["mp_clashscore", "ac_clashscore", "v100_clashscore", "mp_rama_out_pct", "ac_rama_out_pct", "v100_rama_out_pct"]].median())


if __name__ == "__main__":
    main(*sys.argv[1:4])
