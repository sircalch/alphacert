"""
AlphaCert stereochemistry proxies against the wwPDB validation (MolProbity) values.

For a random sample of X-ray entries (one protein entity, 100-600 residues, resolution 1.0-3.5 A):
AlphaCert's heavy-atom clash score and approximate Ramachandran outlier fraction, computed on the
standard amino-acid residues of the first model (first alternate location), are compared with the
MolProbity clashscore and Ramachandran outlier percentage reported by the wwPDB validation pipeline
(PDBe API, validation/global-percentiles). MolProbity adds hydrogens and uses all-atom contacts and
Top8000 contours, so agreement is expected in ranking (Spearman), not in value.

Usage: python benchmark_pdb.py RCSB_IDS_JSON CACHE_DIR OUT_CSV [--n 150] [--seed 7] [--workers 8]
"""
import argparse
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from alphacert.parsers.structure_io import parse_protein_structure     # noqa: E402
from alphacert.core.geometry import evaluate_stereochemistry          # noqa: E402

AA = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
      "PRO", "SER", "THR", "TRP", "TYR", "VAL"}


def get(url, path=None, timeout=120):
    if path and os.path.exists(path) and os.path.getsize(path) > 0:
        return path
    r = subprocess.run(["curl", "-sfL", "--max-time", str(timeout), url] + (["-o", path] if path else []),
                       capture_output=True, text=path is None)
    if r.returncode != 0:
        return None
    return path if path else r.stdout


def protein_only_pdb(cif, out):
    """Writes the standard amino-acid residues of the first model as a PDB file (AlphaCert input)."""
    import gemmi
    st = gemmi.read_structure(cif)
    st.setup_entities()
    st.remove_hydrogens()
    st.remove_alternative_conformations()
    for m in list(st)[1:]:
        st.remove_model(m.num)
    for ch in st[0]:
        for i in reversed(range(len(ch))):
            if ch[i].name not in AA:
                del ch[i]
    st.write_pdb(out)
    return out


def analyse(pid, cache):
    d = os.path.join(cache, pid)
    os.makedirs(d, exist_ok=True)
    cif = get(f"https://files.rcsb.org/download/{pid}.cif", os.path.join(d, f"{pid}.cif"))
    if not cif:
        return {"pdb": pid, "ok": False, "reason": "download failed"}
    vp = os.path.join(d, "validation.json")
    if not os.path.exists(vp):
        txt = get(f"https://www.ebi.ac.uk/pdbe/api/validation/global-percentiles/entry/{pid.lower()}")
        open(vp, "w").write(txt or "{}")
    v = json.load(open(vp)).get(pid.lower(), {})
    if "clashscore" not in v:
        return {"pdb": pid, "ok": False, "reason": "no validation data"}
    pdb = protein_only_pdb(cif, os.path.join(d, "protein.pdb"))
    s = parse_protein_structure(pdb)
    g = evaluate_stereochemistry(s["coordinates"], s["atom_elements"], s["backbone_atoms"])
    import gemmi
    res = gemmi.read_structure(cif).resolution
    return {"pdb": pid, "ok": True, "resolution": res, "n_res": s["n_residues"], "n_atoms": len(s["coordinates"]),
            "ac_clashscore": g.clashscore, "ac_rama_out_pct": 100 * g.frac_rama_outliers,
            "ac_rama_fav_pct": 100 * g.frac_rama_favored,
            "mp_clashscore": v["clashscore"]["rawvalue"],
            "mp_rama_out_pct": v.get("percent-rama-outliers", {}).get("rawvalue")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ids_json")
    ap.add_argument("cache")
    ap.add_argument("out_csv")
    ap.add_argument("--n", type=int, default=150)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    ids = json.load(open(a.ids_json))["result_set"]
    rng = np.random.default_rng(a.seed)
    sample = sorted(rng.choice(ids, size=min(a.n, len(ids)), replace=False))

    def safe(pid):
        try:
            return analyse(pid, a.cache)
        except Exception as e:
            return {"pdb": pid, "ok": False, "reason": f"{type(e).__name__}: {e}"[:150]}

    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        rows = list(ex.map(safe, sample))
    d = pd.DataFrame(rows)
    d.to_csv(a.out_csv, index=False)
    ok = d[d.ok]
    from scipy.stats import spearmanr
    print(f"{len(ok)}/{len(d)} analysed")
    print("clash: Spearman", spearmanr(ok.ac_clashscore, ok.mp_clashscore))
    print("rama outliers: Spearman", spearmanr(ok.ac_rama_out_pct, ok.mp_rama_out_pct))


if __name__ == "__main__":
    main()
