"""
AlphaCert benchmark on AlphaFold DB models of human proteins.

For a random sample of reviewed human UniProt entries (length 100-800) with an AlphaFold DB model:
  * pLDDT: mean and band fractions from AlphaCert (mmCIF) vs the values published by the AFDB API
  * phi/psi: AlphaCert vs gemmi, every residue
  * PAE domains: AlphaCert (Croll clustering + pLDDT filter) vs the TED consensus domains
    (Lau et al., Science 2024; ted.cathdb.info), compared by domain count and by the adjusted Rand
    index of the residue partition (residues outside any domain form one class)
TED was computed on AFDB v4 models. The v4 files are used for the domain comparison when AFDB still
serves them; otherwise the v6 model and PAE (same sequence and residue numbering) are used, and the
version is recorded per protein.

Usage: python benchmark_afdb.py UNIPROT_TSV CACHE_DIR OUT_CSV [--n 150] [--seed 7] [--workers 8]
"""
import argparse
import json
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from alphacert.parsers.structure_io import parse_protein_structure     # noqa: E402
from alphacert.parsers.pae_json import parse_pae_matrix               # noqa: E402
from alphacert.core.plddt import evaluate_plddt_profile               # noqa: E402
from alphacert.core.geometry import evaluate_stereochemistry          # noqa: E402
from alphacert.core.pae import pae_domains                            # noqa: E402


RESOLUTIONS = (0.25, 0.5, 1.0, 2.0)


def get(url, path=None, timeout=120):
    if path and os.path.exists(path) and os.path.getsize(path) > 0:
        return path
    r = subprocess.run(["curl", "-sfL", "--max-time", str(timeout), url] + (["-o", path] if path else []),
                       capture_output=True, text=path is None)
    if r.returncode != 0:
        return None
    return path if path else r.stdout


def fetch(acc, cache):
    d = os.path.join(cache, acc)
    os.makedirs(d, exist_ok=True)
    api_p = os.path.join(d, "api.json")
    if not os.path.exists(api_p):
        txt = get(f"https://alphafold.ebi.ac.uk/api/prediction/{acc}")
        if not txt:
            return None
        open(api_p, "w").write(txt)
    api = json.load(open(api_p))
    if not api:
        return None
    api = api[0]
    files = {"cif6": get(api["cifUrl"], os.path.join(d, "model_v6.cif")),
             "pae6": get(api["paeDocUrl"], os.path.join(d, "pae_v6.json")),
             "cif4": get(api["cifUrl"].replace("_v6", "_v4"), os.path.join(d, "model_v4.cif")),
             "pae4": get(api["paeDocUrl"].replace("_v6", "_v4"), os.path.join(d, "pae_v4.json"))}
    ted_p = os.path.join(d, "ted.json")
    if not os.path.exists(ted_p):
        txt = get(f"https://ted.cathdb.info/api/v1/uniprot/summary/{acc}")
        open(ted_p, "w").write(txt or "{}")
    return api, files, json.load(open(ted_p))


def partition(n, domains):
    lab = np.zeros(n, int)
    for k, res in enumerate(domains, start=1):
        lab[[r for r in res if 0 <= r < n]] = k
    return lab


def ted_domains(ted, levels=("high", "medium")):
    out = []
    for d in ted.get("data", []):
        if d.get("consensus_level") not in levels:
            continue
        res = []
        for seg in d["chopping"].split("_"):
            a, b = seg.split("-")
            res += list(range(int(a) - 1, int(b)))
        out.append(res)
    return out


def analyse(acc, cache):
    from sklearn.metrics import adjusted_rand_score
    import gemmi
    got = fetch(acc, cache)
    if got is None:
        return {"acc": acc, "ok": False, "reason": "no AFDB entry"}
    api, files, ted = got
    if not files["cif6"] or not files["pae6"]:
        return {"acc": acc, "ok": False, "reason": "download failed"}
    s = parse_protein_structure(files["cif6"])
    p = evaluate_plddt_profile(s["per_residue_plddt"])
    g = evaluate_stereochemistry(s["coordinates"], s["atom_elements"], s["backbone_atoms"])
    ch = gemmi.read_structure(files["cif6"])[0][0]
    ref = np.degrees([gemmi.calculate_phi_psi(ch[i - 1], ch[i], ch[i + 1]) for i in range(1, len(ch) - 1)])
    got_pp = np.array(g.phi_psi_angles)
    dpp = float(np.max(np.abs((got_pp - ref + 180) % 360 - 180))) if got_pp.shape == ref.shape else np.nan
    row = {"acc": acc, "ok": True, "n_res": s["n_residues"],
           "mean_plddt": p.mean_plddt, "api_mean_plddt": api.get("globalMetricValue"),
           "vh": p.frac_very_high_90, "api_vh": api.get("fractionPlddtVeryHigh"),
           "conf": p.frac_confident_70_90, "api_conf": api.get("fractionPlddtConfident"),
           "low": p.frac_low_50_70, "api_low": api.get("fractionPlddtLow"),
           "vlow": p.frac_very_low_under_50, "api_vlow": api.get("fractionPlddtVeryLow"),
           "max_phipsi_diff_deg": dpp, "clashscore": g.clashscore, "rama_fav": g.frac_rama_favored,
           "rama_out": g.frac_rama_outliers}
    tds = ted_domains(ted)
    row["ted_domains"] = len(tds)
    row["ted_domains_high"] = len(ted_domains(ted, ("high",)))
    use4 = bool(files["cif4"] and files["pae4"])
    row["domain_model_version"] = "v4" if use4 else "v6"
    if True:
        s4 = parse_protein_structure(files["cif4"] if use4 else files["cif6"])
        pae4 = parse_pae_matrix(files["pae4"] if use4 else files["pae6"])["pae_matrix"]
        n = pae4.shape[0]
        for res in RESOLUTIONS:
            doms = pae_domains(pae4, plddt=s4["per_residue_plddt"], resolution=res)
            tag = f"{res:g}"
            row[f"ac_domains_r{tag}"] = len(doms)
            row[f"ari_r{tag}"] = float(adjusted_rand_score(partition(n, tds), partition(n, doms))) if (tds or doms) else 1.0
        row["ac_domains"] = row["ac_domains_r0.5"]
        row["ari_vs_ted"] = row["ari_r0.5"]
    return row


def safe(acc, cache):
    try:
        return analyse(acc, cache)
    except Exception as e:
        return {"acc": acc, "ok": False, "reason": f"{type(e).__name__}: {e}"[:150]}


def main():
    if os.name == "nt":                       # keep the machine awake while the benchmark runs
        import ctypes
        ctypes.windll.kernel32.SetThreadExecutionState(0x80000000 | 0x00000001)
    ap = argparse.ArgumentParser()
    ap.add_argument("uniprot_tsv")
    ap.add_argument("cache")
    ap.add_argument("out_csv")
    ap.add_argument("--n", type=int, default=150)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    accs = pd.read_csv(a.uniprot_tsv, sep="\t").Entry.tolist()
    rng = np.random.default_rng(a.seed)
    sample = sorted(rng.choice(accs, size=min(a.n, len(accs)), replace=False))

    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        rows = list(ex.map(safe, sample, [a.cache] * len(sample)))
    d = pd.DataFrame(rows)
    d.to_csv(a.out_csv, index=False)
    ok = d[d.ok]
    print(f"{len(ok)}/{len(d)} analysed")
    print("max |mean pLDDT - API|:", (ok.mean_plddt - ok.api_mean_plddt).abs().max())
    for k in ("vh", "conf", "low", "vlow"):
        print(f"max |{k} - API|:", (ok[k] - ok["api_" + k]).abs().max())
    print("max phi/psi diff (deg):", ok.max_phipsi_diff_deg.max())
    if "ac_domains" in ok:
        m = ok.dropna(subset=["ac_domains"])
        print("domain count equal to TED:", float((m.ac_domains == m.ted_domains).mean()), "within 1:",
              float(((m.ac_domains - m.ted_domains).abs() <= 1).mean()), "median ARI:", m.ari_vs_ted.median())


if __name__ == "__main__":
    main()
