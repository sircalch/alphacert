"""
Real AlphaFold DB v6 models: haemoglobin alpha (P69905, globular) and p53 (P04637, two folded
domains plus disordered regions). References are independent of AlphaCert: the mean pLDDT and
band fractions published by the AlphaFold DB API, and gemmi's phi/psi.
"""
import os
import numpy as np
import pytest
from alphacert.parsers.structure_io import parse_protein_structure
from alphacert.parsers.pae_json import parse_pae_matrix
from alphacert.core.plddt import evaluate_plddt_profile
from alphacert.core.geometry import evaluate_stereochemistry
from alphacert.core.pae import evaluate_pae_matrix

D = os.path.join(os.path.dirname(__file__), "data", "afdb")
# AlphaFold DB API: globalMetricValue and fractionPlddt{VeryHigh,Confident,Low,VeryLow}
AFDB = {
    "P69905": (98.06, 0.993, 0.0, 0.007, 0.0),
    "P04637": (75.06, 0.527, 0.071, 0.104, 0.298),
}


@pytest.mark.parametrize("uid", sorted(AFDB))
def test_plddt_matches_afdb(uid):
    s = parse_protein_structure(os.path.join(D, f"AF-{uid}-F1-model_v6.cif"))
    r = evaluate_plddt_profile(s["per_residue_plddt"])
    mean, vh, conf, low, vlow = AFDB[uid]
    assert r.mean_plddt == pytest.approx(mean, abs=0.01)
    assert [round(x, 3) for x in (r.frac_very_high_90, r.frac_confident_70_90, r.frac_low_50_70,
                                  r.frac_very_low_under_50)] == [vh, conf, low, vlow]


def test_pdb_and_mmcif_agree():
    a = parse_protein_structure(os.path.join(D, "AF-P69905-F1-model_v6.pdb"))
    b = parse_protein_structure(os.path.join(D, "AF-P69905-F1-model_v6.cif"))
    assert a["n_residues"] == b["n_residues"] == 142
    np.testing.assert_allclose(a["coordinates"], b["coordinates"])
    np.testing.assert_allclose(a["per_residue_plddt"], b["per_residue_plddt"])
    assert a["atom_elements"] == b["atom_elements"]


def test_phi_psi_match_gemmi():
    gemmi = pytest.importorskip("gemmi")
    f = os.path.join(D, "AF-P69905-F1-model_v6.pdb")
    s = parse_protein_structure(f)
    g = evaluate_stereochemistry(s["coordinates"], s["atom_elements"], s["backbone_atoms"])
    ch = gemmi.read_structure(f)[0][0]
    ref = np.degrees([gemmi.calculate_phi_psi(ch[i - 1], ch[i], ch[i + 1]) for i in range(1, len(ch) - 1)])
    got = np.array(g.phi_psi_angles)
    assert got.shape == ref.shape
    assert np.max(np.abs((got - ref + 180) % 360 - 180)) < 1e-6


def test_high_confidence_model_is_clash_free_and_helical():
    s = parse_protein_structure(os.path.join(D, "AF-P69905-F1-model_v6.pdb"))
    g = evaluate_stereochemistry(s["coordinates"], s["atom_elements"], s["backbone_atoms"])
    assert g.clashscore < 5.0                    # v1.0.0 reported 770 (peptide bonds counted as clashes)
    phi = np.array([p for p, _ in g.phi_psi_angles])
    assert np.median(phi) < -50                  # right-handed helix; v1.0.0 had the sign inverted
    assert g.status == "PASS"


def test_pae_domains():
    hb = parse_protein_structure(os.path.join(D, "AF-P69905-F1-model_v6.cif"))
    r = evaluate_pae_matrix(parse_pae_matrix(os.path.join(D, "AF-P69905-F1-predicted_aligned_error_v6.json"))["pae_matrix"],
                            plddt=hb["per_residue_plddt"])
    assert r.n_compact_domains == 1
    p53 = parse_protein_structure(os.path.join(D, "AF-P04637-F1-model_v6.cif"))
    r = evaluate_pae_matrix(parse_pae_matrix(os.path.join(D, "AF-P04637-F1-predicted_aligned_error_v6.json"))["pae_matrix"],
                            plddt=p53["per_residue_plddt"])
    spans = [(d[0] + 1, d[-1] + 1) for d in r.domains]
    # DNA-binding domain (~94-292) and tetramerisation domain (~325-356); disordered parts excluded
    assert len(spans) == 2
    assert 90 <= spans[0][0] <= 100 and 288 <= spans[0][1] <= 300
    assert 320 <= spans[1][0] <= 330 and 352 <= spans[1][1] <= 362


def test_multimer_iptm_is_primary():
    pae = np.full((100, 100), 4.0)
    pae[:50, 50:] = pae[50:, :50] = 9.0          # inter-chain PAE would pass on its own
    assert evaluate_pae_matrix(pae, chain_boundaries=[50, 50], iptm_score=0.45).status == "FAIL"
    assert evaluate_pae_matrix(pae, chain_boundaries=[50, 50], iptm_score=0.70).status == "WARNING"
    assert evaluate_pae_matrix(pae, chain_boundaries=[50, 50], iptm_score=0.85).status == "PASS"
    assert evaluate_pae_matrix(pae, chain_boundaries=[50, 50]).status == "PASS"
