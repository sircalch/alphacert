"""
Stereochemical quality analysis: steric clashes, Ramachandran dihedral angles, and peptide bonds.

Notes on scope
--------------
* The clash score is a heavy-atom proxy: pairs of non-bonded heavy atoms (more than three covalent
  bonds apart) that overlap by more than 0.4 A of their van der Waals radii, per 1000 heavy atoms;
  N/O pairs count only below 2.5 A (possible hydrogen bonds).
  It is NOT the MolProbity clashscore, which is computed with explicit hydrogens by Probe.
* Ramachandran classification uses the MolProbity Top8000 percentile contours (Richardson Lab,
  CC-BY 4.0; Williams et al., Protein Sci. 27, 293, 2018) for the six MolProbity residue categories
  (general, Gly, cis-Pro, trans-Pro, pre-Pro, Ile/Val) with the MolProbity thresholds. Versions
  before 1.1.0 used approximate rectangular regions that overestimated outliers several-fold.
"""

from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass
import numpy as np
from scipy.spatial import cKDTree


@dataclass
class StereochemistryResult:
    n_clashes: int
    clashscore: float  # heavy-atom clashes per 1000 heavy atoms (not the MolProbity clashscore)
    n_ramachandran_angles: int
    frac_rama_favored: float   # NaN when no backbone was available
    frac_rama_allowed: float
    frac_rama_outliers: float
    n_cis_peptides: int
    n_non_proline_cis_peptides: int
    status: str  # 'PASS', 'WARNING', 'FAIL'
    diagnostic_message: str
    phi_psi_angles: List[Tuple[float, float]]


# Standard Van der Waals radii (Angstroms), Bondi
VDW_RADII = {
    "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80, "H": 1.20, "SE": 1.90, "X": 1.70
}
# Covalent radii (Angstroms), used to infer bonds
COV_RADII = {"C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05, "P": 1.07, "H": 0.31, "SE": 1.20}


def _calc_dihedral(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray, p4: np.ndarray) -> float:
    """Dihedral angle p1-p2-p3-p4 in degrees, IUPAC sign convention (alpha helix phi ~ -60)."""
    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    nb2 = np.linalg.norm(b2)
    if np.linalg.norm(n1) < 1e-6 or np.linalg.norm(n2) < 1e-6 or nb2 < 1e-6:
        return float("nan")
    x = np.dot(n1, n2)
    y = np.dot(np.cross(n1, n2), b2 / nb2)
    return float(np.degrees(np.arctan2(y, x)))


_RAMA_GRIDS = None
# MolProbity thresholds on the Top8000 percentile contours (Williams et al., Protein Sci. 2018):
# favored >= 2%; outlier < 0.05% for the general case and < 0.1% for the other categories
RAMA_FAVORED = 0.02
RAMA_OUTLIER = {"general": 0.0005, "gly": 0.001, "cispro": 0.001, "transpro": 0.001,
                "prepro": 0.001, "ileval": 0.001}


def _rama_grids():
    global _RAMA_GRIDS
    if _RAMA_GRIDS is None:
        import os
        path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "top8000_rama.npz")
        with np.load(path) as z:
            _RAMA_GRIDS = {k: z[k].astype(float) for k in z.files}
    return _RAMA_GRIDS


def rama_category(res_name: str, next_res_name: str, omega_prev: Optional[float]) -> str:
    """MolProbity Ramachandran category of a residue."""
    if res_name == "GLY":
        return "gly"
    if res_name == "PRO":
        return "cispro" if (omega_prev is not None and not np.isnan(omega_prev) and abs(omega_prev) < 30.0) else "transpro"
    if next_res_name == "PRO":
        return "prepro"
    if res_name in ("ILE", "VAL"):
        return "ileval"
    return "general"


def rama_percentile(phi: float, psi: float, category: str = "general") -> float:
    """Bilinear interpolation of the periodic Top8000 percentile grid (2-degree bins centred at -179, -177, ...)."""
    g = _rama_grids()[category]
    x = (phi + 179.0) / 2.0
    y = (psi + 179.0) / 2.0
    i0, j0 = int(np.floor(x)), int(np.floor(y))
    fx, fy = x - i0, y - j0
    i0, j0 = i0 % 180, j0 % 180
    i1, j1 = (i0 + 1) % 180, (j0 + 1) % 180
    return float((1 - fx) * (1 - fy) * g[i0, j0] + fx * (1 - fy) * g[i1, j0]
                 + (1 - fx) * fy * g[i0, j1] + fx * fy * g[i1, j1])


def _classify_ramachandran(phi: float, psi: float, category: str = "general") -> str:
    """FAVORED / ALLOWED / OUTLIER from the Top8000 contours of the residue's category."""
    p = rama_percentile(phi, psi, category)
    if p >= RAMA_FAVORED:
        return "FAVORED"
    if p < RAMA_OUTLIER[category]:
        return "OUTLIER"
    return "ALLOWED"


def _bond_graph_exclusions(coords: np.ndarray, elements: List[str], tree: cKDTree, max_bonds: int = 3) -> set:
    """Pairs (i, j), i < j, separated by at most max_bonds covalent bonds."""
    cov = np.array([COV_RADII.get(e.upper(), 0.77) for e in elements])
    adj: List[List[int]] = [[] for _ in range(len(coords))]
    for i, j in tree.query_pairs(r=2.2, output_type="ndarray"):
        if np.linalg.norm(coords[i] - coords[j]) < cov[i] + cov[j] + 0.4:
            adj[i].append(int(j))
            adj[j].append(int(i))
    excluded = set()
    for s in range(len(coords)):
        frontier, seen = {s}, {s}
        for _ in range(max_bonds):
            nxt = set()
            for a in frontier:
                for b in adj[a]:
                    if b not in seen:
                        seen.add(b)
                        nxt.add(b)
            frontier = nxt
        for t in seen:
            if t > s:
                excluded.add((s, t))
    return excluded


POLAR = {"N", "O"}
HBOND_MIN_DIST = 2.5   # N/O...N/O pairs closer than this are clashes; farther ones may be H-bonds or salt bridges


def count_heavy_atom_clashes(coordinates: np.ndarray, atom_elements: List[str], overlap: float = 0.40) -> int:
    """
    Non-bonded heavy-atom pairs (> 3 bonds apart) with d < r_i + r_j - overlap.

    Pairs of polar atoms (N/O) count only below 2.5 A, since without hydrogens a hydrogen bond or a
    salt bridge (2.6-3.0 A) cannot be told apart from an overlap (Probe likewise allows H-bond overlap).
    """
    elems = [e.upper() for e in atom_elements]
    keep = [i for i, e in enumerate(elems) if e != "H"]
    coords = np.asarray(coordinates, dtype=float)[keep]
    elems = [elems[i] for i in keep]
    if len(coords) < 2:
        return 0
    tree = cKDTree(coords)
    excluded = _bond_graph_exclusions(coords, elems, tree)
    radii = np.array([VDW_RADII.get(e, 1.70) for e in elems])
    n = 0
    for i, j in tree.query_pairs(r=2 * max(radii) - overlap, output_type="ndarray"):
        if (min(i, j), max(i, j)) in excluded:
            continue
        d = np.linalg.norm(coords[i] - coords[j])
        if elems[i] in POLAR and elems[j] in POLAR:
            if d < HBOND_MIN_DIST:
                n += 1
        elif d < radii[i] + radii[j] - overlap:
            n += 1
    return n


def evaluate_stereochemistry(
    coordinates: np.ndarray,
    atom_elements: List[str],
    backbone_atoms: Optional[List[Dict[str, Any]]] = None,
    max_clashscore_pass: float = 10.0,
    min_rama_favored_pass: float = 0.88,
    max_rama_outliers_pass: float = 0.03
) -> StereochemistryResult:
    """
    Evaluates steric clashes, Ramachandran phi/psi distributions, and cis-peptide stereochemistry.

    Parameters
    ----------
    coordinates : np.ndarray
        (N_atoms, 3) Cartesian coordinates.
    atom_elements : list of str
        Element symbols for each atom.
    backbone_atoms : list of dict, optional
        Residue list with backbone atom coordinates {'N', 'CA', 'C', 'O', 'res_name', 'chain'}.
        Without it, the Ramachandran analysis is not performed (fractions are NaN).

    Returns
    -------
    result : StereochemistryResult
        Stereochemical validation metrics and pass/fail score.
    """
    coords = np.asarray(coordinates, dtype=float)
    n_heavy = sum(1 for e in atom_elements if e.upper() != "H")

    # 1. Heavy-atom clashes (all atoms, bonded pairs up to 1-4 excluded via the covalent graph)
    clash_count = count_heavy_atom_clashes(coords, atom_elements) if len(coords) > 1 else 0
    clashscore = clash_count / max(1, n_heavy) * 1000.0

    # 2. Ramachandran angles (only between residues joined by a peptide bond)
    phi_psi_pairs: List[Tuple[float, float]] = []
    n_favored = n_allowed = n_outliers = 0
    n_cis = n_non_pro_cis = 0

    def _linked(a, b):
        return ("C" in a and "N" in b and a.get("chain") == b.get("chain")
                and np.linalg.norm(np.asarray(a["C"]) - np.asarray(b["N"])) < 2.0)

    if backbone_atoms and len(backbone_atoms) >= 3:
        for i in range(1, len(backbone_atoms) - 1):
            prev_res, curr_res, next_res = backbone_atoms[i - 1], backbone_atoms[i], backbone_atoms[i + 1]
            if not all(k in curr_res for k in ("N", "CA", "C")) or "CA" not in next_res:
                continue
            if not (_linked(prev_res, curr_res) and _linked(curr_res, next_res)):
                continue
            c_prev = np.asarray(prev_res["C"])
            n_curr = np.asarray(curr_res["N"])
            ca_curr = np.asarray(curr_res["CA"])
            c_curr = np.asarray(curr_res["C"])
            n_next = np.asarray(next_res["N"])
            ca_next = np.asarray(next_res["CA"])

            phi = _calc_dihedral(c_prev, n_curr, ca_curr, c_curr)
            psi = _calc_dihedral(n_curr, ca_curr, c_curr, n_next)
            omega = _calc_dihedral(ca_curr, c_curr, n_next, ca_next)
            if np.isnan(phi) or np.isnan(psi):
                continue
            phi_psi_pairs.append((phi, psi))

            omega_prev = (_calc_dihedral(np.asarray(prev_res["CA"]), c_prev, n_curr, ca_curr)
                          if "CA" in prev_res else None)
            cat = rama_category(curr_res.get("res_name", ""), next_res.get("res_name", ""), omega_prev)
            reg = _classify_ramachandran(phi, psi, cat)
            if reg == "FAVORED":
                n_favored += 1
            elif reg == "ALLOWED":
                n_allowed += 1
            else:
                n_outliers += 1

            # Cis peptide bond between this residue and the next (|omega| <= 30)
            if not np.isnan(omega) and abs(omega) <= 30.0:
                n_cis += 1
                if next_res.get("res_name", "") != "PRO":
                    n_non_pro_cis += 1

    n_angles = len(phi_psi_pairs)
    if n_angles > 0:
        frac_fav = n_favored / n_angles
        frac_allow = n_allowed / n_angles
        frac_out = n_outliers / n_angles
    else:
        frac_fav = frac_allow = frac_out = float("nan")

    # 3. Overall Stereochemistry Status
    rama_txt = (f"Ramachandran Favored = {frac_fav*100:.1f}%, Outliers = {frac_out*100:.1f}%"
                if n_angles else "Ramachandran not evaluated (no backbone)")
    rama_ok = n_angles == 0 or (frac_fav >= min_rama_favored_pass and frac_out <= max_rama_outliers_pass)
    rama_warn = n_angles == 0 or (frac_fav >= 0.75 and frac_out <= 0.08)
    if clashscore <= max_clashscore_pass and rama_ok:
        status = "PASS"
        diag = f"Good stereochemical quality (heavy-atom clash score = {clashscore:.1f}, {rama_txt})."
    elif clashscore <= 25.0 and rama_warn:
        status = "WARNING"
        diag = f"Moderate stereochemical strain (heavy-atom clash score = {clashscore:.1f}, {rama_txt}). Energy minimization recommended."
    else:
        status = "FAIL"
        diag = f"Severe stereochemical anomalies (heavy-atom clash score = {clashscore:.1f}, {rama_txt})."
    if n_non_pro_cis:
        diag += f" {n_non_pro_cis} non-proline cis peptide bond(s)."

    return StereochemistryResult(
        n_clashes=clash_count,
        clashscore=float(clashscore),
        n_ramachandran_angles=n_angles,
        frac_rama_favored=float(frac_fav),
        frac_rama_allowed=float(frac_allow),
        frac_rama_outliers=float(frac_out),
        n_cis_peptides=n_cis,
        n_non_proline_cis_peptides=n_non_pro_cis,
        status=status,
        diagnostic_message=diag,
        phi_psi_angles=phi_psi_pairs
    )
