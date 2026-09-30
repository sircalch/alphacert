"""
Tests for stereochemical analysis: clashes and Ramachandran angles.
"""

import numpy as np
import pytest
from alphacert.core.geometry import evaluate_stereochemistry, _calc_dihedral, _classify_ramachandran


def test_dihedral_calc():
    p1 = np.array([1.0, 0.0, 0.0])
    p2 = np.array([0.0, 0.0, 0.0])
    p3 = np.array([0.0, 1.0, 0.0])
    p4 = np.array([0.0, 1.0, 1.0])
    
    ang = _calc_dihedral(p1, p2, p3, p4)
    assert np.isclose(ang, -90.0, atol=1e-3)   # IUPAC sign convention (checked against gemmi)


def test_ramachandran_classification():
    # Alpha-helix: phi=-60, psi=-45 -> FAVORED
    assert _classify_ramachandran(-60.0, -45.0) == "FAVORED"
    # Beta-sheet: phi=-120, psi=130 -> FAVORED
    assert _classify_ramachandran(-120.0, 130.0) == "FAVORED"
    # Outlier in the Top8000 general contours: phi=100, psi=-100 (and phi=psi=0)
    assert _classify_ramachandran(100.0, -100.0) == "OUTLIER"
    assert _classify_ramachandran(0.0, 0.0) == "OUTLIER"
    # phi=60, psi=-120 is allowed (not favored) for a general residue in Top8000
    assert _classify_ramachandran(60.0, -120.0) == "ALLOWED"
    # Gly and Pro use their own contours: left-handed region is favored for Gly
    assert _classify_ramachandran(80.0, 10.0, "gly") == "FAVORED"


def test_stereochemistry_evaluation():
    coords = np.array([
        [0.0, 0.0, 0.0],
        [1.5, 0.0, 0.0],
        [2.5, 1.2, 0.0],
        [3.8, 1.2, 0.0]
    ])
    elems = ["N", "CA", "C", "O"]
    
    res = evaluate_stereochemistry(coords, elems)
    assert res.status == "PASS"
    assert res.clashscore <= 10.0
