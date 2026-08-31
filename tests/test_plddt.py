"""
Tests for pLDDT confidence analysis and band distribution.
"""

import numpy as np
import pytest
from alphacert.core.plddt import evaluate_plddt_profile


def test_plddt_high_confidence():
    # 100 residues with pLDDT around 92.0
    vals = [92.0] * 80 + [75.0] * 20
    res = evaluate_plddt_profile(vals)
    
    assert res.status == "PASS"
    assert res.mean_plddt > 85.0
    assert res.frac_very_high_90 == 0.80
    assert res.frac_confident_70_90 == 0.20
    assert res.n_idr_stretches == 0


def test_plddt_disordered_idr_detection():
    # Long stretch of 25 residues with pLDDT = 35.0
    vals = [90.0] * 50 + [35.0] * 25 + [85.0] * 50
    res = evaluate_plddt_profile(vals)
    
    assert res.n_idr_stretches == 1
    assert res.frac_very_low_under_50 > 0.15


def test_plddt_binding_pocket():
    vals = [85.0] * 100
    # Pocket with residues 10, 11, 12 having low pLDDT 55.0
    vals[10] = 55.0
    vals[11] = 58.0
    vals[12] = 52.0
    
    res = evaluate_plddt_profile(vals, binding_site_indices=[10, 11, 12], min_pocket_plddt=80.0)
    assert res.binding_site_plddt < 60.0
    assert res.status == "WARNING"
