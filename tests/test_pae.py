"""
Tests for 2D PAE matrix analysis and complex interfaces.
"""

import numpy as np
import pytest
from alphacert.core.pae import evaluate_pae_matrix


def test_pae_monomer():
    n = 100
    pae = np.ones((n, n)) * 3.5
    res = evaluate_pae_matrix(pae)
    
    assert res.status == "PASS"
    assert res.mean_pae == 3.5
    assert res.frac_high_confidence_pae_under_5 == 1.0


def test_pae_multimer_interface():
    # Chain A: 50 res, Chain B: 50 res
    n = 100
    pae = np.ones((n, n)) * 4.0
    # Cross block high error -> 22.0 A -> non-binding
    pae[:50, 50:] = 22.0
    pae[50:, :50] = 22.0
    
    res = evaluate_pae_matrix(pae, chain_boundaries=[50, 50])
    assert res.inter_chain_ipae == 22.0
    assert res.status == "FAIL"
    assert "Unreliable complex interface" in res.diagnostic_message
