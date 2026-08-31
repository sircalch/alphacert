"""
Tests for structure PDB parser and PAE JSON parser.
"""

import os
import json
import tempfile
import numpy as np
import pytest
from alphacert.parsers.structure_io import parse_protein_structure
from alphacert.parsers.pae_json import parse_pae_matrix


def test_pdb_parser():
    pdb_content = """ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 95.50           N  
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 95.50           C  
ATOM      3  C   ALA A   1       2.000   1.400   0.000  1.00 95.50           C  
ATOM      4  O   ALA A   1       1.300   2.400   0.000  1.00 95.50           O  
ATOM      5  N   VAL A   2       3.300   1.400   0.000  1.00 88.20           N  
ATOM      6  CA  VAL A   2       4.000   2.700   0.000  1.00 88.20           C  
ATOM      7  C   VAL A   2       5.500   2.700   0.000  1.00 88.20           C  
ATOM      8  O   VAL A   2       6.200   1.700   0.000  1.00 88.20           O  
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as f:
        f.write(pdb_content)
        f_path = f.name

    try:
        data = parse_protein_structure(f_path)
        assert data["n_residues"] == 2
        assert len(data["coordinates"]) == 8
        assert data["per_residue_plddt"] == [95.50, 88.20]
        assert data["chains"] == ["A"]
    finally:
        if os.path.exists(f_path):
            os.remove(f_path)


def test_pae_json_parser():
    pae_dict = {
        "predicted_aligned_error": [
            [2.5, 4.0],
            [4.0, 2.5]
        ],
        "iptm": 0.85
    }
    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
        json.dump(pae_dict, f)
        f_path = f.name

    try:
        data = parse_pae_matrix(f_path)
        assert data["pae_matrix"].shape == (2, 2)
        assert data["iptm_score"] == 0.85
    finally:
        if os.path.exists(f_path):
            os.remove(f_path)
