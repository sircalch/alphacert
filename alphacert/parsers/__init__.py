"""
Parsers for protein structures (PDB/mmCIF) and AlphaFold PAE confidence matrices.
"""

from alphacert.parsers.structure_io import parse_protein_structure
from alphacert.parsers.pae_json import parse_pae_matrix

__all__ = [
    "parse_protein_structure",
    "parse_pae_matrix"
]
