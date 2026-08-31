"""
Parser for PDB and mmCIF structure files, extracting coordinates, atom types, and pLDDT B-factors.
"""

from typing import Dict, Any, List, Optional, Tuple
import os
import numpy as np


def parse_protein_structure(filepath: str) -> Dict[str, Any]:
    """
    Parses a PDB or mmCIF protein structure predicted by AlphaFold / ESMFold.

    Parameters
    ----------
    filepath : str
        Path to .pdb or .cif file.

    Returns
    -------
    data : dict
        Coordinates, atom elements, per-residue pLDDT, backbone atoms, chain boundaries.
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
        
    coordinates: List[Tuple[float, float, float]] = []
    atom_elements: List[str] = []
    
    # Residue tracking
    residues_dict: Dict[Tuple[str, int], Dict[str, Any]] = {}
    
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("ATOM  ") or line.startswith("HETATM"):
                # PDB Column definitions
                # Atom name: 12-16
                # Res name: 17-20
                # Chain: 21
                # ResSeq: 22-26
                # X, Y, Z: 30-38, 38-46, 46-54
                # B-factor / pLDDT: 60-66
                # Element: 76-78
                try:
                    atom_name = line[12:16].strip()
                    res_name = line[17:20].strip()
                    chain_id = line[21].strip() or "A"
                    res_seq = int(line[22:26].strip())
                    
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    
                    b_factor = float(line[60:66].strip()) if len(line) >= 66 else 50.0
                    element = line[76:78].strip() if len(line) >= 78 else atom_name[0]
                    if not element:
                        element = atom_name[0]
                        
                    coordinates.append((x, y, z))
                    atom_elements.append(element)
                    
                    res_key = (chain_id, res_seq)
                    if res_key not in residues_dict:
                        residues_dict[res_key] = {
                            "chain_id": chain_id,
                            "res_seq": res_seq,
                            "res_name": res_name,
                            "plddt": b_factor,
                            "atoms": {}
                        }
                    residues_dict[res_key]["atoms"][atom_name] = [x, y, z]
                    
                    # Update CA pLDDT if available
                    if atom_name == "CA":
                        residues_dict[res_key]["plddt"] = b_factor
                        
                except Exception:
                    pass

    # Sort residues
    sorted_res_keys = sorted(residues_dict.keys(), key=lambda k: (k[0], k[1]))
    per_residue_plddt = [residues_dict[k]["plddt"] for k in sorted_res_keys]
    
    # Backbone structures
    backbone_atoms = []
    chains_seen = {}
    for k in sorted_res_keys:
        ch = k[0]
        chains_seen[ch] = chains_seen.get(ch, 0) + 1
        res = residues_dict[k]
        item = {"res_name": res["res_name"], "chain": ch, "seq": k[1]}
        for a_type in ["N", "CA", "C", "O"]:
            if a_type in res["atoms"]:
                item[a_type] = res["atoms"][a_type]
        backbone_atoms.append(item)

    chain_lengths = list(chains_seen.values())

    return {
        "coordinates": np.asarray(coordinates, dtype=float) if coordinates else np.zeros((0, 3)),
        "atom_elements": atom_elements,
        "per_residue_plddt": per_residue_plddt,
        "backbone_atoms": backbone_atoms,
        "n_residues": len(per_residue_plddt),
        "chain_lengths": chain_lengths,
        "chains": list(chains_seen.keys())
    }
