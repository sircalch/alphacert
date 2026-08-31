"""
Generates sample AlphaFold PDB structure and PAE JSON file for testing and tutorials.
"""

import os
import json
import numpy as np


def generate_sample_af_pdb(filepath: str = "sample_alphafold_model.pdb", n_res: int = 150):
    lines = []
    atom_id = 1
    rng = np.random.default_rng(42)
    
    for r in range(1, n_res + 1):
        plddt = 92.5 if 10 <= r <= 130 else (45.0 if r < 10 else 55.0)
        plddt += rng.normal(0, 2.0)
        plddt = max(10.0, min(100.0, plddt))
        
        ca_x = r * 3.8
        ca_y = np.sin(r * 0.2) * 5.0
        ca_z = np.cos(r * 0.2) * 5.0
        
        # N, CA, C, O
        atoms = [
            ("N", ca_x - 1.2, ca_y, ca_z, "N"),
            ("CA", ca_x, ca_y, ca_z, "C"),
            ("C", ca_x + 1.2, ca_y, ca_z, "C"),
            ("O", ca_x + 1.2, ca_y + 1.2, ca_z, "O")
        ]
        
        for name, x, y, z, elem in atoms:
            line = f"ATOM  {atom_id:5d}  {name:<4s}ALA A{r:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 {plddt:5.2f}           {elem:>2s}"
            lines.append(line)
            atom_id += 1

    lines.append("END\n")
    with open(filepath, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"Sample AlphaFold PDB written to: {os.path.abspath(filepath)}")


def generate_sample_pae_json(filepath: str = "sample_alphafold_pae.json", n_res: int = 150):
    pae_mat = np.ones((n_res, n_res)) * 18.0
    rng = np.random.default_rng(42)
    for i in range(n_res):
        for j in range(n_res):
            if 10 <= i < 130 and 10 <= j < 130:
                dist = abs(i - j)
                pae_mat[i, j] = min(12.0, 2.0 + 0.04 * dist + rng.normal(0, 0.3))
    np.fill_diagonal(pae_mat, 0.5)
    
    data = [{
        "predicted_aligned_error": np.round(pae_mat, 2).tolist(),
        "max_predicted_aligned_error": 31.75
    }]
    with open(filepath, "w", encoding="utf-8") as f:
        json.dump(data, f)
    print(f"Sample PAE JSON written to: {os.path.abspath(filepath)}")


if __name__ == "__main__":
    generate_sample_af_pdb()
    generate_sample_pae_json()
