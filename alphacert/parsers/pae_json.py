"""
Parser for AlphaFold and ColabFold Predicted Aligned Error (PAE) files (JSON and NPY).
"""

from typing import Dict, Any, Tuple
import os
import json
import numpy as np


def parse_pae_matrix(filepath: str) -> Dict[str, Any]:
    """
    Parses a PAE matrix from JSON or NPY file.

    Parameters
    ----------
    filepath : str
        Path to JSON or NPY file.

    Returns
    -------
    data : dict
        pae_matrix (np.ndarray), max_predicted_aligned_error, iptm_score.
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
        
    if filepath.endswith(".npy"):
        mat = np.load(filepath)
        return {
            "pae_matrix": mat,
            "max_pae": float(np.max(mat)),
            "iptm_score": None
        }

    with open(filepath, "r", encoding="utf-8") as f:
        data = json.load(f)
        
    pae_arr = None
    iptm_val = None
    
    # 1. Standard AlphaFold JSON: list of dicts or dict with 'predicted_aligned_error'
    if isinstance(data, list) and len(data) > 0:
        item = data[0]
        if isinstance(item, dict):
            if "predicted_aligned_error" in item:
                pae_arr = np.asarray(item["predicted_aligned_error"], dtype=float)
            elif "pae" in item:
                pae_arr = np.asarray(item["pae"], dtype=float)
            if "iptm" in item:
                iptm_val = float(item["iptm"])
    elif isinstance(data, dict):
        if "predicted_aligned_error" in data:
            pae_arr = np.asarray(data["predicted_aligned_error"], dtype=float)
        elif "pae" in data:
            pae_arr = np.asarray(data["pae"], dtype=float)
        elif "distance" in data:
            pae_arr = np.asarray(data["distance"], dtype=float)
            
        if "iptm" in data:
            iptm_val = float(data["iptm"])
        elif "ranking_confidence" in data and "iptm" in str(data):
            try:
                iptm_val = float(data.get("iptm", 0.0))
            except Exception:
                pass

    if pae_arr is None:
        raise ValueError(f"Could not locate PAE matrix array in JSON file: {filepath}")

    return {
        "pae_matrix": pae_arr,
        "max_pae": float(np.max(pae_arr)),
        "iptm_score": iptm_val
    }
