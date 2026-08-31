"""
pLDDT per-residue confidence analysis, band distribution, and IDR detection.
"""

from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass
import numpy as np


@dataclass
class PLDDTAnalysisResult:
    n_residues: int
    mean_plddt: float
    median_plddt: float
    std_plddt: float
    frac_very_high_90: float
    frac_confident_70_90: float
    frac_low_50_70: float
    frac_very_low_under_50: float
    n_idr_stretches: int
    binding_site_plddt: Optional[float]
    status: str  # 'PASS', 'WARNING', 'FAIL'
    diagnostic_message: str
    per_residue_plddt: List[float]


def evaluate_plddt_profile(
    plddt_values: List[float],
    binding_site_indices: Optional[List[int]] = None,
    min_pass_mean: float = 70.0,
    min_pass_frac_70: float = 0.60,
    min_pocket_plddt: float = 80.0
) -> PLDDTAnalysisResult:
    """
    Evaluates residue-level pLDDT confidence scores from AlphaFold / ESMFold.

    Parameters
    ----------
    plddt_values : list of float
        Per-residue pLDDT values on a 0-100 scale.
    binding_site_indices : list of int, optional
        Residue indices (0-indexed or 1-indexed) forming the active/binding site.
    min_pass_mean : float, default 70.0
        Minimum average pLDDT for overall structure pass.
    min_pass_frac_70 : float, default 0.60
        Minimum fraction of confident residues (pLDDT >= 70).
    min_pocket_plddt : float, default 80.0
        Minimum average pLDDT required for binding pocket.

    Returns
    -------
    result : PLDDTAnalysisResult
        Calculated confidence metrics and certification status.
    """
    vals = np.asarray(plddt_values, dtype=float)
    n = len(vals)
    
    if n == 0:
        return PLDDTAnalysisResult(
            n_residues=0,
            mean_plddt=0.0,
            median_plddt=0.0,
            std_plddt=0.0,
            frac_very_high_90=0.0,
            frac_confident_70_90=0.0,
            frac_low_50_70=0.0,
            frac_very_low_under_50=0.0,
            n_idr_stretches=0,
            binding_site_plddt=None,
            status="FAIL",
            diagnostic_message="Empty residue pLDDT array.",
            per_residue_plddt=[]
        )
        
    mean_p = float(np.mean(vals))
    med_p = float(np.median(vals))
    std_p = float(np.std(vals))
    
    n_vh = np.sum(vals >= 90.0)
    n_conf = np.sum((vals >= 70.0) & (vals < 90.0))
    n_low = np.sum((vals >= 50.0) & (vals < 70.0))
    n_vl = np.sum(vals < 50.0)
    
    f_vh = float(n_vh / n)
    f_conf = float(n_conf / n)
    f_low = float(n_low / n)
    f_vl = float(n_vl / n)
    f_good = f_vh + f_conf
    
    # IDR stretches (>= 15 continuous residues with pLDDT < 50)
    is_idr = (vals < 50.0).astype(int)
    idr_count = 0
    current_stretch = 0
    for v in is_idr:
        if v == 1:
            current_stretch += 1
            if current_stretch == 15:
                idr_count += 1
        else:
            current_stretch = 0

    # Pocket specific pLDDT
    pocket_mean = None
    pocket_pass = True
    if binding_site_indices is not None and len(binding_site_indices) > 0:
        valid_indices = [i for i in binding_site_indices if 0 <= i < n]
        if valid_indices:
            pocket_mean = float(np.mean(vals[valid_indices]))
            if pocket_mean < min_pocket_plddt:
                pocket_pass = False

    # Validation Status
    if mean_p >= min_pass_mean and f_good >= min_pass_frac_70 and pocket_pass:
        status = "PASS"
        diag = f"High-quality predicted structure (Mean pLDDT = {mean_p:.1f}, {f_good*100:.1f}% >= 70)."
    elif mean_p >= 60.0 or f_good >= 0.50:
        status = "WARNING"
        reasons = []
        if f_good < min_pass_frac_70:
            reasons.append(f"Only {f_good*100:.1f}% residues have pLDDT >= 70")
        if not pocket_pass and pocket_mean is not None:
            reasons.append(f"Binding pocket pLDDT ({pocket_mean:.1f}) below threshold {min_pocket_plddt}")
        if idr_count > 0:
            reasons.append(f"{idr_count} unstructured IDR loop stretches detected")
        diag = f"Moderate prediction confidence: {'; '.join(reasons)}."
    else:
        status = "FAIL"
        diag = f"Low confidence prediction (Mean pLDDT = {mean_p:.1f} < {min_pass_mean}, {f_vl*100:.1f}% disordered residues). Structure is unsuitable for quantitative modeling."

    return PLDDTAnalysisResult(
        n_residues=n,
        mean_plddt=mean_p,
        median_plddt=med_p,
        std_plddt=std_p,
        frac_very_high_90=f_vh,
        frac_confident_70_90=f_conf,
        frac_low_50_70=f_low,
        frac_very_low_under_50=f_vl,
        n_idr_stretches=idr_count,
        binding_site_plddt=pocket_mean,
        status=status,
        diagnostic_message=diag,
        per_residue_plddt=vals.tolist()
    )
