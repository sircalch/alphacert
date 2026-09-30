"""
Predicted Aligned Error (PAE) matrix analysis, domain boundaries, and complex interface certification.
"""

from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass
import numpy as np


@dataclass
class PAEAnalysisResult:
    matrix_shape: Tuple[int, int]
    mean_pae: float
    median_pae: float
    frac_high_confidence_pae_under_5: float
    n_compact_domains: int
    inter_chain_ipae: Optional[float]
    iptm_score: Optional[float]
    status: str  # 'PASS', 'WARNING', 'FAIL'
    diagnostic_message: str
    pae_matrix: np.ndarray


def evaluate_pae_matrix(
    pae_matrix: np.ndarray,
    chain_boundaries: Optional[List[int]] = None,
    iptm_score: Optional[float] = None,
    warn_ipae_cutoff: float = 10.0,
    fail_ipae_cutoff: float = 15.0
) -> PAEAnalysisResult:
    """
    Evaluates 2D Predicted Aligned Error (PAE) matrix for domain rigidity and multimer interface confidence.

    Parameters
    ----------
    pae_matrix : np.ndarray
        (N, N) matrix of expected position error in Angstroms.
    chain_boundaries : list of int, optional
        Residue lengths of each chain in multimer (e.g. [250, 150]).
    iptm_score : float, optional
        AlphaFold-Multimer interface pTM score.

    Returns
    -------
    result : PAEAnalysisResult
        Rigidity, inter-chain interface error, and complex validation status.
    """
    pae = np.asarray(pae_matrix, dtype=float)
    if pae.ndim != 2 or pae.shape[0] != pae.shape[1]:
        raise ValueError(f"PAE matrix must be a square 2D matrix, got shape {pae.shape}")
        
    n = pae.shape[0]
    mean_p = float(np.mean(pae))
    med_p = float(np.median(pae))
    f_under_5 = float(np.sum(pae < 5.0) / (n * n))
    
    # Estimate compact subdomains (simple block diagonal clustering)
    diag_low = pae < 6.0
    # Number of distinct domain blocks
    block_changes = 1
    window = min(30, max(5, n // 10))
    for i in range(window, n - window, window):
        off_diag_error = np.mean(pae[:i, i:])
        if off_diag_error > 12.0:
            block_changes += 1
    n_domains = min(5, block_changes)

    # Multimer / Inter-chain analysis
    inter_ipae = None
    if chain_boundaries and len(chain_boundaries) > 1:
        c1_len = chain_boundaries[0]
        # Cross-block error between Chain A and Chain B
        cross_block = pae[:c1_len, c1_len:]
        inter_ipae = float(np.mean(cross_block))

    # Determine status
    if inter_ipae is not None:
        if inter_ipae <= warn_ipae_cutoff or (iptm_score and iptm_score >= 0.60):
            status = "PASS"
            diag = f"High-confidence complex interface (Interface PAE = {inter_ipae:.1f} A, ipTM = {iptm_score or 0.0:.2f})."
        elif inter_ipae <= fail_ipae_cutoff or (iptm_score and iptm_score >= 0.40):
            status = "WARNING"
            diag = f"Moderate complex interface confidence (Interface PAE = {inter_ipae:.1f} A). Interface orientation has significant uncertainty."
        else:
            status = "FAIL"
            diag = f"Unreliable complex interface (Interface PAE = {inter_ipae:.1f} A > {fail_ipae_cutoff} A). AlphaFold predicts non-binding or random contact."
    else:
        if f_under_5 >= 0.40 or mean_p <= 8.0:
            status = "PASS"
            diag = f"Well-defined relative domain orientations ({f_under_5*100:.1f}% matrix pairs with PAE < 5 A)."
        elif f_under_5 >= 0.20 or mean_p <= 15.0:
            status = "WARNING"
            diag = f"Multi-domain structure with flexible inter-domain hinges (Mean PAE = {mean_p:.1f} A)."
        else:
            status = "FAIL"
            diag = f"High inter-domain positional uncertainty across full sequence (Mean PAE = {mean_p:.1f} A)."

    return PAEAnalysisResult(
        matrix_shape=(n, n),
        mean_pae=mean_p,
        median_pae=med_p,
        frac_high_confidence_pae_under_5=f_under_5,
        n_compact_domains=n_domains,
        inter_chain_ipae=inter_ipae,
        iptm_score=iptm_score,
        status=status,
        diagnostic_message=diag,
        pae_matrix=pae
    )
