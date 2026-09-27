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
    n_compact_domains: int  # PAE clusters (Croll greedy-modularity) with >= min_domain_size residues
    inter_chain_ipae: Optional[float]  # mean symmetrised PAE over all inter-chain residue pairs
    iptm_score: Optional[float]
    status: str  # 'PASS', 'WARNING', 'FAIL'
    diagnostic_message: str
    pae_matrix: np.ndarray
    domains: Optional[List[List[int]]] = None  # 0-based residue indices of each compact domain


def pae_domains(pae_matrix: np.ndarray, pae_power: float = 1.0, pae_cutoff: float = 5.0,
                resolution: float = 0.5, min_domain_size: int = 30,
                plddt: Optional[List[float]] = None, min_domain_plddt: float = 70.0) -> List[List[int]]:
    """
    Rigid domains from a PAE matrix, following T. Croll's method used by ChimeraX and the AlphaFold
    DB (github.com/tristanic/pae_to_domains): residue pairs with PAE < pae_cutoff are joined by
    edges of weight 1/PAE**pae_power, and the graph is split by greedy modularity maximisation.
    resolution=0.5 is the ChimeraX default. Clusters smaller than min_domain_size residues are
    discarded, and, when per-residue pLDDT is given, so are clusters with mean pLDDT below
    min_domain_plddt (disordered segments have low intra-segment PAE but are not folded domains).
    """
    import networkx as nx
    from networkx.algorithms import community

    pae = np.asarray(pae_matrix, dtype=float)
    weights = 1.0 / np.maximum(pae, 0.2) ** pae_power
    g = nx.Graph()
    g.add_nodes_from(range(pae.shape[0]))
    ii, jj = np.where(pae < pae_cutoff)
    g.add_weighted_edges_from([(int(i), int(j), float(weights[i, j])) for i, j in zip(ii, jj) if i < j])
    clusters = community.greedy_modularity_communities(g, weight="weight", resolution=resolution)
    keep = []
    for c in clusters:
        c = sorted(c)
        if len(c) < min_domain_size:
            continue
        if plddt is not None and len(plddt) == pae.shape[0] and np.mean(np.asarray(plddt)[c]) < min_domain_plddt:
            continue
        keep.append(c)
    return sorted(keep, key=lambda c: c[0])


def evaluate_pae_matrix(
    pae_matrix: np.ndarray,
    chain_boundaries: Optional[List[int]] = None,
    iptm_score: Optional[float] = None,
    warn_ipae_cutoff: float = 10.0,
    fail_ipae_cutoff: float = 15.0,
    plddt: Optional[List[float]] = None
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
        AlphaFold-Multimer / AlphaFold3 interface pTM score (primary criterion for complexes).
    plddt : list of float, optional
        Per-residue pLDDT, used to discard disordered PAE clusters from the domain count.

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
    
    # Compact domains (Croll PAE clustering)
    domains = pae_domains(pae, plddt=plddt)
    n_domains = len(domains)

    # Multimer / inter-chain analysis: mean symmetrised PAE over all residue pairs of different chains
    inter_ipae = None
    if chain_boundaries and len(chain_boundaries) > 1:
        if sum(chain_boundaries) != n:
            raise ValueError(f"Chain lengths {chain_boundaries} do not add up to the PAE size {n}")
        chain_id = np.repeat(np.arange(len(chain_boundaries)), chain_boundaries)
        inter = chain_id[:, None] != chain_id[None, :]
        sym = 0.5 * (pae + pae.T)
        inter_ipae = float(np.mean(sym[inter]))

    # Determine status. For complexes the ipTM, when available, is the primary criterion
    # (AlphaFold: > 0.8 confident, 0.6-0.8 grey zone, < 0.6 likely failed prediction).
    if inter_ipae is not None:
        if iptm_score is not None:
            if iptm_score > 0.80:
                status = "PASS"
                diag = f"High-confidence complex interface (ipTM = {iptm_score:.2f} > 0.80, inter-chain PAE = {inter_ipae:.1f} A)."
            elif iptm_score >= 0.60:
                status = "WARNING"
                diag = f"Grey-zone complex interface (ipTM = {iptm_score:.2f}, 0.60-0.80; inter-chain PAE = {inter_ipae:.1f} A). Verify experimentally or with an orthogonal method."
            else:
                status = "FAIL"
                diag = f"Unreliable complex interface (ipTM = {iptm_score:.2f} < 0.60, inter-chain PAE = {inter_ipae:.1f} A). AlphaFold likely predicts a non-binding or random contact."
        elif inter_ipae <= warn_ipae_cutoff:
            status = "PASS"
            diag = f"High-confidence complex interface (inter-chain PAE = {inter_ipae:.1f} A <= {warn_ipae_cutoff} A; no ipTM supplied)."
        elif inter_ipae <= fail_ipae_cutoff:
            status = "WARNING"
            diag = f"Moderate complex interface confidence (inter-chain PAE = {inter_ipae:.1f} A). Interface orientation has significant uncertainty."
        else:
            status = "FAIL"
            diag = f"Unreliable complex interface (inter-chain PAE = {inter_ipae:.1f} A > {fail_ipae_cutoff} A). AlphaFold likely predicts a non-binding or random contact."
    else:
        if f_under_5 >= 0.40 or mean_p <= 8.0:
            status = "PASS"
            diag = f"Well-defined relative domain orientations ({f_under_5*100:.1f}% matrix pairs with PAE < 5 A, {n_domains} compact domain(s))."
        elif f_under_5 >= 0.20 or mean_p <= 15.0:
            status = "WARNING"
            diag = f"{n_domains} compact domain(s) with uncertain relative placement (mean PAE = {mean_p:.1f} A); treat inter-domain orientations with caution."
        else:
            status = "FAIL"
            diag = f"High positional uncertainty across the sequence (mean PAE = {mean_p:.1f} A, {n_domains} compact domain(s))."

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
        pae_matrix=pae,
        domains=domains
    )
