"""
AlphaFold Model Certification Decision Matrix and Multi-Application Scoring Engine.
"""

from typing import Dict, Any, Optional, List
from dataclasses import dataclass, asdict

from alphacert.core.plddt import PLDDTAnalysisResult, evaluate_plddt_profile
from alphacert.core.pae import PAEAnalysisResult, evaluate_pae_matrix
from alphacert.core.geometry import StereochemistryResult, evaluate_stereochemistry


@dataclass
class AlphaFoldValidationReport:
    overall_status: str  # 'PASS', 'WARNING', 'FAIL'
    validation_score: str
    metadata: Dict[str, Any]
    plddt_result: PLDDTAnalysisResult
    pae_result: Optional[PAEAnalysisResult]
    stereochemistry: Optional[StereochemistryResult]
    docking_readiness: Dict[str, Any]
    md_readiness: Dict[str, Any]
    interface_readiness: Optional[Dict[str, Any]]
    recommendations: List[str]
    provenance: Dict[str, Any]

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def assess_alphafold_quality(
    metadata: Dict[str, Any],
    plddt_values: List[float],
    pae_matrix: Optional[Any] = None,
    coordinates: Optional[Any] = None,
    atom_elements: Optional[List[str]] = None,
    backbone_atoms: Optional[List[Dict[str, Any]]] = None,
    binding_site_indices: Optional[List[int]] = None,
    chain_boundaries: Optional[List[int]] = None,
    iptm_score: Optional[float] = None
) -> AlphaFoldValidationReport:
    """
    Evaluates predicted protein structure confidence and certifies downstream application readiness.

    Parameters
    ----------
    metadata : dict
        Protein name, prediction engine (AlphaFold2, ColabFold, ESMFold), UniProt ID.
    plddt_values : list of float
        Residue pLDDT values.
    pae_matrix : np.ndarray, optional
        2D PAE matrix.
    coordinates : np.ndarray, optional
        Cartesian coordinates.
    atom_elements : list of str, optional
        Element symbols.
    backbone_atoms : list of dict, optional
        Segmented backbone atoms.
    binding_site_indices : list of int, optional
        Residue indices forming target binding pocket.
    chain_boundaries : list of int, optional
        Lengths of chains in complex.
    iptm_score : float, optional
        Interface pTM score.

    Returns
    -------
    report : AlphaFoldValidationReport
        Comprehensive certification report.
    """
    statuses = []
    recommendations = []

    # 1. pLDDT Confidence Evaluation
    plddt_res = evaluate_plddt_profile(plddt_values, binding_site_indices=binding_site_indices)
    statuses.append(plddt_res.status)
    if plddt_res.status != "PASS":
        recommendations.append(plddt_res.diagnostic_message)

    # 2. PAE Matrix Analysis
    pae_res = None
    if pae_matrix is not None:
        pae_res = evaluate_pae_matrix(
            pae_matrix, chain_boundaries=chain_boundaries, iptm_score=iptm_score
        )
        statuses.append(pae_res.status)
        if pae_res.status != "PASS":
            recommendations.append(pae_res.diagnostic_message)

    # 3. Stereochemistry Evaluation
    geom_res = None
    if coordinates is not None and atom_elements is not None:
        geom_res = evaluate_stereochemistry(
            coordinates, atom_elements, backbone_atoms=backbone_atoms
        )
        statuses.append(geom_res.status)
        if geom_res.status != "PASS":
            recommendations.append(geom_res.diagnostic_message)

    # 4. Downstream Application Readiness Certifications
    # A. Docking Readiness
    pocket_plddt = plddt_res.binding_site_plddt or plddt_res.mean_plddt
    if pocket_plddt >= 80.0 and (geom_res is None or geom_res.clashscore <= 15.0):
        dock_st = "PASS"
        dock_msg = f"Binding site is highly confident (pLDDT = {pocket_plddt:.1f} >= 80) and stereochemically sound for ligand docking."
    elif pocket_plddt >= 70.0:
        dock_st = "WARNING"
        dock_msg = f"Moderate binding site confidence (pLDDT = {pocket_plddt:.1f}). Flexible side-chain optimization recommended prior to docking."
    else:
        dock_st = "FAIL"
        dock_msg = f"Binding site has low confidence (pLDDT = {pocket_plddt:.1f} < 70). High risk of false-positive docking poses."

    docking_readiness = {"status": dock_st, "pocket_plddt": pocket_plddt, "message": dock_msg}

    # B. Molecular Dynamics Readiness
    if plddt_res.frac_very_low_under_50 <= 0.15 and (geom_res is None or geom_res.clashscore <= 15.0):
        md_st = "PASS"
        md_msg = f"Globular fold with low disorder ({plddt_res.frac_very_low_under_50*100:.1f}% IDRs). Structure ready for solvated MD simulation."
    elif plddt_res.frac_very_low_under_50 <= 0.35:
        md_st = "WARNING"
        md_msg = f"Structure contains moderate disordered terminal/loop regions ({plddt_res.frac_very_low_under_50*100:.1f}% IDRs). Truncation of flexible termini recommended."
    else:
        md_st = "FAIL"
        md_msg = f"High fraction of intrinsically disordered regions ({plddt_res.frac_very_low_under_50*100:.1f}% IDRs). Simulation box may suffer from large conformational instability."

    md_readiness = {"status": md_st, "disordered_fraction": plddt_res.frac_very_low_under_50, "message": md_msg}

    # C. Interface Readiness (for multimer complexes)
    interface_readiness = None
    if pae_res and pae_res.inter_chain_ipae is not None:
        if pae_res.status == "PASS":
            if_st = "PASS"
            if_msg = f"Interface is certified as authentic and high-confidence (iPAE = {pae_res.inter_chain_ipae:.1f} A, ipTM = {pae_res.iptm_score or 0.0:.2f})."
        elif pae_res.status == "WARNING":
            if_st = "WARNING"
            if_msg = f"Interface is predicted with moderate confidence (iPAE = {pae_res.inter_chain_ipae:.1f} A)."
        else:
            if_st = "FAIL"
            if_msg = f"Interface is predicted with low confidence / non-binding artifact (iPAE = {pae_res.inter_chain_ipae:.1f} A)."
        interface_readiness = {"status": if_st, "ipae": pae_res.inter_chain_ipae, "iptm": pae_res.iptm_score, "message": if_msg}

    # Overall Decision
    if "FAIL" in statuses:
        overall_status = "FAIL"
        validation_score = "STRUCTURE QUALITY = REJECTED (LOW CONFIDENCE / UNPHYSICAL ANOMALIES)"
    elif "WARNING" in statuses:
        overall_status = "WARNING"
        validation_score = "STRUCTURE QUALITY = ACCEPTABLE WITH WARNINGS"
    else:
        overall_status = "PASS"
        validation_score = "STRUCTURE QUALITY = FULLY CERTIFIED (HIGH CONFIDENCE)"

    return AlphaFoldValidationReport(
        overall_status=overall_status,
        validation_score=validation_score,
        metadata=metadata,
        plddt_result=plddt_res,
        pae_result=pae_res,
        stereochemistry=geom_res,
        docking_readiness=docking_readiness,
        md_readiness=md_readiness,
        interface_readiness=interface_readiness,
        recommendations=recommendations,
        provenance={
            "tool": "AlphaCert",
            "version": "1.0.0",
            "citation": "Monreal-Hernández, A. (2026). AlphaCert: An Open-Source Toolkit for Quality-Control, pLDDT/PAE Assessment, and Stereochemical Certification of Predicted Protein Structures."
        }
    )
