"""
AlphaCert: Automated Quality-Control, pLDDT/PAE Assessment, and Stereochemical
Certification for Predicted Protein Structures (AlphaFold2/3, ColabFold, ESMFold).
"""

__version__ = "1.0.0"
__author__ = "Andres Monreal-Hernández"
__license__ = "MIT"

from alphacert.core.plddt import (
    evaluate_plddt_profile,
    PLDDTAnalysisResult
)
from alphacert.core.pae import (
    evaluate_pae_matrix,
    PAEAnalysisResult
)
from alphacert.core.geometry import (
    evaluate_stereochemistry,
    StereochemistryResult
)
from alphacert.core.scoring import assess_alphafold_quality, AlphaFoldValidationReport

__all__ = [
    "__version__",
    "evaluate_plddt_profile",
    "PLDDTAnalysisResult",
    "evaluate_pae_matrix",
    "PAEAnalysisResult",
    "evaluate_stereochemistry",
    "StereochemistryResult",
    "assess_alphafold_quality",
    "AlphaFoldValidationReport"
]
