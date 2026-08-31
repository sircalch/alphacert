"""
Publication-ready vector figure generation for predicted protein structure quality.
"""

from typing import List, Optional
import os
import numpy as np
import matplotlib.pyplot as plt
from alphacert.core.scoring import AlphaFoldValidationReport

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.titlesize': 14,
    'figure.dpi': 300,
    'lines.linewidth': 2.0,
    'grid.alpha': 0.3,
    'grid.linestyle': '--'
})


def generate_alphacert_figures(
    report: AlphaFoldValidationReport,
    output_dir: str,
    formats: List[str] = ("png", "svg", "pdf")
) -> List[str]:
    """
    Generates high-resolution publication charts (pLDDT profile, PAE heatmap, Ramachandran plot).

    Parameters
    ----------
    report : AlphaFoldValidationReport
        Validation report.
    output_dir : str
        Directory to save figures.
    formats : list of str
        Output file formats.

    Returns
    -------
    saved_paths : list of str
        List of generated file paths.
    """
    os.makedirs(output_dir, exist_ok=True)
    saved_files = []

    # 1. pLDDT Profile Chart
    plddt_vals = np.asarray(report.plddt_result.per_residue_plddt, dtype=float)
    if len(plddt_vals) > 0:
        fig, ax = plt.subplots(figsize=(10, 4.5))
        res_idx = np.arange(1, len(plddt_vals) + 1)
        
        # Color bands
        ax.axhspan(90, 100, color="#0053D6", alpha=0.15, label="Very High (>90)")
        ax.axhspan(70, 90, color="#65CBF3", alpha=0.15, label="Confident (70-90)")
        ax.axhspan(50, 70, color="#FFDB13", alpha=0.15, label="Low (50-70)")
        ax.axhspan(0, 50, color="#FF7D45", alpha=0.15, label="Very Low (<50, IDR)")
        
        ax.plot(res_idx, plddt_vals, color="#0f172a", linewidth=1.5)
        ax.axhline(70, color="#0284c7", linestyle="--", linewidth=1.0)
        
        ax.set_xlim(1, len(plddt_vals))
        ax.set_ylim(0, 105)
        ax.set_xlabel("Residue Position")
        ax.set_ylabel("Predicted LDDT (0 - 100)")
        prot = report.metadata.get("name", "Target Protein")
        ax.set_title(f"Per-Residue pLDDT Confidence Profile ({prot}) — Mean: {report.plddt_result.mean_plddt:.1f}")
        ax.grid(True)
        ax.legend(loc="lower right", frameon=True, fontsize=9)
        
        plt.tight_layout()
        for fmt in formats:
            p = os.path.join(output_dir, f"alphacert_plddt_profile.{fmt}")
            plt.savefig(p, dpi=300, bbox_inches="tight")
            saved_files.append(p)
        plt.close()

    # 2. PAE 2D Error Heatmap
    if report.pae_result and report.pae_result.pae_matrix is not None:
        pae_mat = report.pae_result.pae_matrix
        fig, ax = plt.subplots(figsize=(7, 6))
        
        im = ax.imshow(pae_mat, cmap="Greens_r", vmin=0, vmax=30, origin="upper")
        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("Predicted Aligned Error (Å)", rotation=270, labelpad=15)
        
        ax.set_xlabel("Scored Residue")
        ax.set_ylabel("Aligned Residue")
        ax.set_title("Predicted Aligned Error (PAE) Matrix")
        
        plt.tight_layout()
        for fmt in formats:
            p = os.path.join(output_dir, f"alphacert_pae_matrix.{fmt}")
            plt.savefig(p, dpi=300, bbox_inches="tight")
            saved_files.append(p)
        plt.close()

    # 3. Ramachandran Dihedral Plot
    if report.stereochemistry and report.stereochemistry.phi_psi_angles:
        angles = np.asarray(report.stereochemistry.phi_psi_angles)
        fig, ax = plt.subplots(figsize=(6, 6))
        
        # Draw favored region rectangles
        # Alpha-helix: [-160, -35], [-70, 10]
        ax.axvspan(-160, -35, ymin=0.30, ymax=0.52, color="#bbf7d0", alpha=0.5, label="Favored Regions")
        # Beta-sheet: [-180, -50], [80, 180]
        ax.axvspan(-180, -50, ymin=0.72, ymax=1.0, color="#bbf7d0", alpha=0.5)
        
        ax.scatter(angles[:, 0], angles[:, 1], color="#2563eb", s=20, alpha=0.7, edgecolors="none", label="Residue (\u03c6, \u03c8)")
        ax.axhline(0, color="gray", linestyle="--", linewidth=0.6)
        ax.axvline(0, color="gray", linestyle="--", linewidth=0.6)
        
        ax.set_xlim(-180, 180)
        ax.set_ylim(-180, 180)
        ax.set_xlabel(r"$\phi$ Dihedral Angle (degrees)")
        ax.set_ylabel(r"$\psi$ Dihedral Angle (degrees)")
        ax.set_title(f"Ramachandran Plot ({report.stereochemistry.frac_rama_favored*100:.1f}% Favored)")
        ax.grid(True)
        ax.legend(loc="lower left", frameon=True)
        
        plt.tight_layout()
        for fmt in formats:
            p = os.path.join(output_dir, f"alphacert_ramachandran_plot.{fmt}")
            plt.savefig(p, dpi=300, bbox_inches="tight")
            saved_files.append(p)
        plt.close()

    return saved_files
