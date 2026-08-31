# AlphaCert

[![CI](https://github.com/amonreal/alphacert/actions/workflows/test.yml/badge.svg)](https://github.com/amonreal/alphacert/actions)
[![PyPI version](https://img.shields.io/pypi/v/alphacert.svg?color=blue)](https://pypi.org/project/alphacert/)
[![Python versions](https://img.shields.io/pypi/pyversions/alphacert.svg)](https://pypi.org/project/alphacert/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1234571.svg)](https://doi.org/10.5281/zenodo.1234571)

> **Automated Quality-Control, pLDDT/PAE Assessment, and Stereochemical Certification for Predicted Protein Structures (AlphaFold2/3, ColabFold, ESMFold).**

---

## Overview

**AlphaCert** is an open-source scientific toolkit designed to systematically audit, validate, and certify the reliability of 3D protein structures predicted by **AlphaFold2**, **AlphaFold-Multimer**, **AlphaFold3**, **ColabFold**, and **ESMFold**.

Rather than blindly relying on predicted models for molecular docking, molecular dynamics, or drug discovery campaigns, `alphacert` performs an automated multi-criteria structural audit in a single command:

- 📊 **Residue pLDDT Confidence Profile & Classification**:
  - Categorizes residues into the 4 canonical AlphaFold bands:
    - 🔵 **Very High Confidence ($> 90$)**: Suitable for atomic side-chain modeling and docking.
    - 🔷 **Confident ($70 - 90$)**: Reliable backbone topology.
    - 🟡 **Low Confidence ($50 - 70$)**: Loop regions or flexible segments.
    - 🟠 **Very Low Confidence ($< 50$)**: Intrinsically Disordered Regions (IDRs).
- 🧩 **2D Predicted Aligned Error (PAE) Matrix Audit**:
  - Identifies rigid structural subdomains from low intra-block PAE ($< 5$ Å).
  - Evaluates inter-domain hinge flexibility and detects artificial inter-chain contacts in multimer complexes ($\text{ipTM} \ge 0.60$, $\text{iPAE} < 8$ Å).
- 📐 **Stereochemical & Geometry Quality**:
  - Heavy-atom steric clashscore ($d_{ij} < r_{\text{vdw}, i} + r_{\text{vdw}, j} - 0.40$ Å).
  - Ramachandran dihedral distribution ($\phi, \psi$ Favored $\ge 88\%$, Outliers $< 3\%$).
  - Non-proline *cis*-peptide bond identification.
- 🎯 **Target Application Readiness Certification**:
  - 💊 **`DOCKING_READY`**: Verifies that active/binding site residues have $\text{pLDDT} \ge 80$.
  - 🌊 **`MD_READY`**: Verifies that disordered terminal loops do not compromise solvated simulation stability.
  - 🤝 **`INTERFACE_VALID`**: Distinguishes authentic protein-protein interfaces from non-binding artifacts.
- 📑 **Publication Deliverables**:
  - Interactive self-contained `report.html` dashboard.
  - Publication vector plots (pLDDT sequence profile, 2D PAE error heatmap, Ramachandran scatter plot) in SVG, PDF, PNG (300 DPI).
  - Ready-to-compile LaTeX summary tables (`.tex`).
  - Draft **Methods** text snippet and BibTeX citation (`citation.bib`).

```
   AlphaFold PDB / mmCIF + PAE JSON
                  │
                  ▼
  ┌───────────────────────────────────────────────────────────┐
  │                         AlphaCert                         │
  │  ├── Residue-level pLDDT Band Categorization & IDR Audit  │
  │  ├── 2D Predicted Aligned Error (PAE) Domain Decomposition│
  │  ├── Steric Clashscore & Ramachandran Geometry Check      │
  │  └── Downstream Certification (Docking-Ready / MD-Ready)  │
  └───────────────────────────────────────────────────────────┘
                  │
                  ▼
  ┌───────────────────────────────────────────────────────────┐
  │                   Publication Deliverables                │
  │  ├── report.html (Interactive Dashboard & Badges)         │
  │  ├── alphacert_plddt_profile.pdf/svg/png                  │
  │  ├── alphacert_pae_matrix.pdf/svg/png                     │
  │  ├── alphacert_ramachandran_plot.pdf/svg/png              │
  │  ├── alphacert_summary_table.tex / .csv                   │
  │  ├── methods_snippet.txt (Ready for Manuscript)           │
  │  └── citation.bib (BibTeX Reference)                      │
  └───────────────────────────────────────────────────────────┘
```

---

## Installation

### From PyPI
```bash
pip install alphacert
```

### From Source
```bash
git clone https://github.com/amonreal/alphacert.git
cd alphacert
pip install -e .[dev]
```

---

## Quickstart (CLI)

### 1. Run Benchmark Demo (Instant AlphaFold Kinase Model Audit)
```bash
alphacert demo -o my_protein_audit/
```
Open `my_protein_audit/report.html` in any browser to inspect the interactive report!

### 2. Assess Predicted PDB Structure
```bash
alphacert assess -i predicted_structure.pdb --pae predicted_pae.json -o struct_report/
```

### 3. Certify Binding Pocket for Molecular Docking
```bash
alphacert assess -i model.pdb --pocket "45,46,47,78,120,121" -o docking_check/
```

---

## Python API Usage

```python
from alphacert import assess_alphafold_quality
from alphacert.parsers import parse_protein_structure, parse_pae_matrix
from alphacert.reporters import generate_alphacert_figures, generate_alphacert_manuscript_assets, generate_alphacert_html_report

# 1. Parse structure and PAE matrix
struct = parse_protein_structure("alphafold_model.pdb")
pae = parse_pae_matrix("alphafold_pae.json")

# 2. Assess structure quality
report = assess_alphafold_quality(
    metadata={"name": "Target Kinase", "engine": "AlphaFold2"},
    plddt_values=struct["per_residue_plddt"],
    pae_matrix=pae["pae_matrix"],
    coordinates=struct["coordinates"],
    atom_elements=struct["atom_elements"],
    backbone_atoms=struct["backbone_atoms"],
    binding_site_indices=[45, 46, 47, 78, 120]
)

print(f"Overall Certification: {report.overall_status}")
print(f"Mean pLDDT: {report.plddt_result.mean_plddt:.1f}")
print(f"Docking Readiness: {report.docking_readiness['status']}")

# 3. Export publication deliverables
generate_alphacert_figures(report, "output_dir/")
generate_alphacert_manuscript_assets(report, "output_dir/")
generate_alphacert_html_report(report, "output_dir/report.html")
```

---

## Citation

If you use AlphaCert in your research, please cite:

```bibtex
@software{monreal2026alphacert,
  author = {Monreal-Hern{\'a}ndez, Andre},
  title = {{AlphaCert: An Open-Source Toolkit for Quality-Control, pLDDT/PAE Assessment, and Stereochemical Certification of Predicted Protein Structures}},
  year = {2026},
  version = {1.0.0},
  publisher = {Zenodo},
  url = {https://github.com/amonreal/alphacert}
}
```

---

## License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.
