"""
Quickstart API tutorial for AlphaCert.
"""

import os
import sys

# Ensure current script dir is in sys.path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from alphacert import assess_alphafold_quality
from alphacert.parsers import parse_protein_structure, parse_pae_matrix
from alphacert.reporters import (
    generate_alphacert_figures,
    generate_alphacert_manuscript_assets,
    generate_alphacert_html_report
)
from generate_sample_af_models import generate_sample_af_pdb, generate_sample_pae_json


def main():
    print("Running AlphaCert Python API quickstart tutorial...")
    output_dir = "quickstart_alphacert_output"
    os.makedirs(output_dir, exist_ok=True)
    
    pdb_file = os.path.join(output_dir, "sample_alphafold_model.pdb")
    pae_file = os.path.join(output_dir, "sample_alphafold_pae.json")
    generate_sample_af_pdb(pdb_file)
    generate_sample_pae_json(pae_file)
    
    # 1. Parse structure and PAE matrix
    struct_data = parse_protein_structure(pdb_file)
    pae_data = parse_pae_matrix(pae_file)
    
    # 2. Assess structure quality
    report = assess_alphafold_quality(
        metadata={"name": "Sample Target", "engine": "AlphaFold2"},
        plddt_values=struct_data["per_residue_plddt"],
        pae_matrix=pae_data["pae_matrix"],
        coordinates=struct_data["coordinates"],
        atom_elements=struct_data["atom_elements"],
        backbone_atoms=struct_data["backbone_atoms"],
        binding_site_indices=[35, 36, 37, 38, 70, 71]
    )
    
    print(f"\nOverall Certification: {report.overall_status}")
    print(f"Validation Score: {report.validation_score}")
    print(f"Mean pLDDT: {report.plddt_result.mean_plddt:.1f} / 100")
    print(f"Docking Readiness: {report.docking_readiness['status']} ({report.docking_readiness['message']})")
    print(f"MD Readiness: {report.md_readiness['status']}")
    
    # 3. Export all publication assets
    generate_alphacert_figures(report, output_dir)
    assets = generate_alphacert_manuscript_assets(report, output_dir)
    
    with open(assets["methods_text"], "r", encoding="utf-8") as f:
        methods_txt = f.read()
    with open(assets["citation_bib"], "r", encoding="utf-8") as f:
        bib_txt = f.read()
        
    html_p = os.path.join(output_dir, "report.html")
    generate_alphacert_html_report(report, html_p, methods_text=methods_txt, citation_bib=bib_txt)
    
    print(f"\nCompleted! HTML report available at: {os.path.abspath(html_p)}")


if __name__ == "__main__":
    main()
