"""
Command Line Interface (CLI) for AlphaCert.
"""

import sys
import os
import argparse
import numpy as np

from alphacert import __version__
from alphacert.parsers.structure_io import parse_protein_structure
from alphacert.parsers.pae_json import parse_pae_matrix
from alphacert.core.scoring import assess_alphafold_quality
from alphacert.reporters.plot_generator import generate_alphacert_figures
from alphacert.reporters.manuscript_prep import generate_alphacert_manuscript_assets
from alphacert.reporters.html_report import generate_alphacert_html_report


def print_banner():
    banner = rf"""
     /\   | |         | |             / ____|          | |   
    /  \  | |_ __  __ | |__   __ _   | |     ___ _ __| |_  
   / /\ \ | | '_ \/ _` | '_ \ / _` |  | |    / _ \ '__| __| 
  / ____ \| | |_) | (_| | | | | (_| |  | |___|  __/ |  | |_  
 /_/    \_\_|_.__/\__,_|_| |_|\__,_|   \_____\___|_|   \__| v{__version__}

 Predicted Protein Structure Quality-Control & Certification Toolkit
 Monreal-Hernández et al., 2026
"""
    print(banner)


def run_demo(output_dir: str = "alphacert_demo_output"):
    """
    Executes a benchmark demonstration evaluating an AlphaFold2 protein model (320 residues)
    with per-residue pLDDT, 2D PAE error matrix, stereochemical clashes, and Ramachandran angles.
    """
    print(f"\n[AlphaCert] Running demonstration benchmark on AlphaFold2 Kinase Model...")
    os.makedirs(output_dir, exist_ok=True)
    
    metadata = {
        "name": "Target Kinase Domain",
        "engine": "AlphaFold2 v2.3",
        "uniprot_id": "P00533 (EGFR-like)",
        "prediction_date": "2026-08-31"
    }
    
    # 1. Generate realistic pLDDT profile (320 residues)
    rng = np.random.default_rng(42)
    n_res = 320
    
    plddt = np.zeros(n_res)
    # N-term
    plddt[:15] = rng.normal(48.0, 4.0, 15)
    # Core domain (Alpha helices & Beta sheets)
    plddt[15:305] = rng.normal(93.5, 2.5, 290)
    # C-term
    plddt[305:] = rng.normal(55.0, 5.0, 15)
    plddt = np.clip(plddt, 10.0, 100.0)
    
    # 2. Generate 2D PAE matrix (320 x 320)
    pae_mat = np.ones((n_res, n_res)) * 15.0
    for i in range(n_res):
        for j in range(n_res):
            if 15 <= i < 305 and 15 <= j < 305:
                dist_seq = abs(i - j)
                pae_mat[i, j] = min(8.0, 2.0 + 0.02 * dist_seq + rng.normal(0, 0.2))
            else:
                pae_mat[i, j] = 16.0 + rng.normal(0, 1.5)
    pae_mat = np.clip(pae_mat, 0.5, 31.75)
    np.fill_diagonal(pae_mat, 0.5)
    
    # 3. Clean coordinates
    coords = np.array([[i * 3.5, np.sin(i * 0.4) * 5.0, np.cos(i * 0.4) * 5.0] for i in range(n_res * 4)])
    atom_elems = ["N", "C", "C", "O"] * n_res

    # Binding pocket residues: e.g. active site gatekeeper 65 to 85
    binding_pocket = list(range(65, 86))

    print("  -> Performing per-residue pLDDT audit, PAE matrix analysis, and stereochemical check...")
    report = assess_alphafold_quality(
        metadata=metadata,
        plddt_values=plddt.tolist(),
        pae_matrix=pae_mat,
        coordinates=coords,
        atom_elements=atom_elems,
        backbone_atoms=None,
        binding_site_indices=binding_pocket
    )

    print("  -> Generating publication-ready vector figures (pLDDT profile, PAE matrix, and Ramachandran plot)...")
    generate_alphacert_figures(report, output_dir)

    print("  -> Drafting manuscript Methods text snippet, summary LaTeX tables, and BibTeX citations...")
    assets = generate_alphacert_manuscript_assets(report, output_dir)

    with open(assets["methods_text"], "r", encoding="utf-8") as f:
        methods_txt = f.read()
    with open(assets["citation_bib"], "r", encoding="utf-8") as f:
        bib_txt = f.read()

    html_p = os.path.join(output_dir, "report.html")
    print(f"  -> Writing interactive report to {html_p}...")
    generate_alphacert_html_report(report, html_p, methods_text=methods_txt, citation_bib=bib_txt)

    print("\n" + "="*70)
    print(f" [RESULT] Overall Structure Quality Certification: {report.overall_status}")
    print(f" [SCORE]  {report.validation_score}")
    print("="*70)
    print(f" * Target / Engine  : {report.metadata['name']} ({report.metadata['engine']})")
    print(f" * Mean pLDDT       : {report.plddt_result.mean_plddt:.1f} / 100 ({report.plddt_result.frac_very_high_90*100:.1f}% > 90, {report.plddt_result.frac_confident_70_90*100:.1f}% 70-90)")
    if report.pae_result:
        print(f" * PAE Domain Error : {report.pae_result.mean_pae:.1f} A ({report.pae_result.frac_high_confidence_pae_under_5*100:.1f}% pairs < 5 A)")
    if report.stereochemistry:
        print(f" * Stereochemistry  : Clashscore = {report.stereochemistry.clashscore:.1f} | Ramachandran: {report.stereochemistry.frac_rama_favored*100:.1f}% Favored")
    print(f" * Docking Ready    : {report.docking_readiness['status']} (Pocket pLDDT = {report.docking_readiness['pocket_plddt']:.1f})")
    print(f" * MD Ready         : {report.md_readiness['status']} ({report.md_readiness['disordered_fraction']*100:.1f}% IDRs)")
    print("="*70)
    print(f"\nAll outputs successfully saved to: {os.path.abspath(output_dir)}/")
    print(f"Open {os.path.abspath(html_p)} in your browser to inspect the full report.\n")


def run_assess(args):
    """
    Evaluates user-provided structure PDB/CIF and optional PAE file.
    """
    output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)
    
    input_file = args.input
    if not input_file:
        print("[Error] Please specify a structure PDB or mmCIF file with --input.", file=sys.stderr)
        sys.exit(1)
        
    print(f"\n[AlphaCert] Parsing protein structure from {input_file}...")
    struct_data = parse_protein_structure(input_file)
    
    pae_mat = None
    iptm_val = None
    if args.pae:
        print(f"  -> Parsing PAE matrix from {args.pae}...")
        pae_data = parse_pae_matrix(args.pae)
        pae_mat = pae_data["pae_matrix"]
        iptm_val = pae_data.get("iptm_score")
        
    meta = {
        "name": args.name or os.path.splitext(os.path.basename(input_file))[0],
        "engine": args.engine or "AlphaFold / ESMFold",
        "file": input_file
    }
    
    # Parse pocket indices if given as comma-separated string
    pocket_idx = None
    if args.pocket:
        pocket_idx = [int(x.strip()) for x in args.pocket.split(",") if x.strip().isdigit()]

    print("  -> Performing structure validation and downstream application certification...")
    report = assess_alphafold_quality(
        metadata=meta,
        plddt_values=struct_data["per_residue_plddt"],
        pae_matrix=pae_mat,
        coordinates=struct_data["coordinates"],
        atom_elements=struct_data["atom_elements"],
        backbone_atoms=struct_data["backbone_atoms"],
        binding_site_indices=pocket_idx,
        chain_boundaries=struct_data.get("chain_lengths"),
        iptm_score=iptm_val
    )

    print("  -> Generating publication figures...")
    generate_alphacert_figures(report, output_dir)

    print("  -> Generating manuscript text, LaTeX summary table, and BibTeX citations...")
    assets = generate_alphacert_manuscript_assets(report, output_dir)

    with open(assets["methods_text"], "r", encoding="utf-8") as f:
        methods_txt = f.read()
    with open(assets["citation_bib"], "r", encoding="utf-8") as f:
        bib_txt = f.read()

    html_p = os.path.join(output_dir, "report.html")
    print(f"  -> Writing HTML quality report to {html_p}...")
    generate_alphacert_html_report(report, html_p, methods_text=methods_txt, citation_bib=bib_txt)

    print("\n" + "="*70)
    print(f" [RESULT] Overall Structure Quality: {report.overall_status}")
    print(f" [SCORE]  {report.validation_score}")
    print("="*70)
    print(f" * Mean pLDDT       : {report.plddt_result.mean_plddt:.1f} / 100")
    print(f" * Docking Ready    : {report.docking_readiness['status']}")
    print(f" * MD Ready         : {report.md_readiness['status']}")
    print("="*70)
    print(f"\nReport ready at: {os.path.abspath(html_p)}\n")


def print_citation():
    bib = """@software{monreal2026alphacert,
  author = {Monreal-Hern\\'andez, Andre},
  title = {{AlphaCert: An Open-Source Toolkit for Quality-Control, pLDDT/PAE Assessment, and Stereochemical Certification of Predicted Protein Structures}},
  year = {2026},
  version = {1.0.0},
  publisher = {Zenodo},
  url = {https://github.com/amonreal/alphacert}
}"""
    print("\nIf you use AlphaCert in your publications, please cite:\n")
    print("APA Style:")
    print("Monreal-Hernández, A. (2026). AlphaCert: An Open-Source Toolkit for Quality-Control, pLDDT/PAE Assessment, and Stereochemical Certification of Predicted Protein Structures (v1.0.0). Zenodo. https://github.com/amonreal/alphacert\n")
    print("BibTeX:")
    print(bib)
    print()


def main():
    parser = argparse.ArgumentParser(
        prog="alphacert",
        description="AlphaCert: Quality-Control, pLDDT/PAE Assessment, and Stereochemical Certification for Predicted Protein Structures."
    )
    parser.add_argument("-v", "--version", action="version", version=f"alphacert {__version__}")

    subparsers = parser.add_subparsers(dest="command", help="Available subcommands")

    # Assess command
    assess_parser = subparsers.add_parser("assess", help="Assess predicted protein structure PDB/CIF and PAE")
    assess_parser.add_argument("-i", "--input", required=True, help="Path to predicted structure (.pdb, .cif)")
    assess_parser.add_argument("--pae", default=None, help="Path to Predicted Aligned Error file (.json, .npy)")
    assess_parser.add_argument("-o", "--output", default="alphacert_output", help="Directory for output report and assets (default: alphacert_output)")
    assess_parser.add_argument("--name", default=None, help="Protein name / target description")
    assess_parser.add_argument("--engine", default="AlphaFold2", help="Prediction engine (AlphaFold2, AlphaFold3, ColabFold, ESMFold)")
    assess_parser.add_argument("--pocket", default=None, help="Comma-separated residue indices of binding pocket (e.g. '45,46,47,120')")

    # Demo command
    demo_parser = subparsers.add_parser("demo", help="Run AlphaCert on a benchmark AlphaFold2 kinase model")
    demo_parser.add_argument("-o", "--output", default="alphacert_demo_output", help="Output directory (default: alphacert_demo_output)")

    # Cite command
    subparsers.add_parser("cite", help="Display BibTeX and APA citation details")

    if len(sys.argv) == 1:
        print_banner()
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    if args.command == "assess":
        print_banner()
        run_assess(args)
    elif args.command == "demo":
        print_banner()
        run_demo(args.output)
    elif args.command == "cite":
        print_banner()
        print_citation()
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
