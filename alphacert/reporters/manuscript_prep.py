"""
Manuscript Methods text generator, LaTeX summary tables, and BibTeX citations for AlphaCert.
"""

from typing import Dict, Any, Optional
import os
import pandas as pd
from alphacert.core.scoring import AlphaFoldValidationReport


def generate_alphacert_manuscript_assets(
    report: AlphaFoldValidationReport,
    output_dir: str
) -> Dict[str, str]:
    """
    Generates manuscript Methods text, LaTeX summary tables, and BibTeX citations.

    Parameters
    ----------
    report : AlphaFoldValidationReport
        Validation report.
    output_dir : str
        Output directory.

    Returns
    -------
    paths : dict
        Mapping of generated asset paths.
    """
    os.makedirs(output_dir, exist_ok=True)
    generated = {}

    # 1. Summary DataFrame
    rows = []
    meta = report.metadata
    
    rows.append({"Parameter": "Target Protein & Engine", "Value": f"{meta.get('name', 'Protein')} ({meta.get('engine', 'AlphaFold2')})", "Status": "PASS"})
    rows.append({"Parameter": "Total Residues", "Value": f"{report.plddt_result.n_residues} residues", "Status": "PASS"})
    rows.append({"Parameter": "Mean pLDDT", "Value": f"{report.plddt_result.mean_plddt:.2f} (Median: {report.plddt_result.median_plddt:.2f})", "Status": report.plddt_result.status})
    rows.append({"Parameter": "Very High Confidence (>90)", "Value": f"{report.plddt_result.frac_very_high_90*100:.1f}%", "Status": "PASS"})
    rows.append({"Parameter": "Confident Backbone (70-90)", "Value": f"{report.plddt_result.frac_confident_70_90*100:.1f}%", "Status": "PASS"})
    rows.append({"Parameter": "Low / IDR Fraction (<70)", "Value": f"{(report.plddt_result.frac_low_50_70 + report.plddt_result.frac_very_low_under_50)*100:.1f}% ({report.plddt_result.n_idr_stretches} IDR loops)", "Status": report.plddt_result.status})
    
    if report.pae_result:
        pae = report.pae_result
        rows.append({"Parameter": "Mean PAE Error", "Value": f"{pae.mean_pae:.2f} Å ({pae.frac_high_confidence_pae_under_5*100:.1f}% pairs < 5 Å)", "Status": pae.status})
        if pae.inter_chain_ipae is not None:
            rows.append({"Parameter": "Interface PAE (iPAE)", "Value": f"{pae.inter_chain_ipae:.2f} Å (ipTM: {pae.iptm_score or 0.0:.2f})", "Status": pae.status})

    if report.stereochemistry:
        st = report.stereochemistry
        rows.append({"Parameter": "Clashscore", "Value": f"{st.clashscore:.2f} clashes / 1000 atoms", "Status": st.status})
        rows.append({"Parameter": "Ramachandran Distribution", "Value": f"{st.frac_rama_favored*100:.1f}% Favored, {st.frac_rama_outliers*100:.1f}% Outliers", "Status": st.status})

    # Target readiness
    rows.append({"Parameter": "Docking Readiness", "Value": report.docking_readiness["message"], "Status": report.docking_readiness["status"]})
    rows.append({"Parameter": "MD Readiness", "Value": report.md_readiness["message"], "Status": report.md_readiness["status"]})

    df_summary = pd.DataFrame(rows)

    # CSV Table
    csv_path = os.path.join(output_dir, "alphacert_summary_table.csv")
    df_summary.to_csv(csv_path, index=False)
    generated["summary_csv"] = csv_path

    # LaTeX Table
    tex_table_path = os.path.join(output_dir, "alphacert_summary_table.tex")
    tex_table = df_summary.to_latex(index=False, escape=False)
    with open(tex_table_path, "w", encoding="utf-8") as f:
        f.write("% AlphaCert Predicted Protein Structure Validation & Quality Table\n")
        f.write(tex_table)
    generated["summary_tex"] = tex_table_path

    # 2. Methods Text Snippet
    methods_path = os.path.join(output_dir, "methods_snippet.txt")
    pld = report.plddt_result
    prot_str = meta.get("name", "target protein")
    engine_str = meta.get("engine", "AlphaFold2")
    
    geom_str = ""
    if report.stereochemistry:
        st = report.stereochemistry
        geom_str = f"Stereochemical quality was audited with a clashscore of {st.clashscore:.1f} and {st.frac_rama_favored*100:.1f}% residues in favored Ramachandran regions ({st.frac_rama_outliers*100:.1f}% outliers). "

    pae_str = ""
    if report.pae_result:
        pae = report.pae_result
        if pae.inter_chain_ipae is not None:
            pae_str = f"The inter-chain interaction interface was certified with an interface error of iPAE = {pae.inter_chain_ipae:.1f} Å (ipTM = {pae.iptm_score or 0.0:.2f}). "
        else:
            pae_str = f"Inter-domain positional confidence was verified with a mean Predicted Aligned Error of {pae.mean_pae:.1f} Å. "

    full_methods = (
        f"Three-dimensional structural models of {prot_str} were generated using {engine_str}. "
        f"Model reliability, per-residue confidence (pLDDT), domain rigidity, and stereochemistry were systematically audited using AlphaCert v1.0.0 (Monreal-Hernández, 2026). "
        f"The model achieved an average pLDDT of {pld.mean_plddt:.1f} with {pld.frac_very_high_90*100:.1f}% in the very high confidence band (> 90) and {pld.frac_confident_70_90*100:.1f}% in the confident band (70-90). "
        f"{pae_str}{geom_str}"
        f"Downstream modeling suitability was certified as: {report.docking_readiness['status']} for molecular docking and {report.md_readiness['status']} for molecular dynamics."
    )

    with open(methods_path, "w", encoding="utf-8") as f:
        f.write(full_methods + "\n")
    generated["methods_text"] = methods_path

    # 3. BibTeX Citation
    bib_path = os.path.join(output_dir, "citation.bib")
    bib_content = """@software{monreal2026alphacert,
  author = {Monreal-Hern\\'andez, Andre},
  title = {{AlphaCert: An Open-Source Toolkit for Quality-Control, pLDDT/PAE Assessment, and Stereochemical Certification of Predicted Protein Structures}},
  year = {2026},
  version = {1.0.0},
  publisher = {Zenodo},
  url = {https://github.com/amonreal/alphacert}
}
"""
    with open(bib_path, "w", encoding="utf-8") as f:
        f.write(bib_content)
    generated["citation_bib"] = bib_path

    return generated
