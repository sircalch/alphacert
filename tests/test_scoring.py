"""
Tests for scoring, reporting, and CLI demo execution in AlphaCert.
"""

import os
import tempfile
import numpy as np
import pytest
from alphacert.core.scoring import assess_alphafold_quality
from alphacert.reporters.plot_generator import generate_alphacert_figures
from alphacert.reporters.manuscript_prep import generate_alphacert_manuscript_assets
from alphacert.reporters.html_report import generate_alphacert_html_report
from alphacert.cli import run_demo


def test_full_alphafold_validation_pipeline():
    meta = {
        "name": "Kinase Target",
        "engine": "AlphaFold2",
        "uniprot_id": "P12345"
    }
    
    plddt = [92.0] * 80 + [75.0] * 15 + [45.0] * 5
    n = len(plddt)
    pae = np.ones((n, n)) * 3.5
    # Well-separated coordinates along helix
    coords = np.array([[i * 2.5, np.sin(i*0.5)*4.0, np.cos(i*0.5)*4.0] for i in range(n * 4)])
    elems = ["N", "CA", "C", "O"] * n
    
    report = assess_alphafold_quality(
        metadata=meta,
        plddt_values=plddt,
        pae_matrix=pae,
        coordinates=coords,
        atom_elements=elems
    )
    
    assert report.overall_status in ["PASS", "WARNING"]
    assert report.docking_readiness["status"] in ["PASS", "WARNING"]
    assert report.md_readiness["status"] in ["PASS", "WARNING"]
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Plot generation
        plots = generate_alphacert_figures(report, tmpdir, formats=["png", "svg"])
        assert len(plots) > 0
        for p in plots:
            assert os.path.exists(p)
            
        # Manuscript assets
        assets = generate_alphacert_manuscript_assets(report, tmpdir)
        assert os.path.exists(assets["summary_csv"])
        assert os.path.exists(assets["summary_tex"])
        assert os.path.exists(assets["methods_text"])
        assert os.path.exists(assets["citation_bib"])
        
        # HTML report
        html_p = os.path.join(tmpdir, "report.html")
        generate_alphacert_html_report(report, html_p, methods_text="Sample methods", citation_bib="@software{}")
        assert os.path.exists(html_p)
        assert os.path.getsize(html_p) > 500


def test_cli_demo_execution():
    with tempfile.TemporaryDirectory() as tmpdir:
        run_demo(output_dir=tmpdir)
        assert os.path.exists(os.path.join(tmpdir, "report.html"))
        assert os.path.exists(os.path.join(tmpdir, "alphacert_summary_table.csv"))
        assert os.path.exists(os.path.join(tmpdir, "citation.bib"))
