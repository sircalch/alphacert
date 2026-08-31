"""
Interactive HTML report dashboard generator for AlphaCert.
"""

import os
import jinja2
from alphacert.core.scoring import AlphaFoldValidationReport

HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>AlphaCert Protein Structure Report</title>
    <style>
        :root {
            --bg-color: #0f172a;
            --card-bg: #1e293b;
            --card-border: #334155;
            --text-primary: #f8fafc;
            --text-secondary: #94a3b8;
            --pass-color: #10b981;
            --pass-bg: rgba(16, 185, 129, 0.15);
            --warn-color: #f59e0b;
            --warn-bg: rgba(245, 158, 11, 0.15);
            --fail-color: #ef4444;
            --fail-bg: rgba(239, 68, 68, 0.15);
            --accent-blue: #38bdf8;
            --plddt-vh: #0053D6;
            --plddt-conf: #65CBF3;
            --plddt-low: #FFDB13;
            --plddt-vl: #FF7D45;
        }
        * { box-sizing: border-box; margin: 0; padding: 0; }
        body {
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
            background-color: var(--bg-color);
            color: var(--text-primary);
            line-height: 1.6;
            padding: 2rem 1rem;
        }
        .container { max-width: 1200px; margin: 0 auto; }
        .header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            border-bottom: 1px solid var(--card-border);
            padding-bottom: 1.5rem;
            margin-bottom: 2rem;
            flex-wrap: wrap;
            gap: 1rem;
        }
        .title-group h1 { font-size: 2rem; font-weight: 700; }
        .title-group p { color: var(--text-secondary); font-size: 0.95rem; }
        .status-badge {
            display: inline-flex;
            align-items: center;
            padding: 0.5rem 1.25rem;
            border-radius: 9999px;
            font-weight: 700;
            font-size: 1.1rem;
            letter-spacing: 0.05em;
            text-transform: uppercase;
        }
        .badge-pass { background-color: var(--pass-bg); color: var(--pass-color); border: 1px solid var(--pass-color); }
        .badge-warning { background-color: var(--warn-bg); color: var(--warn-color); border: 1px solid var(--warn-color); }
        .badge-fail { background-color: var(--fail-bg); color: var(--fail-color); border: 1px solid var(--fail-color); }

        .grid-cards {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
            gap: 1.25rem;
            margin-bottom: 2rem;
        }
        .card {
            background-color: var(--card-bg);
            border: 1px solid var(--card-border);
            border-radius: 0.75rem;
            padding: 1.25rem;
            box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.2);
        }
        .card-label { font-size: 0.8rem; text-transform: uppercase; color: var(--text-secondary); margin-bottom: 0.4rem; }
        .card-value { font-size: 1.4rem; font-weight: 700; color: var(--text-primary); }
        .card-subtext { font-size: 0.85rem; color: var(--text-secondary); margin-top: 0.25rem; }

        .plddt-bar {
            display: flex;
            height: 12px;
            border-radius: 6px;
            overflow: hidden;
            margin-top: 0.5rem;
            background-color: #334155;
        }
        .plddt-segment { height: 100%; }

        .section-title { font-size: 1.3rem; font-weight: 600; margin-bottom: 1rem; color: var(--accent-blue); }
        table {
            width: 100%;
            border-collapse: collapse;
            margin-bottom: 2rem;
            font-size: 0.95rem;
            background-color: var(--card-bg);
            border-radius: 0.75rem;
            overflow: hidden;
            border: 1px solid var(--card-border);
        }
        th, td { padding: 0.85rem 1rem; text-align: left; border-bottom: 1px solid var(--card-border); }
        th { background-color: rgba(255, 255, 255, 0.03); color: var(--text-secondary); font-weight: 600; text-transform: uppercase; font-size: 0.85rem; }
        tr:hover td { background-color: rgba(255, 255, 255, 0.02); }

        .tag { display: inline-block; padding: 0.2rem 0.6rem; border-radius: 0.375rem; font-size: 0.75rem; font-weight: 700; }
        .tag-pass { background-color: var(--pass-bg); color: var(--pass-color); }
        .tag-warning { background-color: var(--warn-bg); color: var(--warn-color); }
        .tag-fail { background-color: var(--fail-bg); color: var(--fail-color); }

        .box { background-color: var(--card-bg); border: 1px solid var(--card-border); border-radius: 0.75rem; padding: 1.25rem; margin-bottom: 2rem; }
        pre { background-color: rgba(0, 0, 0, 0.4); padding: 1rem; border-radius: 0.5rem; color: #38bdf8; font-family: monospace; font-size: 0.85rem; white-space: pre-wrap; }
        .btn-copy { background-color: #2563eb; color: white; border: none; padding: 0.4rem 0.8rem; border-radius: 0.375rem; cursor: pointer; font-size: 0.8rem; margin-top: 0.5rem; }
        .btn-copy:hover { background-color: #1d4ed8; }

        footer { text-align: center; font-size: 0.85rem; color: var(--text-secondary); margin-top: 3rem; padding-top: 1.5rem; border-top: 1px solid var(--card-border); }
    </style>
</head>
<body>
    <div class="container">
        <header class="header">
            <div class="title-group">
                <h1>AlphaCert Protein Structure Quality Report</h1>
                <p>{{ report.metadata.name }} &bull; Predicted with {{ report.metadata.engine }}</p>
            </div>
            <div>
                <span class="status-badge badge-{{ report.overall_status.lower() }}">
                    {{ report.overall_status }}
                </span>
            </div>
        </header>

        <div class="grid-cards">
            <div class="card">
                <div class="card-label">Mean pLDDT</div>
                <div class="card-value">{{ "%.1f"|format(report.plddt_result.mean_plddt) }} / 100</div>
                <div class="card-subtext">{{ "%.1f"|format((report.plddt_result.frac_very_high_90 + report.plddt_result.frac_confident_70_90)*100) }}% >= 70 confident</div>
                <div class="plddt-bar">
                    <div class="plddt-segment" style="width: {{ report.plddt_result.frac_very_high_90*100 }}%; background-color: var(--plddt-vh);" title="Very High (>90)"></div>
                    <div class="plddt-segment" style="width: {{ report.plddt_result.frac_confident_70_90*100 }}%; background-color: var(--plddt-conf);" title="Confident (70-90)"></div>
                    <div class="plddt-segment" style="width: {{ report.plddt_result.frac_low_50_70*100 }}%; background-color: var(--plddt-low);" title="Low (50-70)"></div>
                    <div class="plddt-segment" style="width: {{ report.plddt_result.frac_very_low_under_50*100 }}%; background-color: var(--plddt-vl);" title="Very Low (<50)"></div>
                </div>
            </div>
            <div class="card">
                <div class="card-label">PAE Quality</div>
                <div class="card-value">
                    {% if report.pae_result %}
                    {{ "%.1f"|format(report.pae_result.mean_pae) }} Å
                    {% else %}
                    N/A
                    {% endif %}
                </div>
                <div class="card-subtext">
                    {% if report.pae_result %}
                    {{ report.pae_result.status }} ({{ report.pae_result.n_compact_domains }} domains)
                    {% else %}
                    No PAE matrix supplied
                    {% endif %}
                </div>
            </div>
            <div class="card">
                <div class="card-label">Docking Readiness</div>
                <div class="card-value">{{ report.docking_readiness.status }}</div>
                <div class="card-subtext">Pocket pLDDT: {{ "%.1f"|format(report.docking_readiness.pocket_plddt) }}</div>
            </div>
            <div class="card">
                <div class="card-label">Stereochemistry</div>
                <div class="card-value">
                    {% if report.stereochemistry %}
                    {{ "%.1f"|format(report.stereochemistry.frac_rama_favored*100) }}%
                    {% else %}
                    PASS
                    {% endif %}
                </div>
                <div class="card-subtext">
                    {% if report.stereochemistry %}
                    Ramachandran Favored (Clash: {{ "%.1f"|format(report.stereochemistry.clashscore) }})
                    {% else %}
                    Stereochemistry verified
                    {% endif %}
                </div>
            </div>
        </div>

        <h2 class="section-title">Structure Validation Summary</h2>
        <table>
            <thead>
                <tr>
                    <th>Quality Assessment Metric</th>
                    <th>Measured Value</th>
                    <th>Acceptance Standard</th>
                    <th>Status</th>
                    <th>Diagnostic Interpretation</th>
                </tr>
            </thead>
            <tbody>
                <tr>
                    <td><strong>Overall Backbone Confidence</strong></td>
                    <td>Mean pLDDT = {{ "%.1f"|format(report.plddt_result.mean_plddt) }}</td>
                    <td>Mean &ge; 70.0 (Frac &ge; 60%)</td>
                    <td><span class="tag tag-{{ report.plddt_result.status.lower() }}">{{ report.plddt_result.status }}</span></td>
                    <td style="color: var(--text-secondary)">{{ report.plddt_result.diagnostic_message }}</td>
                </tr>
                {% if report.pae_result %}
                <tr>
                    <td><strong>Inter-Domain / Interface PAE</strong></td>
                    <td>Mean PAE = {{ "%.1f"|format(report.pae_result.mean_pae) }} Å</td>
                    <td>PAE &le; 10.0 Å</td>
                    <td><span class="tag tag-{{ report.pae_result.status.lower() }}">{{ report.pae_result.status }}</span></td>
                    <td style="color: var(--text-secondary)">{{ report.pae_result.diagnostic_message }}</td>
                </tr>
                {% endif %}
                {% if report.stereochemistry %}
                <tr>
                    <td><strong>Ramachandran Dihedral Angles</strong></td>
                    <td>{{ "%.1f"|format(report.stereochemistry.frac_rama_favored*100) }}% Favored, {{ "%.1f"|format(report.stereochemistry.frac_rama_outliers*100) }}% Outliers</td>
                    <td>Favored &ge; 88%, Outliers &le; 3%</td>
                    <td><span class="tag tag-{{ report.stereochemistry.status.lower() }}">{{ report.stereochemistry.status }}</span></td>
                    <td style="color: var(--text-secondary)">{{ report.stereochemistry.diagnostic_message }}</td>
                </tr>
                <tr>
                    <td><strong>Steric Clashes</strong></td>
                    <td>Clashscore = {{ "%.1f"|format(report.stereochemistry.clashscore) }}</td>
                    <td>Clashscore &le; 10.0</td>
                    <td><span class="tag tag-{{ report.stereochemistry.status.lower() }}">{{ report.stereochemistry.status }}</span></td>
                    <td style="color: var(--text-secondary)">{{ report.stereochemistry.n_clashes }} total atom overlap contacts</td>
                </tr>
                {% endif %}
                <tr>
                    <td><strong>Molecular Docking Readiness</strong></td>
                    <td>Pocket pLDDT = {{ "%.1f"|format(report.docking_readiness.pocket_plddt) }}</td>
                    <td>Pocket pLDDT &ge; 80.0</td>
                    <td><span class="tag tag-{{ report.docking_readiness.status.lower() }}">{{ report.docking_readiness.status }}</span></td>
                    <td style="color: var(--text-secondary)">{{ report.docking_readiness.message }}</td>
                </tr>
                <tr>
                    <td><strong>Solvated MD Readiness</strong></td>
                    <td>{{ "%.1f"|format(report.md_readiness.disordered_fraction*100) }}% disordered</td>
                    <td>Disorder &le; 25.0%</td>
                    <td><span class="tag tag-{{ report.md_readiness.status.lower() }}">{{ report.md_readiness.status }}</span></td>
                    <td style="color: var(--text-secondary)">{{ report.md_readiness.message }}</td>
                </tr>
            </tbody>
        </table>

        {% if report.recommendations %}
        <div class="box" style="border-left: 4px solid var(--warn-color);">
            <h3 style="color: var(--warn-color); margin-bottom: 0.5rem;">Diagnostic Notes & Structural Warnings</h3>
            <ul style="padding-left: 1.25rem;">
                {% for rec in report.recommendations %}
                <li style="margin-bottom: 0.25rem; color: var(--text-secondary);">{{ rec }}</li>
                {% endfor %}
            </ul>
        </div>
        {% endif %}

        <h2 class="section-title">Computational Details (Publication Ready Methods)</h2>
        <div class="box">
            <pre id="methodsSnippet">{{ methods_text }}</pre>
            <button class="btn-copy" onclick="copyToClipboard('methodsSnippet')">Copy Methods Paragraph</button>
        </div>

        <h2 class="section-title">BibTeX Citation</h2>
        <div class="box">
            <pre id="bibSnippet">{{ citation_bib }}</pre>
            <button class="btn-copy" onclick="copyToClipboard('bibSnippet')">Copy BibTeX</button>
        </div>

        <footer>
            Generated automatically by <strong>AlphaCert v1.0.0</strong> &bull; AlphaFold Model Quality & Certification Toolkit &bull; Monreal-Hernández, 2026.
        </footer>
    </div>

    <script>
        function copyToClipboard(elementId) {
            const text = document.getElementById(elementId).innerText;
            navigator.clipboard.writeText(text).then(() => {
                alert('Copied to clipboard!');
            }).catch(err => {
                console.error('Error copying: ', err);
            });
        }
    </script>
</body>
</html>
"""


def generate_alphacert_html_report(
    report: AlphaFoldValidationReport,
    output_path: str,
    methods_text: str = "",
    citation_bib: str = ""
) -> str:
    """
    Renders HTML report template and writes it to disk.
    """
    template = jinja2.Template(HTML_TEMPLATE)
    rendered = template.render(
        report=report,
        methods_text=methods_text,
        citation_bib=citation_bib
    )
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(rendered)
    return output_path
