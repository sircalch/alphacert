"""
Parser for PDB and mmCIF structure files, extracting coordinates, atom types, and pLDDT B-factors.

AlphaFold, ColabFold and ESMFold write the per-residue pLDDT in the B-factor column. Only the
first model and the first alternate location are read; water molecules are ignored. Residues keep
the order of the file (so insertion codes are handled).
"""

from typing import Dict, Any, List, Tuple
import os
import shlex
import numpy as np

_WATER = {"HOH", "WAT", "DOD", "H2O"}


def _pdb_atoms(filepath: str):
    """Yields atom dicts from a PDB file (first MODEL only)."""
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("ENDMDL"):
                break
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            try:
                atom_name = line[12:16].strip()
                yield {
                    "atom": atom_name,
                    "altloc": line[16].strip(),
                    "res_name": line[17:20].strip(),
                    "chain": line[21].strip() or "A",
                    "res_seq": int(line[22:26]),
                    "icode": line[26].strip(),
                    "xyz": (float(line[30:38]), float(line[38:46]), float(line[46:54])),
                    "b": float(line[60:66]) if len(line) >= 66 and line[60:66].strip() else float("nan"),
                    "element": (line[76:78].strip() if len(line) >= 78 else "") or atom_name.lstrip("0123456789")[:1],
                }
            except ValueError:
                continue


def _cif_atoms(filepath: str):
    """Yields atom dicts from the _atom_site loop of an mmCIF file (first model only)."""
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.read().splitlines()
    i, n = 0, len(lines)
    while i < n:
        if lines[i].strip() == "loop_" and i + 1 < n and lines[i + 1].strip().startswith("_atom_site."):
            break
        i += 1
    if i >= n:
        return
    i += 1
    cols: List[str] = []
    while i < n and lines[i].strip().startswith("_atom_site."):
        cols.append(lines[i].strip().split(".", 1)[1].split()[0])
        i += 1
    idx = {c: k for k, c in enumerate(cols)}

    def col(row, *names, default=None):
        for nm in names:
            if nm in idx and row[idx[nm]] not in ("?", "."):
                return row[idx[nm]]
        return default

    first_model = None
    while i < n:
        s = lines[i].strip()
        i += 1
        if not s or s.startswith("#") or s.startswith("loop_") or s.startswith("_"):
            if s.startswith(("loop_", "_")):
                break
            continue
        row = shlex.split(s, posix=True) if ('"' in s or "'" in s) else s.split()
        if len(row) != len(cols):
            continue
        model = col(row, "pdbx_PDB_model_num", default="1")
        if first_model is None:
            first_model = model
        if model != first_model:
            break
        try:
            yield {
                "atom": col(row, "auth_atom_id", "label_atom_id"),
                "altloc": col(row, "label_alt_id", default=""),
                "res_name": col(row, "auth_comp_id", "label_comp_id"),
                "chain": col(row, "auth_asym_id", "label_asym_id", default="A"),
                "res_seq": int(col(row, "auth_seq_id", "label_seq_id")),
                "icode": col(row, "pdbx_PDB_ins_code", default=""),
                "xyz": (float(col(row, "Cartn_x")), float(col(row, "Cartn_y")), float(col(row, "Cartn_z"))),
                "b": float(col(row, "B_iso_or_equiv", default="nan")),
                "element": col(row, "type_symbol", default=""),
            }
        except (TypeError, ValueError):
            continue


def parse_protein_structure(filepath: str) -> Dict[str, Any]:
    """
    Parses a PDB or mmCIF protein structure predicted by AlphaFold / ColabFold / ESMFold.

    Parameters
    ----------
    filepath : str
        Path to .pdb/.ent or .cif/.mmcif file.

    Returns
    -------
    data : dict
        Coordinates, atom elements, per-residue pLDDT (CA B-factor), backbone atoms,
        chain lengths and chain identifiers.
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
    ext = os.path.splitext(filepath)[1].lower()
    atoms = _cif_atoms(filepath) if ext in (".cif", ".mmcif") else _pdb_atoms(filepath)

    coordinates: List[Tuple[float, float, float]] = []
    atom_elements: List[str] = []
    residues: Dict[Tuple[str, int, str], Dict[str, Any]] = {}   # insertion-ordered

    for a in atoms:
        if a["res_name"] in _WATER or a["altloc"] not in ("", "A", "1"):
            continue
        key = (a["chain"], a["res_seq"], a["icode"])
        res = residues.setdefault(key, {"res_name": a["res_name"], "plddt": a["b"], "atoms": {}})
        if a["atom"] in res["atoms"]:
            continue
        coordinates.append(a["xyz"])
        atom_elements.append(a["element"].upper())
        res["atoms"][a["atom"]] = list(a["xyz"])
        if a["atom"] == "CA":
            res["plddt"] = a["b"]

    if not residues:
        raise ValueError(f"No atoms could be read from {filepath}")

    per_residue_plddt = [r["plddt"] for r in residues.values()]
    backbone_atoms = []
    chains_seen: Dict[str, int] = {}
    for (ch, seq, icode), res in residues.items():
        chains_seen[ch] = chains_seen.get(ch, 0) + 1
        item = {"res_name": res["res_name"], "chain": ch, "seq": seq, "icode": icode}
        for a_type in ("N", "CA", "C", "O"):
            if a_type in res["atoms"]:
                item[a_type] = res["atoms"][a_type]
        backbone_atoms.append(item)

    return {
        "coordinates": np.asarray(coordinates, dtype=float),
        "atom_elements": atom_elements,
        "per_residue_plddt": per_residue_plddt,
        "backbone_atoms": backbone_atoms,
        "n_residues": len(per_residue_plddt),
        "chain_lengths": list(chains_seen.values()),
        "chains": list(chains_seen.keys())
    }
