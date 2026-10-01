#!/usr/bin/env python3
"""
build_browser.py  --  Turn a MOE index CSV into a standalone, searchable
chemical-structure browser (one self-contained .html file).

Filters to small molecules (n_res == 1, MW in [MW_MIN, MW_MAX]), renders each
SMILES to an embedded SVG, and writes an HTML page with:
  * live text search over compound name AND smiles string
  * substructure search: paste a SMARTS/SMILES query, matching atoms highlight
  * the source .moe file shown for every hit

Usage:
    python build_browser.py INPUT.csv [-o OUTPUT.html]
                            [--mw-min 200] [--mw-max 1000] [--nres 1]

Requires: rdkit, pandas   (pip install rdkit pandas)

v6

Updates:  
v6 checkbox & export functionality added including filter for common crystallization compds
v6 also added easier copy option for filenames and preset queries for common substructure searches
v5 stereo chem rendering added, and some cosmetic changes to the HTML page (larger header, larger subheader, larger lobster logo)
"""


import argparse
import base64
import html
import json
import sys
import os

import pandas as pd
from rdkit import Chem
from rdkit import DataStructs
from rdkit import RDLogger
from rdkit.Chem import Draw
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.ML.Cluster import Butina

DEFAULT_PRESETS = [
    ("amide", "C(=O)N"),
    ("nitrile", "C#N"),
    ("halogen", "[F,Cl,Br,I]"),
    ("acid", "C(=O)[OH]"),
    ("aromatic ring", "c1ccccc1"),
]


def load_presets(path):
    """Read 'label = SMARTS' lines. Splits on the FIRST '=' only, so SMARTS
    containing '=' survive; '#' only starts a comment at the start of a line,
    so 'C#C' is safe. Returns DEFAULT_PRESETS if the file is absent."""
    if not path or not os.path.exists(path):
        return list(DEFAULT_PRESETS)
    presets = []
    with open(path, encoding="utf-8") as f:
        for n, raw in enumerate(f, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if "=" not in line:
                print(f"  presets: line {n} has no '=', skipped: {line}")
                continue
            label, smarts = (p.strip() for p in line.split("=", 1))
            if not label or not smarts:
                continue
            if Chem.MolFromSmarts(smarts) is None:
                print(f"  presets: line {n} invalid SMARTS, skipped: {smarts}")
                continue
            presets.append((label, smarts))
    return presets


def presets_html(presets):
    rows = [
        '    <span class="presetlabel">quick queries:</span>',
    ]
    for label, smarts in presets:
        rows.append('    <button data-q="%s">%s</button>'
                    % (html.escape(smarts, quote=True), html.escape(label)))
    rows.append('    <button data-q="" class="clearq">clear</button>')
    return "\n".join(rows)

def cluster_families(mols_by_idx, cutoff=0.35):
    """Butina-cluster a dict {record_index: mol}. Returns:
       family_of: {record_index: family_id}, ordered so family 0 is largest.
       Also returns on-bits per record for JS-side Tanimoto."""
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    idxs = list(mols_by_idx.keys())
    fps = {i: gen.GetFingerprint(mols_by_idx[i]) for i in idxs}
    onbits = {i: list(fps[i].GetOnBits()) for i in idxs}

    n = len(idxs)
    family_of = {}
    if n == 0:
        return family_of, onbits
    if n == 1:
        family_of[idxs[0]] = 0
        return family_of, onbits

    # lower-triangle distance list for Butina
    fp_list = [fps[i] for i in idxs]
    dists = []
    for row in range(1, n):
        sims = DataStructs.BulkTanimotoSimilarity(fp_list[row], fp_list[:row])
        dists.extend(1.0 - s for s in sims)

    clusters = Butina.ClusterData(dists, n, cutoff, isDistData=True)
    # clusters are tuples of positions into idxs, already largest-first
    for fam_id, members in enumerate(clusters):
        for pos in members:
            family_of[idxs[pos]] = fam_id
    return family_of, onbits


def mol_to_svg(mol, size=(260, 200), highlight_atoms=None):
    """Render an RDKit mol to an inline SVG string (no external files).
    Preserves and annotates stereochemistry (wedge/dash bonds, R/S/E/Z labels).
    If highlight_atoms is given, those atoms are drawn highlighted."""
    # make sure stereo is perceived from whatever the SMILES encoded
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    d = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
    opts = d.drawOptions()
    opts.clearBackground = False
    opts.addStereoAnnotation = True   # show (R)/(S)/(E)/(Z) labels on the drawing
    if highlight_atoms:
        rdMolDraw2D.PrepareAndDrawMolecule(d, mol, highlightAtoms=highlight_atoms)
    else:
        rdMolDraw2D.PrepareAndDrawMolecule(d, mol)
    d.FinishDrawing()
    svg = d.GetDrawingText()
    return svg.replace("<?xml version='1.0' encoding='iso-8859-1'?>\n", "")


def stereo_count(mol):
    """Return (n_chiral_centers, n_stereo_double_bonds) actually defined on mol."""
    try:
        centers = Chem.FindMolChiralCenters(
            mol, includeUnassigned=False, useLegacyImplementation=False)
        n_chiral = len(centers)
    except Exception:
        n_chiral = 0
    n_ez = sum(1 for b in mol.GetBonds()
               if b.GetStereo() != Chem.BondStereo.STEREONONE)
    return n_chiral, n_ez


# Common crystallization / buffer additives to filter out.
# Matched two ways: canonical SMILES (exact structure) and name substrings.
_ADDITIVE_SMILES_RAW = {
    'HEPES':          'OCCN1CCN(CCS(=O)(=O)O)CC1',
    'MES':            'OCC[NH+]1CCOCC1',
    'Tris':           'OCC(N)(CO)CO',
    'glycerol':       'OCC(O)CO',
    'ethylene glycol':'OCCO',
    'PEG fragment':   'OCCOCCO',
    'sulfate':        'OS(=O)(=O)O',
    'phosphate':      'OP(=O)(O)O',
    'acetate':        'CC(=O)O',
    'formate':        'OC=O',
    'citrate':        'OC(=O)CC(O)(CC(=O)O)C(=O)O',
    'EDTA':           'OC(=O)CN(CCN(CC(=O)O)CC(=O)O)CC(=O)O',
    'DMSO':           'CS(C)=O',
    'imidazole':      'c1c[nH]cn1',
    'malonate':       'OC(=O)CC(=O)O',
    'tartrate':       'OC(=O)C(O)C(O)C(=O)O',
    'MPD':            'CC(O)CC(C)(C)O',
    'BME':            'OCCS',
    'DTT':            'OCC(S)C(S)CO',
    'glucose':        'OCC1OC(O)C(O)C(O)C1O',
    'benzoate':       'OC(=O)c1ccccc1',
}

# Name substrings (lowercased) that flag a chain as an additive regardless of structure.
_ADDITIVE_NAME_PATTERNS = [
    'hepes', 'mes', 'tris', 'glycerol', 'peg', 'sulfate', 'sulphate',
    'phosphate', 'acetate', 'formate', 'citrate', 'edta', 'dmso',
    'imidazole', 'malonate', 'tartrate', 'mpd', 'glucose', 'benzoate',
    'water', 'wat', 'hoh', 'buffer', 'cryo', 'ion',
]


def repair_moe_smiles(smi):
    """Second-chance parse for SMILES that RDKit rejects on valence grounds.

    MOE sometimes exports hypervalent sulfur/selenium groups (sulfoximines,
    sulfonimidoyls) in a charge-separated form such as [S@@+2]([OH0])(=N...),
    which RDKit reads as an over-valent atom and refuses to parse. This finds
    any S/Se carrying a >=+2 formal charge next to a bare (single-bonded,
    uncharged, H-free) oxygen, promotes that S-O bond to a double bond, and
    neutralizes the charge -- i.e. rewrites [S+2]-[O] as neutral S=O. Returns
    an RDKit mol on success, or None if it still can't be sanitized.
    """
    m = Chem.MolFromSmiles(smi, sanitize=False)
    if m is None:
        return None
    changed = False
    for atom in m.GetAtoms():
        if atom.GetSymbol() in ("S", "Se") and atom.GetFormalCharge() >= 2:
            for nbr in atom.GetNeighbors():
                if (nbr.GetSymbol() == "O" and nbr.GetTotalNumHs() == 0
                        and nbr.GetFormalCharge() == 0):
                    b = m.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                    if b is not None and b.GetBondType() == Chem.BondType.SINGLE:
                        b.SetBondType(Chem.BondType.DOUBLE)
                        changed = True
            atom.SetFormalCharge(0)
            changed = True
    if not changed:
        return None
    try:
        Chem.SanitizeMol(m)
        return m
    except Exception:
        return None


def file_date_of(row):
    """Return the row's file_date as a string, or '' if the column is absent
    or empty. Dates are written by the SVL indexer as 'YYYY-MM-DD HH:MM:SS',
    which sorts correctly as plain text -- no parsing needed."""
    try:
        v = row["file_date"]
    except (KeyError, IndexError):
        return ""
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return ""
    s = str(v).strip()
    return "" if s.lower() in ("nan", "none") else s


def mark_versions(records):
    """Identify duplicate copies of the same structure across files and mark
    which copy is the most recent.

    Molecules get copied forward into newer .moe files, so the same structure
    can appear many times. Records are grouped by canonical SMILES (a rename
    doesn't fool it, and a real structural change correctly splits the group).
    Within each group the newest file_date wins; exactly one record is marked
    is_latest=1 so a "latest only" view shows one card per structure.

    Returns the number of structures that exist in more than one copy.
    """
    groups = {}
    for i, r in enumerate(records):
        key = r["canonical"] or ("__unparsed__%d" % i)  # never merge bad mols
        groups.setdefault(key, []).append(i)

    n_dup = 0
    for key, idxs in groups.items():
        # newest first; blank dates sort last; stable tiebreak on original order
        idxs.sort(key=lambda i: (records[i]["date"] or "", -i), reverse=True)
        for rank, i in enumerate(idxs):
            records[i]["is_latest"] = 1 if rank == 0 else 0
            records[i]["n_versions"] = len(idxs)
        if len(idxs) > 1:
            n_dup += 1
    return n_dup


def build_additive_set():
    """Canonicalize the additive SMILES once; return a set of canonical strings."""
    s = set()
    for smi in _ADDITIVE_SMILES_RAW.values():
        m = Chem.MolFromSmiles(smi)
        if m is not None:
            s.add(Chem.MolToSmiles(m))
    return s


def is_additive(mol, canon_smiles, name, additive_set):
    """True if this compound looks like a crystallographic/buffer additive."""
    nm = (name or '').lower()
    for pat in _ADDITIVE_NAME_PATTERNS:
        if pat in nm:
            return True
    if canon_smiles and canon_smiles in additive_set:
        return True
    return False


def build(input_csv, output_html, mw_min, mw_max, nres, smarts=None,
          filter_additives=True, presets_file=None):
    # Silence RDKit's own parse/valence chatter -- expected first-attempt
    # failures are handled by repair_moe_smiles and reported in our own summary.
    RDLogger.DisableLog("rdApp.*")

    presets = load_presets(presets_file)
    print(f"  presets: {len(presets)} quick-query button(s)")

    df = pd.read_csv(input_csv)

    required = {"file", "compound", "n_res", "n_atoms", "mw", "smiles"}
    missing = required - set(df.columns)
    if missing:
        sys.exit(f"ERROR: CSV missing columns: {', '.join(sorted(missing))}")

    # --- Optional SMARTS query for substructure filter (done in Python/RDKit) ---
    query = None
    if smarts:
        query = Chem.MolFromSmarts(smarts)
        if query is None:
            sys.exit(f"ERROR: could not parse SMARTS query: {smarts}")

    # --- Filter to small molecules ---
    df["mw"] = pd.to_numeric(df["mw"], errors="coerce")
    df["n_res"] = pd.to_numeric(df["n_res"], errors="coerce")

    n_before = len(df)
    mask = (
        (df["n_res"] == nres)
        & (df["mw"] >= mw_min)
        & (df["mw"] <= mw_max)
        & df["smiles"].notna()
        & (df["smiles"].astype(str).str.len() > 0)
    )
    df = df[mask].copy().reset_index(drop=True)
    print(f"Filtered {n_before} -> {len(df)} rows "
          f"(n_res=={nres}, {mw_min} <= MW <= {mw_max}, non-empty SMILES)")

    # --- Render each molecule; keep canonical SMILES for JS substructure search ---
    records = []
    mols_by_idx = {}   # record index -> RDKit mol, for clustering
    n_bad = 0
    n_repaired = 0   # SMILES rescued by repair_moe_smiles
    n_smarts_kept = 0
    n_with_stereo = 0   # how many compounds carry defined stereochemistry
    n_additives = 0     # how many buffer/crystallography additives filtered out

    additive_set = build_additive_set() if filter_additives else set()

    for _, row in df.iterrows():
        smi = str(row["smiles"]).strip()
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            # second-chance: repair MOE's charge-separated hypervalent S/Se
            mol = repair_moe_smiles(smi)
            if mol is not None:
                n_repaired += 1
        name = str(row["compound"])

        # Drop common crystallization / buffer additives (HEPES, glycerol, etc.)
        if filter_additives:
            canon_for_test = Chem.MolToSmiles(mol) if mol is not None else ""
            if is_additive(mol, canon_for_test, name, additive_set):
                n_additives += 1
                continue

        # SMARTS substructure filter (skip non-matching compounds)
        hl = None
        if query is not None:
            if mol is None:
                continue  # can't match a broken mol
            match = mol.GetSubstructMatch(query)
            if not match:
                continue  # doesn't contain the query motif
            hl = list(match)
            n_smarts_kept += 1

        if mol is None:
            n_bad += 1
            svg = ("<div class='badmol'>SMILES failed to parse:<br>"
                   + html.escape(smi) + "</div>")
            canon = ""
            n_chiral, n_ez = 0, 0
        else:
            svg = mol_to_svg(mol, highlight_atoms=hl)
            canon = Chem.MolToSmiles(mol)   # isomeric by default: keeps stereo
            n_chiral, n_ez = stereo_count(mol)
            if n_chiral or n_ez:
                n_with_stereo += 1
        rec = {
            "file": str(row["file"]),
            "date": file_date_of(row),   # 'YYYY-MM-DD HH:MM:SS' or '' if absent
            "compound": name,
            "n_atoms": int(row["n_atoms"]) if pd.notna(row["n_atoms"]) else "",
            "mw": round(float(row["mw"]), 1),
            "smiles": smi,
            "canonical": canon,
            "svg": svg,
            "stereo": f"{n_chiral}R/S, {n_ez}E/Z" if (n_chiral or n_ez) else "",
            "family": -1,      # filled in after clustering
            "onbits": [],      # Morgan on-bits, for JS-side Tanimoto
            "is_latest": 1,    # filled in by mark_versions()
            "n_versions": 1,   # how many copies of this structure exist overall
        }
        if mol is not None:
            mols_by_idx[len(records)] = mol
        records.append(rec)

    if filter_additives and n_additives:
        print(f"  filtered out {n_additives} buffer/crystallography additive(s)")

    # --- Cluster into chemical families and attach fingerprints ---
    family_of, onbits = cluster_families(mols_by_idx)
    for i, fam in family_of.items():
        records[i]["family"] = fam
    for i, bits in onbits.items():
        records[i]["onbits"] = bits
    n_families = (max(family_of.values()) + 1) if family_of else 0
    print(f"  clustered into {n_families} chemical famil"
          f"{'y' if n_families == 1 else 'ies'}")

    # --- Identify copies of the same structure across files; mark the newest ---
    n_dup = mark_versions(records)
    n_dated = sum(1 for r in records if r["date"])
    has_dates = n_dated > 0
    if has_dates:
        print(f"  file dates: {n_dated}/{len(records)} rows dated")
    else:
        print("  file dates: none found (no 'file_date' column) -- "
              "date sort and latest-version view disabled")
    if n_dup:
        n_hidden = sum(1 for r in records if not r["is_latest"])
        print(f"  {n_dup} structure(s) appear in multiple copies "
              f"({n_hidden} older copies can be hidden via 'latest only')")

    # Order records so families are contiguous (family 0 = largest first).
    # Unparseable mols (family -1) sink to the end.
    records.sort(key=lambda r: (r["family"] if r["family"] >= 0 else 10**9,
                                r["compound"]))

    if query is not None:
        print(f"  SMARTS '{smarts}' matched {n_smarts_kept} compounds")
    if n_repaired:
        print(f"  repaired {n_repaired} SMILES (charge-separated hypervalent S/Se)")
    if n_bad:
        print(f"  ({n_bad} SMILES could not be parsed and are flagged in the page)")
    n_total = len(records)
    print(f"  stereochemistry: {n_with_stereo}/{n_total} compounds carry defined "
          f"stereo (wedge/dash or E/Z)")
    if n_total and n_with_stereo == 0:
        print("  NOTE: no stereo found in any SMILES. If your compounds ARE chiral, "
              "MOE likely exported flat SMILES -- check the sm_ExtractUnique call.")

    smarts_note = ""
    if smarts:
        smarts_note = f" &middot; substructure: {html.escape(smarts)}"

    data_json = json.dumps(records)
    page = HTML_TEMPLATE.replace("/*DATA*/", data_json) \
                        .replace("/*PRESETS*/", presets_html(presets)) \
                        .replace("/*HASDATES*/", "true" if has_dates else "false") \
                        .replace("{{COUNT}}", str(len(records))) \
                        .replace("{{MWMIN}}", str(mw_min)) \
                        .replace("{{MWMAX}}", str(mw_max)) \
                        .replace("{{NRES}}", str(nres)) \
                        .replace("{{SMARTSNOTE}}", smarts_note)

    with open(output_html, "w", encoding="utf-8") as f:
        f.write(page)
    print(f"Wrote {output_html}  ({len(page)//1024} KB)")


HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>MOE Compound Browser</title>
<script src="RDKit_minimal.js"></script>
<style>
  :root {
    --bg: #0f1115; --panel: #171a21; --edge: #262b36;
    --ink: #e6e9ef; --muted: #9aa4b2; --accent: #6ea8fe; --accent2: #7ee0c0;
    --bad: #ff6b6b;
  }
  * { box-sizing: border-box; }
  body {
    margin: 0; background: var(--bg); color: var(--ink);
    font: 15px/1.45 -apple-system, Segoe UI, Roboto, sans-serif;
  }
  header {
    position: sticky; top: 0; z-index: 10; background: var(--panel);
    border-bottom: 1px solid var(--edge); padding: 14px 20px;
  }
  h1 { margin: 0 0 4px; font-size: 40px; font-weight: 700; }
  .sub { color: var(--muted); font-size: 13px; }
  .controls { display: flex; flex-wrap: wrap; gap: 10px; margin-top: 12px; align-items: center; }
  input[type=text] {
    background: var(--bg); border: 1px solid var(--edge); color: var(--ink);
    padding: 8px 11px; border-radius: 8px; font-size: 14px; min-width: 240px;
  }
  input[type=text]:focus { outline: none; border-color: var(--accent); }
  .hint { color: var(--muted); font-size: 12px; }
  .presets { display:flex; align-items:center; gap:6px; flex-wrap:wrap; margin-top:10px; }
  .presetlabel { color: var(--muted); font-size:12px; margin-right:2px; }
  .presets button { border:1px solid var(--edge); background:var(--bg); color:var(--ink);
                    border-radius:14px; padding:3px 11px; font-size:12px; cursor:pointer; }
  .presets button:hover { border-color: var(--accent); color: var(--accent); }
  .presets button.clearq { color: var(--muted); }
  .presets button.clearq:hover { border-color: var(--muted); color: var(--ink); }
  #status { color: var(--accent2); font-size: 13px; margin-left: auto; }
  .grid {
    display: grid; grid-template-columns: repeat(auto-fill, minmax(280px, 1fr));
    gap: 14px; padding: 18px 20px;
  }
  .card {
    background: var(--panel); border: 1px solid var(--edge); border-radius: 12px;
    padding: 12px; display: flex; flex-direction: column; gap: 8px;
  }
  .card svg { background: #fff; border-radius: 8px; width: 100%; height: auto; }
  .badmol { background:#2a1414; color:var(--bad); padding:20px; border-radius:8px;
            font-size:12px; text-align:center; word-break:break-all; }
  .name { font-weight: 600; font-size: 14px; word-break: break-all; }
  .meta { color: var(--muted); font-size: 12px; }
  .file { color: var(--accent); font-size: 12px; word-break: break-all; }
  .smi  { color: var(--muted); font-size: 11px; word-break: break-all;
          font-family: ui-monospace, Consolas, monospace; }
  .copy { cursor: pointer; border: 1px solid var(--edge); background: var(--bg);
          color: var(--muted); border-radius: 6px; font-size: 11px; padding: 2px 6px; }
  .copy:hover { color: var(--ink); border-color: var(--accent); }
  mark { background: var(--accent); color: #06121f; border-radius: 3px; padding: 0 2px; }
  .empty { padding: 40px; text-align: center; color: var(--muted); }
  .card.clickable { cursor: pointer; }
  .card.clickable:hover { border-color: var(--accent); }
  .card.refcard { border-color: var(--accent2); box-shadow: 0 0 0 1px var(--accent2); }
  .fam { display:inline-block; font-size:11px; padding:1px 7px; border-radius:10px;
         color:#06121f; font-weight:600; }
  .sim { color: var(--accent2); font-size: 12px; font-weight:600; }
  .stereo { color: #ffd166; font-weight: 600; }
  .date { color: var(--muted); font-size: 11px; font-family: ui-monospace, Consolas, monospace; }
  .datelabel { color: var(--muted); opacity: 0.75; }
  .older { display:inline-block; font-size:10px; padding:1px 6px; border-radius:9px;
           background:#4a3a12; color:#ffd166; font-weight:600; margin-left:4px; }
  .newest { display:inline-block; font-size:10px; padding:1px 6px; border-radius:9px;
            background:#123c30; color:var(--accent2); font-weight:600; margin-left:4px; }
  .card.stale { opacity: 0.62; }
  .card.stale:hover { opacity: 1; }
  .datectl { display:flex; align-items:center; gap:10px; }
  .datectl button.active { border-color: var(--accent2); color: var(--accent2); }
  .chk { display:flex; align-items:center; gap:5px; color:var(--muted);
         font-size:12px; cursor:pointer; }
  .chk input { accent-color: var(--accent2); cursor:pointer; }
  .banner { display:flex; align-items:center; gap:12px; flex-wrap:wrap;
            padding:8px 20px; background:var(--bg); border-bottom:1px solid var(--edge);
            color:var(--muted); font-size:13px; }
  .banner button { border:1px solid var(--edge); background:var(--panel); color:var(--ink);
                   border-radius:8px; padding:5px 12px; font-size:13px; cursor:pointer; }
  .banner button:hover { border-color: var(--accent); }
  .selbar { display:flex; align-items:center; gap:10px; flex-wrap:wrap;
            padding:8px 20px; background:var(--panel); border-bottom:1px solid var(--edge);
            font-size:13px; }
  .selbar button { border:1px solid var(--edge); background:var(--bg); color:var(--ink);
                   border-radius:8px; padding:5px 12px; font-size:13px; cursor:pointer; }
  .selbar button:hover { border-color: var(--accent); }
  .selbar button.exp { background:var(--accent2); color:#06121f; font-weight:600; border-color:var(--accent2); }
  .selbar button.exp:hover { filter:brightness(1.1); }
  #selcount { color: var(--accent2); font-weight:600; }
  .card { position: relative; }
  .selbox { position:absolute; top:10px; left:10px; width:20px; height:20px;
            cursor:pointer; z-index:3; accent-color: var(--accent2); }
  .card.selected { border-color: var(--accent2); box-shadow: 0 0 0 1px var(--accent2); }
</style>
</head>
<body>
<header>
  <div style="display:flex; align-items:center; gap:12px;">
    <svg width="120" height="135" viewBox="0 0 72 78" xmlns="http://www.w3.org/2000/svg" aria-label="MOE-CLAW cartoon lobster logo" style="flex-shrink:0;">
      <g fill="#ff5a4d" stroke="#7a1a10" stroke-width="2.4" stroke-linejoin="round" stroke-linecap="round">
        <path d="M30 16 Q30 27 34 36" fill="none"/>
        <path d="M42 16 Q42 27 38 36" fill="none"/>
        <path d="M36 34 Q47 34 47 46 Q47 52 36 52 Q25 52 25 46 Q25 34 36 34 Z"/>
        <circle cx="32" cy="40" r="2.4" fill="#7a1a10" stroke="none"/>
        <circle cx="40" cy="40" r="2.4" fill="#7a1a10" stroke="none"/>
        <path d="M29 51 Q36 55 43 51 Q42 57 41 58 Q36 61 31 58 Q30 57 29 51 Z"/>
        <path d="M31 58 Q36 61 41 58 Q40 63 39 64 Q36 66 33 64 Q32 63 31 58 Z"/>
        <path d="M33 64 Q36 66 39 64 Q41 71 36 73 Q31 71 33 64 Z"/>
        <path d="M27 48 L21 51" fill="none"/>
        <path d="M28 51 L23 55" fill="none"/>
        <path d="M45 48 L51 51" fill="none"/>
        <path d="M44 51 L49 55" fill="none"/>
        <path d="M26 43 Q19 44 15 39" fill="none"/>
        <path d="M15 39 Q8 36 8 29 Q8 25 12 25 Q16 25 17 30 L17 34 Q16 39 15 39 Z"/>
        <path d="M12 25 Q11 20 14 19 Q16 21 15 25 Z"/>
        <path d="M17 30 Q20 27 19 24 Q16 25 15 28 Z" fill="#ff6f5e"/>
        <path d="M46 43 Q53 44 57 39" fill="none"/>
        <path d="M57 39 Q64 36 64 29 Q64 25 60 25 Q56 25 55 30 L55 34 Q56 39 57 39 Z"/>
        <path d="M60 25 Q61 20 58 19 Q56 21 57 25 Z"/>
        <path d="M55 30 Q52 27 53 24 Q56 25 57 28 Z" fill="#ff6f5e"/>
      </g>
    </svg>
    <div>
      <h1 style="margin:0;">MOE&#8288;-&#8288;<span style="color:#ff5a4d;">CLAW</span></h1>
      <div class="sub">{{COUNT}} small molecules &middot; filtered n_res=={{NRES}}, {{MWMIN}}&ndash;{{MWMAX}} MW{{SMARTSNOTE}}</div>
    </div>
  </div>
  <div class="controls">
    <input id="q" type="text" placeholder="Search name or SMILES text…" autocomplete="off">
    <input id="sub" type="text" placeholder="Substructure query e.g. C#C or C=CC(=O)N" autocomplete="off">
    <span class="hint" id="subhint">substructure: loading…</span>
    <span id="status"></span>
  </div>
  <div class="presets" id="presets">
/*PRESETS*/
  </div>
</header>
<div class="banner" id="banner">
  <span id="sortlabel">Sorted by chemical family. Click any compound to sort by similarity to it.</span>
  <button id="resetsort" style="display:none;">Back to family view</button>
  <span style="flex:1;"></span>
  <span class="datectl" id="datectl">
    <button id="sortdate">Newest first</button>
    <label class="chk"><input type="checkbox" id="latestonly"> latest version only</label>
  </span>
</div>
<div class="selbar" id="selbar">
  <button id="selall">Select all shown</button>
  <button id="selnone">Clear selection</button>
  <span id="selcount">0 selected</span>
  <span style="flex:1;"></span>
  <button id="expcsv" class="exp">Export CSV</button>
  <button id="expsdf" class="exp">Export SDF</button>
</div>
<div class="grid" id="grid"></div>

<script>
const DATA = /*DATA*/;
const HAS_DATES = /*HASDATES*/;
let rdkit = null;

// Date-driven view state
let dateSort = false;    // true = newest file first
let latestOnly = false;  // true = one card per structure (the most recent copy)

// Load local RDKit WASM. locateFile points the loader at the .wasm sitting
// next to this HTML file, so no CDN and no internet are needed.
if (typeof initRDKitModule === "function") {
  initRDKitModule({ locateFile: (f) => f })   // f === 'RDKit_minimal.wasm', same folder
    .then(m => { rdkit = m;
                 document.getElementById('subhint').textContent = 'substructure: ready'; })
    .catch(e => { document.getElementById('subhint').textContent =
                    'substructure: failed to load (' + e + ')'; });
} else {
  document.getElementById('subhint').textContent =
    'substructure: RDKit_minimal.js not found — keep it beside this HTML';
}

const grid = document.getElementById('grid');
const qEl = document.getElementById('q');
const subEl = document.getElementById('sub');
const statusEl = document.getElementById('status');
const bannerEl = document.getElementById('banner');
const sortLabelEl = document.getElementById('sortlabel');
const resetBtn = document.getElementById('resetsort');
const selAllBtn = document.getElementById('selall');
const selNoneBtn = document.getElementById('selnone');
const selCountEl = document.getElementById('selcount');
const expCsvBtn = document.getElementById('expcsv');
const expSdfBtn = document.getElementById('expsdf');
const sortDateBtn = document.getElementById('sortdate');
const latestOnlyEl = document.getElementById('latestonly');
const dateCtlEl = document.getElementById('datectl');

// If the CSV carried no file_date column, hide the date controls entirely
// rather than offering buttons that can't do anything.
if (!HAS_DATES) dateCtlEl.style.display = 'none';

// Selection is keyed by each record's index into DATA, so it survives
// filtering, searching, and re-sorting.
const selected = new Set();
function updateSelCount(){ selCountEl.textContent = selected.size + ' selected'; }

// Distinct-ish palette for family badges (cycles if more families than colors)
const FAM_COLORS = ['#7ee0c0','#6ea8fe','#ffd166','#ff9d92','#c39bd3','#f7a072',
                    '#90be6d','#f4978e','#8ecae6','#e0aaff','#bde0fe','#ffc8dd'];
function famColor(f){ return f < 0 ? '#555' : FAM_COLORS[f % FAM_COLORS.length]; }

// Tanimoto over sorted on-bit index lists (Jaccard of two sets)
function tanimoto(a, b){
  if(!a.length || !b.length) return 0;
  let i=0, j=0, inter=0;
  while(i < a.length && j < b.length){
    if(a[i] === b[j]){ inter++; i++; j++; }
    else if(a[i] < b[j]) i++;
    else j++;
  }
  return inter / (a.length + b.length - inter);
}

// Pre-sort each record's on-bits once so tanimoto's merge works
DATA.forEach(r => { if(r.onbits) r.onbits.sort((x,y)=>x-y); });

let refIndex = null;   // when set, sort by similarity to DATA[refIndex]

function esc(s){ return s.replace(/[&<>"]/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c])); }

function highlight(text, term){
  if(!term) return esc(text);
  const i = text.toLowerCase().indexOf(term.toLowerCase());
  if(i < 0) return esc(text);
  return esc(text.slice(0,i)) + '<mark>' + esc(text.slice(i,i+term.length)) + '</mark>' + esc(text.slice(i+term.length));
}

function subMatch(rec, qmol){
  if(!rec.canonical) return false;
  let m = null;
  try {
    m = rdkit.get_mol(rec.canonical);
    if(!m || !m.is_valid()) return false;
    const res = JSON.parse(m.get_substruct_match(qmol));
    return res && res.atoms && res.atoms.length > 0;
  } catch(e){ return false; }
  finally { if(m) m.delete(); }
}

function render(){
  const q = qEl.value.trim();
  const subq = subEl.value.trim();

  let qmol = null;
  if(subq && rdkit){
    try {
      qmol = rdkit.get_qmol(subq);
      if(!qmol || !qmol.is_valid()){ if(qmol) qmol.delete(); qmol = null; }
    } catch(e){ qmol = null; }
  }

  // attach original index so we can reference DATA after filtering
  let rows = DATA.map((rec, idx) => ({rec, idx})).filter(({rec}) => {
    if(latestOnly && !rec.is_latest) return false;
    if(q){
      const hay = (rec.compound + ' ' + rec.smiles).toLowerCase();
      if(!hay.includes(q.toLowerCase())) return false;
    }
    if(qmol && !subMatch(rec, qmol)) return false;
    return true;
  });

  if(qmol) qmol.delete();

  // Sort. Similarity (when a reference is picked) wins over date, which wins
  // over the default family order baked in by the builder.
  let refBits = null;
  if(refIndex !== null){
    refBits = DATA[refIndex].onbits;
    rows.forEach(o => { o.sim = tanimoto(refBits, o.rec.onbits); });
    rows.sort((a,b) => b.sim - a.sim);
    sortLabelEl.innerHTML = 'Sorted by similarity to <span class="sim">' +
      esc(DATA[refIndex].compound) + '</span> (Tanimoto, Morgan r2).';
    resetBtn.style.display = '';
  } else if(dateSort){
    // ISO 'YYYY-MM-DD HH:MM:SS' strings sort correctly as plain text.
    // Undated rows sink to the bottom.
    rows.sort((a,b) => {
      const da = a.rec.date || '', db = b.rec.date || '';
      if(da === db) return a.rec.compound.localeCompare(b.rec.compound);
      if(!da) return 1;
      if(!db) return -1;
      return db.localeCompare(da);           // newest first
    });
    sortLabelEl.textContent =
      'Sorted by file date, newest first. Click any compound to sort by similarity to it.';
    resetBtn.style.display = 'none';
  } else {
    sortLabelEl.textContent =
      'Sorted by chemical family. Click any compound to sort by similarity to it.';
    resetBtn.style.display = 'none';
  }

  sortDateBtn.classList.toggle('active', dateSort && refIndex === null);

  statusEl.textContent = rows.length + ' / ' + DATA.length + ' shown';

  if(rows.length === 0){
    grid.innerHTML = '<div class="empty">No compounds match.</div>';
    return;
  }

  grid.innerHTML = rows.map(({rec, idx, sim}) => {
    const isRef = (idx === refIndex);
    const isSel = selected.has(idx);
    const simLine = (refIndex !== null)
      ? `<div class="meta">similarity <span class="sim">${sim.toFixed(3)}</span></div>`
      : '';
    const famBadge = (rec.family >= 0)
      ? `<span class="fam" style="background:${famColor(rec.family)};">family ${rec.family}</span>`
      : `<span class="fam" style="background:#555;color:#ddd;">unparsed</span>`;
    // Version badge: only meaningful when the same structure exists in >1 copy
    let verBadge = '';
    if(HAS_DATES && rec.n_versions > 1){
      verBadge = rec.is_latest
        ? `<span class="newest" title="most recent of ${rec.n_versions} copies">newest of ${rec.n_versions}</span>`
        : `<span class="older" title="a newer copy of this structure exists">older</span>`;
    }
    const dateLine = rec.date
      ? `<div class="date"><span class="datelabel">file date:</span> ${esc(rec.date)}${verBadge}</div>`
      : '';
    const isStale = HAS_DATES && rec.n_versions > 1 && !rec.is_latest;
    return `
    <div class="card clickable${isRef ? ' refcard' : ''}${isSel ? ' selected' : ''}${isStale ? ' stale' : ''}" data-idx="${idx}">
      <input type="checkbox" class="selbox" data-idx="${idx}" ${isSel ? 'checked' : ''} title="select for export">
      ${rec.svg}
      <div class="name">${highlight(rec.compound, q)} ${famBadge}</div>
      <div class="meta">MW ${rec.mw} &middot; ${rec.n_atoms} atoms${rec.stereo ? ' &middot; <span class="stereo">' + rec.stereo + '</span>' : ''}</div>
      ${simLine}
      <div class="file">${esc(rec.file)}
        <button class="copy" data-copy="${esc(rec.file)}">copy</button>
      </div>
      ${dateLine}
      <div class="smi">${highlight(rec.smiles, q)}
        <button class="copy" data-copy="${esc(rec.smiles)}">copy</button>
      </div>
    </div>`;
  }).join('');

  // checkbox toggles selection (stop it bubbling to the card's similarity-sort click)
  grid.querySelectorAll('.selbox').forEach(box => {
    box.onclick = (e) => {
      e.stopPropagation();
      const i = parseInt(box.dataset.idx, 10);
      if(box.checked) selected.add(i); else selected.delete(i);
      box.closest('.card').classList.toggle('selected', box.checked);
      updateSelCount();
    };
  });

  // click a card -> sort by similarity to it (ignore clicks on copy button or checkbox)
  grid.querySelectorAll('.card').forEach(card => {
    card.onclick = (e) => {
      if(e.target.classList.contains('copy')) return;
      if(e.target.classList.contains('selbox')) return;
      refIndex = parseInt(card.dataset.idx, 10);
      render();
      window.scrollTo({top:0, behavior:'smooth'});
    };
  });
  grid.querySelectorAll('.copy').forEach(b => {
    b.onclick = (e) => { e.stopPropagation();
                        navigator.clipboard.writeText(b.dataset.copy);
                        b.textContent = 'copied'; setTimeout(()=>b.textContent='copy',900); };
  });

  // remember what's currently shown, for "select all shown"
  shownIndices = rows.map(o => o.idx);
  updateSelCount();
}

let shownIndices = [];   // indices currently visible after filter/search

resetBtn.onclick = () => { refIndex = null; render(); };

sortDateBtn.onclick = () => {
  dateSort = !dateSort;
  refIndex = null;      // leaving similarity mode so the date order is visible
  render();
};
latestOnlyEl.onchange = () => { latestOnly = latestOnlyEl.checked; render(); };

selAllBtn.onclick = () => {
  shownIndices.forEach(i => selected.add(i));
  render();
};
selNoneBtn.onclick = () => {
  selected.clear();
  render();
};

// ---- Export helpers ----
function triggerDownload(filename, text, mime){
  const blob = new Blob([text], {type: mime});
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url; a.download = filename;
  document.body.appendChild(a); a.click();
  document.body.removeChild(a);
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

function selectedRecords(){
  // preserve DATA order for a stable export
  const out = [];
  DATA.forEach((rec, i) => { if(selected.has(i)) out.push(rec); });
  return out;
}

function csvField(s){
  s = String(s == null ? '' : s);
  return '"' + s.replace(/"/g, '""') + '"';   // CSV-quote, escaping internal quotes
}

function exportCSV(){
  const recs = selectedRecords();
  if(recs.length === 0){ alert('Nothing selected to export.'); return; }
  const header = ['compound','file','file_date','mw','n_atoms','stereo','family','smiles'];
  const lines = [header.join(',')];
  recs.forEach(r => {
    lines.push([
      csvField(r.compound), csvField(r.file), csvField(r.date), csvField(r.mw),
      csvField(r.n_atoms), csvField(r.stereo),
      csvField(r.family), csvField(r.canonical || r.smiles)
    ].join(','));
  });
  triggerDownload('moeclaw_selection.csv', lines.join('\n'), 'text/csv');
}

function exportSDF(){
  const recs = selectedRecords();
  if(recs.length === 0){ alert('Nothing selected to export.'); return; }
  // SMILES-only SDF: empty atom/bond block, SMILES + metadata as data fields.
  // The isomeric canonical SMILES preserves stereochemistry.
  const blocks = recs.map(r => {
    const smi = r.canonical || r.smiles;
    const name = r.compound || '';
    // minimal V2000 header with zero atoms/bonds
    let b = name + '\n';
    b += '  MOE-CLAW\n';
    b += '\n';
    b += '  0  0  0  0  0  0  0  0  0  0999 V2000\n';
    b += 'M  END\n';
    b += '>  <SMILES>\n' + smi + '\n\n';
    b += '>  <compound>\n' + name + '\n\n';
    b += '>  <file>\n' + (r.file || '') + '\n\n';
    b += '>  <file_date>\n' + (r.date || '') + '\n\n';
    b += '>  <MW>\n' + r.mw + '\n\n';
    b += '>  <stereo>\n' + (r.stereo || '') + '\n\n';
    b += '>  <family>\n' + r.family + '\n\n';
    b += '$$$$\n';
    return b;
  });
  triggerDownload('moeclaw_selection.sdf', blocks.join(''), 'chemical/x-mdl-sdfile');
}

expCsvBtn.onclick = exportCSV;
expSdfBtn.onclick = exportSDF;

let t;
function debounced(){ clearTimeout(t); t = setTimeout(render, 120); }
qEl.addEventListener('input', debounced);
subEl.addEventListener('input', debounced);

// Preset quick-query chips: populate the substructure box and search immediately
document.querySelectorAll('#presets button').forEach(btn => {
  btn.onclick = () => {
    subEl.value = btn.dataset.q;   // empty string for the "clear" chip
    render();
    subEl.focus();
  };
});

render();
</script>
</body>
</html>
"""


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Build a searchable HTML structure browser from a MOE index CSV.")
    ap.add_argument("input_csv", help="Path to the MOE index CSV")
    ap.add_argument("-o", "--output", default="moe_browser.html", help="Output HTML file")
    ap.add_argument("--mw-min", type=float, default=200.0)
    ap.add_argument("--mw-max", type=float, default=1000.0)
    ap.add_argument("--nres", type=int, default=1)
    ap.add_argument("--smarts", default=None,
                    help="Optional SMARTS query; only compounds containing it are shown, "
                         "with matching atoms highlighted (e.g. 'C#C' or 'C=CC(=O)N')")
    ap.add_argument("--keep-additives", action="store_true",
                    help="Do NOT filter out common buffer/crystallography additives "
                         "(HEPES, glycerol, sulfate, etc.). By default they are removed.")
    ap.add_argument("--presets",
                    default=os.path.join(
                        os.path.dirname(os.path.abspath(__file__)), "presets.txt"),
                    help="File of 'label = SMARTS' quick-query presets "
                         "(default: presets.txt beside this script; "
                         "built-in generics if absent)")
    args = ap.parse_args()
    build(args.input_csv, args.output, args.mw_min, args.mw_max, args.nres,
          args.smarts, filter_additives=not args.keep_additives, presets_file=args.presets)