# MOE-CLAW

A tool for indexing compounds out of MOE `.moe` files and browsing them as an
interactive, searchable, offline chemical-structure gallery.

MOE-CLAW has two halves:

1. **An SVL indexer** (`moe_claw.svl`) that runs inside MOE (headless via
   `moebatch`), walks a folder of `.moe` files, and writes a CSV of every
   compound it finds — name, residue/atom counts, molecular weight, and an
   **isomeric SMILES** (stereochemistry perceived from the 3D pose).

2. **A browser builder** (`MOEClaw_buildbrowser.py`) that turns that CSV into a
   single self-contained HTML page. The page renders every structure and lets
   you search by name, search by substructure, sort by chemical similarity,
   filter out buffer/crystallography additives, select compounds, and export
   the selection as CSV or SDF.

The page runs entirely offline in a browser once built.

---

## What you need

- **MOE** (Molecular Operating Environment) — for the SVL indexing step only.
  MOE is licensed software from CCG and is **not** included or installed here.
  If you already have a CSV of compounds (columns: `file, compound, n_res,
  n_atoms, mw, smiles`), you can skip the MOE half entirely and use the browser
  builder on its own.
- **Python 3.10+** with **RDKit** and **pandas** — installed below.
- A **web browser** (Chrome, Edge, or Firefox).

---

## Setup

### 1. Clone the repo

```
git clone https://github.com/<your-username>/MOE_Claw.git
cd MOE_Claw
```

### 2. Create the Python environment

**With conda (recommended — RDKit installs most reliably this way):**

```
conda env create -f environment.yml
conda activate moeclaw
```

**Or with pip:**

```
python -m venv .venv
# Windows:  .venv\Scripts\activate
# mac/Linux: source .venv/bin/activate
pip install -r requirements.txt
```

### 3. Fetch the RDKit-JS runtime

The in-page substructure search needs two RDKit-JS files
(`RDKit_minimal.js` + `RDKit_minimal.wasm`, ~7 MB). They are included in
git but may need updating depending on your local environment. You can download them with:

```
python fetch_rdkit_js.py
```

This places both files in the repo folder. You only need to run it once. (The
launcher scripts below also do this automatically on first run.)

---

## Usage

### Step 1 — Index your `.moe` files (MOE required)

Edit `run_index.bat` (or run `moebatch` directly) to point at your MOE install
and your folder of `.moe` files. The indexer writes `moe_index.csv`:

```
"C:\path\to\moe\bin\moebatch.exe" -licwait -run moe_claw.svl ^
    -dir "C:/path/to/your/moe_files" ^
    -csv "C:/path/to/your/moe_files/moe_index.csv"
```

The SVL script perceives R/S stereochemistry from each 3D pose and forces it
into the SMILES, so chirality survives into the browser.

> **Google Drive note:** if your `.moe` files live on Google Drive File Stream,
> mark the folder "Available offline" first so `ReadAuto` isn't waiting on the
> cloud. Reading many large files can take several minutes — the indexer prints
> `[n/total] reading ...` progress so you can see it's working.

### Step 1b — Incremental updates (optional & under development - mileage may vary!)

Reindexing a whole folder is slow. `moeclaw_update.py` reads only `.moe` files
that are new or have changed since the last index, then merges the results into
the existing CSV:

```
python moeclaw_update.py --dir "C:/path/to/moe_files" ^
                         --moebatch "C:/Program Files/moe2024/bin/moebatch.exe"
```

Useful flags: `--dry-run` (report what's stale, never starts MOE), `--force`
(reindex everything — **required after any edit to `moe_claw.svl`**, since the
files haven't changed but every existing row is stale), `--keep-deleted`.

> **Known issue:** the staleness check currently flags every file as changed, so
> it reindexes everything. See [NOTES.md](NOTES.md) for the diagnosis in
> progress. The merge logic itself is verified.

### Step 2 — Build and open the browser

**Windows:** double-click `run_moeclaw.bat`, or:

```
run_moeclaw.bat path\to\moe_index.csv
```

**macOS / Linux:**

```
./run_moeclaw.sh path/to/moe_index.csv
```
### Step 2B - Already have an index file?

'''
Run MOEClaw - browseropener.bat (this will look for an existing html file and use it)
'''

The launcher activates the conda env, fetches RDKit-JS if needed, builds
`moe_browser.html`, starts a local server on port 8000, and opens the page.
**Leave the terminal window open** while you use it; press **Ctrl+C** to stop
the server when done.

**Or run the builder directly:**

```
python MOEClaw_buildbrowser.py moe_index.csv -o moe_browser.html
python -m http.server 8000
# then open http://localhost:8000/moe_browser.html
```

> **Why a local server?** Chrome refuses to load the RDKit WASM from a
> `file://` page (opening the HTML by double-click), so substructure search
> won't work that way. Serving over `http://localhost` fixes it. Firefox is
> more permissive if you prefer double-click, but the server is the reliable
> route.

---

## Builder options

```
python MOEClaw_buildbrowser.py INPUT.csv [options]

  -o, --output FILE     Output HTML (default: moe_browser.html)
  --mw-min N            Minimum molecular weight (default: 200)
  --mw-max N            Maximum molecular weight (default: 1000)
  --nres N              Keep only chains with this residue count (default: 1,
                        i.e. small molecules; drops protein/peptide chains)
  --smarts "PATTERN"    Pre-filter to compounds containing a SMARTS motif,
                        highlighting the matching atoms (e.g. "C=CC(=O)N")
  --keep-additives      Do NOT filter common buffers/cryoprotectants
                        (HEPES, glycerol, sulfate, etc.). They are removed
                        by default.
```

---

## What the browser page does

- **Renders** every compound as a 2D structure with stereochemistry
  (wedge/dash bonds and R/S/E/Z labels).
- **Text search** over compound name and SMILES.
- **Substructure search** — type a SMARTS/SMILES, or click a preset chip
  (alkyne, warhead, sulfinamide, indole, etc.). Matching atoms highlight.
- **Chemical similarity** — compounds are grouped into families by Butina
  clustering; click any structure to re-sort everything by Tanimoto similarity
  to it.
- **Additive filtering** — common crystallization/buffer molecules are removed
  by name and by structure.
- **File dates** — each card shows the date of the `.moe` file it came from.
  Sort **newest first**, or tick **latest version only** to collapse copies of
  the same structure down to the most recently modelled one (older copies are
  dimmed and badged when shown).
- **Selection + export** — tick compounds (or "select all shown"), then export
  the selection as CSV or SDF.
- **Copy buttons** on both the filename and the SMILES of each card.

---

## Files in this repo

| File | Purpose |
|------|---------|
| `moe_claw.svl` | MOE/SVL indexer — walks `.moe` files, writes the CSV (incl. file dates) |
| `MOEClaw_buildbrowser.py` | Builds the HTML browser from the CSV |
| `moeclaw_update.py` | Incremental indexing — reindex only new/changed `.moe` files |
| `fetch_rdkit_js.py` | Downloads the RDKit-JS runtime (run once) |
| `run_moeclaw.bat` / `.sh` | One-click build + serve launchers |
| `environment.yml` | Conda environment spec |
| `requirements.txt` | pip alternative |

Generated files (`moe_index.csv`, `moe_browser.html`) and the fetched RDKit-JS
binaries are git-ignored.

---

## Troubleshooting

**"substructure: failed to load" in the page** — you opened the HTML as a
`file://` page in Chrome. Use the launcher (which serves over http), or run
`python -m http.server 8000` and open `http://localhost:8000/moe_browser.html`.

**Some structures show a red "failed to parse" card** — the SMILES for those
couldn't be read by RDKit. The builder already auto-repairs one common MOE
quirk (charge-separated hypervalent sulfur in sulfoximines); anything still
flagged is a different issue — check the raw SMILES shown on the card.

**"conda env not found"** — create it: `conda env create -f environment.yml`.

**Port 8000 already in use** — an earlier server is still running, or Windows
is holding the port. Wait a minute, close the old terminal, or change the port
in the launcher.

**No stereochemistry in the structures** — the SMILES were exported flat. Make
sure the SVL indexer's `aSetForceRS [atoms, aRSChirality atoms]` line ran
(it perceives R/S from the 3D pose before extracting SMILES).

---

## License

MIT — see [LICENSE](LICENSE).
