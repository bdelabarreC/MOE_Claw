#!/usr/bin/env bash
# MOE-CLAW: build the browser from a CSV and serve it (macOS / Linux).
#
#   1. Activates the conda "moeclaw" environment
#   2. Fetches RDKit-JS runtime the first time (if not already present)
#   3. Builds moe_browser.html from your CSV
#   4. Starts a local server and opens the page
#
# Usage:
#   ./run_moeclaw.sh                 # uses moe_index.csv in this folder
#   ./run_moeclaw.sh path/to/index.csv

set -euo pipefail

APPDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT="MOEClaw_buildbrowser.py"
CONDA_ENV="moeclaw"
PORT=8000

INPUT_CSV="${1:-$APPDIR/moe_index.csv}"
OUTPUT_HTML="$APPDIR/moe_browser.html"

echo
echo "==== MOE-CLAW ===="
echo "  Folder    : $APPDIR"
echo "  Input CSV : $INPUT_CSV"
echo "  Conda env : $CONDA_ENV"
echo

[ -f "$APPDIR/$SCRIPT" ] || { echo "ERROR: $SCRIPT not found here."; exit 1; }
[ -f "$INPUT_CSV" ] || { echo "ERROR: input CSV not found: $INPUT_CSV"; exit 1; }

# Activate conda
if command -v conda >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$CONDA_ENV" || {
        echo "ERROR: could not activate '$CONDA_ENV'."
        echo "Create it first:  conda env create -f environment.yml"
        exit 1
    }
else
    echo "WARNING: conda not found on PATH; assuming rdkit/pandas are available."
fi

# Fetch RDKit-JS if missing
if [ ! -f "$APPDIR/RDKit_minimal.wasm" ]; then
    echo "First run: fetching RDKit-JS runtime (~7 MB)..."
    python "$APPDIR/fetch_rdkit_js.py" --dir "$APPDIR"
fi

# Build
echo "Building HTML from CSV..."
python "$APPDIR/$SCRIPT" "$INPUT_CSV" -o "$OUTPUT_HTML"
echo "Build complete: $OUTPUT_HTML"

# Serve + open
echo
echo "Starting local server on http://localhost:$PORT/"
echo "  Leave this terminal open; press Ctrl+C to stop the server."
echo
cd "$APPDIR"

# open the browser (macOS: open, Linux: xdg-open)
URL="http://localhost:$PORT/moe_browser.html"
( sleep 1
  if command -v open >/dev/null 2>&1; then open "$URL"
  elif command -v xdg-open >/dev/null 2>&1; then xdg-open "$URL"
  fi ) &

python -m http.server "$PORT"
