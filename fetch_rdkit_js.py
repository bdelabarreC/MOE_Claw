#!/usr/bin/env python3
"""
fetch_rdkit_js.py -- download the RDKit-JS runtime (RDKit_minimal.js + .wasm)
that MOE-CLAW's browser needs for in-page substructure search.

These two files (~7 MB total) are NOT stored in the git repo because the .wasm
is a large binary. This script fetches them from the public npm registry and
places them next to your HTML so the offline browser works.

Usage:
    python fetch_rdkit_js.py            # downloads into the current folder
    python fetch_rdkit_js.py --dir DIR  # downloads into DIR

No external dependencies -- uses only the Python standard library.
"""

import argparse
import io
import json
import os
import sys
import tarfile
import urllib.request

REGISTRY = "https://registry.npmjs.org/@rdkit/rdkit"
FILES = ("RDKit_minimal.js", "RDKit_minimal.wasm")


def main():
    ap = argparse.ArgumentParser(description="Fetch RDKit-JS runtime for MOE-CLAW.")
    ap.add_argument("--dir", default=".", help="Destination folder (default: current)")
    ap.add_argument("--version", default=None,
                    help="Specific @rdkit/rdkit version (default: latest)")
    args = ap.parse_args()

    os.makedirs(args.dir, exist_ok=True)

    # Skip if both already present
    if all(os.path.exists(os.path.join(args.dir, f)) for f in FILES):
        print("RDKit-JS files already present in", os.path.abspath(args.dir))
        print("(delete them and re-run if you want to refresh.)")
        return

    print("Looking up @rdkit/rdkit on npm ...")
    with urllib.request.urlopen(REGISTRY, timeout=60) as r:
        meta = json.load(r)
    version = args.version or meta["dist-tags"]["latest"]
    tarball = meta["versions"][version]["dist"]["tarball"]
    print(f"  version {version}")
    print(f"  tarball {tarball}")

    print("Downloading package tarball ...")
    with urllib.request.urlopen(tarball, timeout=180) as r:
        data = r.read()
    print(f"  {len(data) // 1024} KB downloaded")

    print("Extracting runtime files ...")
    found = {}
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as tf:
        for member in tf.getmembers():
            base = os.path.basename(member.name)
            # take them from the package's dist/ folder, once each
            if (base in FILES and base not in found
                    and "/dist/" in member.name.replace("\\", "/")):
                fobj = tf.extractfile(member)
                if fobj is None:
                    continue
                out = os.path.join(args.dir, base)
                with open(out, "wb") as w:
                    w.write(fobj.read())
                found[base] = out
                print(f"  wrote {out}")

    missing = [f for f in FILES if f not in found]
    if missing:
        sys.exit(f"ERROR: could not find {missing} in the package. "
                 "The package layout may have changed.")

    # sanity: the wasm should start with the WebAssembly magic bytes
    wasm = os.path.join(args.dir, "RDKit_minimal.wasm")
    with open(wasm, "rb") as f:
        magic = f.read(4)
    if magic != b"\x00asm":
        sys.exit("ERROR: downloaded .wasm does not look valid (bad magic bytes).")

    print("\nDone. RDKit-JS runtime is ready in", os.path.abspath(args.dir))


if __name__ == "__main__":
    main()
