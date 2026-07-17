#!/usr/bin/env python3
"""
moeclaw_update.py -- incrementally update a MOE-CLAW index CSV.

Only reads .moe files that are NEW or have CHANGED since the last index, then
merges the results into the master CSV. Reading .moe files is the slow part of
MOE-CLAW, so skipping unchanged ones makes routine updates fast.

How "changed" is decided: each row in the master CSV carries the file_date of
the .moe it came from. This compares that recorded date against the file's
current modification time on disk. That is deliberately NOT a "newer than the
last run" cutoff -- files copied or synced in often keep their original
timestamps, so a cutoff would silently skip them.

Usage:
    python moeclaw_update.py --dir DIR [--csv master.csv] [--moebatch PATH]
    python moeclaw_update.py --dir DIR --force        # reindex everything
    python moeclaw_update.py --dir DIR --dry-run      # just report what's stale

Requires: MOE (moebatch) for the indexing step. No RDKit needed.
"""

import argparse
import csv
import datetime
import os
import shutil
import subprocess
import sys
import tempfile

SVL_SCRIPT = "moe_claw.svl"
FIELDS = ["file", "file_date", "compound", "n_res", "n_atoms", "mw", "smiles"]

# SVL writes whole seconds; os.path.getmtime is a float. Allow a small slop so
# rounding alone never triggers a pointless reindex.
DATE_TOLERANCE_SEC = 2.0


def iso_from_mtime(ts):
    """Format a POSIX mtime the same way the SVL indexer does."""
    return datetime.datetime.fromtimestamp(ts).strftime("%Y-%m-%d %H:%M:%S")


def parse_iso(s):
    try:
        return datetime.datetime.strptime(s.strip(), "%Y-%m-%d %H:%M:%S")
    except (ValueError, AttributeError):
        return None


def read_master(path):
    """Return (rows, recorded_dates) from the master CSV. Missing file -> empty."""
    if not os.path.exists(path):
        return [], {}
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    recorded = {}
    for r in rows:
        fn = (r.get("file") or "").strip()
        if fn:
            # if a file somehow has rows with differing dates, keep the newest
            d = (r.get("file_date") or "").strip()
            if fn not in recorded or d > recorded[fn]:
                recorded[fn] = d
    return rows, recorded


def scan_dir(directory):
    """Return {basename: (fullpath, mtime)} for every .moe in directory."""
    out = {}
    for name in os.listdir(directory):
        if name.lower().endswith(".moe"):
            full = os.path.join(directory, name)
            if os.path.isfile(full):
                out[name] = (full, os.path.getmtime(full))
    return out


def classify(on_disk, recorded, force=False):
    """Split files into new / changed / unchanged (basenames)."""
    new, changed, unchanged = [], [], []
    for name, (full, mtime) in sorted(on_disk.items()):
        if force:
            changed.append(name)
            continue
        if name not in recorded:
            new.append(name)
            continue
        rec = parse_iso(recorded[name])
        cur = datetime.datetime.fromtimestamp(mtime)
        if rec is None:
            changed.append(name)
        elif abs((cur - rec).total_seconds()) > DATE_TOLERANCE_SEC:
            changed.append(name)
        else:
            unchanged.append(name)
    return new, changed, unchanged


def run_moebatch(moebatch, svl, files, out_csv):
    """Invoke MOE once for the whole stale set (one startup, not one per file)."""
    cmd = [moebatch, "-licwait", "-run", svl]
    for f in files:
        cmd += ["-file", f]
    cmd += ["-csv", out_csv]
    print("\nRunning MOE on %d file(s)..." % len(files))
    print("  (this is the slow part -- each .moe must be fully loaded)\n")
    proc = subprocess.run(cmd)
    return proc.returncode


def merge(master_rows, partial_rows, reindexed, deleted, keep_deleted):
    """Replace rows for reindexed files, drop rows for deleted files, append new.

    Replacing (not appending) matters: a file that gained or lost compounds
    would otherwise leave stale rows behind alongside the fresh ones.
    """
    drop = set(reindexed)
    if not keep_deleted:
        drop |= set(deleted)
    kept = [r for r in master_rows if (r.get("file") or "").strip() not in drop]
    return kept + partial_rows


def write_csv(path, rows):
    """Write atomically: build a temp file, then replace, so an interrupted
    write can't leave you with a truncated master index."""
    d = os.path.dirname(os.path.abspath(path)) or "."
    fd, tmp = tempfile.mkstemp(suffix=".csv", dir=d)
    os.close(fd)
    try:
        with open(tmp, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=FIELDS, extrasaction="ignore")
            w.writeheader()
            for r in rows:
                w.writerow(r)
        shutil.move(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def main():
    ap = argparse.ArgumentParser(
        description="Incrementally update a MOE-CLAW index CSV.")
    ap.add_argument("--dir", required=True, help="Folder of .moe files")
    ap.add_argument("--csv", default=None,
                    help="Master index CSV (default: <dir>/moe_index.csv)")
    ap.add_argument("--moebatch", default="moebatch",
                    help="Path to moebatch executable")
    ap.add_argument("--svl", default=None,
                    help="Path to moe_claw.svl (default: next to this script)")
    ap.add_argument("--force", action="store_true",
                    help="Reindex every file, ignoring dates. Use this after "
                         "changing the SVL logic -- existing rows are then "
                         "stale even though the .moe files have not changed.")
    ap.add_argument("--dry-run", action="store_true",
                    help="Report what would be indexed; do not run MOE")
    ap.add_argument("--keep-deleted", action="store_true",
                    help="Keep rows for .moe files that no longer exist "
                         "(default: purge them)")
    args = ap.parse_args()

    directory = os.path.abspath(args.dir)
    if not os.path.isdir(directory):
        sys.exit("ERROR: not a directory: %s" % directory)

    master_csv = args.csv or os.path.join(directory, "moe_index.csv")
    svl = args.svl or os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   SVL_SCRIPT)
    if not args.dry_run and not os.path.exists(svl):
        sys.exit("ERROR: SVL script not found: %s" % svl)

    master_rows, recorded = read_master(master_csv)
    on_disk = scan_dir(directory)
    if not on_disk:
        sys.exit("No .moe files found in %s" % directory)

    new, changed, unchanged = classify(on_disk, recorded, force=args.force)
    deleted = sorted(set(recorded) - set(on_disk))
    stale = new + changed

    print("Index : %s" % master_csv)
    print("Folder: %s" % directory)
    print("  %4d .moe file(s) on disk" % len(on_disk))
    print("  %4d already indexed and unchanged (skipped)" % len(unchanged))
    print("  %4d new" % len(new))
    print("  %4d changed since last index" % len(changed))
    if deleted:
        what = "kept" if args.keep_deleted else "will be purged from the index"
        print("  %4d indexed file(s) no longer on disk (%s)" % (len(deleted), what))

    if new:
        print("\nNew files:")
        for n in new[:20]:
            print("   +", n)
        if len(new) > 20:
            print("   ... and %d more" % (len(new) - 20))
    if changed:
        print("\nChanged files:")
        for n in changed[:20]:
            print("   *", n)
        if len(changed) > 20:
            print("   ... and %d more" % (len(changed) - 20))

    if not stale and not (deleted and not args.keep_deleted):
        print("\nNothing to do -- index is up to date. MOE was not started.")
        return

    if args.dry_run:
        print("\n[dry run] would index %d file(s); no changes made." % len(stale))
        return

    partial_rows = []
    if stale:
        fd, partial = tempfile.mkstemp(suffix=".csv")
        os.close(fd)
        try:
            paths = [on_disk[n][0].replace("\\", "/") for n in stale]
            rc = run_moebatch(args.moebatch, svl, paths, partial.replace("\\", "/"))
            if rc != 0:
                sys.exit("ERROR: moebatch exited with code %d. "
                         "Master index left unchanged." % rc)
            if not os.path.exists(partial) or os.path.getsize(partial) == 0:
                sys.exit("ERROR: MOE produced no output. "
                         "Master index left unchanged.")
            with open(partial, newline="", encoding="utf-8") as f:
                partial_rows = list(csv.DictReader(f))
        finally:
            if os.path.exists(partial):
                os.remove(partial)

        indexed_files = {(r.get("file") or "").strip() for r in partial_rows}
        missing = [n for n in stale if n not in indexed_files]
        if missing:
            print("\nWARNING: MOE returned no rows for %d file(s):" % len(missing))
            for n in missing[:10]:
                print("   !", n)
            print("  Their old rows (if any) are left in place.")
            # don't drop rows for files MOE didn't actually reindex
            stale = [n for n in stale if n in indexed_files]

    merged = merge(master_rows, partial_rows, stale, deleted, args.keep_deleted)
    write_csv(master_csv, merged)

    print("\nIndex updated: %s" % master_csv)
    print("  %d row(s) total (%d added/refreshed from %d file(s))"
          % (len(merged), len(partial_rows), len(stale)))
    print("\nNow rebuild the browser:")
    print("  python MOEClaw_buildbrowser.py \"%s\" -o moe_browser.html" % master_csv)


if __name__ == "__main__":
    main()
