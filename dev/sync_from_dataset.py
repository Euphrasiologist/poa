#!/usr/bin/env python3
"""Re-sync this repo's analysis/ code from the dataset repo it was split out of.

    python3 dev/sync_from_dataset.py ../plant_organellar_database [<commit>]

Exports <commit> (default HEAD) of the dataset repo with `git archive`, so
only committed code is taken, never a working tree mid-change. Copies
every module's src/, README.md and trans_splicing/reference/ (the bundled
GenBank references, exon sets and templates), leaves out per-dataset
material (results/, work/, logs/, gff_export/ground_truth/, the dataset's
own analysis/README.md and run_all.sh) - except editing/results/profiles/,
which ships as editing/reference/profiles/ (the fallback for datasets too
small to build their own) - and applies the rewrites below so
code and data can live in different trees. Files under analysis/ that no
longer exist upstream are deleted. The source commit is written to
dev/SYNCED_FROM.

Rewrites (each asserts it matched, so an upstream refactor fails loudly
instead of syncing half-converted code):
  - Python: CODE_DIR (this checkout's analysis/) for imports and bundled
    reference data; ROOT_DIR = $PLANT_ORGANELLE_DATA_ROOT (default: this
    checkout) for data/ and every module's results/ and work/.
  - Shell: species_discovery.sh, args.sh, linearize's resolve scripts
    honour $PLANT_ORGANELLE_DATA_ROOT the same way.
  - Cluster-only install paths (/software/team301/...) become env-var
    overridable and are only used when they exist.
"""
from __future__ import annotations

import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

EXCLUDE = [
    re.compile(r"^analysis/[^/]+/(results|work|logs)(/|$)"),
    re.compile(r"^analysis/gff_export/ground_truth/"),
    re.compile(r"^analysis/(README\.md|run_all\.sh)$"),
    re.compile(r"(^|/)(__pycache__|\.cache)/"),
]

PY_HEADER = (
    "# CODE_DIR is this checkout's analysis/ (shared modules, bundled reference\n"
    "# data); data/ and every module's results/ and work/ live under the data\n"
    "# root - PLANT_ORGANELLE_DATA_ROOT, or this checkout if unset.\n"
    "CODE_DIR = Path(__file__).resolve().parents[2]\n"
    'ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT", str(CODE_DIR.parent)))\n'
    'ANALYSIS_DIR = ROOT_DIR / "analysis"\n'
)
ARCHIVE_PATHSPECS = [
    ":(glob)analysis/*/src/**",
    ":(glob)analysis/*/README.md",
    "analysis/common",
    "analysis/trans_splicing/reference",
    "analysis/editing/results/profiles",
]
EDITING_PROFILES_SRC ="analysis/editing/results/profiles/"
EDITING_PROFILES_DEST = "analysis/editing/reference/profiles/"
DATA_ROOT_SH ='"${PLANT_ORGANELLE_DATA_ROOT:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"'


def sub(text: str, pattern: str, repl: str, path: str, count: int = 0, flags: int = 0) -> str:
    new, n = re.subn(pattern, repl, text, count=count, flags=flags)
    if n == 0:
        sys.exit(f"[err] sync: rewrite {pattern!r} did not match {path} - upstream changed, update this script")
    return new


def ensure_import_os(text: str) -> str:
    if re.search(r"^import os$", text, re.M):
        return text
    # after the last top-of-file plain `import x` line
    return re.sub(r"^(import (?!os\b)\w+)$", r"import os\n\1", text, count=1, flags=re.M)


def rewrite_python(rel: str, text: str) -> str:
    if "ANALYSIS_DIR = Path(__file__).resolve().parents[2]" in text:
        text = sub(text, r"ANALYSIS_DIR = Path\(__file__\)\.resolve\(\)\.parents\[2\]\n(ROOT_DIR = ANALYSIS_DIR\.parent\n)?",
                   PY_HEADER.replace("\\", "\\\\"), rel, count=1)
        text = ensure_import_os(text)
        # code-side paths: shared modules and bundled reference data
        text = text.replace('ANALYSIS_DIR / "common"', 'CODE_DIR / "common"')
        text = text.replace('ANALYSIS_DIR / "trans_splicing" / "src"', 'CODE_DIR / "trans_splicing" / "src"')
        text = text.replace('ANALYSIS_DIR / "trans_splicing" / "reference"', 'CODE_DIR / "trans_splicing" / "reference"')
        text = text.replace("ANALYSIS_DIR.parent", "ROOT_DIR")

    if rel == "analysis/common/species_discovery.py":
        text = sub(text, r'Path\(args\.root_dir\) if args\.root_dir else Path\(__file__\)\.resolve\(\)\.parents\[2\]',
                   'Path(args.root_dir or os.environ.get("PLANT_ORGANELLE_DATA_ROOT") or Path(__file__).resolve().parents[2])',
                   rel)
        text = sub(text, r'help="repo root \(default: 3 levels up from this file\)"',
                   'help="data root (default: $PLANT_ORGANELLE_DATA_ROOT, else this checkout)"', rel)
        text = ensure_import_os(text)

    if rel == "analysis/denovo_annotation/src/00_run_nhmmscan.py":
        text = sub(text, r'"mito": "(/software/team301/OatkDB/[^"]+)",\n(\s*)"pltd": "(/software/team301/OatkDB/[^"]+)",',
                   r'"mito": os.environ.get("OATKDB_MITO_FAM", "\1"),\n\2"pltd": os.environ.get("OATKDB_PLTD_FAM", "\3"),',
                   rel)

    if rel == "analysis/denovo_annotation/src/01_run_trnascan.py":
        # a conda tRNAscan-SE finds its own config; the cluster build needs -c
        text = sub(text, r'TRNASCAN_CONF = "(/software/team301/[^"]+)"',
                   r'TRNASCAN_CONF = os.environ.get("TRNASCAN_CONF", "\1")', rel)
        text = sub(text, r'\[TRNASCAN, "-c", TRNASCAN_CONF, ',
                   '[TRNASCAN, *(["-c", TRNASCAN_CONF] if Path(TRNASCAN_CONF).exists() else []), ', rel)
        text = ensure_import_os(text)

    if rel == "analysis/orf_scan/src/02_scan_pfam.py":
        text = sub(text, r'PFAM_DB = "(/software/team301/[^"]+)"', r'PFAM_DB = os.environ.get("PFAM_DB", "\1")', rel)
        text = ensure_import_os(text)

    if rel == "analysis/editing/src/01_scan_editing.py":
        text = sub(text, r'( *)profiles_dir = ANALYSIS_DIR / "editing" / "results" / "profiles"\n',
                   r'\1profiles_dir = ANALYSIS_DIR / "editing" / "results" / "profiles"\n'
                   r'\1if not any(profiles_dir.glob("*/*.pssm")):\n'
                   r'\1    # a dataset too small for 00_build_reference_profiles.py (>=5 species\n'
                   r'\1    # per gene) uses the bundled ones, built on the 1250-species dataset\n'
                   r'\1    profiles_dir = CODE_DIR / "editing" / "reference" / "profiles"\n'
                   r'\1    print(f"[info] scan_editing: no dataset-built profiles, using bundled {profiles_dir}", file=sys.stderr)\n',
                   rel)

    if 'env = {"MAFFT_BINARIES": MAFFT_BINARIES_DIR, "PATH": "/usr/bin:/bin"}' in text:
        # the minimal env only suits the cluster's hand-built mafft; a conda
        # mafft needs its own PATH and no MAFFT_BINARIES override
        text = sub(text, r'MAFFT_BINARIES_DIR = "(/software/team301/[^"]+)"',
                   r'MAFFT_BINARIES_DIR = os.environ.get("MAFFT_BINARIES", "\1")', rel)
        text = sub(text, r'( *)env = \{"MAFFT_BINARIES": MAFFT_BINARIES_DIR, "PATH": "/usr/bin:/bin"\}\n',
                   r'\1env = ({"MAFFT_BINARIES": MAFFT_BINARIES_DIR, "PATH": "/usr/bin:/bin"}\n'
                   r'\1       if Path(MAFFT_BINARIES_DIR).is_dir() else dict(os.environ))\n', rel)
        text = ensure_import_os(text)
    return text


def rewrite_shell(rel: str, text: str) -> str:
    if rel == "analysis/common/species_discovery.sh":
        text = sub(text, r'ROOT_DIR="\$\(cd "\$\{SCRIPT_DIR\}/\.\./\.\." && pwd\)"', f"ROOT_DIR={DATA_ROOT_SH}", rel)
        # per data root, not per checkout: one install can serve several datasets
        text = sub(text, r'CACHE_DIR="\$\{SCRIPT_DIR\}/\.cache"', 'CACHE_DIR="${ROOT_DIR}/analysis/.cache"', rel)
    if rel == "analysis/common/args.sh":
        text = sub(text, r'  ANALYSIS_DIR="\$\(cd "\$\{script_dir\}/\.\./\.\." && pwd\)"\n  ROOT_DIR="\$\(cd "\$\{ANALYSIS_DIR\}/\.\." && pwd\)"\n  export ANALYSIS_DIR ROOT_DIR',
                   '  CODE_DIR="$(cd "${script_dir}/../.." && pwd)"\n'
                   '  ROOT_DIR="${PLANT_ORGANELLE_DATA_ROOT:-$(cd "${CODE_DIR}/.." && pwd)}"\n'
                   '  ANALYSIS_DIR="${ROOT_DIR}/analysis"\n'
                   '  export CODE_DIR ANALYSIS_DIR ROOT_DIR', rel)
    if rel.startswith("analysis/linearize/src/") and 'ANALYSIS_DIR="$(cd "${SRC_DIR}/../.." && pwd)"' in text:
        text = sub(text, r'ANALYSIS_DIR="\$\(cd "\$\{SRC_DIR\}/\.\./\.\." && pwd\)"\nROOT_DIR="\$\(cd "\$\{ANALYSIS_DIR\}/\.\." && pwd\)"',
                   'ROOT_DIR="${PLANT_ORGANELLE_DATA_ROOT:-$(cd "${SRC_DIR}/../../.." && pwd)}"\n'
                   'ANALYSIS_DIR="${ROOT_DIR}/analysis"', rel)
    return text


def main() -> None:
    if len(sys.argv) not in (2, 3):
        sys.exit(__doc__)
    src = Path(sys.argv[1]).resolve()
    rev = sys.argv[2] if len(sys.argv) == 3 else "HEAD"
    commit = subprocess.check_output(["git", "-C", str(src), "rev-parse", rev], text=True).strip()
    subject = subprocess.check_output(["git", "-C", str(src), "log", "-1", "--format=%h %ad %s", "--date=short", commit],
                                      text=True).strip()

    with tempfile.TemporaryDirectory() as tmp:
        tar = Path(tmp) / "src.tar"
        # only what can be synced - the committed results/ are far larger than the code
        subprocess.run(["git", "-C", str(src), "archive", "-o", str(tar), commit, *ARCHIVE_PATHSPECS], check=True)
        with tarfile.open(tar) as tf:
            tf.extractall(tmp, filter="data")
        exported = sorted(p.relative_to(tmp).as_posix() for p in (Path(tmp) / "analysis").rglob("*") if p.is_file())
        wanted = {rel: rel for rel in exported if not any(x.search(rel) for x in EXCLUDE)}
        # dataset-built editing profiles ship as the fallback for small datasets
        for rel in exported:
            if rel.startswith(EDITING_PROFILES_SRC):
                wanted[EDITING_PROFILES_DEST + rel[len(EDITING_PROFILES_SRC):]] = rel

        n_changed = 0
        for rel, src_rel in sorted(wanted.items()):
            data = (Path(tmp) / src_rel).read_bytes()
            if rel.endswith((".py", ".sh")):
                # this repo's conda env name (environment.yml)
                data = data.replace(b"conda activate plant_organellar_database", b"conda activate plant_organelle_annotator")
            if rel.endswith(".py"):
                data = rewrite_python(rel, data.decode()).encode()
            elif rel.endswith(".sh"):
                data = rewrite_shell(rel, data.decode()).encode()
            dest = REPO / rel
            if not dest.exists() or dest.read_bytes() != data:
                dest.parent.mkdir(parents=True, exist_ok=True)
                dest.write_bytes(data)
                shutil.copymode(Path(tmp) / src_rel, dest)
                n_changed += 1

    keep = set(wanted)
    tracked = subprocess.check_output(["git", "-C", str(REPO), "ls-files", "analysis"], text=True).split()
    removed = [rel for rel in tracked if rel not in keep]
    for rel in removed:
        (REPO / rel).unlink(missing_ok=True)
    (REPO / "dev" / "SYNCED_FROM").write_text(f"{commit}\n{subject}\n")
    print(f"[info] sync: {len(wanted)} files from {subject}; {n_changed} written, {len(removed)} removed",
          file=sys.stderr)


if __name__ == "__main__":
    main()
