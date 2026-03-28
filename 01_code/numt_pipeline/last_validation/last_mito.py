#!/usr/bin/env python3
"""
last_mito.py — LAST alignment vs rCRS (NC_012920.1) for NUMT validation.

Runs last-train + lastal | last-split against the mitochondrial reference,
parses the MAF output, and produces a TSV with the same schema as
excel_blast_mito.tsv (+ LAST-specific columns: mismap, is_split, alignment_method).

This is a *validation* tool — it runs in parallel with BLAST, not as a replacement.
LAST's adaptive scoring (last-train) and split alignment (last-split) can catch
divergent or chimeric NUMTs that BLAST misses.

Usage:
  python last_validation/last_mito.py \
    --query   02_data/processed/excel_pipeline_query.fasta \
    --db-fasta 02_data/raw/NC_012920.1.fasta \
    --workdir 02_data/processed/last_work \
    --output  02_data/processed/excel_last_mito.tsv
"""

import argparse
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

# Import helpers from run_pipeline_excel.py
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from run_pipeline_excel import (
    MTDNA_FEATURES, BLAST_COLS, classify_hit, annotate_region, log,
)


# ─── LAST binary checks ─────────────────────────────────────────────────────

LAST_BINARIES = ["lastdb", "last-train", "lastal", "last-split"]


def ensure_last_binaries() -> bool:
    """Check that all LAST binaries are on PATH."""
    missing = [b for b in LAST_BINARIES if shutil.which(b) is None]
    if missing:
        log(f"FATAL — LAST binaries missing: {', '.join(missing)}")
        log("Install with: micromamba install -n NUMT -c bioconda last")
        return False
    return True


# ─── LAST workflow steps ─────────────────────────────────────────────────────

def build_lastdb(fasta: Path, db_prefix: Path) -> bool:
    """Build a LAST database from a FASTA file. Skip if already built."""
    if db_prefix.with_suffix(".suf").exists():
        log(f"LAST DB already exists: {db_prefix}")
        return True

    db_prefix.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["lastdb", "-Q", "0", "-R", "01", str(db_prefix), str(fasta)]
    log(f"Building LAST DB: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log(f"FATAL — lastdb failed: {result.stderr}")
        return False
    log("LAST DB built.")
    return True


def train_scoring(db_prefix: Path, query: Path, out_params: Path) -> bool:
    """Train substitution/gap parameters on the actual data."""
    if out_params.exists() and out_params.stat().st_mtime > query.stat().st_mtime:
        log(f"LAST params already trained: {out_params}")
        return True

    cmd = ["last-train", "--revsym", "-C2", str(db_prefix), str(query)]
    log(f"Training LAST scoring: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log(f"FATAL — last-train failed: {result.stderr}")
        return False

    out_params.write_text(result.stdout)
    log(f"LAST params saved: {out_params}")
    return True


def run_lastal(db_prefix: Path, query: Path, params: Path, out_maf: Path) -> bool:
    """Run lastal with trained params, pipe through last-split for chimera detection."""
    if out_maf.exists() and out_maf.stat().st_mtime > query.stat().st_mtime:
        log(f"LAST MAF already exists: {out_maf}")
        return True

    # lastal -p {params} -D1e6 {db} {query} | last-split -m1 > {maf}
    cmd_lastal = [
        "lastal", "-p", str(params), "-D1e6", str(db_prefix), str(query),
    ]
    cmd_split = ["last-split", "-m1"]
    log(f"Running lastal | last-split ...")

    p_lastal = subprocess.Popen(cmd_lastal, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p_split = subprocess.Popen(
        cmd_split, stdin=p_lastal.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    p_lastal.stdout.close()

    stdout, stderr_split = p_split.communicate()
    _, stderr_lastal = p_lastal.communicate()

    if p_lastal.returncode != 0:
        log(f"FATAL — lastal failed: {stderr_lastal.decode()}")
        return False
    if p_split.returncode != 0:
        log(f"FATAL — last-split failed: {stderr_split.decode()}")
        return False

    out_maf.write_bytes(stdout)
    log(f"LAST MAF written: {out_maf}")
    return True


# ─── MAF parsing ────────────────────────────────────────────────────────────

def _parse_lambda_from_params(params_path: Path) -> float | None:
    """Extract lambda value from last-train params file header."""
    try:
        for line in params_path.read_text().splitlines():
            m = re.search(r'lambda[=:\s]+([\d.eE+-]+)', line, re.IGNORECASE)
            if m:
                return float(m.group(1))
    except Exception:
        pass
    return None


def parse_maf(maf_path: Path, params_path: Path | None = None) -> pd.DataFrame:
    """
    Parse a MAF file produced by lastal | last-split into a DataFrame.

    Each alignment block (starting with 'a') contains:
      - 'a' line: score=NNN, mismap=X.XX (from last-split)
      - 's' lines: subject first, query second (when query is second arg to lastal)

    MAF coordinates are 0-based, start + alnSize = end. We convert to 1-based
    inclusive to match BLAST convention.
    """
    lam = _parse_lambda_from_params(params_path) if params_path else None

    rows = []
    current_score = None
    current_mismap = None
    s_lines = []

    with open(maf_path) as fh:
        for line in fh:
            line = line.rstrip()

            if line.startswith("a "):
                # Flush previous block
                if s_lines and len(s_lines) >= 2:
                    rows.append(_build_row(s_lines, current_score, current_mismap, lam))
                s_lines = []

                # Parse score and mismap
                score_m = re.search(r'score=(\d+)', line)
                current_score = int(score_m.group(1)) if score_m else 0
                mismap_m = re.search(r'mismap=([\d.eE+-]+)', line)
                current_mismap = float(mismap_m.group(1)) if mismap_m else None

            elif line.startswith("s "):
                s_lines.append(line)

    # Flush last block
    if s_lines and len(s_lines) >= 2:
        rows.append(_build_row(s_lines, current_score, current_mismap, lam))

    if not rows:
        return pd.DataFrame()

    return pd.DataFrame(rows)


def _build_row(s_lines: list[str], score: int, mismap: float | None, lam: float | None) -> dict:
    """Build a result dict from two MAF 's' lines (subject, query)."""
    # s line format: s name start alnSize strand srcSize alignment
    s_subj = s_lines[0].split()
    s_query = s_lines[1].split()

    sseqid = s_subj[1]
    s_start_0 = int(s_subj[2])
    s_aln_size = int(s_subj[3])
    s_strand = s_subj[4]
    s_src_size = int(s_subj[5])
    s_aln_text = s_subj[6]

    qseqid = s_query[1]
    q_start_0 = int(s_query[2])
    q_aln_size = int(s_query[3])
    q_strand = s_query[4]
    q_src_size = int(s_query[5])
    q_aln_text = s_query[6]

    # Convert to 1-based inclusive coords (BLAST convention)
    # Subject
    if s_strand == "+":
        sstart = s_start_0 + 1
        send = s_start_0 + s_aln_size
    else:
        sstart = s_src_size - s_start_0
        send = s_src_size - s_start_0 - s_aln_size + 1

    # Query
    if q_strand == "+":
        qstart = q_start_0 + 1
        qend = q_start_0 + q_aln_size
    else:
        qstart = q_src_size - q_start_0
        qend = q_src_size - q_start_0 - q_aln_size + 1

    # Compute pident, mismatch, gapopen from aligned strings
    matches = 0
    mismatches = 0
    gap_opens = 0
    aligned_cols = len(s_aln_text)
    in_gap_s = False
    in_gap_q = False

    for sc, qc in zip(s_aln_text, q_aln_text):
        if sc == "-":
            if not in_gap_s:
                gap_opens += 1
                in_gap_s = True
            in_gap_q = False
        elif qc == "-":
            if not in_gap_q:
                gap_opens += 1
                in_gap_q = True
            in_gap_s = False
        else:
            in_gap_s = False
            in_gap_q = False
            if sc.upper() == qc.upper():
                matches += 1
            else:
                mismatches += 1

    pident = (matches / aligned_cols * 100.0) if aligned_cols > 0 else 0.0
    length = q_aln_size  # aligned length on query side

    # E-value approximation from LAST score
    evalue = math.exp(-score * lam) if lam and score else float("nan")

    return {
        "qseqid": qseqid,
        "sseqid": sseqid,
        "pident": round(pident, 2),
        "length": length,
        "mismatch": mismatches,
        "gapopen": gap_opens,
        "qstart": qstart,
        "qend": qend,
        "sstart": sstart,
        "send": send,
        "evalue": evalue,
        "bitscore": score,  # LAST raw score, not numerically comparable to BLAST bitscore
        "qlen": q_src_size,
        "slen": s_src_size,
        "mismap": mismap,
    }


# ─── Enrichment ─────────────────────────────────────────────────────────────

def enrich_results(df: pd.DataFrame) -> pd.DataFrame:
    """Add allele_id, gene_symbol, hit_tier, mtdna_region, and LAST-specific columns."""
    if df.empty:
        return df

    # Extract AlleleID and gene from FASTA header (format: AlleleID|Gene|ClinSig|Length)
    df["allele_id"] = df["qseqid"].apply(
        lambda x: int(x.split("|")[0]) if "|" in str(x) else None
    )
    df["gene_symbol"] = df["qseqid"].apply(
        lambda x: x.split("|")[1] if "|" in str(x) and len(x.split("|")) > 1 else None
    )

    # Hit tier classification (same logic as BLAST)
    df["hit_tier"] = df.apply(classify_hit, axis=1)

    # mtDNA region annotation
    regions = df.apply(lambda r: annotate_region(int(r["sstart"]), int(r["send"])), axis=1)
    df["mtdna_region"] = regions.apply(lambda x: x[0])
    df["mtdna_region_type"] = regions.apply(lambda x: x[1])

    # Coverage and combined score
    df["query_coverage"] = (df["length"] / df["qlen"]) * 100.0
    df["combined_score"] = df["pident"] * (df["length"] / df["qlen"]) / 100.0

    # LAST-specific: detect split alignments (>1 block per query)
    split_counts = df.groupby("qseqid").size()
    df["is_split"] = df["qseqid"].map(split_counts) > 1

    # Alignment method tag
    df["alignment_method"] = "LAST"

    return df


# ─── Main ────────────────────────────────────────────────────────────────────

def run_last_mito(query: Path, db_fasta: Path, workdir: Path, output: Path) -> pd.DataFrame:
    """Full LAST workflow: lastdb → last-train → lastal|last-split → parse → enrich."""
    if not ensure_last_binaries():
        return pd.DataFrame()

    workdir.mkdir(parents=True, exist_ok=True)
    db_prefix = workdir / "mito_lastdb"
    params_file = workdir / "mito.par"
    maf_file = workdir / "mito.maf"

    if not build_lastdb(db_fasta, db_prefix):
        return pd.DataFrame()

    if not train_scoring(db_prefix, query, params_file):
        return pd.DataFrame()

    if not run_lastal(db_prefix, query, params_file, maf_file):
        return pd.DataFrame()

    log("Parsing MAF output ...")
    df = parse_maf(maf_file, params_path=params_file)

    if df.empty:
        log("LAST mito: no hits found.")
        return pd.DataFrame()

    log(f"LAST mito: {len(df):,} alignment blocks parsed")
    df = enrich_results(df)

    # Save TSV
    output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output, sep="\t", index=False)
    log(f"LAST mito results saved: {output}")

    return df


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="LAST alignment vs rCRS for NUMT validation.",
    )
    parser.add_argument(
        "--query", type=Path, required=True,
        help="Query FASTA (e.g. excel_pipeline_query.fasta)",
    )
    parser.add_argument(
        "--db-fasta", type=Path, required=True,
        help="rCRS FASTA (e.g. NC_012920.1.fasta)",
    )
    parser.add_argument(
        "--workdir", type=Path, required=True,
        help="Working directory for LAST intermediate files",
    )
    parser.add_argument(
        "--output", type=Path, required=True,
        help="Output TSV path",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_last_mito(args.query, args.db_fasta, args.workdir, args.output)
