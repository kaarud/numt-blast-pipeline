#!/usr/bin/env python3
"""
01_clinvar_fetch.py — Download and filter ClinVar variants for candidate NUMT detection.

BIOLOGICAL RATIONALE
====================
NUMTs (Nuclear Mitochondrial DNA segments) are fragments of mitochondrial DNA that
have been transferred and integrated into the nuclear genome over evolutionary time.
The human nuclear genome harbours hundreds of such insertions, ranging from small
fragments to near-complete mitochondrial genome copies (~16.5 kb).

When a NUMT insertion occurs within a coding exon it can disrupt gene function and
cause Mendelian disease.  A confirmed example is a 138 bp mitochondrial fragment
inserted into exon 17 of TSC2, causing Tuberous Sclerosis Complex (Bhatt et al.).

WHY ClinVar INSERTIONS > 50 bp?
================================
ClinVar (https://www.ncbi.nlm.nih.gov/clinvar/) is NCBI's public archive linking
human genetic variants to clinical phenotypes.  Large insertions (>50 bp) that are
classified as pathogenic or of uncertain significance may represent unrecognised
NUMTs because:
  1. Short-read sequencing pipelines often mischaracterise NUMTs as novel insertions.
  2. NUMTs in the human reference range from ~30 bp to ~16 kb; a 50 bp floor
     captures the majority while excluding noisy short indels.
  3. Insertions of mitochondrial origin can be identified downstream by aligning the
     inserted sequence to the revised Cambridge Reference Sequence (rCRS).

This script performs Step 1 of the pipeline: bulk download of ClinVar's
variant_summary.txt.gz flat-file, followed by filtering for large insertion/indel
variants that are plausible NUMT candidates.
"""

# ─── Standard-library imports ────────────────────────────────────────────────
import argparse
import gzip
import re
import sys
from datetime import datetime
from pathlib import Path

# ─── Third-party imports ─────────────────────────────────────────────────────
import pandas as pd
import requests
from tqdm import tqdm


# ─── Constants ───────────────────────────────────────────────────────────────

# NCBI distributes ClinVar flat files via HTTPS (the FTP mirror is also
# reachable over HTTPS, which avoids firewall issues with passive FTP).
CLINVAR_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
)

# Project root is two levels above this script (01_code/numt_pipeline/../../)
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

# Paths relative to project root, using pathlib to stay OS-agnostic
RAW_DIR = PROJECT_ROOT / "02_data" / "raw"
PROCESSED_DIR = PROJECT_ROOT / "02_data" / "processed"
RAW_FILE = RAW_DIR / "variant_summary.txt.gz"
OUTPUT_FILE = PROCESSED_DIR / "clinvar_insertions_gt50bp.tsv"

# Minimum insertion length (in bp) to consider as a NUMT candidate
MIN_LENGTH_BP = 50

# When Length is missing, we fall back to inspecting the HGVS Name column.
# This regex matches "ins" followed by at least MIN_LENGTH_BP characters of
# inserted sequence (nucleotides or other annotation characters).
HGVS_INS_PATTERN = re.compile(
    rf"ins[A-Za-z0-9]{{  {MIN_LENGTH_BP},}}",  # note: spaces inside {} are stripped below
    re.IGNORECASE,
)
# Clean up the pattern — re.compile does not allow spaces inside quantifiers,
# so we build it properly:
HGVS_INS_PATTERN = re.compile(
    rf"ins[A-Za-z0-9]{{{MIN_LENGTH_BP},}}", re.IGNORECASE
)


# ─── Helper: timestamped logging ────────────────────────────────────────────

def log(message: str) -> None:
    """Print a message prefixed with the current timestamp."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")


# ─── Step 1: Download variant_summary.txt.gz ────────────────────────────────

def download_clinvar(force: bool = False) -> Path:
    """
    Stream-download ClinVar's variant_summary.txt.gz with a progress bar.

    Parameters
    ----------
    force : bool
        If True, re-download even when the local file already exists.

    Returns
    -------
    Path
        Path to the downloaded gzip file.

    Raises
    ------
    SystemExit
        On unrecoverable network errors (DNS failure, HTTP 4xx/5xx, timeout).
    """
    RAW_DIR.mkdir(parents=True, exist_ok=True)  # ensure output directory exists

    if RAW_FILE.exists() and not force:
        log(f"File already exists: {RAW_FILE}  (use --force to re-download)")
        return RAW_FILE

    log(f"Downloading ClinVar variant summary from:\n         {CLINVAR_URL}")

    try:
        # stream=True defers body download so we can iterate in chunks
        response = requests.get(CLINVAR_URL, stream=True, timeout=60)
        response.raise_for_status()  # raises HTTPError for 4xx/5xx
    except requests.exceptions.ConnectionError as exc:
        log(f"FATAL — Connection error (check your network / DNS): {exc}")
        sys.exit(1)
    except requests.exceptions.Timeout:
        log("FATAL — Request timed out after 60 s.  Retry later.")
        sys.exit(1)
    except requests.exceptions.HTTPError as exc:
        log(f"FATAL — HTTP error: {exc}")
        sys.exit(1)
    except requests.exceptions.RequestException as exc:
        log(f"FATAL — Unexpected request error: {exc}")
        sys.exit(1)

    # Content-Length may be absent; default to 0 so tqdm shows a spinner
    total_bytes = int(response.headers.get("Content-Length", 0))
    chunk_size = 1024 * 256  # 256 KB chunks — balance between progress updates and I/O overhead

    # Write to a temporary path first, then rename — avoids leaving a partial
    # file on disk if the download is interrupted.
    tmp_file = RAW_FILE.with_suffix(".tmp")

    with (
        open(tmp_file, "wb") as fh,
        tqdm(
            total=total_bytes,
            unit="B",
            unit_scale=True,
            desc="variant_summary.txt.gz",
            disable=False,
        ) as pbar,
    ):
        for chunk in response.iter_content(chunk_size=chunk_size):
            fh.write(chunk)
            pbar.update(len(chunk))

    tmp_file.rename(RAW_FILE)  # atomic on same filesystem
    log(f"Download complete: {RAW_FILE}  ({RAW_FILE.stat().st_size / 1e6:.1f} MB)")
    return RAW_FILE


# ─── Step 2: Parse the gzipped TSV into a DataFrame ─────────────────────────

def load_clinvar(path: Path) -> pd.DataFrame:
    """
    Read variant_summary.txt.gz into a pandas DataFrame.

    The file is tab-delimited with a single header row.  We use the
    built-in gzip support of pandas and set dtype=str initially so that
    no column is silently coerced (e.g. chromosome "X" being parsed as NaN).

    Returns
    -------
    pd.DataFrame
        Raw, unfiltered ClinVar data.
    """
    log(f"Parsing gzipped TSV: {path}")

    df = pd.read_csv(
        path,
        sep="\t",
        dtype=str,            # read everything as string first — we cast numerics below
        low_memory=False,     # avoid mixed-type inference across chunks
        on_bad_lines="warn",  # log malformed rows instead of crashing
    )

    log(f"Loaded {len(df):,} variant records with {len(df.columns)} columns")
    return df


# ─── Step 3: Filter for large insertions / indels ───────────────────────────

def filter_large_insertions(df: pd.DataFrame) -> pd.DataFrame:
    """
    Retain variants that are plausible NUMT candidates.

    Primary filter (Length-based):
        Type ∈ {"Insertion", "Indel"}  AND  Length >= 50

    Secondary filter (HGVS-based, catches missing Length values):
        Type ∈ {"Insertion", "Indel"}  AND  Name matches  /ins[ACGT…]{50,}/i

    Parameters
    ----------
    df : pd.DataFrame
        Raw ClinVar DataFrame (all columns as str).

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame with a new numeric ``Length_bp`` column.
    """
    log("Applying filters: Type ∈ {Insertion, Indel} AND Length ≥ 50 bp …")

    # --- Normalise the Type column for case-insensitive matching ----
    # Some ClinVar records use mixed case; lowering avoids missed matches.
    type_lower = df["Type"].str.lower()
    is_ins_or_indel = type_lower.isin(["insertion", "indel"])

    # --- Derive insertion length from available columns ----
    # ClinVar dropped the dedicated "Length" column. We compute it from the
    # VCF alleles: len(ALT) - len(REF) gives the net inserted bases.
    # "na" strings and missing values are coerced to NaN.
    if "Length" in df.columns:
        length_numeric = pd.to_numeric(df["Length"], errors="coerce")
    elif "AlternateAlleleVCF" in df.columns and "ReferenceAlleleVCF" in df.columns:
        alt_len = df["AlternateAlleleVCF"].str.len().where(
            df["AlternateAlleleVCF"].notna() & (df["AlternateAlleleVCF"] != "na"), other=pd.NA
        )
        ref_len = df["ReferenceAlleleVCF"].str.len().where(
            df["ReferenceAlleleVCF"].notna() & (df["ReferenceAlleleVCF"] != "na"), other=pd.NA
        )
        length_numeric = pd.to_numeric(alt_len - ref_len, errors="coerce")
    else:
        log("WARNING: No Length or VCF allele columns found — relying on HGVS filter only.")
        length_numeric = pd.Series(pd.NA, index=df.index, dtype="float64")

    # --- Primary filter: numeric length available and >= threshold ----
    primary_mask = is_ins_or_indel & (length_numeric >= MIN_LENGTH_BP)

    # --- Secondary filter: Length missing but HGVS name implies long insertion ----
    # This rescues variants where Length was not annotated but the HGVS
    # nomenclature encodes the inserted sequence literally (e.g.,
    # "NM_000548.5:c.1832_1833ins<72-nt sequence>").
    length_missing = length_numeric.isna()
    hgvs_long_ins = df["Name"].str.contains(HGVS_INS_PATTERN, na=False)
    secondary_mask = is_ins_or_indel & length_missing & hgvs_long_ins

    # --- Combine both masks ----
    combined_mask = primary_mask | secondary_mask

    filtered = df.loc[combined_mask].copy()

    # Add a clean numeric length column for downstream analysis
    filtered["Length_bp"] = length_numeric[filtered.index]

    n_primary = primary_mask.sum()
    n_secondary = (secondary_mask & ~primary_mask).sum()  # exclusively from secondary
    log(
        f"Primary filter (Length ≥ {MIN_LENGTH_BP}): {n_primary:,} variants  |  "
        f"Secondary filter (HGVS pattern): {n_secondary:,} additional variants"
    )
    log(f"Total passing: {len(filtered):,} variants")

    return filtered


# ─── Step 4: Save filtered results ──────────────────────────────────────────

def save_results(df: pd.DataFrame) -> Path:
    """Write the filtered DataFrame to a tab-separated file."""
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTPUT_FILE, sep="\t", index=False)
    log(f"Saved {len(df):,} variants to: {OUTPUT_FILE}")
    return OUTPUT_FILE


# ─── Step 5: Print summary statistics ───────────────────────────────────────

def print_summary(total_variants: int, filtered_df: pd.DataFrame) -> None:
    """
    Display a human-readable summary of the filtering results.

    Includes total counts, breakdown by Type, and breakdown by
    ClinicalSignificance — the latter is especially relevant because
    pathogenic/likely-pathogenic NUMTs are the most clinically actionable.
    """
    print("\n" + "=" * 72)
    print("  ClinVar Large-Insertion Filter — Summary")
    print("=" * 72)

    print(f"\n  Total variants in ClinVar download : {total_variants:>10,}")
    print(f"  Variants passing filter (≥{MIN_LENGTH_BP} bp ins/indel): {len(filtered_df):>10,}")
    pct = len(filtered_df) / total_variants * 100 if total_variants else 0
    print(f"  Fraction retained                  : {pct:>10.3f} %")

    # --- Breakdown by Type ----
    print(f"\n  {'Breakdown by Type':─<50}")
    if "Type" in filtered_df.columns:
        type_counts = filtered_df["Type"].value_counts()
        for variant_type, count in type_counts.items():
            print(f"    {variant_type:<35} {count:>8,}")

    # --- Breakdown by ClinicalSignificance ----
    # This column uses free-text values; the most common are:
    # Pathogenic, Likely pathogenic, Uncertain significance, Benign, etc.
    sig_col = "ClinicalSignificance"
    if sig_col in filtered_df.columns:
        print(f"\n  {'Breakdown by ClinicalSignificance':─<50}")
        sig_counts = filtered_df[sig_col].value_counts().head(15)  # top 15 to avoid clutter
        for sig, count in sig_counts.items():
            # Truncate very long significance strings for display
            sig_display = sig[:45] + "…" if len(str(sig)) > 45 else sig
            print(f"    {sig_display:<47} {count:>8,}")
        remaining = len(filtered_df) - sig_counts.sum()
        if remaining > 0:
            print(f"    {'(other categories)':<47} {remaining:>8,}")

    print("\n" + "=" * 72 + "\n")


# ─── CLI entry point ────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    """Define and parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Step 1 of the NUMT pipeline: download ClinVar variant_summary.txt.gz "
            "and filter for large insertions/indels (≥50 bp) that may represent "
            "unrecognised nuclear-mitochondrial DNA insertions."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-download variant_summary.txt.gz even if it already exists locally.",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=MIN_LENGTH_BP,
        help=f"Minimum insertion length in bp (default: {MIN_LENGTH_BP}).",
    )
    return parser.parse_args()


def main() -> None:
    """Orchestrate the full download-filter-save workflow."""
    args = parse_args()

    # Allow the threshold to be overridden from the CLI
    global MIN_LENGTH_BP, HGVS_INS_PATTERN  # noqa: PLW0603 — intentional override
    MIN_LENGTH_BP = args.min_length
    HGVS_INS_PATTERN = re.compile(
        rf"ins[A-Za-z0-9]{{{MIN_LENGTH_BP},}}", re.IGNORECASE
    )

    log("=" * 60)
    log("NUMT Pipeline — Step 1: ClinVar Fetch & Filter")
    log("=" * 60)

    # Step 1 — Download
    gz_path = download_clinvar(force=args.force)

    # Step 2 — Parse
    df_raw = load_clinvar(gz_path)
    total_variants = len(df_raw)

    # Step 3 — Filter
    df_filtered = filter_large_insertions(df_raw)

    # Step 4 — Save
    save_results(df_filtered)

    # Step 5 — Summarise
    print_summary(total_variants, df_filtered)

    log("Done.")


if __name__ == "__main__":
    main()
