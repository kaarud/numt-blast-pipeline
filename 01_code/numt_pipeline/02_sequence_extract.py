#!/usr/bin/env python3
"""
02_sequence_extract.py — Extract inserted sequences from ClinVar HGVS notation.

BIOLOGICAL RATIONALE
====================
In the previous step (01_clinvar_fetch.py) we downloaded ClinVar's variant summary
and filtered for large insertions (>50 bp) that could represent NUMTs — fragments
of mitochondrial DNA that have been transferred into the nuclear genome.

To confirm whether an insertion is mitochondrial in origin, we need the actual
nucleotide sequence of the inserted fragment so we can BLAST it against the human
mitochondrial genome (the revised Cambridge Reference Sequence, NC_012920.1).

WHAT IS HGVS NOTATION?
=======================
HGVS (Human Genome Variation Society) nomenclature is the international standard
for describing sequence variants in human DNA.  The general form for an insertion is:

    <reference>:<coordinate_type>.<position>ins<description>

where <description> can be:

  • A literal nucleotide sequence:  NM_000548.5:c.1832_1833insATTGCCAGTAATGC
    → The inserted bases are spelled out directly.  This is the simplest case and
      gives us the sequence we need for BLAST alignment.

  • A numeric length only:  NM_000548.5:c.1832_1833ins50
    → The variant is known to be a 50-bp insertion, but the submitter did not
      provide the actual sequence.  We flag these as NEEDS_FETCH because the
      sequence might be retrievable from NCBI's Entrez API using the AlleleID.

  • Bracket notation with multiple segments:  NM_X:c.100_101ins[ATCG;GCTA]
    → HGVS allows complex insertions to be described as a series of segments
      separated by semicolons inside square brackets.  Each segment can be a
      nucleotide string, a repeat, or a reference to another sequence.

  • Inversion within an insertion:  ...ins[...inv...]
    → When part of the inserted material is an inverted copy of another region,
      HGVS uses the "inv" keyword.  These are too complex to extract a simple
      linear sequence from, so we flag them as COMPLEX.

  • A reference to another genomic locus:  ins(NC_012920.1:m.1-72)
    → This is the most exciting case for NUMT research!  It means the insertion
      is explicitly annotated as originating from the mitochondrial genome
      (NC_012920.1 is the human mtDNA RefSeq accession).  The coordinate range
      (e.g., m.1-72) tells us exactly which mitochondrial region was inserted.
      We flag these as DIRECT_MITO_REF — they are already confirmed NUMTs.

This script parses each of these formats, extracts what it can, and categorises
every variant so that downstream steps know how to handle it.
"""

# ─── Standard-library imports ────────────────────────────────────────────────
import argparse
import re
import sys
from datetime import datetime
from pathlib import Path

# ─── Third-party imports ─────────────────────────────────────────────────────
import pandas as pd


# ─── Constants: Project paths ────────────────────────────────────────────────

# Project root is two levels above this script (01_code/numt_pipeline/../../)
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

# Default input: the filtered TSV produced by Step 1
DEFAULT_INPUT = PROJECT_ROOT / "02_data" / "processed" / "clinvar_insertions_gt50bp.tsv"

# Default outputs
DEFAULT_OUTPUT = PROJECT_ROOT / "02_data" / "processed" / "clinvar_insertions_with_sequences.tsv"
DEFAULT_FASTA = PROJECT_ROOT / "02_data" / "processed" / "extracted_sequences.fasta"


# ─── Constants: HGVS regex patterns ─────────────────────────────────────────
#
# Each pattern is annotated with the HGVS format it targets and an example.
# Patterns are applied in order of specificity: more specific patterns first,
# with the broadest fallback last.  This avoids, e.g., a numeric-only pattern
# consuming what is actually a direct sequence string.

# Pattern 6 — Mitochondrial reference insertion
# Matches: ins(NC_012920.1:m.1-72)  or  ins(NC_012920.1:m.577-647)
# The parenthesised block references the human mitochondrial genome (RefSeq
# NC_012920.1) with a coordinate range on the mitochondrial sequence.
# This is the "smoking gun" for a confirmed NUMT — the submitter explicitly
# states that the inserted material comes from mtDNA.
#
# Regex breakdown:
#   ins\(               — literal "ins(" (the parenthesis is part of HGVS syntax)
#   (NC_012920\.\d+)    — capture group 1: the mitochondrial RefSeq accession
#                          (NC_012920 followed by a version number like .1)
#   :m\.                — literal ":m." — mitochondrial coordinate prefix
#   (\d+)[-_](\d+)      — capture groups 2 & 3: start and end positions of the
#                          mitochondrial region (separated by hyphen or underscore)
#   \)                  — closing parenthesis
RE_MITO_REF = re.compile(
    r"ins\((NC_012920\.\d+):m\.(\d+)[-_](\d+)\)"
)

# Pattern 5 — Inversion within an insertion
# Matches: any HGVS name containing "ins" followed at some point by "inv"
# Example: NM_X:c.100_101ins[ATCG;200_250inv]
# Inversions create a reversed-complement segment, making it impossible to
# extract a simple linear sequence without resolving the inversion coordinates
# against the reference.  We flag these as COMPLEX.
#
# Regex breakdown:
#   ins       — the insertion keyword
#   .*        — any characters in between (greedy)
#   inv       — the inversion keyword somewhere after ins
RE_INV_IN_INS = re.compile(
    r"ins.*inv"
)

# Pattern 4 — Bracket notation (complex/multi-segment insertion)
# Matches: NM_X:c.100_101ins[ATCG;GCTA]  or  ins[ATCGATCG]
# Bracket notation can contain multiple semicolon-separated segments.
# Each segment that looks like a nucleotide string (only A/C/G/T/N) is
# extracted; segments that are coordinates or references are skipped.
#
# Regex breakdown:
#   ins\[       — literal "ins[" — opening bracket notation
#   ([^\]]+)    — capture group 1: everything inside the brackets (up to the
#                 closing bracket), which may contain semicolons separating
#                 multiple sequence segments
#   \]          — closing bracket
RE_BRACKET_INS = re.compile(
    r"ins\[([^\]]+)\]"
)

# Pattern 1 & 2 — Direct nucleotide sequence (coding or genomic coordinates)
# Matches: NM_000548.5:c.1832_1833insATTGCCAGTAATGC
#          NM_000071.3(CBS):c.845_846insGAAGGGT... (p.His123fs)  ← protein suffix
# The inserted sequence is spelled out after the "ins" keyword.
# The original pattern used $ (end of string) which failed when ClinVar appends
# a protein-change annotation like " (p.His1244fs)" after the nucleotide sequence.
# Fix: allow optional trailing whitespace + parenthesised annotation.
#
# Regex breakdown:
#   ins             — the insertion keyword in HGVS
#   ([ACGTNacgtn]+) — capture group 1: one or more nucleotide characters
#   (?:\s*\(p\..*)?$ — optional protein annotation suffix, then end of string
RE_DIRECT_SEQ = re.compile(
    r"ins([ACGTNacgtn]+)(?:\s*\(p\..*)?$"
)

# Pattern 1b — Deletion-insertion (delins) with nucleotide sequence
# Matches: NM_005502.4(ABCA1):c.1584_1597delinsCGGGCGTGGTGGCAGGAGCTG...
#          NM_000251.3(MSH2):c.243_273delinsCTGACAA...  (p.Xxx)
# ClinVar uses "delins" (deletion-insertion) when existing bases are replaced
# by a new sequence.  The inserted sequence follows "delins" and may also be
# followed by a protein annotation suffix.
# These were previously all UNPARSEABLE because only "ins" was handled.
#
# Regex breakdown:
#   delins          — deletion-insertion keyword
#   ([ACGTNacgtn]+) — capture group 1: the replacing nucleotide sequence
#   (?:\s*\(p\..*)?$ — optional protein annotation suffix
RE_DELINS = re.compile(
    r"delins([ACGTNacgtn]+)(?:\s*\(p\..*)?$",
    re.IGNORECASE,
)

# Pattern 3 — Length-only insertion (no sequence provided)
# Matches: NM_000548.5:c.1832_1833ins50
#          NM_X:c.100_101ins1234
# The number after "ins" indicates the length of the insertion in base pairs,
# but the actual nucleotide sequence was not submitted to ClinVar.  These
# variants need their sequences fetched from NCBI Entrez (a later pipeline step).
#
# Regex breakdown:
#   ins       — the insertion keyword
#   (\d+)     — capture group 1: one or more digits (the insertion length)
#   $         — end of string (ensures it's purely a number, not followed by
#               nucleotide characters which would make it Pattern 1/2)
RE_LENGTH_ONLY = re.compile(
    r"ins(\d+)$"
)


# ─── Helper: timestamped logging ────────────────────────────────────────────

def log(message: str) -> None:
    """Print a message prefixed with the current timestamp."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")


# ─── Core parser: classify and extract from a single HGVS name ──────────────

def parse_hgvs_insertion(name: str) -> dict:
    """
    Parse an HGVS variant name and attempt to extract the inserted sequence.

    The function applies regex patterns in order of specificity:
      1. Mitochondrial reference  → DIRECT_MITO_REF
      2. Inversion in insertion   → COMPLEX
      3. Bracket notation         → DIRECT_SEQUENCE (if all segments are nucleotides)
                                    or COMPLEX (if some segments are not)
      4. Direct nucleotide string → DIRECT_SEQUENCE
      5. Length-only              → LENGTH_ONLY
      6. Fallback                 → UNPARSEABLE

    Parameters
    ----------
    name : str
        The HGVS variant name from ClinVar's Name column.

    Returns
    -------
    dict
        Keys: extracted_sequence, extraction_status, sequence_length_extracted,
              mito_ref_range.
    """
    # Default return structure — will be updated by whichever pattern matches
    result = {
        "extracted_sequence": None,
        "extraction_status": "UNPARSEABLE",
        "sequence_length_extracted": None,
        "mito_ref_range": None,
    }

    # Guard against missing or non-string values (NaN from pandas)
    if not isinstance(name, str) or not name.strip():
        return result

    # ── Check 1: Mitochondrial reference insertion ───────────────────────
    # This is the most specific and most interesting case for NUMT research.
    # If the HGVS name explicitly references the mitochondrial genome, this
    # variant is already annotated as a NUMT by the ClinVar submitter.
    match = RE_MITO_REF.search(name)
    if match:
        accession = match.group(1)          # e.g. NC_012920.1
        start = match.group(2)              # e.g. 1
        end = match.group(3)                # e.g. 72
        result["extraction_status"] = "DIRECT_MITO_REF"
        result["mito_ref_range"] = f"{start}-{end}"
        # We cannot extract the sequence directly from the name — it would
        # need to be looked up from the mitochondrial reference.  But the
        # coordinate range tells us exactly where to look.
        result["sequence_length_extracted"] = abs(int(end) - int(start)) + 1
        return result

    # ── Check 2: Inversion within an insertion ───────────────────────────
    # If "inv" appears after "ins", the insertion involves an inverted segment.
    # Resolving this requires knowing the reference sequence and applying the
    # reverse-complement operation — beyond simple regex extraction.
    if RE_INV_IN_INS.search(name):
        result["extraction_status"] = "COMPLEX"
        return result

    # ── Check 3: Bracket notation ────────────────────────────────────────
    # Bracket notation can encode multiple segments: ins[ATCG;GCTA]
    # We attempt to extract and concatenate all nucleotide-only segments.
    match = RE_BRACKET_INS.search(name)
    if match:
        inner = match.group(1)  # everything inside the brackets
        segments = inner.split(";")

        # Check each segment: is it a pure nucleotide string?
        nucleotide_segments = []
        all_are_nucleotides = True
        for seg in segments:
            seg = seg.strip()
            if re.fullmatch(r"[ACGTNacgtn]+", seg):
                nucleotide_segments.append(seg.upper())
            else:
                # This segment is a coordinate reference, repeat notation, etc.
                all_are_nucleotides = False

        if nucleotide_segments and all_are_nucleotides:
            # All segments are nucleotide strings — concatenate them
            seq = "".join(nucleotide_segments)
            result["extracted_sequence"] = seq
            result["extraction_status"] = "DIRECT_SEQUENCE"
            result["sequence_length_extracted"] = len(seq)
        else:
            # Some segments could not be parsed as nucleotides — flag as COMPLEX
            # rather than silently discarding information.
            result["extraction_status"] = "COMPLEX"
        return result

    # ── Check 4a: Deletion-insertion with nucleotide sequence ────────────
    # delinsATCG... — same logic as DIRECT_SEQUENCE, different HGVS keyword.
    # Must be checked BEFORE RE_DIRECT_SEQ because "delins" contains "ins"
    # and would otherwise partially match the ins pattern.
    match = RE_DELINS.search(name)
    if match:
        seq = match.group(1).upper()
        result["extracted_sequence"] = seq
        result["extraction_status"] = "DIRECT_SEQUENCE"
        result["sequence_length_extracted"] = len(seq)
        return result

    # ── Check 4b: Direct nucleotide sequence ─────────────────────────────
    # ins[ACGTNacgtn]+ optionally followed by a protein annotation (p.Xxx).
    match = RE_DIRECT_SEQ.search(name)
    if match:
        seq = match.group(1).upper()
        result["extracted_sequence"] = seq
        result["extraction_status"] = "DIRECT_SEQUENCE"
        result["sequence_length_extracted"] = len(seq)
        return result

    # ── Check 5: Length-only insertion ────────────────────────────────────
    # The HGVS name ends with "ins" followed by a number — the insertion
    # length is known but the sequence was not submitted.
    match = RE_LENGTH_ONLY.search(name)
    if match:
        length = int(match.group(1))
        result["extraction_status"] = "LENGTH_ONLY"
        result["sequence_length_extracted"] = length
        return result

    # ── Fallback: UNPARSEABLE ────────────────────────────────────────────
    # The HGVS name does not match any known insertion pattern.  This could
    # be a deletion, substitution, duplication, or a non-standard notation.
    return result


# ─── Step 1: Load the filtered ClinVar TSV ──────────────────────────────────

def load_input(input_path: Path) -> pd.DataFrame:
    """
    Read the filtered ClinVar insertions TSV from Step 1.

    Parameters
    ----------
    input_path : Path
        Path to the tab-separated file produced by 01_clinvar_fetch.py.

    Returns
    -------
    pd.DataFrame
        DataFrame with all columns from Step 1.
    """
    log(f"Reading input file: {input_path}")

    if not input_path.exists():
        log(f"FATAL — Input file not found: {input_path}")
        log("Have you run 01_clinvar_fetch.py first?")
        sys.exit(1)

    df = pd.read_csv(input_path, sep="\t", dtype=str, low_memory=False)
    log(f"Loaded {len(df):,} variant records with {len(df.columns)} columns")
    return df


# ─── Step 2: Apply the HGVS parser across all rows ──────────────────────────

def extract_sequences(df: pd.DataFrame) -> pd.DataFrame:
    """
    Parse the HGVS Name column for every variant and add extraction columns.

    This applies `parse_hgvs_insertion()` row by row, then merges the results
    back into the DataFrame as four new columns:
      - extracted_sequence
      - extraction_status
      - sequence_length_extracted
      - mito_ref_range

    Parameters
    ----------
    df : pd.DataFrame
        The filtered ClinVar DataFrame (must contain a 'Name' column).

    Returns
    -------
    pd.DataFrame
        The input DataFrame augmented with extraction columns.
    """
    log("Parsing HGVS names to extract inserted sequences …")

    # Verify the required column exists
    if "Name" not in df.columns:
        log("FATAL — Input DataFrame is missing the 'Name' column.")
        log("Available columns: " + ", ".join(df.columns[:20]))
        sys.exit(1)

    # Apply the parser to each Name value and expand the resulting dicts
    # into a DataFrame with one column per dict key.
    parsed = df["Name"].apply(parse_hgvs_insertion).apply(pd.Series)

    # Concatenate the new columns with the original DataFrame
    df_out = pd.concat([df, parsed], axis=1)

    # Convert sequence_length_extracted to numeric (it may be int or None)
    df_out["sequence_length_extracted"] = pd.to_numeric(
        df_out["sequence_length_extracted"], errors="coerce"
    )

    n_extracted = (df_out["extraction_status"] == "DIRECT_SEQUENCE").sum()
    log(f"Sequence extraction complete: {n_extracted:,} direct sequences obtained")

    return df_out


# ─── Step 3: Save results ───────────────────────────────────────────────────

def save_tsv(df: pd.DataFrame, output_path: Path) -> None:
    """Write the augmented DataFrame to a tab-separated file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep="\t", index=False)
    log(f"Saved {len(df):,} variants to: {output_path}")


def save_fasta(df: pd.DataFrame, fasta_path: Path) -> None:
    """
    Write a FASTA file containing only DIRECT_SEQUENCE entries.

    Each entry uses the ClinVar AlleleID as the sequence identifier, which
    provides a stable, unique key for linking back to ClinVar records.
    FASTA format:
        >AlleleID_XXXXX
        ATCGATCG...

    Parameters
    ----------
    df : pd.DataFrame
        The augmented DataFrame (must contain 'AlleleID', 'extracted_sequence',
        and 'extraction_status' columns).
    fasta_path : Path
        Output path for the FASTA file.
    """
    fasta_path.parent.mkdir(parents=True, exist_ok=True)

    # Filter to only rows where we have a directly extracted sequence
    direct = df[df["extraction_status"] == "DIRECT_SEQUENCE"].copy()

    if direct.empty:
        log("WARNING — No DIRECT_SEQUENCE entries found; FASTA file will be empty.")

    # Determine the ID column — ClinVar uses "AlleleID" or "#AlleleID"
    # (the variant_summary.txt header sometimes includes a '#' prefix)
    id_col = None
    for candidate in ["AlleleID", "#AlleleID"]:
        if candidate in df.columns:
            id_col = candidate
            break

    if id_col is None:
        log("WARNING — No AlleleID column found. Using row index as FASTA IDs.")

    n_written = 0
    with open(fasta_path, "w") as fh:
        for idx, row in direct.iterrows():
            seq = row["extracted_sequence"]
            if seq is None or (isinstance(seq, float)):
                continue

            # Build the FASTA header
            if id_col is not None and pd.notna(row.get(id_col)):
                seq_id = f"AlleleID_{row[id_col]}"
            else:
                seq_id = f"row_{idx}"

            # Write FASTA entry — sequence wrapped at 80 characters per line
            # (standard FASTA convention for readability)
            fh.write(f">{seq_id}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")

            n_written += 1

    log(f"Wrote {n_written:,} sequences to FASTA: {fasta_path}")


# ─── Step 4: Print summary breakdown ────────────────────────────────────────

def print_summary(df: pd.DataFrame) -> None:
    """
    Display a human-readable breakdown of extraction results by status.

    This helps the researcher understand how many variants have usable
    sequences for the next BLAST step versus how many need additional work
    (e.g., fetching from Entrez, manual review).
    """
    print("\n" + "=" * 72)
    print("  Sequence Extraction — Summary")
    print("=" * 72)

    print(f"\n  Total variants processed: {len(df):>10,}")

    # ── Breakdown by extraction_status ───────────────────────────────────
    print(f"\n  {'Breakdown by extraction_status':─<50}")

    status_counts = df["extraction_status"].value_counts()
    for status, count in status_counts.items():
        pct = count / len(df) * 100 if len(df) > 0 else 0
        print(f"    {status:<30} {count:>8,}  ({pct:5.1f}%)")

    # ── Sequence length statistics for DIRECT_SEQUENCE entries ───────────
    direct = df[df["extraction_status"] == "DIRECT_SEQUENCE"]
    if not direct.empty:
        lengths = pd.to_numeric(direct["sequence_length_extracted"], errors="coerce")
        lengths = lengths.dropna()
        if not lengths.empty:
            print(f"\n  {'DIRECT_SEQUENCE length statistics':─<50}")
            print(f"    Min length (bp)     : {lengths.min():>10,.0f}")
            print(f"    Max length (bp)     : {lengths.max():>10,.0f}")
            print(f"    Mean length (bp)    : {lengths.mean():>10,.1f}")
            print(f"    Median length (bp)  : {lengths.median():>10,.1f}")

    # ── DIRECT_MITO_REF entries — the confirmed NUMTs ────────────────────
    mito = df[df["extraction_status"] == "DIRECT_MITO_REF"]
    if not mito.empty:
        print(f"\n  {'DIRECT_MITO_REF entries (confirmed mtDNA references)':─<50}")
        for _, row in mito.iterrows():
            name_short = str(row.get("Name", "?"))[:60]
            mito_range = row.get("mito_ref_range", "?")
            print(f"    {name_short:<62} mtDNA range: {mito_range}")

    print("\n" + "=" * 72 + "\n")


# ─── CLI entry point ────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    """Define and parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Step 2 of the NUMT pipeline: parse HGVS insertion names to extract "
            "inserted nucleotide sequences.  Sequences are needed for downstream "
            "BLAST alignment against the human mitochondrial genome (rCRS)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help=(
            "Path to the filtered ClinVar TSV from Step 1 "
            f"(default: {DEFAULT_INPUT.relative_to(PROJECT_ROOT)})"
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=(
            "Path for the output TSV with extracted sequences "
            f"(default: {DEFAULT_OUTPUT.relative_to(PROJECT_ROOT)})"
        ),
    )
    parser.add_argument(
        "--fasta",
        type=Path,
        default=DEFAULT_FASTA,
        help=(
            "Path for the output FASTA file of extracted sequences "
            f"(default: {DEFAULT_FASTA.relative_to(PROJECT_ROOT)})"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Orchestrate the sequence extraction workflow."""
    args = parse_args()

    log("=" * 60)
    log("NUMT Pipeline — Step 2: Sequence Extraction from HGVS")
    log("=" * 60)

    # Step 1 — Load the filtered ClinVar TSV
    df = load_input(args.input)

    # Step 2 — Parse HGVS names and extract sequences
    df = extract_sequences(df)

    # Step 3 — Save the augmented TSV and the FASTA file
    save_tsv(df, args.output)
    save_fasta(df, args.fasta)

    # Step 4 — Print human-readable summary
    print_summary(df)

    log("Done.")


if __name__ == "__main__":
    main()
