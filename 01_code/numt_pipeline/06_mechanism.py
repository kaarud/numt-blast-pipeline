#!/usr/bin/env python3
"""
06_mechanism.py — Investigate the molecular mechanism of NUMT insertion.

WHY THIS STEP MATTERS
=====================
In Steps 4–5 we identified insertion variants that are of mitochondrial origin
(confirmed NUMTs and novel NUMT candidates).  Knowing THAT a fragment of mtDNA
ended up in the nuclear genome is only half the story.  The other half is HOW.

Understanding the insertion mechanism has both scientific and clinical value:
  - It strengthens the NUMT hypothesis: if the flanking sequences show hallmarks
    of a known mechanism (e.g., target-site duplications typical of L1-mediated
    retrotransposition), the insertion is far more likely to be genuine.
  - It has implications for recurrence risk: de novo L1-mediated insertions are
    ongoing in the human population, while NAHR events are tied to specific
    genomic architectures (e.g., segmental duplications).
  - It informs detection strategy: TPRT-mediated NUMTs create target-site
    duplications that are detectable by paired-end and long-read sequencing.

THE THREE MAIN NUMT INSERTION MECHANISMS
========================================

1. NON-ALLELIC HOMOLOGOUS RECOMBINATION (NAHR)
   -------------------------------------------
   NAHR occurs when two non-allelic DNA sequences that share high similarity
   ("direct repeats" or "segmental duplications") misalign during meiosis or
   mitotic DNA repair, leading to illegitimate recombination.  The result is a
   deletion, duplication, or — in our case — insertion of foreign DNA between
   the two repeat copies.

   Sequence signature:
   - Long direct repeats (>100 bp) flanking the insertion site
   - The repeats are in the same orientation (direct, not inverted)
   - The inserted DNA replaces the sequence between the repeat copies

   Example at the DNA level:
     NORMAL:   5'---[REPEAT_A]---100bp---[REPEAT_A']---3'
     AFTER:    5'---[REPEAT_A]---mtDNA_FRAGMENT---[REPEAT_A']---3'

   For computational purposes, we search for exact repeats >= 10 bp within
   the flanking region as a relaxed proxy for this mechanism.  True NAHR
   typically requires much longer repeats (hundreds to thousands of bp), but
   shorter repeats can facilitate microhomology-mediated recombination events.

2. NON-HOMOLOGOUS END JOINING (NHEJ) / MICROHOMOLOGY-MEDIATED END JOINING (MMEJ)
   -------------------------------------------------------------------------------
   NHEJ is the primary DNA double-strand break (DSB) repair pathway in human
   cells.  It ligates broken DNA ends with little or no regard for sequence
   homology.  A subtype, MMEJ (also called "alternative end joining"), uses
   short stretches of microhomology (2–25 bp) at the junction to align the
   broken ends before ligation.

   When a DSB occurs in a nuclear gene, and a fragment of mitochondrial DNA
   is captured and inserted during the repair process, the junctions between
   nuclear and mitochondrial DNA typically show microhomology: a few base
   pairs that are identical between the end of the nuclear sequence and the
   beginning of the mtDNA fragment (or vice versa).

   Sequence signature:
   - Microhomology of 2–25 bp at one or both insertion junctions
   - No large flanking repeats
   - No target-site duplications
   - May have small deletions at the insertion site

   Example at the DNA level (junction):
     LEFT FLANK:  ...ATCGATCG|ATCG        <-- last 4 bp of left flank
     MT INSERT:              ATCG|AATTCC.. <-- first 4 bp of insertion
                              ^^^^
                         microhomology (4 bp)

   NHEJ: 0–1 bp of overlap (blunt or near-blunt ligation)
   MMEJ: 2–25 bp of microhomology (most commonly 4–12 bp)

3. TARGET-PRIMED REVERSE TRANSCRIPTION (TPRT) VIA L1 RETROTRANSPOSONS
   -------------------------------------------------------------------
   This is the most common mechanism for RECENT NUMT insertions in the human
   genome.  The L1 (LINE-1) retrotransposon encodes two proteins: ORF1p (an
   RNA-binding protein) and ORF2p (which has both endonuclease and reverse
   transcriptase activity).  Although L1 normally mobilises its own RNA, the
   ORF2p endonuclease can nick chromosomal DNA at a consensus site, and the
   reverse transcriptase can then copy nearby RNA — including mitochondrial
   RNA that has leaked into the cytoplasm — into the nick.  This "hijacking"
   of the L1 machinery is called TPRT.

   Sequence signatures (ALL should be present for a confident call):
   - L1 endonuclease cleavage consensus: TTTT/A (or degenerate [TC]TTT[TA])
     within ~20 bp upstream of the insertion site.  The enzyme nicks the
     bottom strand at a T-rich sequence.
   - Target-site duplications (TSDs): identical sequences of 7–20 bp that
     appear both immediately upstream AND immediately downstream of the
     insertion.  TSDs are created when the L1 endonuclease makes staggered
     cuts on the two DNA strands, and the resulting single-strand gaps are
     filled in after insertion.
   - Antisense orientation: TPRT-mediated NUMTs are frequently inserted in
     the antisense orientation relative to the original mtDNA sequence,
     because reverse transcription initiates from the 3' end of the RNA.
   - A poly-A tail at the 3' end of the insertion (from polyadenylation of
     the mitochondrial transcript).

   Example at the DNA level:
     5'---TTTTATTTT[A]---[TSD]---[mtDNA_antisense]---[TSD]---3'
          ^^^^^^^^^^            target-site duplications
          L1 EN site

LIMITATIONS OF COMPUTATIONAL MECHANISM PREDICTION
==================================================
This script performs SEQUENCE-BASED PREDICTIONS ONLY.  Important caveats:

  1. Flanking sequences are fetched from the REFERENCE genome (hg38).  The
     actual patient genome may differ — especially at the insertion site.
     Only long-read sequencing of the patient DNA can reveal the true
     junction sequences.

  2. Microhomology of a few base pairs can occur by chance.  Short
     microhomologies (2–3 bp) are not strong evidence for MMEJ.  We require
     >= 4 bp for a meaningful call.

  3. TSD detection is approximate: we look for short identical sequences in
     the immediate flanking region, but the exact breakpoint is often not
     precisely defined in ClinVar records.

  4. L1 endonuclease site consensus is degenerate and relatively common in
     AT-rich regions.  Its presence alone is suggestive but not conclusive.

  5. Experimental validation (PCR across the junction, Sanger sequencing,
     or long-read sequencing) is required to confirm any mechanism prediction.
"""

# ─── Standard-library imports ────────────────────────────────────────────────
import argparse
import json
import re
import sys
import time
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path

# ─── Third-party imports ─────────────────────────────────────────────────────
try:
    import requests
    REQUESTS_AVAILABLE = True
except ImportError:
    REQUESTS_AVAILABLE = False

try:
    from Bio import Entrez, SeqIO
    from Bio.Seq import Seq
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

try:
    import pymongo
    from pymongo import UpdateOne
    from pymongo.errors import ConnectionFailure, ServerSelectionTimeoutError
    PYMONGO_AVAILABLE = True
except ImportError:
    PYMONGO_AVAILABLE = False


# ─── Constants: Project paths ────────────────────────────────────────────────

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

# MongoDB connection
MONGO_HOST = "localhost"
MONGO_PORT = 27017
MONGO_DB = "numt_db"
MONGO_COLLECTION = "clinvar_variants"

# Output paths
FLANKING_DIR = PROJECT_ROOT / "02_data" / "processed" / "flanking_sequences"
MECHANISM_JSON = PROJECT_ROOT / "02_data" / "processed" / "mechanism_analysis.json"
MECHANISM_TSV = PROJECT_ROOT / "02_data" / "processed" / "mechanism_summary.tsv"

# ─── NCBI chromosome accession map ──────────────────────────────────────────
# Used to build the correct NCBI accession for Entrez efetch calls.
# ClinVar records report chromosome names like "1", "2", ..., "X", "Y".
# Entrez efetch requires the RefSeq accession (e.g., NC_000001.11).

CHR_TO_ACCESSION = {
    "1": "NC_000001.11",   "2": "NC_000002.12",   "3": "NC_000003.12",
    "4": "NC_000004.12",   "5": "NC_000005.10",   "6": "NC_000006.12",
    "7": "NC_000007.14",   "8": "NC_000008.11",   "9": "NC_000009.12",
    "10": "NC_000010.11",  "11": "NC_000011.10",  "12": "NC_000012.12",
    "13": "NC_000013.11",  "14": "NC_000014.9",   "15": "NC_000015.10",
    "16": "NC_000016.10",  "17": "NC_000017.11",  "18": "NC_000018.10",
    "19": "NC_000019.10",  "20": "NC_000020.11",  "21": "NC_000021.9",
    "22": "NC_000022.11",  "X": "NC_000023.11",   "Y": "NC_000024.10",
}

# ─── L1 endonuclease consensus ──────────────────────────────────────────────
# The L1 ORF2p endonuclease recognises and nicks the bottom (antisense) strand
# at a T-rich consensus.  On the top (sense) strand this appears as:
#   TTTT/A  (canonical)  or more broadly  [TC]TTT[TA]
# The degenerate regex below captures both canonical and near-canonical sites.

L1_EN_CONSENSUS_REGEX = re.compile(r"[TC]TTT[TA]")

# Hardcoded test AlleleID (included for pipeline testing)
TEST_ALLELE_ID = 1234567

# ─── UCSC REST API ───────────────────────────────────────────────────────────
UCSC_API = "https://api.genome.ucsc.edu/getData/track"

# ─── Bolshoy dinucleotide roll angles (degrees) ──────────────────────────────
# Used to compute a local DNA curvature proxy over flanking sequences.
# Source: Bolshoy et al. (1993) PNAS 90:2312–2316.
# Interpretation: higher mean roll angle ≈ more curved DNA.
# Threshold for "high curvature": > 1.5 °/bp.
BOLSHOY_ROLL: dict[str, float] = {
    "AA": 1.1, "TT": 1.1, "AT": 4.0, "TA": 0.5,
    "CA": 1.4, "TG": 1.4, "GT": 1.5, "AC": 1.5,
    "AG": 1.1, "CT": 1.1, "GA": 1.1, "TC": 1.1,
    "GG": 1.9, "CC": 1.9, "GC": 0.5, "CG": 3.4,
}

# ─── Mechanism interpretation text (for statistics sheet) ───────────────────
MECHANISM_INTERPRETATIONS: dict[str, str] = {
    "NHEJ": (
        "Non-Homologous End Joining (NHEJ) ligates double-strand break (DSB) ends with "
        "minimal sequence homology. mtDNA fragments captured during nuclear DSB repair "
        "show blunt or near-blunt junctions (0–3 bp). NHEJ is the dominant DSB repair "
        "pathway in G1/G0 cells and at open chromatin loci (Tsuji et al. 2012; "
        "Xue et al. 2023)."
    ),
    "MMEJ/MMBIR": (
        "Microhomology-Mediated End Joining (MMEJ, alt-NHEJ) or Microhomology-Mediated "
        "Break-Induced Replication (MMBIR) — both mechanisms use short overlapping sequences "
        "to mediate DSB repair or template switching during replication. Junctions with "
        "microhomology ≥4 bp suggest capture of mtDNA during MMEJ- or MMBIR-mediated repair. "
        "MMEJ and MMBIR are computationally indistinguishable from reference sequences alone: "
        "their distinction requires patient-derived long-read sequencing to resolve templated "
        "insertions or complex junction architecture (Hastings et al. 2009; Bíró et al. 2024). "
        "Genomic context: open chromatin, AT-rich, high DNA curvature (Tsuji et al. 2012)."
    ),
    "NAHR": (
        "Non-Allelic Homologous Recombination (NAHR) requires direct repeats flanking "
        "the insertion site. Misalignment of repeat copies during meiosis or DSB repair "
        "allows mtDNA integration between them. Typically associated with larger NUMTs "
        "(>5 kb) near segmental duplications. Short direct repeats (10–15 bp) are "
        "computationally detected but may occur by chance "
        "(Hazkani-Covo et al. 2003; Tsuji et al. 2012)."
    ),
    "TPRT": (
        "Target-Primed Reverse Transcription (TPRT) uses LINE-1 (L1) ORF2p endonuclease "
        "and reverse transcriptase to insert mtRNA into nicked chromosomal DNA. "
        "Hallmarks: target-site duplications (7–20 bp), poly-A tail, L1 EN site upstream, "
        "antisense orientation relative to rCRS. TPRT is the most common mechanism for "
        "recent NUMT insertions in the human genome (Xue et al. 2023)."
    ),
    "UNKNOWN": (
        "No definitive sequence signatures (microhomology, TSDs, direct repeats, L1 EN "
        "site) were detected in the 500 bp flanking window. The mechanism could not be "
        "determined computationally. Possible explanations: imprecise ClinVar breakpoints, "
        "complex rearrangement, or signatures outside the analysed window."
    ),
}


# ─── Helper: timestamped logging ────────────────────────────────────────────

def log(message: str) -> None:
    """Print a message prefixed with the current timestamp."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")


# ─── Helper: reverse complement ─────────────────────────────────────────────

def reverse_complement(seq: str) -> str:
    """
    Return the reverse complement of a DNA sequence string.

    Uses Biopython if available, otherwise a manual complement table.
    This is needed to determine the orientation of the NUMT insertion
    relative to the mitochondrial genome.
    """
    if BIOPYTHON_AVAILABLE:
        return str(Seq(seq).reverse_complement())
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(complement)[::-1]


# ─── 1. Fetch NUMT candidates from MongoDB ──────────────────────────────────

def fetch_candidates(collection) -> list:
    """
    Query MongoDB for confirmed NUMTs and novel NUMT candidates.

    These are variants classified in Step 5 (05_blast_genome.py) as either:
      - CONFIRMED_NUMT:       mtDNA hit + nuclear hit at the ClinVar position
      - NOVEL_NUMT_CANDIDATE: mtDNA hit but no nuclear hit (not in reference)

    We also include AlleleID 1234567 as a hardcoded test case (if present)
    to ensure the pipeline can be exercised even when no real candidates exist.

    Returns
    -------
    list of dict
        MongoDB documents for each candidate variant.
    """
    query = {
        "$or": [
            {"numt_classification": {"$in": ["CONFIRMED_NUMT", "NOVEL_NUMT_CANDIDATE"]}},
            {"AlleleID": TEST_ALLELE_ID},
        ]
    }
    candidates = list(collection.find(query))
    log(f"Fetched {len(candidates)} NUMT candidate(s) from MongoDB.")

    # Deduplicate in case the test AlleleID is also a confirmed/novel candidate
    seen_ids = set()
    unique = []
    for doc in candidates:
        aid = doc.get("AlleleID")
        if aid not in seen_ids:
            seen_ids.add(aid)
            unique.append(doc)

    if len(unique) < len(candidates):
        log(f"  (Deduplicated to {len(unique)} unique candidates.)")

    return unique


# ─── 2. Fetch flanking sequences from NCBI ──────────────────────────────────

def fetch_flanking_sequence(
    chrom: str,
    start: int,
    stop: int,
    allele_id: int,
    flank_size: int = 500,
    email: str = "numt_pipeline@example.com",
) -> dict | None:
    """
    Fetch flanking genomic sequence around an insertion site via NCBI Entrez.

    For a variant at chr16:2000000-2000000, this fetches:
      - LEFT  flank: chr16:1999500-1999999  (500 bp upstream)
      - RIGHT flank: chr16:2000001-2000500  (500 bp downstream)

    Results are cached to disk so re-runs do not re-query NCBI.

    Parameters
    ----------
    chrom : str
        Chromosome name as reported in ClinVar (e.g., "16", "X").
    start : int
        1-based start coordinate of the insertion site.
    stop : int
        1-based stop coordinate of the insertion site.
    allele_id : int
        ClinVar AlleleID, used for caching and logging.
    flank_size : int
        Number of base pairs to fetch on each side (default 500).
    email : str
        Email for NCBI Entrez (required by NCBI usage policy).

    Returns
    -------
    dict or None
        {"left_flank": str, "right_flank": str, "full_flank": str} or None
        if the fetch failed.
    """
    if not BIOPYTHON_AVAILABLE:
        log("ERROR: Biopython is required for Entrez efetch. pip install biopython")
        return None

    # Check for cached result
    cache_file = FLANKING_DIR / f"{allele_id}_flanking.fasta"
    if cache_file.exists():
        log(f"  AlleleID {allele_id}: using cached flanking sequence.")
        return _parse_cached_flanking(cache_file, flank_size)

    # Resolve chromosome to NCBI accession
    chrom_str = str(chrom).replace("chr", "")
    accession = CHR_TO_ACCESSION.get(chrom_str)
    if not accession:
        log(f"  WARNING: AlleleID {allele_id}: unknown chromosome '{chrom}', skipping.")
        return None

    Entrez.email = email

    # Calculate coordinates for left and right flanks.
    # NCBI Entrez efetch uses 1-based coordinates with seq_start and seq_stop.
    left_start = max(1, start - flank_size)
    left_stop = max(1, start - 1)
    right_start = stop + 1
    right_stop = stop + flank_size

    left_seq = _entrez_fetch_region(accession, left_start, left_stop, allele_id, "left")
    if left_seq is None:
        return None

    # Respect NCBI rate limits: no more than 3 requests per second
    time.sleep(0.4)

    right_seq = _entrez_fetch_region(accession, right_start, right_stop, allele_id, "right")
    if right_seq is None:
        return None

    # Cache to FASTA file
    FLANKING_DIR.mkdir(parents=True, exist_ok=True)
    with open(cache_file, "w") as fh:
        fh.write(f">AlleleID_{allele_id}_left_flank {accession}:{left_start}-{left_stop}\n")
        fh.write(left_seq + "\n")
        fh.write(f">AlleleID_{allele_id}_right_flank {accession}:{right_start}-{right_stop}\n")
        fh.write(right_seq + "\n")

    log(f"  AlleleID {allele_id}: fetched {len(left_seq)} bp left + {len(right_seq)} bp right flank.")

    return {
        "left_flank": left_seq.upper(),
        "right_flank": right_seq.upper(),
        "full_flank": (left_seq + right_seq).upper(),
    }


def _entrez_fetch_region(
    accession: str, start: int, stop: int, allele_id: int, label: str
) -> str | None:
    """
    Fetch a specific genomic region via Entrez efetch.

    Parameters
    ----------
    accession : str
        NCBI RefSeq accession (e.g., NC_000016.10).
    start : int
        1-based start coordinate.
    stop : int
        1-based stop coordinate.
    allele_id : int
        For logging purposes.
    label : str
        "left" or "right", for logging.

    Returns
    -------
    str or None
        The nucleotide sequence, or None on failure.
    """
    try:
        handle = Entrez.efetch(
            db="nuccore",
            id=accession,
            rettype="fasta",
            retmode="text",
            seq_start=start,
            seq_stop=stop,
        )
        record = SeqIO.read(handle, "fasta")
        handle.close()
        return str(record.seq).upper()
    except Exception as exc:
        log(f"  WARNING: AlleleID {allele_id}: Entrez efetch failed for {label} flank: {exc}")
        return None


def _parse_cached_flanking(cache_file: Path, flank_size: int) -> dict | None:
    """
    Parse a cached flanking-sequence FASTA file.

    The file contains two records: left_flank and right_flank.

    Returns
    -------
    dict or None
        {"left_flank": str, "right_flank": str, "full_flank": str}
    """
    try:
        records = list(SeqIO.parse(str(cache_file), "fasta"))
        if len(records) < 2:
            log(f"  WARNING: cached file {cache_file} has < 2 records, re-fetching.")
            cache_file.unlink()
            return None
        left_seq = str(records[0].seq).upper()
        right_seq = str(records[1].seq).upper()
        return {
            "left_flank": left_seq,
            "right_flank": right_seq,
            "full_flank": left_seq + right_seq,
        }
    except Exception as exc:
        log(f"  WARNING: failed to parse cached file {cache_file}: {exc}")
        return None


# ─── 3. Microhomology analysis ──────────────────────────────────────────────
#
# WHAT IS MICROHOMOLOGY?
# ----------------------
# When a double-strand break (DSB) is repaired by MMEJ (microhomology-mediated
# end joining), the two broken ends are aligned using a short stretch of
# identical sequence (the "microhomology") before ligation.  If a mitochondrial
# DNA fragment is captured during this repair, the junction between the nuclear
# DNA and the mtDNA will show microhomology.
#
# Concretely, we compare the last N bp of the LEFT flank (nuclear DNA upstream
# of the insertion) to the first N bp of the RIGHT flank (nuclear DNA downstream
# of the insertion).  If these overlap by >= 4 bp, it suggests MMEJ/NHEJ.
#
# Example:
#   Left flank ends with:   ...GCATGCATCG
#   Right flank starts with:     GCATCGATTT...
#                                 ^^^^
#                            4 bp microhomology ("CATC")
#
# Wait — actually, the comparison is between the END of the left flank and
# the START of the right flank.  We slide a window of size 2..25 and check
# if the last W bp of the left flank exactly match the first W bp of the
# right flank.  The longest exact match is the microhomology.

def analyse_microhomology(
    left_flank: str, right_flank: str, min_length: int = 4, max_length: int = 25
) -> dict:
    """
    Detect microhomology at the insertion junction.

    Compares the last N bp of the left flank to the first N bp of the right
    flank for N = max_length down to min_length.  Returns the longest exact
    match found.

    Parameters
    ----------
    left_flank : str
        Upstream flanking sequence (e.g., 500 bp).
    right_flank : str
        Downstream flanking sequence (e.g., 500 bp).
    min_length : int
        Minimum microhomology length to report (default 4).
    max_length : int
        Maximum microhomology length to scan (default 25).

    Returns
    -------
    dict
        {"found": bool, "sequence": str, "length": int, "position": str}
        position describes where the microhomology was found relative to
        the insertion breakpoint.
    """
    result = {"found": False, "sequence": "", "length": 0, "position": ""}

    # Scan from longest to shortest — return the first (longest) match
    for n in range(min(max_length, len(left_flank), len(right_flank)), min_length - 1, -1):
        left_end = left_flank[-n:]
        right_start = right_flank[:n]
        if left_end == right_start:
            result = {
                "found": True,
                "sequence": left_end,
                "length": n,
                "position": f"last {n} bp of left flank = first {n} bp of right flank",
            }
            return result

    return result


# ─── 4. Direct repeat analysis ──────────────────────────────────────────────
#
# WHAT ARE DIRECT REPEATS IN THE CONTEXT OF NAHR?
# ------------------------------------------------
# Non-Allelic Homologous Recombination (NAHR) requires two copies of a similar
# sequence (direct repeats) that flank the insertion site.  During meiosis or
# DNA repair, these repeats can misalign and undergo illegitimate recombination,
# resulting in the insertion (or deletion) of the intervening DNA.
#
# We search for exact repeats of >= min_repeat_len bp within the combined
# flanking sequence.  A repeat that has one copy in the left flank and one
# copy in the right flank is especially interesting, because it suggests the
# insertion occurred between the two copies.
#
# This is a computationally simplified search: we use a sliding-window
# comparison rather than a full Smith–Waterman alignment.  This means we
# only find EXACT repeats, missing repeats with mismatches.  For a
# comprehensive analysis, tools like RepeatMasker or BLAST of the flanking
# region against itself would be more appropriate.
#
# Example:
#   Left flank:   ...ATCGATCGATCG[insertion]
#   Right flank:               [insertion]ATCGATCGATCG...
#                                          ^^^^^^^^^^^^
#                                    12 bp direct repeat

def analyse_direct_repeats(
    full_flank: str,
    left_flank_len: int,
    min_repeat_len: int = 10,
) -> dict:
    """
    Search for exact direct repeats within the flanking sequence.

    Uses a sliding window to find identical subsequences of at least
    min_repeat_len bp.  Prioritises repeats that span the insertion
    breakpoint (one copy in left flank, one in right flank).

    Parameters
    ----------
    full_flank : str
        Concatenated left_flank + right_flank sequence.
    left_flank_len : int
        Length of the left flank (to determine which side a repeat falls on).
    min_repeat_len : int
        Minimum length of repeat to report (default 10).

    Returns
    -------
    dict
        {"found": bool, "sequence": str, "positions": [int, int],
         "repeat_length": int, "spans_breakpoint": bool}
    """
    result = {
        "found": False,
        "sequence": "",
        "positions": [],
        "repeat_length": 0,
        "spans_breakpoint": False,
    }

    seq_len = len(full_flank)
    if seq_len < min_repeat_len * 2:
        return result

    best_repeat = ""
    best_positions = []
    best_spans = False

    # To keep computation tractable, limit the window sizes we scan.
    # We scan from the largest feasible window down to min_repeat_len,
    # and stop at the first hit that spans the breakpoint, or keep the
    # longest overall hit.
    max_window = min(100, seq_len // 2)

    for window_size in range(max_window, min_repeat_len - 1, -1):
        found_spanning = False
        for i in range(seq_len - window_size * 2 + 1):
            kmer = full_flank[i : i + window_size]

            # Search for the same kmer downstream (non-overlapping)
            search_start = i + window_size
            j = full_flank.find(kmer, search_start)

            if j != -1:
                # Determine if the repeat spans the breakpoint
                # (one copy in left flank, one in right flank)
                copy1_in_left = i < left_flank_len
                copy2_in_right = j >= left_flank_len
                spans = copy1_in_left and copy2_in_right

                if spans:
                    # This is the most biologically meaningful configuration
                    result = {
                        "found": True,
                        "sequence": kmer,
                        "positions": [i, j],
                        "repeat_length": window_size,
                        "spans_breakpoint": True,
                    }
                    return result  # Best possible: spanning + longest

                if window_size > len(best_repeat):
                    best_repeat = kmer
                    best_positions = [i, j]
                    best_spans = spans

        # If we found a spanning repeat at this window size, return it
        if found_spanning:
            break

    if best_repeat:
        result = {
            "found": True,
            "sequence": best_repeat,
            "positions": best_positions,
            "repeat_length": len(best_repeat),
            "spans_breakpoint": best_spans,
        }

    return result


# ─── 5. Target-site duplication (TSD) detection ─────────────────────────────
#
# WHAT ARE TARGET-SITE DUPLICATIONS?
# -----------------------------------
# When a retrotransposon (like L1) inserts into the genome, the endonuclease
# makes staggered cuts on the two strands of the target DNA.  After the new
# DNA is inserted and the gaps are filled in, the short sequence between the
# two cut sites is DUPLICATED — appearing once before and once after the
# insertion.  This duplicated sequence is called a "target-site duplication"
# or TSD.
#
# For L1-mediated insertions, TSDs are typically 7–20 bp long.  Their presence
# is a hallmark of retrotransposon-mediated insertion and strong evidence for
# the TPRT mechanism.
#
# To detect TSDs, we look for identical sequences of 7–20 bp that appear at
# the END of the left flank AND at the START of the right flank.  This is
# subtly different from microhomology: microhomology means the same sequence
# is shared at the junction, while a TSD means the sequence is DUPLICATED
# on both sides.
#
# Example:
#   Before insertion:  5'---GCATCGATCGATCG---3'
#                               ^^^^^^^^^^^
#                               future TSD (11 bp)
#
#   After L1-mediated insertion:
#     5'---GCATC[GATCGATCGATCG]---mtDNA_insert---[GATCGATCGATCG]---3'
#                ^^^^^^^^^^^^^                     ^^^^^^^^^^^^^
#                TSD (copy 1)                      TSD (copy 2)
#
# In our flanking-sequence framework:
#   Left flank ends with:     ...GCATCGATCGATCGATCG
#   Right flank starts with:  GATCGATCGATCGATTT...
#   The last 13 bp of left == first 13 bp of right?  No — that would be
#   microhomology, not TSD.
#
# IMPORTANT DISTINCTION from microhomology:
# TSDs are ADJACENT to the insertion on BOTH sides.  In our data, because
# the insertion is removed (we only have flanking genomic DNA), a TSD would
# appear as the last K bp of the left flank being identical to the first K bp
# of the right flank.  This is formally the same test as microhomology but
# with a DIFFERENT SIZE RANGE (7–20 bp for TSDs vs 2–6 bp for microhomology).
#
# In practice, we detect TSDs by looking for matches in the 7–20 bp range.

def detect_tsd(
    left_flank: str, right_flank: str, min_tsd: int = 7, max_tsd: int = 20
) -> dict:
    """
    Detect target-site duplications (TSDs) at the insertion breakpoint.

    TSDs are identical sequences of 7–20 bp found at both ends of the
    insertion site.  They are the hallmark of retrotransposon-mediated
    insertion (TPRT).

    Parameters
    ----------
    left_flank : str
        Upstream flanking sequence.
    right_flank : str
        Downstream flanking sequence.
    min_tsd : int
        Minimum TSD length (default 7 bp).
    max_tsd : int
        Maximum TSD length (default 20 bp).

    Returns
    -------
    dict
        {"found": bool, "sequence": str, "length": int}
    """
    result = {"found": False, "sequence": "", "length": 0}

    # Scan from longest to shortest — return the first (longest) match
    for n in range(min(max_tsd, len(left_flank), len(right_flank)), min_tsd - 1, -1):
        left_end = left_flank[-n:]
        right_start = right_flank[:n]
        if left_end == right_start:
            result = {
                "found": True,
                "sequence": left_end,
                "length": n,
            }
            return result

    return result


# ─── 6. L1 endonuclease site detection ──────────────────────────────────────
#
# The L1 ORF2p endonuclease cleaves the bottom strand of genomic DNA at a
# T-rich consensus sequence.  On the sense (top) strand, this appears as:
#
#   Canonical:    5'---TTTT/A---3'   (nick between T and A)
#   Degenerate:   5'---[TC]TTT[TA]---3'
#
# We search the last 20 bp of the left flank (i.e., immediately upstream of
# the insertion) for this motif.  If found, it is consistent with L1-mediated
# insertion.
#
# NOTE: T-rich sequences are common in the genome, especially in AT-rich
# regions and in Alu elements.  The presence of this motif alone is not
# sufficient evidence for TPRT — it must be combined with TSD detection
# and/or antisense orientation for a confident call.

def detect_l1_site(left_flank: str, window: int = 20) -> dict:
    """
    Check for L1 endonuclease cleavage consensus near the insertion site.

    Searches the last `window` bp of the left flank for the degenerate
    L1 EN consensus: [TC]TTT[TA].

    Parameters
    ----------
    left_flank : str
        Upstream flanking sequence.
    window : int
        Number of bp upstream of the insertion to search (default 20).

    Returns
    -------
    dict
        {"found": bool, "motif": str, "position": int}
        position is relative to the end of the left flank (negative = upstream).
    """
    result = {"found": False, "motif": "", "position": 0}

    # Extract the region immediately upstream of the insertion
    search_region = left_flank[-window:] if len(left_flank) >= window else left_flank

    match = L1_EN_CONSENSUS_REGEX.search(search_region)
    if match:
        # Position relative to end of left flank
        pos_from_end = -(len(search_region) - match.start())
        result = {
            "found": True,
            "motif": match.group(),
            "position": pos_from_end,
        }

    return result


# ─── 7. Insertion orientation analysis ───────────────────────────────────────
#
# When a mitochondrial DNA fragment is inserted into the nuclear genome, it
# can be in either orientation:
#
#   SENSE:     The inserted sequence reads the same as the mtDNA reference
#              (rCRS, NC_012920.1).  5'→3' of the insert matches 5'→3' of
#              the mitochondrial heavy strand.
#
#   ANTISENSE: The inserted sequence is the reverse complement of the mtDNA
#              reference.  This is expected for TPRT-mediated insertions
#              because reverse transcription proceeds 3'→5' on the RNA
#              template, producing a cDNA that is antisense to the original.
#
# To determine orientation, we compare the extracted_sequence (from the ClinVar
# record) to the mtDNA BLAST alignment.  If the mito_blast results include
# strand information (sstart > send indicates minus strand), we use that.
# Otherwise, we compare identity to forward vs reverse complement.

def determine_orientation(doc: dict) -> str:
    """
    Determine the orientation of the NUMT insertion relative to the mtDNA.

    Uses the mitochondrial BLAST results from Step 4 if available.  When
    sstart > send in the BLAST alignment, the match is to the minus strand
    (antisense).

    Parameters
    ----------
    doc : dict
        MongoDB document for the variant.  Expected fields:
        - mito_blast.best_hit.sstart, .send (from Step 4)
        - extracted_sequence (from Step 2)

    Returns
    -------
    str
        "SENSE", "ANTISENSE", or "UNKNOWN"
    """
    # Try to use BLAST strand information from Step 4
    mito_blast = doc.get("mito_blast", {})
    best_hit = mito_blast.get("best_hit", {})

    sstart = best_hit.get("sstart")
    send = best_hit.get("send")

    if sstart is not None and send is not None:
        try:
            sstart = int(sstart)
            send = int(send)
            # In BLAST output, sstart > send means the query matched the
            # minus strand of the subject (mtDNA).  This means the insertion
            # is in the ANTISENSE orientation relative to the rCRS.
            if sstart <= send:
                return "SENSE"
            else:
                return "ANTISENSE"
        except (ValueError, TypeError):
            pass

    # Fallback: check if the document has explicit strand annotation
    strand = best_hit.get("strand", "")
    if strand:
        strand_lower = str(strand).lower()
        if "minus" in strand_lower or "antisense" in strand_lower:
            return "ANTISENSE"
        elif "plus" in strand_lower or "sense" in strand_lower:
            return "SENSE"

    return "UNKNOWN"


# ─── 8. Predict mechanism and generate summary ──────────────────────────────

def predict_mechanism(
    microhomology: dict,
    direct_repeats: dict,
    tsd: dict,
    l1_site: dict,
    orientation: str,
    at_combined: float = 0.0,
    l1_within_5kb: bool = False,
    curvature: float = 0.0,
) -> tuple[str, str]:
    """
    Predict the most likely insertion mechanism based on sequence evidence.

    Decision logic (in order of priority):

    1. TPRT (Target-Primed Reverse Transcription via L1):
       - TSD + L1 site + ANTISENSE orientation → HIGH
       - TSD + L1 site (any orientation) → MEDIUM
       - TSD + L1_within_5kb (UCSC, no local EN site) → MEDIUM
       - TSD alone → LOW

    2. NAHR (Non-Allelic Homologous Recombination):
       - Spanning direct repeats ≥200 bp (approaching MEPS) → MEDIUM
       - Shorter repeats: fall through to junction-based classification

    3. MMEJ/MMBIR (Microhomology-Mediated End Joining / Break-Induced Replication):
       - Microhomology ≥8 bp → MEDIUM
       - Microhomology 4–7 bp + high curvature (>1.5 °/bp) → MEDIUM
       - Microhomology ≥4 bp → LOW

    4. NHEJ (Non-Homologous End Joining):
       - Microhomology 0–3 bp (including blunt-end) → LOW
       - This is the default when no TSD, no long repeats, and junction
         microhomology ≤ 3 bp (Wei et al. 2022)

    Parameters
    ----------
    microhomology : dict
        Result from analyse_microhomology().
    direct_repeats : dict
        Result from analyse_direct_repeats().
    tsd : dict
        Result from detect_tsd().
    l1_site : dict
        Result from detect_l1_site() (regex on 20 bp upstream).
    orientation : str
        "SENSE", "ANTISENSE", or "UNKNOWN".
    at_combined : float
        AT% of combined flanking sequences (from analyse_at_content()).
    l1_within_5kb : bool
        True if a LINE/L1 element was found within 5 kb via UCSC (from
        get_ucsc_repeats()).
    curvature : float
        Mean dinucleotide roll angle in °/bp (from compute_curvature()).

    Returns
    -------
    tuple of (str, str)
        (predicted_mechanism, confidence)
        mechanism is one of: "TPRT", "NAHR", "NHEJ", "MMEJ/MMBIR", "UNKNOWN"
        confidence is one of: "HIGH", "MEDIUM", "LOW"
    """
    # ── TPRT: requires TSD ──────────────────────────────────────────────────
    if tsd["found"] and l1_site["found"]:
        if orientation == "ANTISENSE":
            return "TPRT", "HIGH"      # All three canonical hallmarks
        else:
            return "TPRT", "MEDIUM"    # TSD + L1 site, orientation non-canonical

    if tsd["found"] and not l1_site["found"]:
        if l1_within_5kb:
            return "TPRT", "MEDIUM"    # TSD + L1 element nearby (UCSC)
        return "TPRT", "LOW"           # TSD alone

    # L1 site alone is not called TPRT (motif too common in AT-rich regions)

    # ── NAHR: long direct repeats spanning the breakpoint ───────────────────
    # True NAHR requires substantial homology approaching the Minimal Efficient
    # Processing Segment (MEPS ≥ 200–337 bp, Chen & Bhargava 2005). Short
    # direct repeats (< 200 bp) are insufficient for homologous recombination
    # and should not override junction-based evidence (NHEJ/MMEJ). L1 proximity
    # is irrelevant to NAHR — it is contextual for DSB repair mechanisms.
    if direct_repeats["found"] and direct_repeats.get("spans_breakpoint", False):
        if direct_repeats["repeat_length"] >= 200:
            return "NAHR", "MEDIUM"    # Approaching MEPS — genuine NAHR candidate
        # Short repeats (< 200 bp): noted in output, fall through to
        # junction-based classification (NHEJ/MMEJ) below.

    # ── MMEJ/MMBIR: microhomology ≥4 bp ─────────────────────────────────────
    # MMEJ and MMBIR share the same junction microhomology signature; their
    # distinction requires patient-derived long-read sequencing (Hastings 2009).
    if microhomology["found"] and microhomology["length"] >= 4:
        if microhomology["length"] >= 8:
            return "MMEJ/MMBIR", "MEDIUM"    # Substantial microhomology
        if curvature > 1.5:
            return "MMEJ/MMBIR", "MEDIUM"    # Short microhomology + high curvature
        return "MMEJ/MMBIR", "LOW"           # Short microhomology only

    # ── NHEJ: microhomology 0–3 bp ──────────────────────────────────────────
    # Junction microhomology ≤ 3 bp is consistent with classical NHEJ
    # (Wei et al. 2022). Blunt-end ligation (0 bp) is also NHEJ. Short direct
    # repeats and L1 proximity are noted as genomic context but do not change
    # the mechanism classification.
    return "NHEJ", "LOW"


def generate_interpretation(
    allele_id: int,
    gene: str,
    mechanism: str,
    confidence: str,
    microhomology: dict,
    tsd: dict,
    l1_site: dict,
    orientation: str,
    direct_repeats: dict,
    at_content: dict | None = None,
    repeat_context: dict | None = None,
) -> str:
    """
    Generate a deterministic scientific interpretation of the mechanism.

    This produces 2–3 sentences describing the predicted mechanism in
    scientific language suitable for a research report.  No AI generation —
    these are templated strings filled with analysis results.

    Parameters
    ----------
    allele_id : int
        ClinVar AlleleID.
    gene : str
        Gene symbol.
    mechanism : str
        Predicted mechanism (TPRT, NAHR, NHEJ, MMEJ, UNKNOWN).
    confidence : str
        Confidence level (HIGH, MEDIUM, LOW).
    microhomology, tsd, l1_site, direct_repeats : dict
        Analysis results.
    orientation : str
        Insertion orientation.

    Returns
    -------
    str
        2–3 sentence scientific interpretation.
    """
    intro = (
        f"AlleleID {allele_id} ({gene}): "
        f"The predicted insertion mechanism is {mechanism} "
        f"(confidence: {confidence})."
    )

    if mechanism == "TPRT":
        details = (
            f" A target-site duplication of {tsd['length']} bp "
            f"(\"{tsd['sequence']}\") was detected flanking the insertion, "
            f"consistent with retrotransposon-mediated integration."
        )
        if l1_site["found"]:
            details += (
                f" An L1 endonuclease consensus motif (\"{l1_site['motif']}\") "
                f"was identified {abs(l1_site['position'])} bp upstream of "
                f"the insertion site."
            )
        if orientation == "ANTISENSE":
            details += (
                " The insertion is in the antisense orientation relative to "
                "the mitochondrial reference, as expected for L1-mediated "
                "TPRT events."
            )

    elif mechanism == "NAHR":
        details = (
            f" Direct repeats of {direct_repeats['repeat_length']} bp "
            f"(\"{direct_repeats['sequence'][:20]}{'...' if direct_repeats['repeat_length'] > 20 else ''}\") "
            f"were identified flanking the insertion site, suggesting "
            f"non-allelic homologous recombination as the integration mechanism."
        )
        if direct_repeats.get("spans_breakpoint"):
            details += (
                " The repeat copies span the insertion breakpoint, with one "
                "copy upstream and one downstream."
            )

    elif mechanism == "MMEJ/MMBIR":
        details = (
            f" Microhomology of {microhomology['length']} bp "
            f"(\"{microhomology['sequence']}\") was detected at the "
            f"insertion junction, consistent with microhomology-mediated "
            f"end joining (MMEJ) or break-induced replication (MMBIR). "
            f"These mechanisms are computationally indistinguishable; "
            f"distinction requires patient-derived long-read sequencing."
        )
        details += (
            " No target-site duplications or flanking direct repeats were "
            "detected, arguing against TPRT or NAHR mechanisms."
        )

    elif mechanism == "NHEJ":
        if microhomology["found"]:
            details = (
                f" Only {microhomology['length']} bp of microhomology was "
                f"detected at the junction, consistent with classical "
                f"non-homologous end joining (NHEJ)."
            )
        else:
            details = (
                " No microhomology, target-site duplications, or direct "
                "repeats were detected at the junction. The blunt-end "
                "signature is consistent with classical NHEJ repair."
            )

    else:  # UNKNOWN
        details = (
            " No definitive sequence signatures (microhomology, TSDs, "
            "direct repeats, or L1 endonuclease sites) were detected. "
            "The insertion mechanism could not be determined computationally."
        )

    # Genomic context sentence (AT%, repeat element)
    context = ""
    if at_content:
        at_val = at_content.get("at_combined", 0.0)
        if at_val > 60.0:
            context += (
                f" The insertion flanks are AT-rich ({at_val:.1f}% AT), "
                "consistent with the AT-biased genomic context of known NUMT "
                "hotspots (Tsuji et al. 2012)."
            )
    if repeat_context and repeat_context.get("nearest_repeat_class", "N/A") != "N/A":
        rep_class  = repeat_context["nearest_repeat_class"]
        rep_family = repeat_context.get("nearest_repeat_family", "")
        rep_dist   = repeat_context.get("nearest_repeat_dist_bp", -1)
        if rep_dist >= 0:
            context += (
                f" Nearest genomic repeat element: {rep_class}"
                f"/{rep_family} at {rep_dist:,} bp from the insertion site."
            )
        if repeat_context.get("l1_within_5kb"):
            context += " A LINE-1 (L1) element is present within 5 kb."

    caveat = (
        " Note: this is a computational prediction based on reference "
        "genome flanking sequences; experimental validation (e.g., PCR "
        "across the junction, long-read sequencing) is required to "
        "confirm the mechanism."
    )

    return intro + details + context + caveat


# ─── 8b. AT content and DNA curvature ───────────────────────────────────────

def analyse_at_content(left_flank: str, right_flank: str) -> dict:
    """
    Compute AT% for the left flank, right flank, and their combination.

    AT-rich flanking sequences are a known feature of NUMT insertion sites
    (Tsuji et al. 2012: "NUMT insertion sites occur immediately adjacent to
    A+T oligomers").  High AT% supports a NHEJ/MMEJ insertion context.

    Parameters
    ----------
    left_flank : str
        Upstream flanking sequence (e.g. 500 bp).
    right_flank : str
        Downstream flanking sequence.

    Returns
    -------
    dict
        {"at_left": float, "at_right": float, "at_combined": float}
        Values are percentages rounded to 1 decimal place.
    """
    def _at_pct(seq: str) -> float:
        if not seq:
            return 0.0
        return round(100.0 * sum(1 for b in seq.upper() if b in "AT") / len(seq), 1)

    return {
        "at_left":     _at_pct(left_flank),
        "at_right":    _at_pct(right_flank),
        "at_combined": _at_pct(left_flank + right_flank),
    }


def compute_curvature(sequence: str) -> float:
    """
    Compute mean local DNA curvature using the Bolshoy dinucleotide roll model.

    For each consecutive dinucleotide step, look up its roll angle (degrees)
    from BOLSHOY_ROLL.  The mean roll angle over the full sequence is returned
    as a scalar proxy for intrinsic DNA curvature.

    Threshold for "high curvature": > 1.5 °/bp (genome-wide mean is ~1.3 °/bp).

    Reference: Bolshoy et al. (1993) PNAS 90:2312–2316.

    Parameters
    ----------
    sequence : str
        DNA sequence (typically left_flank + right_flank, ~1000 bp).

    Returns
    -------
    float
        Mean roll angle in degrees/bp.  0.0 if sequence is too short.
    """
    seq = sequence.upper()
    if len(seq) < 2:
        return 0.0
    rolls = [BOLSHOY_ROLL.get(seq[i:i + 2], 1.5) for i in range(len(seq) - 1)]
    return round(sum(rolls) / len(rolls), 3)


# ─── 8c. UCSC REST API queries ───────────────────────────────────────────────

@lru_cache(maxsize=256)
def _ucsc_get(url: str, timeout: int = 15) -> dict:
    """
    Cached GET request to the UCSC REST API.

    Results are cached in memory so that repeated calls with the same URL
    (e.g. the same genomic window for two candidates with nearby coordinates)
    do not trigger duplicate HTTP requests.

    Parameters
    ----------
    url : str
        Full UCSC REST API URL.
    timeout : int
        Request timeout in seconds.

    Returns
    -------
    dict
        Parsed JSON response, or {} on any failure (network error, HTTP error,
        JSON parse error).  Failures are silent to keep the pipeline running.
    """
    if not REQUESTS_AVAILABLE:
        return {}
    try:
        resp = requests.get(url, timeout=timeout)
        resp.raise_for_status()
        return resp.json()
    except Exception:
        return {}


def get_ucsc_repeats(chrom: str, pos: int, window: int = 5000) -> dict:
    """
    Query the UCSC RepeatMasker track (rmsk) for repeat elements near pos.

    Searches a ±window bp region around the insertion site.  Returns the
    nearest repeat element and a flag for any L1 (LINE/L1 family) within
    the window.  The 5 kb window mirrors the search radius used by
    Tsuji et al. (2012) to associate NUMTs with retrotransposon-rich loci.

    Parameters
    ----------
    chrom : str
        Chromosome name (ClinVar format, e.g. "16" or "chr16").
    pos : int
        1-based ClinVar Start coordinate.
    window : int
        Search radius in bp (default 5000).

    Returns
    -------
    dict
        {
          "nearest_repeat_class":  str,   # LINE, SINE, LTR, DNA, Simple_repeat …
          "nearest_repeat_family": str,   # L1, Alu, MIR …
          "nearest_repeat_name":   str,   # e.g. L1HS
          "nearest_repeat_dist_bp": int,  # distance from pos to element midpoint
          "l1_within_5kb":         bool,  # True if any L1 element in the window
        }
    """
    chrom_ucsc = "chr" + str(chrom).replace("chr", "")
    ucsc_start = max(0, pos - window - 1)   # 0-based, inclusive
    ucsc_end   = pos + window                # 0-based, exclusive

    url = (f"{UCSC_API}?genome=hg38;track=rmsk;"
           f"chrom={chrom_ucsc};start={ucsc_start};end={ucsc_end}")
    data = _ucsc_get(url)
    elements = data.get("rmsk", [])

    default = {
        "nearest_repeat_class":  "N/A",
        "nearest_repeat_family": "N/A",
        "nearest_repeat_name":   "N/A",
        "nearest_repeat_dist_bp": -1,
        "l1_within_5kb":         False,
    }

    if not elements:
        return default

    # Find nearest element (by midpoint distance to pos)
    pos0 = pos - 1  # convert to 0-based for distance calc
    best = None
    best_dist = float("inf")
    for el in elements:
        el_start = int(el.get("genoStart", el.get("chromStart", 0)))
        el_end   = int(el.get("genoEnd",   el.get("chromEnd",   0)))
        mid  = (el_start + el_end) // 2
        dist = abs(mid - pos0)
        if dist < best_dist:
            best_dist = dist
            best = el

    if best is None:
        return default

    rep_class  = best.get("repClass",  best.get("class",  "N/A"))
    rep_family = best.get("repFamily", best.get("family", "N/A"))
    rep_name   = best.get("repName",   best.get("name",   "N/A"))

    # Flag L1 elements: repClass == "LINE" AND repFamily contains "L1"
    l1_within = any(
        str(el.get("repClass",  el.get("class",  ""))).upper() == "LINE"
        and "L1" in str(el.get("repFamily", el.get("family", "")))
        for el in elements
    )

    return {
        "nearest_repeat_class":   rep_class,
        "nearest_repeat_family":  rep_family,
        "nearest_repeat_name":    rep_name,
        "nearest_repeat_dist_bp": int(best_dist),
        "l1_within_5kb":          l1_within,
    }


def get_ucsc_recomb_rate(chrom: str, pos: int) -> dict:
    """
    Query UCSC recombination rate track for a genomic position.

    Uses the UCSC 'recombAvg' track (deCODE average recombination rate, hg38).
    Recombination rate >10 cM/Mb is used as a proxy for a meiotic hotspot;
    true PRDM9 binding site ChIP-seq data is not available as a UCSC track.

    Parameters
    ----------
    chrom : str
        Chromosome (e.g. "16").
    pos : int
        1-based ClinVar Start coordinate.

    Returns
    -------
    dict
        {
          "recombination_rate_cM_Mb": float,
          "in_recombination_hotspot": bool,   # rate > 10 cM/Mb
        }
    """
    chrom_ucsc = "chr" + str(chrom).replace("chr", "")
    # Query a ±50 kb window — recombAvg bigWig segments can be sparse in coding regions
    ucsc_start = max(0, pos - 50_000)
    ucsc_end   = pos + 50_000

    url = (f"{UCSC_API}?genome=hg38;track=recombAvg;"
           f"chrom={chrom_ucsc};start={ucsc_start};end={ucsc_end}")
    data = _ucsc_get(url)
    records = data.get("recombAvg", [])

    default = {
        "recombination_rate_cM_Mb": 0.0,
        "in_recombination_hotspot": False,
    }

    if not records:
        return default

    # Try several possible field names for the rate value
    rate: float | None = None
    for rec in records:
        for field in ("value", "decodeAvg", "deCODEAvg", "rate", "score",
                      "decodeAvgRate", "deCODE_avg"):
            if field in rec:
                try:
                    rate = float(rec[field])
                    break
                except (ValueError, TypeError):
                    continue
        if rate is not None:
            break

    if rate is None:
        return default

    return {
        "recombination_rate_cM_Mb": round(rate, 3),
        "in_recombination_hotspot": rate > 10.0,
    }


# ─── 9. Analyse a single candidate ──────────────────────────────────────────

def analyse_candidate(
    doc: dict,
    flank_size: int = 500,
    min_microhomology: int = 4,
    email: str = "numt_pipeline@example.com",
) -> dict | None:
    """
    Run the full mechanism analysis for a single NUMT candidate.

    Steps:
      1. Extract coordinates from the MongoDB document
      2. Fetch flanking sequences from NCBI (or cache)
      3. Run microhomology analysis
      4. Run direct repeat analysis
      5. Detect target-site duplications
      6. Check for L1 endonuclease site
      7. Determine insertion orientation
      8. Predict mechanism and confidence
      9. Generate scientific interpretation

    Parameters
    ----------
    doc : dict
        MongoDB document for the variant.
    flank_size : int
        Number of bp to fetch on each side of the insertion.
    min_microhomology : int
        Minimum microhomology length to report.
    email : str
        Email for NCBI Entrez.

    Returns
    -------
    dict or None
        The mechanism_analysis dict, or None if analysis could not be performed.
    """
    allele_id = doc.get("AlleleID", "UNKNOWN")
    gene = doc.get("GeneSymbol", doc.get("gene", "UNKNOWN"))
    chrom = doc.get("Chromosome", doc.get("chromosome"))
    start = doc.get("Start", doc.get("start"))
    stop = doc.get("Stop", doc.get("stop"))

    log(f"\n{'=' * 70}")
    log(f"Analysing AlleleID {allele_id} ({gene}) — chr{chrom}:{start}-{stop}")
    log(f"{'=' * 70}")

    # ── Validate coordinates ─────────────────────────────────────────────
    if not chrom or not start or not stop:
        log(f"  WARNING: AlleleID {allele_id}: missing coordinates (Chromosome={chrom}, "
            f"Start={start}, Stop={stop}). Skipping.")
        return None

    try:
        start = int(start)
        stop = int(stop)
    except (ValueError, TypeError):
        log(f"  WARNING: AlleleID {allele_id}: non-numeric coordinates. Skipping.")
        return None

    # ── Fetch flanking sequences ─────────────────────────────────────────
    flanks = fetch_flanking_sequence(chrom, start, stop, allele_id, flank_size, email)
    if flanks is None:
        log(f"  WARNING: AlleleID {allele_id}: could not fetch flanking sequences. Skipping.")
        return None

    left_flank = flanks["left_flank"]
    right_flank = flanks["right_flank"]
    full_flank = flanks["full_flank"]

    # ── NEW: AT content ──────────────────────────────────────────────────
    at_content = analyse_at_content(left_flank, right_flank)
    log(f"  AT content: left={at_content['at_left']}%  "
        f"right={at_content['at_right']}%  "
        f"combined={at_content['at_combined']}%")

    # ── NEW: DNA curvature (Bolshoy model) ────────────────────────────────
    curvature = compute_curvature(full_flank)
    log(f"  DNA curvature (Bolshoy): {curvature:.3f} °/bp"
        f"{'  [HIGH]' if curvature > 1.5 else ''}")

    # ── NEW: UCSC RepeatMasker (±5 kb) ───────────────────────────────────
    log(f"  Querying UCSC RepeatMasker (rmsk) within 5 kb ...")
    repeat_context = get_ucsc_repeats(chrom, start)
    log(f"  >> Nearest repeat: {repeat_context['nearest_repeat_class']}/"
        f"{repeat_context['nearest_repeat_family']} "
        f"({repeat_context['nearest_repeat_name']}) "
        f"at {repeat_context['nearest_repeat_dist_bp']} bp"
        f"  L1_within_5kb={repeat_context['l1_within_5kb']}")

    # ── NEW: UCSC recombination rate ──────────────────────────────────────
    log(f"  Querying UCSC recombination rate ...")
    recomb_data = get_ucsc_recomb_rate(chrom, start)
    log(f"  >> Recombination rate: {recomb_data['recombination_rate_cM_Mb']} cM/Mb"
        f"{'  [HOTSPOT]' if recomb_data['in_recombination_hotspot'] else ''}")

    # ── 3. Microhomology analysis ────────────────────────────────────────
    log(f"  Checking for microhomology at the insertion junction...")
    microhomology = analyse_microhomology(left_flank, right_flank, min_length=min_microhomology)
    if microhomology["found"]:
        log(f"  >> Microhomology detected: {microhomology['length']} bp "
            f"(\"{microhomology['sequence']}\")")
    else:
        log(f"  >> No microhomology >= {min_microhomology} bp detected.")

    # ── 4. Direct repeat analysis ────────────────────────────────────────
    log(f"  Searching for direct repeats in flanking region...")
    direct_repeats = analyse_direct_repeats(full_flank, len(left_flank))
    if direct_repeats["found"]:
        log(f"  >> Direct repeat detected: {direct_repeats['repeat_length']} bp "
            f"at positions {direct_repeats['positions']}"
            f"{' (SPANS BREAKPOINT)' if direct_repeats['spans_breakpoint'] else ''}")
    else:
        log(f"  >> No direct repeats >= 10 bp detected.")

    # ── 5. TSD detection ─────────────────────────────────────────────────
    log(f"  Checking for target-site duplications (TSDs)...")
    tsd = detect_tsd(left_flank, right_flank)
    if tsd["found"]:
        log(f"  >> TSD detected: {tsd['length']} bp (\"{tsd['sequence']}\")")
    else:
        log(f"  >> No TSDs (7–20 bp) detected.")

    # ── 6. L1 endonuclease site ──────────────────────────────────────────
    log(f"  Checking for L1 endonuclease consensus upstream of insertion...")
    l1_site = detect_l1_site(left_flank)
    if l1_site["found"]:
        log(f"  >> L1 EN site detected: \"{l1_site['motif']}\" "
            f"at position {l1_site['position']} relative to insertion")
    else:
        log(f"  >> No L1 EN consensus detected in last 20 bp of left flank.")

    # ── 7. Insertion orientation ─────────────────────────────────────────
    log(f"  Determining insertion orientation...")
    orientation = determine_orientation(doc)
    log(f"  >> Orientation: {orientation}")

    # ── 8. Predict mechanism ─────────────────────────────────────────────
    mechanism, confidence = predict_mechanism(
        microhomology, direct_repeats, tsd, l1_site, orientation,
        at_combined=at_content["at_combined"],
        l1_within_5kb=repeat_context["l1_within_5kb"],
        curvature=curvature,
    )
    log(f"  >> PREDICTED MECHANISM: {mechanism} (confidence: {confidence})")

    # ── 9. Generate interpretation ───────────────────────────────────────
    interpretation = generate_interpretation(
        allele_id, gene, mechanism, confidence,
        microhomology, tsd, l1_site, orientation, direct_repeats,
        at_content=at_content,
        repeat_context=repeat_context,
    )
    log(f"  >> Interpretation: {interpretation}")

    # ── Build summary dict ───────────────────────────────────────────────
    mechanism_analysis = {
        "allele_id": allele_id,
        "gene": gene,
        "microhomology": {
            "found": microhomology["found"],
            "sequence": microhomology["sequence"],
            "length": microhomology["length"],
        },
        "direct_repeats": {
            "found": direct_repeats["found"],
            "sequence": direct_repeats["sequence"],
            "positions": direct_repeats["positions"],
        },
        "tsd": {
            "found": tsd["found"],
            "sequence": tsd["sequence"],
        },
        "l1_site": {
            "found": l1_site["found"],
            "motif": l1_site.get("motif", ""),
        },
        "insertion_orientation": orientation,
        "predicted_mechanism": mechanism,
        "mechanism_confidence": confidence,
        "interpretation": interpretation,
        # ── New fields ──────────────────────────────────────────────────
        "at_content": at_content,
        "dna_curvature_mean": curvature,
        "repeat_context": repeat_context,
        "recombination": recomb_data,
        "analysis_timestamp": datetime.now(timezone.utc).isoformat(),
    }

    return mechanism_analysis


# ─── 10. MongoDB update ─────────────────────────────────────────────────────

def update_mongodb(collection, results: list[dict]) -> int:
    """
    Write mechanism_analysis results back to MongoDB.

    Uses bulk $set operations to update each variant document with its
    mechanism analysis.

    Parameters
    ----------
    collection : pymongo.collection.Collection
        The clinvar_variants collection.
    results : list of dict
        List of mechanism_analysis dicts (each must have 'allele_id').

    Returns
    -------
    int
        Number of documents successfully updated.
    """
    if not results:
        log("No results to update in MongoDB.")
        return 0

    operations = []
    for result in results:
        allele_id = result["allele_id"]
        operations.append(
            UpdateOne(
                {"AlleleID": allele_id},
                {"$set": {"mechanism_analysis": result}},
            )
        )

    bulk_result = collection.bulk_write(operations)
    updated = bulk_result.modified_count
    log(f"Updated {updated} document(s) in MongoDB with mechanism analysis.")
    return updated


# ─── 11. Output: summary table, JSON report, TSV ────────────────────────────

def print_summary_table(results: list[dict]) -> None:
    """
    Print a formatted summary table of mechanism analysis results.

    Parameters
    ----------
    results : list of dict
        List of mechanism_analysis dicts.
    """
    if not results:
        log("No results to display.")
        return

    # Header
    header = (
        f"{'AlleleID':>12}  {'Gene':<12}  {'Mechanism':<10}  {'Conf':<8}  "
        f"{'Microhom':<10}  {'TSD':<10}  {'L1_site':<8}  {'Orient':<10}"
    )
    separator = "-" * len(header)

    print(f"\n{separator}")
    print("MECHANISM ANALYSIS SUMMARY")
    print(separator)
    print(header)
    print(separator)

    for r in results:
        mh_str = f"{r['microhomology']['length']}bp" if r["microhomology"]["found"] else "none"
        tsd_str = f"{r['tsd'].get('length', len(r['tsd']['sequence']))}bp" if r["tsd"]["found"] else "none"
        l1_str = r["l1_site"]["motif"] if r["l1_site"]["found"] else "none"

        print(
            f"{r['allele_id']:>12}  {r['gene']:<12}  {r['predicted_mechanism']:<10}  "
            f"{r['mechanism_confidence']:<8}  {mh_str:<10}  {tsd_str:<10}  "
            f"{l1_str:<8}  {r['insertion_orientation']:<10}"
        )

    print(separator)
    print(f"Total candidates analysed: {len(results)}")

    # Mechanism breakdown
    mechanism_counts = {}
    for r in results:
        m = r["predicted_mechanism"]
        mechanism_counts[m] = mechanism_counts.get(m, 0) + 1
    print("Mechanism breakdown:")
    for mech, count in sorted(mechanism_counts.items()):
        print(f"  {mech}: {count}")
    print(separator + "\n")


def save_json_report(results: list[dict], output_path: Path) -> None:
    """
    Save the full mechanism analysis results as a JSON file.

    Parameters
    ----------
    results : list of dict
        List of mechanism_analysis dicts.
    output_path : Path
        Output file path.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    report = {
        "pipeline_step": "06_mechanism",
        "description": "NUMT insertion mechanism analysis",
        "analysis_date": datetime.now(timezone.utc).isoformat(),
        "total_candidates": len(results),
        "results": results,
    }

    with open(output_path, "w") as fh:
        json.dump(report, fh, indent=2, default=str)

    log(f"JSON report saved to {output_path}")


def save_tsv_summary(results: list[dict], output_path: Path) -> None:
    """
    Save a tab-separated summary of mechanism analysis results.

    Parameters
    ----------
    results : list of dict
        List of mechanism_analysis dicts.
    output_path : Path
        Output file path.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    columns = [
        "allele_id", "gene", "predicted_mechanism", "mechanism_confidence",
        "microhomology_found", "microhomology_length", "microhomology_sequence",
        "tsd_found", "tsd_sequence",
        "l1_site_found", "l1_site_motif",
        "direct_repeats_found", "direct_repeats_length",
        "insertion_orientation",
        "at_content_combined", "dna_curvature_mean",
        "nearest_repeat_class", "nearest_repeat_family", "nearest_repeat_dist_bp",
        "l1_within_5kb",
        "recombination_rate_cM_Mb", "in_recombination_hotspot",
        "interpretation",
    ]

    with open(output_path, "w") as fh:
        fh.write("\t".join(columns) + "\n")

        for r in results:
            at = r.get("at_content", {})
            rep = r.get("repeat_context", {})
            rec = r.get("recombination", {})
            row = [
                str(r["allele_id"]),
                r["gene"],
                r["predicted_mechanism"],
                r["mechanism_confidence"],
                str(r["microhomology"]["found"]),
                str(r["microhomology"]["length"]),
                r["microhomology"]["sequence"],
                str(r["tsd"]["found"]),
                r["tsd"]["sequence"],
                str(r["l1_site"]["found"]),
                r["l1_site"].get("motif", ""),
                str(r["direct_repeats"]["found"]),
                str(r["direct_repeats"].get("repeat_length",
                    len(r["direct_repeats"]["sequence"]))),
                r["insertion_orientation"],
                str(at.get("at_combined", "")),
                str(r.get("dna_curvature_mean", "")),
                rep.get("nearest_repeat_class", "N/A"),
                rep.get("nearest_repeat_family", "N/A"),
                str(rep.get("nearest_repeat_dist_bp", -1)),
                str(rep.get("l1_within_5kb", False)),
                str(rec.get("rate_cM_Mb", rec.get("recombination_rate_cM_Mb", ""))),
                str(rec.get("in_hotspot", rec.get("in_recombination_hotspot", False))),
                r.get("interpretation", "").replace("\t", " ").replace("\n", " "),
            ]
            fh.write("\t".join(row) + "\n")

    log(f"TSV summary saved to {output_path}")


# ─── Main pipeline ──────────────────────────────────────────────────────────

def main(args: argparse.Namespace) -> None:
    """
    Orchestrate the mechanism analysis pipeline.

    1. Connect to MongoDB and fetch NUMT candidates
    2. For each candidate, run the full mechanism analysis
    3. Update MongoDB with results
    4. Save outputs (JSON, TSV, summary table)
    """
    log("=" * 70)
    log("Step 6: NUMT insertion mechanism analysis")
    log("=" * 70)

    # ── Dependency checks ────────────────────────────────────────────────
    if not BIOPYTHON_AVAILABLE:
        log("FATAL: Biopython is required.  Install with: pip install biopython")
        sys.exit(1)

    if not PYMONGO_AVAILABLE:
        log("FATAL: pymongo is required.  Install with: pip install pymongo")
        sys.exit(1)

    # ── Connect to MongoDB ───────────────────────────────────────────────
    log(f"Connecting to MongoDB at {MONGO_HOST}:{MONGO_PORT}...")
    try:
        client = pymongo.MongoClient(
            MONGO_HOST, MONGO_PORT, serverSelectionTimeoutMS=5000
        )
        client.admin.command("ping")
        db = client[MONGO_DB]
        collection = db[MONGO_COLLECTION]
        log(f"Connected to {MONGO_DB}.{MONGO_COLLECTION} "
            f"({collection.estimated_document_count()} documents).")
    except (ConnectionFailure, ServerSelectionTimeoutError) as exc:
        log(f"FATAL: Could not connect to MongoDB: {exc}")
        sys.exit(1)

    # ── Fetch candidates ─────────────────────────────────────────────────
    candidates = fetch_candidates(collection)
    if not candidates:
        log("No NUMT candidates found. Nothing to analyse.")
        log("  Hint: run Steps 4–5 first, or ensure AlleleID 1234567 exists for testing.")
        return

    # ── Analyse each candidate ───────────────────────────────────────────
    results = []
    for i, doc in enumerate(candidates, 1):
        log(f"\n--- Candidate {i}/{len(candidates)} ---")
        analysis = analyse_candidate(
            doc,
            flank_size=args.flank_size,
            min_microhomology=args.min_microhomology,
            email=args.email,
        )
        if analysis is not None:
            results.append(analysis)

        # Respect NCBI rate limits between candidates
        if i < len(candidates):
            time.sleep(0.5)

    # ── Summary ──────────────────────────────────────────────────────────
    log(f"\nAnalysis complete. {len(results)}/{len(candidates)} candidates analysed successfully.")

    if not results:
        log("No results to report (all candidates failed analysis).")
        return

    # ── Print summary table ──────────────────────────────────────────────
    print_summary_table(results)

    # ── Save outputs ─────────────────────────────────────────────────────
    save_json_report(results, MECHANISM_JSON)
    save_tsv_summary(results, MECHANISM_TSV)

    # ── Update MongoDB ───────────────────────────────────────────────────
    if not args.dry_run:
        update_mongodb(collection, results)
    else:
        log("Dry-run mode: MongoDB NOT updated.")

    # ── Print interpretations ────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("MECHANISM INTERPRETATIONS")
    print("=" * 70)
    for r in results:
        print(f"\n{r.get('interpretation', 'No interpretation available.')}")
    print("\n" + "=" * 70)

    log("Step 6 complete.")


# ─── Argument parsing ───────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments with the following fields:
        - min_microhomology (int): minimum microhomology length to report
        - flank_size (int): bp to fetch on each side of the insertion
        - email (str): email for NCBI Entrez
        - dry_run (bool): if True, do not update MongoDB
    """
    parser = argparse.ArgumentParser(
        description=(
            "Step 6: Investigate the molecular mechanism of NUMT insertion. "
            "Analyses flanking sequences for microhomology, direct repeats, "
            "target-site duplications, and L1 endonuclease sites to predict "
            "whether the NUMT was inserted via TPRT, NHEJ, MMEJ/MMBIR, or NAHR."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python 06_mechanism.py\n"
            "  python 06_mechanism.py --flank-size 1000 --min-microhomology 3\n"
            "  python 06_mechanism.py --dry-run\n"
            "\n"
            "This script requires:\n"
            "  - MongoDB running locally with numt_db.clinvar_variants populated\n"
            "  - Biopython (pip install biopython)\n"
            "  - pymongo (pip install pymongo)\n"
            "  - Internet access for NCBI Entrez efetch\n"
        ),
    )

    parser.add_argument(
        "--min-microhomology",
        type=int,
        default=4,
        metavar="BP",
        help=(
            "Minimum microhomology length (bp) to report as evidence for "
            "NHEJ/MMEJ.  Shorter matches (2-3 bp) often occur by chance.  "
            "Default: 4."
        ),
    )

    parser.add_argument(
        "--flank-size",
        type=int,
        default=500,
        metavar="BP",
        help=(
            "Number of base pairs to fetch on each side (upstream/downstream) "
            "of the insertion site.  Larger flanks allow detection of more "
            "distant repeats but increase NCBI query time.  Default: 500."
        ),
    )

    parser.add_argument(
        "--email",
        type=str,
        default="numt_pipeline@example.com",
        metavar="EMAIL",
        help=(
            "Email address for NCBI Entrez (required by NCBI usage policy).  "
            "Default: numt_pipeline@example.com."
        ),
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help=(
            "Run the analysis but do NOT update MongoDB.  Useful for testing "
            "or when you only want the output files (JSON, TSV)."
        ),
    )

    return parser.parse_args()


# ─── Entry point ─────────────────────────────────────────────────────────────

if __name__ == "__main__":
    args = parse_args()
    main(args)
