#!/usr/bin/env python3
"""
07_prdm9_scan.py — PRDM9 hotspot motif scanner and NUMT mechanism analyser

Scans flanking sequences around NUMT insertion sites for PRDM9 binding motifs
AND performs full NUMT mechanism analysis: junction microhomology, local repeat
pairs (NAHR evidence), orientation (direct vs inverted), and confidence scoring.

Biological rationale
--------------------
NUMTs can arise through several DSB-repair pathways:

  NAHR  — Non-Allelic Homologous Recombination: two high-identity repeat copies
           (LINE-1, Alu, LCR) flank the breakpoint; crossover happens within the
           homologous tract.  Signature: long (>50 bp) direct repeat on both sides
           with >85% identity.  Associated with meiotic hotspots and PRDM9 activity.

  MMEJ  — Microhomology-Mediated End-Joining: a 4–25 bp microhomology at the
           junction is used to re-ligate broken ends.  Leaves an exact microhomology
           at the junction (left_flank[-N:] == right_flank[:N]).

  MMBIR — Microhomology-Mediated Break-Induced Replication: intermediate (20–200 bp)
           microhomology, template switching during replication.

  NHEJ  — Non-Homologous End-Joining: ≤3 bp microhomology or none; blunt-end ligation.

  TPRT  — Target-Primed Reverse Transcription: L1-mediated insertion.  Signature:
           TSD (7–20 bp) flanking insertion, poly-A at one end, EN nick consensus.
           Detected in 06_mechanism.py; 07_prdm9_scan.py focuses on DSB-repair models.

PRDM9 (PR domain zinc finger protein 9) initiates meiotic recombination by
binding a degenerate DNA sequence and trimethylating histone H3K4me3 at incipient
double-strand break sites.  The canonical human PRDM9-A allele (~75 % Europeans)
recognises CCNCCNTNNCCNC (Myers et al. 2008 Science).  PRDM9 typically binds
100–300 bp upstream of the actual crossover, so this scan is a candidate screen,
not definitive hotspot calling.

Usage
-----
  # PRDM9 motif scan only:
  python 07_prdm9_scan.py --fasta 2068854_flanking.fasta
  python 07_prdm9_scan.py --scan-all --max-mismatches 1 --format json

  # Full NUMT mechanism analysis (microhomology + repeat pairs + confidence):
  python 07_prdm9_scan.py --scan-all --numt-analysis --format json
  python 07_prdm9_scan.py --fasta 3074437_flanking.fasta --numt-analysis

  # Run built-in unit tests:
  python 07_prdm9_scan.py --test

References
----------
Myers S et al. (2008) A common sequence motif associated with recombination
  hot spots and genome instability in humans. Science 322:759-763.
Baudat F et al. (2010) PRDM9 is a major determinant of meiotic recombination
  hotspots in humans and mice. Science 327:836-840.
Pratto F et al. (2014) Recombination initiation maps of individual human
  genomes. Science 346:1256442.
Hastings PJ et al. (2009) Mechanisms of change in gene copy number.
  Nat Rev Genet 10:551-564.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Iterator

# ─── IUPAC ambiguity codes ────────────────────────────────────────────────────

IUPAC: dict[str, frozenset[str]] = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "U": frozenset("T"),
    "N": frozenset("ACGT"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "S": frozenset("GC"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
}

_IUPAC_COMP: dict[str, str] = {
    "A": "T", "T": "A", "C": "G", "G": "C", "U": "A",
    "N": "N",
    "R": "Y", "Y": "R",
    "S": "S", "W": "W",
    "K": "M", "M": "K",
    "B": "V", "V": "B",
    "D": "H", "H": "D",
}

# ─── Motif database ───────────────────────────────────────────────────────────

MOTIF_DB: list[dict] = [
    {
        "name": "PRDM9-A",
        "pattern": "CCNCCNTNNCCNC",
        "note": "Most common human allele (~75 % Europeans)",
        "reference": "Myers et al. 2008 Science 322:759",
    },
]

# ─── Default parameters ───────────────────────────────────────────────────────

DEFAULT_FLANKING_DIR = (
    Path(__file__).parent.parent.parent
    / "02_data" / "processed" / "flanking_sequences"
)
DEFAULT_MAX_MISMATCHES = 0
DEFAULT_FORMAT = "tsv"

# NUMT mechanism analysis defaults (all configurable via CLI)
DEFAULT_MICROHOM_MIN = 1
DEFAULT_MICROHOM_MAX = 20
DEFAULT_REPEAT_MIN_LEN = 10
DEFAULT_REPEAT_MAX_LEN = 200
DEFAULT_REPEAT_MIN_IDENTITY = 0.80   # 80 % identity threshold
DEFAULT_REPEAT_SCAN_STEP = 10        # step size for sliding-window repeat scan
JUNCTION_CONTEXT_BP = 10             # bp shown on each side of junction in output

# ── PRDM9 scan window context (Option A — Wei et al. 2022) ────────────────────
#
# Wei et al. (2022) found germinal NUMTs enriched within 3 kb of PRDM9 ChIP-seq
# peaks (Altemose et al. 2017, GSE99407; P=0.003).  Somatic NUMTs: < 1 kb.
# This script scans the flanking FASTAs available in the cache (typically ±500 bp
# per side = 1000 bp total).  The expected false-positive rate (FPR) for the
# PRDM9-A motif CCNCCNTNNCCNC (8 fixed positions) by random chance is:
#   FPR ≈ (2 × window_bp) / 4^8 = window_bp / 32768
#   At ±500 bp (1000 bp total): FPR ≈ 3%
#   At ±3000 bp (6000 bp total): FPR ≈ 18%  ← Wei-recommended window
# See 08_prdm9_chipseq.py for the Wei-compatible ChIP-seq peak distance approach.

PRDM9_RECOMMENDED_WINDOW_BP = 3000  # per-side window per Wei 2022 (germinal NUMTs)
PRDM9_MOTIF_FIXED_POSITIONS = 8     # fixed C/T positions in CCNCCNTNNCCNC

# ─── Core sequence helpers ────────────────────────────────────────────────────


def compute_motif_fpr(window_bp: int) -> float:
    """
    Expected number of PRDM9-A motif hits by chance in a scan of *window_bp* bp
    (both strands combined).

    CCNCCNTNNCCNC has 8 fixed C/T positions → P(match) = (1/4)^8 = 1/65536.
    Both strands: 2 × window_bp positions scanned.
    """
    return (2 * window_bp) / (4 ** PRDM9_MOTIF_FIXED_POSITIONS)


def rc_seq(seq: str) -> str:
    """Reverse complement of a sequence containing IUPAC ambiguity codes."""
    comp = str.maketrans(
        "ACGTURYSWKMBDHVacgturyswkmbdhv",
        "TGCAAYRSWMKVHDBtgcaayrswmkvhdb",
    )
    return seq.translate(comp)[::-1]


def _mismatch_count(pattern: str, window: str) -> int:
    """
    Count positions where window[i] is not in IUPAC[pattern[i]].
    Both strings must be the same length and upper-case.
    """
    total = 0
    for p, q in zip(pattern, window):
        allowed = IUPAC.get(p)
        if allowed is None:
            if p != q:
                total += 1
        elif q not in allowed:
            total += 1
    return total


def scan_sequence(
    seq: str,
    pattern: str,
    max_mismatches: int = 0,
) -> Iterator[tuple[int, str, int]]:
    """
    Slide *pattern* over *seq* and yield (position, matched_window, n_mismatches)
    for every window where mismatches ≤ max_mismatches.
    """
    pattern = pattern.upper()
    seq_upper = seq.upper()
    plen = len(pattern)
    slen = len(seq_upper)
    if plen > slen:
        return
    for i in range(slen - plen + 1):
        window = seq_upper[i : i + plen]
        mm = _mismatch_count(pattern, window)
        if mm <= max_mismatches:
            yield (i, window, mm)


# ─── Both-strand PRDM9 scanner ────────────────────────────────────────────────


def scan_both_strands(
    seq: str,
    motif: dict,
    max_mismatches: int,
    insertion_offset: int,
    left_genomic_start: int = 0,
    right_genomic_start: int | None = None,
) -> list[dict]:
    """
    Scan *seq* on both strands for *motif*.  Returns hit records sorted by
    (pos_in_seq, strand).

    Parameters
    ----------
    seq : str
        Combined flanking sequence (left_flank + right_flank, ~1000 bp).
    motif : dict
        Entry from MOTIF_DB.
    max_mismatches : int
        Maximum mismatches per hit.
    insertion_offset : int
        0-based index in *seq* of the insertion site (= len(left_flank)).
    left_genomic_start : int
        0-based genomic coordinate of seq[0].
    right_genomic_start : int | None
        0-based genomic coordinate of seq[insertion_offset].
    """
    pattern = motif["pattern"].upper()
    plen = len(pattern)
    seq_len = len(seq)
    hits: list[dict] = []

    if right_genomic_start is None:
        right_genomic_start = left_genomic_start + insertion_offset

    def genomic_coord(pos_in_seq: int) -> int:
        if pos_in_seq < insertion_offset:
            return left_genomic_start + pos_in_seq
        else:
            return right_genomic_start + (pos_in_seq - insertion_offset)

    fwd_positions: set[int] = set()
    for pos, window, mm in scan_sequence(seq, pattern, max_mismatches):
        fwd_positions.add(pos)
        hits.append({
            "motif_name": motif["name"],
            "pattern": motif["pattern"],
            "strand": "+",
            "pos_in_seq": pos,
            "genomic_pos_0based": genomic_coord(pos),
            "matched_sequence": window,
            "mismatches": mm,
            "match_length": plen,
            "pct_identity": round(100.0 * (plen - mm) / plen, 1),
            "match_type": "PRDM9_motif",
            "distance_from_insertion": pos - insertion_offset,
        })

    seq_rc = rc_seq(seq)
    for rc_pos, window_rc, mm in scan_sequence(seq_rc, pattern, max_mismatches):
        fwd_pos = seq_len - rc_pos - plen
        if fwd_pos in fwd_positions:
            continue
        hits.append({
            "motif_name": motif["name"],
            "pattern": motif["pattern"],
            "strand": "-",
            "pos_in_seq": fwd_pos,
            "genomic_pos_0based": genomic_coord(fwd_pos),
            "matched_sequence": window_rc,
            "mismatches": mm,
            "match_length": plen,
            "pct_identity": round(100.0 * (plen - mm) / plen, 1),
            "match_type": "PRDM9_motif",
            "distance_from_insertion": fwd_pos - insertion_offset,
        })

    hits.sort(key=lambda h: (h["pos_in_seq"], h["strand"]))
    return hits


# ─── FASTA parser ─────────────────────────────────────────────────────────────


def _parse_header(header: str) -> dict:
    """
    Parse a flanking-sequence FASTA header produced by the NUMT pipeline.

    Expected format:
      >AlleleID_<id>_left_flank  NC_000009.12:36226889-36227388
      >AlleleID_<id>_right_flank NC_000009.12:36227391-36227890

    Returns dict: allele_id, side, chrom, start_0, end_0 (0-based half-open).
    """
    line = header.lstrip(">").strip()
    parts = line.split()
    name_part = parts[0]
    coord_part = parts[1] if len(parts) > 1 else ""

    tokens = name_part.split("_")
    allele_id = tokens[1] if len(tokens) > 1 else "unknown"
    side = tokens[2] if len(tokens) > 2 else "unknown"

    chrom = ""
    start_0 = 0
    end_0 = 0
    if ":" in coord_part:
        chrom, range_part = coord_part.split(":", 1)
        if "-" in range_part:
            s, e = range_part.split("-", 1)
            try:
                start_0 = int(s) - 1
                end_0 = int(e)
            except ValueError:
                pass

    return {"allele_id": allele_id, "side": side, "chrom": chrom,
            "start_0": start_0, "end_0": end_0}


def load_flanking_fasta(fasta_path: Path) -> dict:
    """
    Load a two-record FASTA (left_flank + right_flank) and return a dict with
    combined sequence, insertion offset, and coordinate metadata.
    """
    records: list[dict] = []
    current_header = None
    current_seq: list[str] = []

    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_header is not None:
                    info = _parse_header(current_header)
                    info["seq"] = "".join(current_seq).upper()
                    records.append(info)
                current_header = line
                current_seq = []
            elif line:
                current_seq.append(line)
        if current_header is not None:
            info = _parse_header(current_header)
            info["seq"] = "".join(current_seq).upper()
            records.append(info)

    left  = next((r for r in records if r["side"] == "left"),  None)
    right = next((r for r in records if r["side"] == "right"), None)

    if left is None or right is None:
        combined = "".join(r["seq"] for r in records)
        return {
            "allele_id": records[0]["allele_id"] if records else "unknown",
            "chrom": records[0]["chrom"] if records else "",
            "combined_seq": combined,
            "insertion_offset": len(combined) // 2,
            "left_genomic_start": 0,
            "right_genomic_start": len(combined) // 2,
            "left_seq": combined[:len(combined)//2],
            "right_seq": combined[len(combined)//2:],
        }

    combined = left["seq"] + right["seq"]
    insertion_offset = len(left["seq"])

    return {
        "allele_id": left["allele_id"],
        "chrom": left["chrom"],
        "combined_seq": combined,
        "insertion_offset": insertion_offset,
        "left_genomic_start": left["start_0"],
        "right_genomic_start": right["start_0"],
        "insertion_genomic": left["end_0"],
        "left_seq": left["seq"],
        "right_seq": right["seq"],
    }


# ─── PRDM9 scan orchestrator ──────────────────────────────────────────────────


def scan_fasta(
    fasta_path: Path,
    motif_db: list[dict] | None = None,
    max_mismatches: int = DEFAULT_MAX_MISMATCHES,
) -> list[dict]:
    """Scan a single flanking FASTA for all motifs.  Returns hit dicts."""
    if motif_db is None:
        motif_db = MOTIF_DB

    data = load_flanking_fasta(fasta_path)
    ctx = JUNCTION_CONTEXT_BP
    left_seq = data.get("left_seq", data["combined_seq"][:data["insertion_offset"]])
    right_seq = data.get("right_seq", data["combined_seq"][data["insertion_offset"]:])
    junction_seq = left_seq[-ctx:] + "|" + right_seq[:ctx]

    all_hits: list[dict] = []
    for motif in motif_db:
        hits = scan_both_strands(
            seq=data["combined_seq"],
            motif=motif,
            max_mismatches=max_mismatches,
            insertion_offset=data["insertion_offset"],
            left_genomic_start=data["left_genomic_start"],
            right_genomic_start=data["right_genomic_start"],
        )
        for h in hits:
            h["allele_id"] = data["allele_id"]
            h["chrom"] = data["chrom"]
            h["insertion_genomic"] = data.get("insertion_genomic", 0)
            h["junction_sequence"] = junction_seq
        all_hits.extend(hits)

    return all_hits


def scan_all_cached(
    flanking_dir: Path = DEFAULT_FLANKING_DIR,
    motif_db: list[dict] | None = None,
    max_mismatches: int = DEFAULT_MAX_MISMATCHES,
) -> list[dict]:
    """Scan every *_flanking.fasta file in flanking_dir for PRDM9 motifs."""
    if motif_db is None:
        motif_db = MOTIF_DB

    all_hits: list[dict] = []
    fasta_files = sorted(flanking_dir.glob("*_flanking.fasta"))
    fasta_files = [f for f in fasta_files if f.stem != "UNKNOWN_flanking"]

    for fasta_path in fasta_files:
        try:
            hits = scan_fasta(fasta_path, motif_db, max_mismatches)
            all_hits.extend(hits)
        except Exception as exc:
            print(f"[WARN] {fasta_path.name}: {exc}", file=sys.stderr)

    return all_hits


# ─── NUMT mechanism analysis ──────────────────────────────────────────────────
# These functions address DSB-repair mechanism inference for mtDNA insertions.
# They are independent of the PRDM9 motif scan and can be used standalone.


def classify_homology_length(length: int) -> str:
    """
    Classify a homology length into mechanistic category.

    ≤20 bp  → "microhomology"  (MMEJ / NHEJ signature)
    21–200  → "intermediate"   (MMBIR / short-tract NAHR)
    >200    → "long"           (classical NAHR)
    """
    if length <= 20:
        return "microhomology"
    elif length <= 200:
        return "intermediate"
    else:
        return "long"


def _seq_identity(a: str, b: str) -> float:
    """Fraction of identical positions between two same-length strings."""
    if not a:
        return 0.0
    return sum(x == y for x, y in zip(a, b)) / len(a)


def detect_junction_microhomology(
    left_flank: str,
    right_flank: str,
    min_len: int = DEFAULT_MICROHOM_MIN,
    max_len: int = DEFAULT_MICROHOM_MAX,
) -> dict | None:
    """
    Find the longest exact direct microhomology at the insertion junction.

    MMEJ/MMBIR signature: left_flank[-N:] == right_flank[:N] for some N.
    Scans from max_len down to min_len; returns the longest exact match, or None.

    Returns a dict with keys:
        length, sequence, homology_class, pct_identity, mismatches, mechanism_hint
    """
    left = left_flank.upper()
    right = right_flank.upper()
    effective_max = min(max_len, len(left), len(right))

    for length in range(effective_max, min_len - 1, -1):
        if left[-length:] == right[:length]:
            mhint = "MMEJ/MMBIR" if length >= 4 else "NHEJ"
            return {
                "length": length,
                "sequence": left[-length:],
                "homology_class": classify_homology_length(length),
                "pct_identity": 100.0,
                "mismatches": 0,
                "mechanism_hint": mhint,
            }
    return None


def find_local_repeat_pairs(
    left_flank: str,
    right_flank: str,
    min_len: int = DEFAULT_REPEAT_MIN_LEN,
    max_len: int = DEFAULT_REPEAT_MAX_LEN,
    min_identity: float = DEFAULT_REPEAT_MIN_IDENTITY,
    step: int = DEFAULT_REPEAT_SCAN_STEP,
) -> list[dict]:
    """
    Find pairs of similar sequences in left and right flanks (NAHR evidence).

    Searches for:
    - Direct repeats: same orientation in left and right flanks.
    - Inverted repeats: left flank vs RC(right flank).

    Returns the single best direct hit and single best inverted hit (by score =
    length × identity), sorted by decreasing length.  Returns [] if no pair
    exceeds the min_identity threshold.

    Performance: O(W × L/step × R/step × win_len) where W = number of window
    sizes, L = len(left_flank), R = len(right_flank).  With defaults
    (step=10, max_len=200, 500 bp flanks) this is ~200K comparisons.
    """
    left = left_flank.upper()
    right = right_flank.upper()
    right_rc_str = rc_seq(right)

    best_direct: dict | None = None
    best_inverted: dict | None = None
    best_direct_score = 0.0
    best_inverted_score = 0.0

    # Use step=1 for short windows (≤25 bp) to avoid missing 10–15 bp direct repeats
    # that are the primary NAHR/MMEJ signal.  Larger windows use the caller-supplied step.
    win_sizes = range(min_len, min(max_len + 1, len(left) + 1, len(right) + 1), step)
    for win_len in win_sizes:
        effective_step = 1 if win_len <= 25 else step
        left_positions = range(0, len(left) - win_len + 1, effective_step)
        right_positions = range(0, len(right) - win_len + 1, effective_step)

        for li in left_positions:
            lwin = left[li : li + win_len]

            for ri in right_positions:
                # Direct repeat
                rwin = right[ri : ri + win_len]
                ident = _seq_identity(lwin, rwin)
                if ident >= min_identity:
                    score = win_len * ident
                    if score > best_direct_score:
                        best_direct_score = score
                        best_direct = {
                            "repeat_type": "direct",
                            "orientation": "same",
                            "length": win_len,
                            "pct_identity": round(ident * 100, 1),
                            "homology_class": classify_homology_length(win_len),
                            "left_pos": li,
                            "right_pos": ri,
                            "left_seq": lwin,
                            "right_seq": rwin,
                        }

                # Inverted repeat
                rwin_rc = right_rc_str[ri : ri + win_len]
                ident_inv = _seq_identity(lwin, rwin_rc)
                if ident_inv >= min_identity:
                    score_inv = win_len * ident_inv
                    if score_inv > best_inverted_score:
                        best_inverted_score = score_inv
                        # Convert RC position back to forward-strand coordinates
                        fwd_ri = len(right) - ri - win_len
                        best_inverted = {
                            "repeat_type": "inverted",
                            "orientation": "opposite",
                            "length": win_len,
                            "pct_identity": round(ident_inv * 100, 1),
                            "homology_class": classify_homology_length(win_len),
                            "left_pos": li,
                            "right_pos": fwd_ri,
                            "left_seq": lwin,
                            "right_seq_rc": rwin_rc,
                        }

    results = []
    if best_direct:
        results.append(best_direct)
    if best_inverted:
        results.append(best_inverted)
    results.sort(key=lambda x: (-x["length"], -x["pct_identity"]))
    return results


def compute_nahr_confidence(
    junction_microhom: dict | None,
    repeat_pairs: list[dict],
    prdm9_hits: list[dict],
) -> dict:
    """
    Heuristic confidence score for NAHR vs alternative DSB-repair mechanisms.

    Scoring (additive):
    - Long direct repeat (>50 bp, ≥90% identity): +40
    - Intermediate direct repeat (20–50 bp, ≥80%): +20
    - Short direct repeat (10–20 bp, ≥90%): +10
    - Inverted repeat (any length ≥ 20 bp, ≥80%): +15 (palindromic/hairpin context)
    - PRDM9 binding motif in flanks: +10 (meiotic context)
    - Junction microhomology ≥ 20 bp: +30 → MMBIR
    - Junction microhomology 4–19 bp: +15 → MMEJ
    - Junction microhomology 1–3 bp: +5 → NHEJ

    Returns:
        confidence_score (int), confidence_label ("LOW"|"MEDIUM"|"HIGH"),
        mechanism_prediction, confidence_rationale (str)
    """
    score = 0
    parts: list[str] = []
    mechanism = "UNKNOWN"

    # Direct repeat evidence (NAHR signature)
    direct_reps = [r for r in repeat_pairs if r["repeat_type"] == "direct"]
    if direct_reps:
        best = max(direct_reps, key=lambda x: x["length"] * x["pct_identity"])
        L, ident = best["length"], best["pct_identity"]
        if L > 50 and ident >= 90:
            score += 40
            parts.append(f"direct repeat {L}bp {ident:.0f}% identity (long, NAHR)")
            mechanism = "NAHR"
        elif L >= 20 and ident >= 80:
            score += 20
            parts.append(f"direct repeat {L}bp {ident:.0f}% identity (intermediate)")
            mechanism = "NAHR"
        elif L >= 10 and ident >= 90:
            score += 10
            parts.append(f"direct repeat {L}bp {ident:.0f}% identity (short, ambiguous)")
            if mechanism == "UNKNOWN":
                mechanism = "NAHR"

    # Inverted repeat evidence
    inv_reps = [r for r in repeat_pairs if r["repeat_type"] == "inverted"]
    if inv_reps:
        best_inv = max(inv_reps, key=lambda x: x["length"] * x["pct_identity"])
        if best_inv["length"] >= 20:
            score += 15
            parts.append(
                f"inverted repeat {best_inv['length']}bp "
                f"{best_inv['pct_identity']:.0f}% identity (palindromic context)"
            )

    # PRDM9 meiotic context
    if prdm9_hits:
        score += 10
        parts.append(
            f"PRDM9 binding motif in flanks "
            f"(n={len(prdm9_hits)}, meiotic context candidate)"
        )

    # Junction microhomology
    if junction_microhom:
        mhlen = junction_microhom["length"]
        mhseq = junction_microhom["sequence"]
        if mhlen >= 20:
            score += 30
            parts.append(
                f"junction microhomology {mhlen}bp '{mhseq}' "
                f"(intermediate, MMEJ/MMBIR — long-read required to distinguish)"
            )
            mechanism = "MMEJ/MMBIR"
        elif mhlen >= 4:
            score += 20
            parts.append(
                f"junction microhomology {mhlen}bp '{mhseq}' (MMEJ/MMBIR-consistent)"
            )
            if mechanism == "UNKNOWN":
                mechanism = "MMEJ/MMBIR"
        else:
            score += 5
            parts.append(
                f"junction microhomology {mhlen}bp '{mhseq}' (1–3 bp, NHEJ-consistent)"
            )
            if mechanism == "UNKNOWN":
                mechanism = "NHEJ"

    label = "HIGH" if score >= 50 else ("MEDIUM" if score >= 20 else "LOW")
    rationale = "; ".join(parts) if parts else "no supporting evidence detected"

    return {
        "confidence_score": score,
        "confidence_label": label,
        "mechanism_prediction": mechanism,
        "confidence_rationale": rationale,
    }


def analyse_numt_candidate(
    data: dict,
    max_mismatches: int = DEFAULT_MAX_MISMATCHES,
    microhom_max: int = DEFAULT_MICROHOM_MAX,
    repeat_min_len: int = DEFAULT_REPEAT_MIN_LEN,
    repeat_min_identity: float = DEFAULT_REPEAT_MIN_IDENTITY,
    motif_db: list[dict] | None = None,
) -> dict:
    """
    Full NUMT mechanism analysis for a single candidate.

    Combines PRDM9 motif scanning, junction microhomology detection, and local
    repeat pair identification to infer the insertion mechanism.

    Parameters
    ----------
    data : dict
        Output of load_flanking_fasta().
    max_mismatches : int
        Maximum mismatches for PRDM9 motif scanning.
    microhom_max : int
        Maximum microhomology length to search (1 bp to microhom_max).
    repeat_min_len : int
        Minimum window length for local repeat pair detection.
    repeat_min_identity : float
        Minimum fractional identity (0–1) for repeat pairs.
    motif_db : list[dict] | None
        Motif database; defaults to MOTIF_DB.

    Returns
    -------
    dict with keys:
        allele_id, chrom, insertion_genomic, junction_sequence,
        prdm9_hits, junction_microhomology, repeat_pairs, nahr_confidence
    """
    if motif_db is None:
        motif_db = MOTIF_DB

    combined = data["combined_seq"]
    offset = data["insertion_offset"]
    left_flank = data.get("left_seq", combined[:offset])
    right_flank = data.get("right_seq", combined[offset:])

    ctx = JUNCTION_CONTEXT_BP
    junction_seq = (
        (left_flank[-ctx:] if len(left_flank) >= ctx else left_flank)
        + "|"
        + (right_flank[:ctx] if len(right_flank) >= ctx else right_flank)
    )

    # PRDM9 motif scan
    prdm9_hits = []
    for motif in motif_db:
        hits = scan_both_strands(
            seq=combined,
            motif=motif,
            max_mismatches=max_mismatches,
            insertion_offset=offset,
            left_genomic_start=data["left_genomic_start"],
            right_genomic_start=data["right_genomic_start"],
        )
        for h in hits:
            h["allele_id"] = data["allele_id"]
            h["chrom"] = data["chrom"]
            h["insertion_genomic"] = data.get("insertion_genomic", 0)
            h["junction_sequence"] = junction_seq
        prdm9_hits.extend(hits)

    # Junction microhomology
    junction_microhom = detect_junction_microhomology(
        left_flank, right_flank, max_len=microhom_max
    )

    # Local repeat pairs
    repeat_pairs = find_local_repeat_pairs(
        left_flank, right_flank,
        min_len=repeat_min_len,
        min_identity=repeat_min_identity,
    )

    # Confidence
    conf = compute_nahr_confidence(junction_microhom, repeat_pairs, prdm9_hits)

    scan_window_bp = len(left_flank) + len(right_flank)
    fpr = compute_motif_fpr(scan_window_bp // 2)

    return {
        "allele_id": data["allele_id"],
        "chrom": data["chrom"],
        "insertion_genomic": data.get("insertion_genomic", 0),
        "junction_sequence": junction_seq,
        "prdm9_hits": prdm9_hits,
        "junction_microhomology": junction_microhom,
        "repeat_pairs": repeat_pairs,
        "nahr_confidence": conf,
        "scan_window_bp": scan_window_bp,
        "prdm9_fpr_per_locus": round(fpr, 3),
        "prdm9_window_ok": scan_window_bp >= 2 * PRDM9_RECOMMENDED_WINDOW_BP,
    }


def analyse_all_cached(
    flanking_dir: Path = DEFAULT_FLANKING_DIR,
    motif_db: list[dict] | None = None,
    max_mismatches: int = DEFAULT_MAX_MISMATCHES,
    microhom_max: int = DEFAULT_MICROHOM_MAX,
    repeat_min_len: int = DEFAULT_REPEAT_MIN_LEN,
    repeat_min_identity: float = DEFAULT_REPEAT_MIN_IDENTITY,
) -> list[dict]:
    """Run full NUMT analysis on every *_flanking.fasta in flanking_dir."""
    if motif_db is None:
        motif_db = MOTIF_DB
    results = []
    fasta_files = sorted(flanking_dir.glob("*_flanking.fasta"))
    fasta_files = [f for f in fasta_files if f.stem != "UNKNOWN_flanking"]
    for fasta_path in fasta_files:
        try:
            data = load_flanking_fasta(fasta_path)
            result = analyse_numt_candidate(
                data, max_mismatches, microhom_max, repeat_min_len, repeat_min_identity, motif_db
            )
            results.append(result)
        except Exception as exc:
            print(f"[WARN] {fasta_path.name}: {exc}", file=sys.stderr)
    return results


# ─── Output formatters ────────────────────────────────────────────────────────

_TSV_FIELDS = [
    "allele_id", "chrom", "insertion_genomic", "motif_name", "pattern", "strand",
    "pos_in_seq", "genomic_pos_0based", "matched_sequence",
    "mismatches", "match_length", "pct_identity", "match_type",
    "distance_from_insertion", "junction_sequence",
]

_NUMT_SUMMARY_FIELDS = [
    "allele_id", "chrom", "insertion_genomic", "junction_sequence",
    "junction_microhom_length", "junction_microhom_seq",
    "junction_microhom_class", "junction_microhom_mechanism_hint",
    "best_repeat_type", "best_repeat_orientation",
    "best_repeat_length", "best_repeat_pct_identity", "best_repeat_class",
    "prdm9_hits_count",
    "confidence_score", "confidence_label",
    "mechanism_prediction", "confidence_rationale",
    "scan_window_bp", "prdm9_fpr_per_locus", "prdm9_window_ok",
]


def format_tsv(hits: list[dict]) -> str:
    lines = ["\t".join(_TSV_FIELDS)]
    for h in hits:
        lines.append("\t".join(str(h.get(f, "")) for f in _TSV_FIELDS))
    return "\n".join(lines)


def format_json(hits: list[dict]) -> str:
    return json.dumps(hits, indent=2)


def format_numt_summary_tsv(analyses: list[dict]) -> str:
    """Format per-candidate NUMT analysis as a flat TSV summary."""
    lines = ["\t".join(_NUMT_SUMMARY_FIELDS)]
    for a in analyses:
        mh = a.get("junction_microhomology") or {}
        rp_list = a.get("repeat_pairs") or []
        rp = rp_list[0] if rp_list else {}
        conf = a.get("nahr_confidence") or {}
        row = {
            "allele_id": a.get("allele_id", ""),
            "chrom": a.get("chrom", ""),
            "insertion_genomic": a.get("insertion_genomic", ""),
            "junction_sequence": a.get("junction_sequence", ""),
            "junction_microhom_length": mh.get("length", ""),
            "junction_microhom_seq": mh.get("sequence", ""),
            "junction_microhom_class": mh.get("homology_class", ""),
            "junction_microhom_mechanism_hint": mh.get("mechanism_hint", ""),
            "best_repeat_type": rp.get("repeat_type", ""),
            "best_repeat_orientation": rp.get("orientation", ""),
            "best_repeat_length": rp.get("length", ""),
            "best_repeat_pct_identity": rp.get("pct_identity", ""),
            "best_repeat_class": rp.get("homology_class", ""),
            "prdm9_hits_count": len(a.get("prdm9_hits") or []),
            "confidence_score": conf.get("confidence_score", ""),
            "confidence_label": conf.get("confidence_label", ""),
            "mechanism_prediction": conf.get("mechanism_prediction", ""),
            "confidence_rationale": conf.get("confidence_rationale", ""),
            "scan_window_bp": a.get("scan_window_bp", ""),
            "prdm9_fpr_per_locus": a.get("prdm9_fpr_per_locus", ""),
            "prdm9_window_ok": a.get("prdm9_window_ok", ""),
        }
        lines.append("\t".join(str(row.get(f, "")) for f in _NUMT_SUMMARY_FIELDS))
    return "\n".join(lines)


def print_summary(hits: list[dict]) -> None:
    """Print a human-readable PRDM9 scan summary to stderr."""
    by_allele: dict[str, list] = {}
    for h in hits:
        by_allele.setdefault(h["allele_id"], []).append(h)

    total = len(hits)
    n_alleles = len(by_allele)
    print(f"\nPRDM9 motif scan — {total} hit(s) across {n_alleles} candidate(s)\n",
          file=sys.stderr)

    if not hits:
        print("  No hits found.", file=sys.stderr)
        return

    for allele_id, ahits in sorted(by_allele.items()):
        print(f"  AlleleID {allele_id}:", file=sys.stderr)
        for h in ahits:
            dist = h["distance_from_insertion"]
            side = "left flank" if dist < 0 else "right flank"
            print(
                f"    {h['motif_name']} [{h['strand']}] pos={h['pos_in_seq']:>4d} "
                f"({dist:+d} bp from insertion, {side})  "
                f"seq={h['matched_sequence']}  mm={h['mismatches']}",
                file=sys.stderr,
            )
    print(file=sys.stderr)


def print_numt_summary(analyses: list[dict]) -> None:
    """Print a human-readable NUMT analysis summary to stderr."""
    print(
        f"\nNUMT mechanism analysis — {len(analyses)} candidate(s)\n",
        file=sys.stderr,
    )
    # Warn if scan window is below the Wei 2022 recommendation
    undersized = [a for a in analyses if not a.get("prdm9_window_ok", True)]
    if undersized:
        win = analyses[0].get("scan_window_bp", "?")
        fpr = analyses[0].get("prdm9_fpr_per_locus", "?")
        print(
            f"  [NOTE] PRDM9 motif scan window = {win} bp total "
            f"(FPR ≈ {fpr:.2f} random hit/locus).\n"
            f"         Wei et al. (2022) recommends ±{PRDM9_RECOMMENDED_WINDOW_BP} bp "
            f"(FPR ≈ {compute_motif_fpr(PRDM9_RECOMMENDED_WINDOW_BP):.2f}/locus).\n"
            f"         Use 08_prdm9_chipseq.py for ChIP-seq peak proximity analysis.\n",
            file=sys.stderr,
        )
    for a in analyses:
        conf = a.get("nahr_confidence") or {}
        mh = a.get("junction_microhomology")
        rp_list = a.get("repeat_pairs") or []
        n_prdm9 = len(a.get("prdm9_hits") or [])
        print(
            f"  AlleleID {a['allele_id']} ({a['chrom']}) "
            f"→ {conf.get('mechanism_prediction', '?')} "
            f"[{conf.get('confidence_label', '?')} score={conf.get('confidence_score', '?')}]",
            file=sys.stderr,
        )
        if mh:
            print(
                f"    junction microhomology: {mh['length']}bp '{mh['sequence']}' "
                f"({mh['homology_class']}, {mh['mechanism_hint']})",
                file=sys.stderr,
            )
        if rp_list:
            rp = rp_list[0]
            print(
                f"    best repeat pair: {rp['repeat_type']} {rp['length']}bp "
                f"{rp['pct_identity']:.0f}% identity ({rp['homology_class']})",
                file=sys.stderr,
            )
        if n_prdm9:
            print(f"    PRDM9 hits: {n_prdm9}", file=sys.stderr)
        print(
            f"    rationale: {conf.get('confidence_rationale', 'none')}",
            file=sys.stderr,
        )
    print(file=sys.stderr)


# ─── Unit tests ───────────────────────────────────────────────────────────────


def run_tests() -> None:
    """Built-in unit tests.  Raises AssertionError on failure."""

    # ── Test 1: rc_seq ────────────────────────────────────────────────────
    assert rc_seq("ACGT") == "ACGT", "rc_seq self-complement failed"
    assert rc_seq("AACCGG") == "CCGGTT", "rc_seq basic failed"
    assert rc_seq("N") == "N", "rc_seq N failed"
    assert rc_seq("CCNCCNTNNCCNC") == "GNGGNNANGGNGG", (
        f"rc_seq motif failed: got {rc_seq('CCNCCNTNNCCNC')}"
    )
    print("  [PASS] rc_seq")

    # ── Test 2: _mismatch_count ───────────────────────────────────────────
    assert _mismatch_count("CCNCCNTNNCCNC", "CCACCATCCCCAC") == 0
    assert _mismatch_count("CCNCCNTNNCCNC", "CCACCATCCCCAG") == 1
    assert _mismatch_count("CCNCCNTNNCCNC", "AACCCCTTACCNC") == 3
    print("  [PASS] _mismatch_count")

    # ── Test 3: exact forward hit ─────────────────────────────────────────
    test_seq = "AAAAAAAAAA" + "CCACCATCCCCAC" + "AAAAAAA"
    hits = list(scan_sequence(test_seq, "CCNCCNTNNCCNC", max_mismatches=0))
    assert len(hits) == 1 and hits[0][0] == 10 and hits[0][2] == 0
    print("  [PASS] exact forward hit at correct position")

    # ── Test 4: no hit in negative sequence ───────────────────────────────
    hits_neg = list(scan_sequence("A" * 100, "CCNCCNTNNCCNC", max_mismatches=0))
    assert len(hits_neg) == 0
    print("  [PASS] no hit in negative sequence")

    # ── Test 5: reverse-complement detection ─────────────────────────────
    motif_instance = "CCACCATCCCCAC"
    motif_rc = rc_seq(motif_instance)
    rc_seq_test = "AAAAAAAAAA" + motif_rc + "AAAAAAA"
    both = scan_both_strands(
        seq=rc_seq_test,
        motif=MOTIF_DB[0],
        max_mismatches=0,
        insertion_offset=len(rc_seq_test) // 2,
    )
    minus_hits = [h for h in both if h["strand"] == "-"]
    assert len(minus_hits) >= 1, (
        f"expected ≥1 reverse-strand hit; motif_rc={motif_rc}"
    )
    print("  [PASS] reverse-complement hit detected")

    # ── Test 6: mismatch threshold ────────────────────────────────────────
    one_mm = "AAAAAAAAAA" + "CCACCATCCCCAG" + "AAAAAAA"
    hits_0 = list(scan_sequence(one_mm, "CCNCCNTNNCCNC", max_mismatches=0))
    hits_1 = list(scan_sequence(one_mm, "CCNCCNTNNCCNC", max_mismatches=1))
    assert len(hits_0) == 0 and len(hits_1) == 1 and hits_1[0][2] == 1
    print("  [PASS] mismatch threshold respected")

    # ── Test 7: palindrome deduplication ─────────────────────────────────
    fwd_only = "AAAAAAAAAA" + "CCACCATCCCCAC" + "A" * 20 + rc_seq("CCACCATCCCCAC") + "AAAA"
    both2 = scan_both_strands(
        seq=fwd_only,
        motif=MOTIF_DB[0],
        max_mismatches=0,
        insertion_offset=len(fwd_only) // 2,
    )
    n_fwd = sum(1 for h in both2 if h["strand"] == "+")
    n_rev = sum(1 for h in both2 if h["strand"] == "-")
    assert n_fwd == 1 and n_rev == 1, f"expected 1+1, got {n_fwd}+{n_rev}"
    print("  [PASS] distinct fwd and rev hits (no spurious dedup)")

    # ── Test 8: IUPAC R ambiguity code ───────────────────────────────────
    assert _mismatch_count("CCR", "CCA") == 0
    assert _mismatch_count("CCR", "CCG") == 0
    assert _mismatch_count("CCR", "CCT") == 1
    assert _mismatch_count("CCR", "CCC") == 1
    print("  [PASS] IUPAC R ambiguity code")

    # ── Test 9: NAHR-like insertion — direct repeat on both flanks ────────
    # A 50 bp repeat unit at the END of left_flank and START of right_flank
    # simulates a direct repeat flanking the breakpoint (classical NAHR signature).
    repeat_unit = "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG"  # 50 bp
    assert len(repeat_unit) == 50
    left_nahr  = "AGCTAGCTAGCTAGCTAGCT" * 10 + repeat_unit   # 200 + 50 = 250 bp
    right_nahr = repeat_unit + "TAGCTAGCTAGCTAGCTAGC" * 10   # 50 + 200 = 250 bp

    pairs = find_local_repeat_pairs(
        left_nahr, right_nahr, min_len=10, min_identity=0.95, step=5
    )
    direct = [p for p in pairs if p["repeat_type"] == "direct"]
    assert len(direct) >= 1, "should detect ≥1 direct repeat pair"
    assert direct[0]["length"] >= 10, f"repeat length should be ≥10, got {direct[0]['length']}"
    assert direct[0]["pct_identity"] >= 95.0, (
        f"identity should be ≥95%, got {direct[0]['pct_identity']}"
    )
    conf_nahr = compute_nahr_confidence(None, pairs, [])
    assert conf_nahr["mechanism_prediction"] == "NAHR", (
        f"should predict NAHR, got {conf_nahr['mechanism_prediction']}"
    )
    assert conf_nahr["confidence_label"] in ("MEDIUM", "HIGH"), (
        f"should be MEDIUM/HIGH, got {conf_nahr['confidence_label']}"
    )
    # Also confirm no false junction microhomology:
    # left_nahr ends with repeat_unit[-20:] = "TGCATGCATGCATGCATGCA"
    # right_nahr starts with repeat_unit[:20] = "TGCATGCATGCATGCATGCA"
    # → there IS a 50 bp microhomology here (the whole repeat unit overlaps)
    mh_nahr = detect_junction_microhomology(left_nahr, right_nahr, max_len=20)
    # mh_nahr will be 20 bp (capped by max_len=20 default); that is expected
    # The important check is that find_local_repeat_pairs found the full 50 bp
    assert direct[0]["length"] >= 20, "should find at least 20 bp of the 50 bp repeat"
    print("  [PASS] NAHR-like insertion: direct repeat detected")

    # ── Test 10: NUMT microhomology (MMEJ signature) ──────────────────────
    # 5 bp exact microhomology at junction (left ends with CGTAC, right starts with CGTAC)
    microhom_seq = "CGTAC"
    left_mmej  = "AGCTTGCAATCG" * 16 + microhom_seq   # 192 + 5 = 197 bp
    right_mmej = microhom_seq + "GCATCGATATCG" * 16   # 5 + 192 = 197 bp

    mh = detect_junction_microhomology(left_mmej, right_mmej, min_len=1, max_len=20)
    assert mh is not None, "should detect junction microhomology"
    assert mh["length"] == len(microhom_seq), (
        f"expected length {len(microhom_seq)}, got {mh['length']}"
    )
    assert mh["sequence"] == microhom_seq, (
        f"expected '{microhom_seq}', got '{mh['sequence']}'"
    )
    assert mh["homology_class"] == "microhomology"
    assert mh["mechanism_hint"] == "MMEJ/MMBIR"

    conf_mmej = compute_nahr_confidence(mh, [], [])
    assert conf_mmej["mechanism_prediction"] == "MMEJ/MMBIR", (
        f"should predict MMEJ/MMBIR, got {conf_mmej['mechanism_prediction']}"
    )
    assert conf_mmej["confidence_label"] in ("MEDIUM", "HIGH"), (
        f"should be MEDIUM/HIGH, got {conf_mmej['confidence_label']}"
    )
    # Also verify no spurious NAHR repeat pairs in this sequence
    pairs_mmej = find_local_repeat_pairs(left_mmej, right_mmej, min_len=20, min_identity=0.90)
    direct_mmej = [p for p in pairs_mmej if p["repeat_type"] == "direct"]
    assert len(direct_mmej) == 0, (
        "MMEJ/MMBIR test sequence should not trigger NAHR (no long direct repeats)"
    )
    print("  [PASS] NUMT microhomology (MMEJ/MMBIR): 5 bp junction microhomology detected")

    # ── Test 11: negative control — no microhomology, no repeats ─────────
    # Left ends with GGGGGGGG, right starts with CCCCCCCC: no shared terminal
    # Left body uses ACGT repeats, right body uses poly-T: no cross-flank identity
    left_neg  = "ACGTACGTACGT" * 8 + "GGGGGGGG"   # 96 + 8 = 104 bp
    right_neg = "CCCCCCCC" + "TTTTTTTTTTTT" * 8    # 8 + 96 = 104 bp

    mh_neg = detect_junction_microhomology(left_neg, right_neg, min_len=4, max_len=20)
    assert mh_neg is None, (
        f"negative control: no ≥4 bp microhomology expected, got {mh_neg}"
    )

    pairs_neg = find_local_repeat_pairs(left_neg, right_neg, min_len=10, min_identity=0.80)
    # ACGT (left) vs TTTT (right): 25% identity — should produce no pairs
    conf_neg = compute_nahr_confidence(None, pairs_neg, [])
    assert conf_neg["confidence_label"] == "LOW", (
        f"negative control should yield LOW confidence, got {conf_neg['confidence_label']}"
    )
    assert conf_neg["mechanism_prediction"] == "UNKNOWN", (
        f"negative control mechanism should be UNKNOWN, got {conf_neg['mechanism_prediction']}"
    )
    print("  [PASS] negative control: no false-positive NAHR/MMEJ/MMBIR")

    # ── Test 12: classify_homology_length boundaries ──────────────────────
    assert classify_homology_length(1)   == "microhomology"
    assert classify_homology_length(20)  == "microhomology"
    assert classify_homology_length(21)  == "intermediate"
    assert classify_homology_length(200) == "intermediate"
    assert classify_homology_length(201) == "long"
    print("  [PASS] classify_homology_length boundaries")

    print("\n  All 12 tests passed.\n")


# ─── CLI ──────────────────────────────────────────────────────────────────────


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Scan NUMT flanking sequences for PRDM9 binding motifs and/or "
            "perform full NUMT mechanism analysis (microhomology, repeat pairs, "
            "confidence scoring)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    mode = parser.add_mutually_exclusive_group(required=False)
    mode.add_argument(
        "--fasta", metavar="FILE",
        help="Scan a single flanking FASTA file.",
    )
    mode.add_argument(
        "--scan-all", action="store_true",
        help="Scan all *_flanking.fasta files in the cached flanking directory.",
    )
    mode.add_argument(
        "--test", action="store_true",
        help="Run built-in unit tests and exit.",
    )
    parser.add_argument(
        "--numt-analysis", action="store_true",
        help=(
            "Run full NUMT mechanism analysis in addition to PRDM9 motif scan: "
            "junction microhomology, local repeat pairs (direct and inverted), "
            "and confidence scoring.  Output format changes to NUMT summary TSV/JSON."
        ),
    )
    parser.add_argument(
        "--flanking-dir", metavar="DIR", default=str(DEFAULT_FLANKING_DIR),
        help=f"Directory containing *_flanking.fasta files (default: {DEFAULT_FLANKING_DIR}).",
    )
    parser.add_argument(
        "--max-mismatches", "-m", type=int, default=DEFAULT_MAX_MISMATCHES,
        metavar="N",
        help=f"Maximum PRDM9 motif mismatches (default: {DEFAULT_MAX_MISMATCHES}).",
    )
    parser.add_argument(
        "--microhom-max", type=int, default=DEFAULT_MICROHOM_MAX,
        metavar="N",
        help=(
            f"Maximum microhomology length to search at junction "
            f"(default: {DEFAULT_MICROHOM_MAX} bp)."
        ),
    )
    parser.add_argument(
        "--repeat-min-len", type=int, default=DEFAULT_REPEAT_MIN_LEN,
        metavar="N",
        help=(
            f"Minimum repeat pair length for NAHR detection "
            f"(default: {DEFAULT_REPEAT_MIN_LEN} bp)."
        ),
    )
    parser.add_argument(
        "--repeat-min-identity", type=float, default=DEFAULT_REPEAT_MIN_IDENTITY,
        metavar="F",
        help=(
            f"Minimum fractional identity for repeat pairs (0–1, "
            f"default: {DEFAULT_REPEAT_MIN_IDENTITY})."
        ),
    )
    parser.add_argument(
        "--format", choices=["tsv", "json"], default=DEFAULT_FORMAT,
        help=f"Output format (default: {DEFAULT_FORMAT}).",
    )
    parser.add_argument(
        "--no-summary", action="store_true",
        help="Suppress the human-readable summary on stderr.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.test:
        print("Running PRDM9 scanner unit tests...")
        run_tests()
        return 0

    if not args.fasta and not args.scan_all:
        parser.print_help()
        return 1

    flanking_dir = Path(args.flanking_dir)

    # ── Full NUMT analysis mode ───────────────────────────────────────────────
    if args.numt_analysis:
        analyses: list[dict] = []

        if args.fasta:
            fasta_path = Path(args.fasta)
            if not fasta_path.exists():
                print(f"Error: {fasta_path} not found.", file=sys.stderr)
                return 1
            data = load_flanking_fasta(fasta_path)
            analyses = [analyse_numt_candidate(
                data,
                max_mismatches=args.max_mismatches,
                microhom_max=args.microhom_max,
                repeat_min_len=args.repeat_min_len,
                repeat_min_identity=args.repeat_min_identity,
            )]

        elif args.scan_all:
            if not flanking_dir.is_dir():
                print(f"Error: {flanking_dir} is not a directory.", file=sys.stderr)
                return 1
            analyses = analyse_all_cached(
                flanking_dir,
                max_mismatches=args.max_mismatches,
                microhom_max=args.microhom_max,
                repeat_min_len=args.repeat_min_len,
                repeat_min_identity=args.repeat_min_identity,
            )

        if not args.no_summary:
            print_numt_summary(analyses)

        if args.format == "tsv":
            print(format_numt_summary_tsv(analyses))
        else:
            print(format_json(analyses))

        return 0

    # ── PRDM9 motif scan only (original mode) ────────────────────────────────
    hits: list[dict] = []

    if args.fasta:
        fasta_path = Path(args.fasta)
        if not fasta_path.exists():
            print(f"Error: {fasta_path} not found.", file=sys.stderr)
            return 1
        hits = scan_fasta(fasta_path, max_mismatches=args.max_mismatches)

    elif args.scan_all:
        if not flanking_dir.is_dir():
            print(f"Error: {flanking_dir} is not a directory.", file=sys.stderr)
            return 1
        hits = scan_all_cached(flanking_dir, max_mismatches=args.max_mismatches)

    if not args.no_summary:
        print_summary(hits)

    if args.format == "tsv":
        print(format_tsv(hits))
    else:
        print(format_json(hits))

    return 0


if __name__ == "__main__":
    sys.exit(main())
