#!/usr/bin/env python3
"""
08_prdm9_chipseq.py — PRDM9 ChIP-seq peak proximity analysis for NUMT candidates

Implements the methodology of Wei et al. (2022, Nature) using experimental PRDM9
binding sites from Altemose et al. (2017) ChIP-seq data (GSE99407, 170 198 peaks).

For each NUMT candidate, computes the distance to the nearest PRDM9 ChIP-seq peak
on the same chromosome and classifies:
  PRDM9_SOMATIC  : distance < 1 kb  — enriched for somatic NUMTs (P=0.003, Wei 2022)
  PRDM9_GERMINAL : distance 1–3 kb  — enriched for germinal NUMTs (P=0.003, Wei 2022)
  NO_PRDM9       : distance > 3 kb  — no significant enrichment

Validates significance via a 1000-iteration permutation test with size-matched
random genomic positions as controls (Wei 2022 methodology).

Coordinate system
-----------------
GSE99407 peaks are aligned to hg19/GRCh37 (confirmed: chr19 max pos 59,087,683 >
hg38 chr19 size 58,617,616).  NUMT flanking FASTAs use hg38/GRCh38 NC accessions
(.10/.11/.12 suffixes).  This script lifts NUMT insertion positions from hg38 → hg19
using a chain file (download with --download-chain or supply via --chain).

Compare with 07_prdm9_scan.py (motif-sequence scan, ±500 bp) to assess
concordance between sequence-based and ChIP-seq-based approaches.

Usage
-----
  # One-time: download chain file
  python 08_prdm9_chipseq.py --download-chain 02_data/raw/hg38ToHg19.over.chain.gz

  python 08_prdm9_chipseq.py \\
      --peaks 02_data/raw/GSE99407_ChIPseq_Peaks.*.txt.gz \\
      --chain 02_data/raw/hg38ToHg19.over.chain.gz

  python 08_prdm9_chipseq.py \\
      --peaks 02_data/raw/GSE99407_*.txt.gz \\
      --chain 02_data/raw/hg38ToHg19.over.chain.gz \\
      --permutations 1000 --format json

References
----------
Altemose N et al. (2017) — GEO dataset GSE99407.
  PRDM9 ChIP-seq peaks in YFP-tagged human PRDM9 B cells (hg19).
Wei W et al. (2022) Nuclear-embedded mitochondrial DNA sequences in 66,083
  human genomes. Nature 611:105-114. https://doi.org/10.1038/s41586-022-05288-7
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import random
import sys
import urllib.request
from pathlib import Path

# pandas/numpy are required; they are part of the NUMT conda environment
try:
    import numpy as np
    import pandas as pd
except ImportError:
    sys.exit(
        "Error: pandas and numpy are required.\n"
        "  micromamba activate NUMT && pip install pandas numpy"
    )

# pyliftover is required for hg38 → hg19 coordinate conversion
try:
    from pyliftover import LiftOver
    _PYLIFTOVER_AVAILABLE = True
except ImportError:
    _PYLIFTOVER_AVAILABLE = False

# ─── Genomic constants ────────────────────────────────────────────────────────

# hg19 chromosome sizes (GRCh37 — used by GSE99407 peaks)
HG19_CHR_SIZES: dict[str, int] = {
    "chr1": 249250621, "chr2": 243199373, "chr3": 198022430,
    "chr4": 191154276, "chr5": 180915260, "chr6": 171115067,
    "chr7": 159138663, "chr8": 146364022, "chr9": 141213431,
    "chr10": 135534747, "chr11": 135006516, "chr12": 133851895,
    "chr13": 115169878, "chr14": 107349540, "chr15": 102531392,
    "chr16": 90354753, "chr17": 81195210, "chr18": 78077248,
    "chr19": 59128983, "chr20": 63025520, "chr21": 48129895,
    "chr22": 51304566, "chrX": 155270560, "chrY": 59373566,
}

# hg38 chromosome sizes (GRCh38.p14 / Ensembl 108) — used by NUMT FASTAs
HG38_CHR_SIZES: dict[str, int] = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559,
    "chr4": 190214555, "chr5": 181538259, "chr6": 170805979,
    "chr7": 159345973, "chr8": 145138636, "chr9": 138394717,
    "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
    "chr16": 90338345, "chr17": 83257441, "chr18": 80373285,
    "chr19": 58617616, "chr20": 64444167, "chr21": 46709983,
    "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
}

CHAIN_URL = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"
)

# NC RefSeq accession → UCSC chromosome name (hg38)
NC_TO_CHR: dict[str, str] = {
    "NC_000001.11": "chr1",  "NC_000002.12": "chr2",  "NC_000003.12": "chr3",
    "NC_000004.12": "chr4",  "NC_000005.10": "chr5",  "NC_000006.12": "chr6",
    "NC_000007.14": "chr7",  "NC_000008.11": "chr8",  "NC_000009.12": "chr9",
    "NC_000010.11": "chr10", "NC_000011.10": "chr11", "NC_000012.12": "chr12",
    "NC_000013.11": "chr13", "NC_000014.9":  "chr14", "NC_000015.10": "chr15",
    "NC_000016.10": "chr16", "NC_000017.11": "chr17", "NC_000018.10": "chr18",
    "NC_000019.10": "chr19", "NC_000020.11": "chr20", "NC_000021.9":  "chr21",
    "NC_000022.11": "chr22", "NC_000023.11": "chrX",  "NC_000024.10": "chrY",
}

# Wei et al. (2022) distance thresholds (Fig. 5e)
PRDM9_SOMATIC_THRESHOLD_BP  = 1000   # < 1 kb → somatic NUMT enrichment
PRDM9_GERMINAL_THRESHOLD_BP = 3000   # < 3 kb → germinal NUMT enrichment

DEFAULT_FLANKING_DIR = (
    Path(__file__).parent.parent.parent
    / "02_data" / "processed" / "flanking_sequences"
)
DEFAULT_CHAIN = (
    Path(__file__).parent.parent.parent
    / "02_data" / "raw" / "hg38ToHg19.over.chain.gz"
)
DEFAULT_PERMUTATIONS = 1000
DEFAULT_FORMAT = "tsv"

# ─── Chain file download ──────────────────────────────────────────────────────


def download_chain(dest: Path) -> None:
    """Download the hg38→hg19 chain file from UCSC if not already present."""
    if dest.exists():
        print(f"Chain file already present: {dest}", file=sys.stderr)
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading chain file from UCSC → {dest} …", file=sys.stderr)
    urllib.request.urlretrieve(CHAIN_URL, dest)
    print(f"  Done. Size: {dest.stat().st_size / 1e6:.1f} MB", file=sys.stderr)


# ─── Genome build validator ───────────────────────────────────────────────────


def detect_peaks_genome(peaks_path: Path) -> str:
    """
    Heuristic: read the max chr19 position in the peaks file.
    If > hg38 chr19 size (58,617,616) → hg19.  Otherwise → unknown (assume hg38).
    """
    opener = gzip.open if str(peaks_path).endswith(".gz") else open
    max_pos = 0
    with opener(peaks_path, "rt") as fh:
        next(fh)  # skip header
        for line in fh:
            parts = line.split("\t", 3)
            if len(parts) >= 3 and parts[0] == "chr19":
                try:
                    max_pos = max(max_pos, int(parts[2]))
                except ValueError:
                    pass
    if max_pos > HG38_CHR_SIZES["chr19"]:
        return "hg19"
    return "hg38"


# ─── Liftover helper ─────────────────────────────────────────────────────────


def build_liftover(chain_path: Path) -> "LiftOver | None":
    """Load a pyliftover LiftOver object from *chain_path*.  Returns None if unavailable."""
    if not _PYLIFTOVER_AVAILABLE:
        return None
    if not chain_path.exists():
        return None
    return LiftOver(str(chain_path))


def liftover_position(lo: "LiftOver", chrom: str, pos_1based: int) -> tuple[str, int] | None:
    """
    Convert a 1-based position (hg38) to hg19 using pyliftover.
    Returns (chrom_hg19, pos_1based_hg19) or None if unmapped.

    pyliftover uses 0-based coordinates internally.
    """
    result = lo.convert_coordinate(chrom, pos_1based - 1)
    if not result:
        return None
    chrom_out, pos_0, strand, score = result[0]
    return chrom_out, pos_0 + 1   # back to 1-based


# ─── ChIP-seq peak loader ────────────────────────────────────────────────────


def load_prdm9_peaks(peaks_path: Path) -> pd.DataFrame:
    """
    Load GSE99407 ChIP-seq peak file → DataFrame(chr, center, ci_start, ci_stop).

    The file is tab-separated, first row is header.  Key columns:
      chr          — UCSC chromosome name (chr1 … chrY)
      center_start — peak centre start (0-based half-open)
      center_stop  — peak centre stop
      enrichment   — ChIP enrichment ratio (informational)
      pvalue       — -log10 p-value (informational)
    """
    opener = gzip.open if str(peaks_path).endswith(".gz") else open
    with opener(peaks_path, "rt") as fh:
        df = pd.read_csv(fh, sep="\t", low_memory=False,
                         usecols=["chr", "center_start", "center_stop",
                                  "enrichment", "pvalue"])

    df = df[df["chr"].isin(HG38_CHR_SIZES)].copy()
    df["center"] = (df["center_start"] + df["center_stop"]) // 2
    df = df.sort_values(["chr", "center"]).reset_index(drop=True)
    return df


# ─── FASTA header parser ──────────────────────────────────────────────────────


def parse_fasta_coords(fasta_path: Path) -> dict | None:
    """
    Parse the two-record flanking FASTA produced by the NUMT pipeline and return:
      allele_id, chrom (UCSC), insertion_pos (0-based, = end of left flank),
      left_len, right_len (flanking sequence lengths in bp).

    Returns None if the file cannot be parsed.
    """
    records: list[dict] = []
    current_name = None
    current_len = 0

    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current_name is not None:
                    records.append({"header": current_name, "seq_len": current_len})
                current_name = line
                current_len = 0
            else:
                current_len += len(line)
        if current_name is not None:
            records.append({"header": current_name, "seq_len": current_len})

    left = next((r for r in records if "_left_flank" in r["header"]), None)
    right = next((r for r in records if "_right_flank" in r["header"]), None)
    if left is None or right is None:
        return None

    # Parse: >AlleleID_3074437_left_flank NC_000019.10:7524665-7525164
    parts = left["header"].lstrip(">").split()
    name_tokens = parts[0].split("_")
    allele_id = name_tokens[1] if len(name_tokens) > 1 else "unknown"

    coord_part = parts[1] if len(parts) > 1 else ""
    chrom_ucsc = ""
    insertion_pos = 0
    if ":" in coord_part:
        nc_acc, range_part = coord_part.split(":", 1)
        chrom_ucsc = NC_TO_CHR.get(nc_acc, nc_acc)
        if "-" in range_part:
            _, end_str = range_part.split("-", 1)
            try:
                insertion_pos = int(end_str)   # 1-based end of left flank = insertion point
            except ValueError:
                pass

    return {
        "allele_id": allele_id,
        "chrom": chrom_ucsc,
        "insertion_pos": insertion_pos,   # 1-based coordinate
        "left_len": left["seq_len"],
        "right_len": right["seq_len"],
    }


def load_all_numt_coords(flanking_dir: Path) -> list[dict]:
    """Load coordinates for all cached flanking FASTAs (excludes UNKNOWN)."""
    coords = []
    for fasta_path in sorted(flanking_dir.glob("*_flanking.fasta")):
        if "UNKNOWN" in fasta_path.stem:
            continue
        rec = parse_fasta_coords(fasta_path)
        if rec and rec["chrom"]:
            coords.append(rec)
    return coords


# ─── Distance calculation ─────────────────────────────────────────────────────


def nearest_peak_distance(chrom: str, pos: int, peaks_df: pd.DataFrame) -> dict:
    """
    Return the distance (bp) and metadata of the nearest PRDM9 ChIP-seq peak
    on *chrom* relative to *pos* (1-based).

    Returns a dict with:
      distance_bp   — absolute distance (0 if overlapping the peak centre)
      peak_center   — genomic position of the nearest peak centre
      peak_enrichment, peak_pvalue — from the ChIP-seq data
    """
    subset = peaks_df[peaks_df["chr"] == chrom]
    if subset.empty:
        return {"distance_bp": None, "peak_center": None,
                "peak_enrichment": None, "peak_pvalue": None}

    centers = subset["center"].values
    idx = np.searchsorted(centers, pos)

    candidates = []
    if idx > 0:
        candidates.append(idx - 1)
    if idx < len(centers):
        candidates.append(idx)

    best_idx = min(candidates, key=lambda i: abs(centers[i] - pos))
    dist = abs(int(centers[best_idx]) - pos)
    row = subset.iloc[best_idx]

    return {
        "distance_bp": dist,
        "peak_center": int(centers[best_idx]),
        "peak_enrichment": float(row["enrichment"]),
        "peak_pvalue": float(row["pvalue"]),
    }


def classify_prdm9(distance_bp: int | None) -> str:
    """Classify PRDM9 association by distance (Wei et al. 2022 thresholds)."""
    if distance_bp is None:
        return "NO_PRDM9"
    if distance_bp < PRDM9_SOMATIC_THRESHOLD_BP:
        return "PRDM9_SOMATIC"
    if distance_bp < PRDM9_GERMINAL_THRESHOLD_BP:
        return "PRDM9_GERMINAL"
    return "NO_PRDM9"


# ─── Permutation test ─────────────────────────────────────────────────────────


def run_permutation_test(
    numt_records: list[dict],
    peaks_df: pd.DataFrame,
    n_iter: int = DEFAULT_PERMUTATIONS,
    seed: int = 42,
) -> dict:
    """
    Permutation test: compare the observed mean distance to PRDM9 peaks against
    n_iter sets of random control positions with matched chromosome distribution.

    Methodology follows Wei et al. (2022):
    - Control positions are drawn uniformly from the same chromosomes as the NUMTs
    - The test statistic is the mean distance to the nearest PRDM9 peak
    - P-value = fraction of permutations with mean distance ≤ observed

    Returns
    -------
    dict with: observed_mean_dist, permutation_mean, permutation_sd,
               p_value, n_iter, n_numt_with_peaks
    """
    rng = random.Random(seed)

    # Observed distances (exclude NUMTs with no peak on their chromosome)
    obs_dists = []
    for rec in numt_records:
        d = rec.get("distance_bp")
        if d is not None:
            obs_dists.append(d)

    if not obs_dists:
        return {"observed_mean_dist": None, "p_value": None,
                "n_iter": n_iter, "n_numt_with_peaks": 0}

    obs_mean = float(np.mean(obs_dists))

    # Build chromosome pool weighted by NUMT chromosome distribution
    # Use hg19 sizes since peaks are in hg19 and NUMTs are lifted to hg19
    chrom_pool = [rec["chrom_hg19"] for rec in numt_records
                  if rec.get("distance_bp") is not None and rec.get("chrom_hg19")]

    perm_means = []
    for _ in range(n_iter):
        perm_dists = []
        for chrom in chrom_pool:
            chrom_size = HG19_CHR_SIZES.get(chrom, 0)
            if chrom_size == 0:
                continue
            rand_pos = rng.randint(1, chrom_size)
            peak_info = nearest_peak_distance(chrom, rand_pos, peaks_df)
            d = peak_info.get("distance_bp")
            if d is not None:
                perm_dists.append(d)
        if perm_dists:
            perm_means.append(float(np.mean(perm_dists)))

    if not perm_means:
        return {"observed_mean_dist": obs_mean, "p_value": None,
                "n_iter": n_iter, "n_numt_with_peaks": len(obs_dists)}

    perm_arr = np.array(perm_means)
    p_value = float(np.mean(perm_arr <= obs_mean))

    return {
        "observed_mean_dist": round(obs_mean, 1),
        "permutation_mean": round(float(np.mean(perm_arr)), 1),
        "permutation_sd": round(float(np.std(perm_arr)), 1),
        "p_value": round(p_value, 4),
        "n_iter": n_iter,
        "n_numt_with_peaks": len(obs_dists),
    }


# ─── Main analysis ────────────────────────────────────────────────────────────


def analyse_all(
    peaks_path: Path,
    flanking_dir: Path = DEFAULT_FLANKING_DIR,
    chain_path: Path | None = None,
    n_permutations: int = DEFAULT_PERMUTATIONS,
) -> tuple[list[dict], dict]:
    """
    Full analysis: load peaks, liftover NUMT positions hg38→hg19, compute distances,
    run permutation test.

    Returns (per_numt_records, permutation_result).
    """
    print(f"Loading PRDM9 ChIP-seq peaks: {peaks_path.name} …", file=sys.stderr)
    peaks_genome = detect_peaks_genome(peaks_path)
    peaks_df = load_prdm9_peaks(peaks_path)
    print(f"  {len(peaks_df):,} peaks loaded across {peaks_df['chr'].nunique()} chromosomes "
          f"(detected build: {peaks_genome}).", file=sys.stderr)

    # Liftover setup
    lo = None
    if peaks_genome == "hg19":
        if chain_path and chain_path.exists():
            lo = build_liftover(chain_path)
            if lo:
                print(f"  Liftover loaded: hg38 → hg19 (chain: {chain_path.name})",
                      file=sys.stderr)
            else:
                print("  [WARN] pyliftover unavailable — run: pip install pyliftover",
                      file=sys.stderr)
        else:
            print(
                f"  [WARN] Peaks are hg19 but no chain file provided or found.\n"
                f"         Supply --chain <hg38ToHg19.over.chain.gz> for accurate distances.\n"
                f"         Without liftover, coordinates are INCOMPATIBLE — distances are wrong.",
                file=sys.stderr,
            )

    print(f"Loading NUMT flanking sequences: {flanking_dir} …", file=sys.stderr)
    numt_coords = load_all_numt_coords(flanking_dir)
    print(f"  {len(numt_coords)} NUMT candidates found (hg38 coordinates from FASTA headers).",
          file=sys.stderr)

    print("Lifting NUMT positions hg38 → hg19 and computing distances …", file=sys.stderr)
    records = []
    n_liftover_fail = 0
    for rec in numt_coords:
        chrom_hg38 = rec["chrom"]
        pos_hg38 = rec["insertion_pos"]

        # Liftover hg38 → hg19 if peaks are hg19
        chrom_hg19, pos_hg19 = chrom_hg38, pos_hg38
        liftover_status = "not_needed"
        if peaks_genome == "hg19" and lo is not None:
            lifted = liftover_position(lo, chrom_hg38, pos_hg38)
            if lifted:
                chrom_hg19, pos_hg19 = lifted
                liftover_status = "ok"
            else:
                liftover_status = "failed"
                n_liftover_fail += 1
        elif peaks_genome == "hg19" and lo is None:
            liftover_status = "skipped_no_chain"

        peak_info = nearest_peak_distance(chrom_hg19, pos_hg19, peaks_df)
        dist = peak_info["distance_bp"]
        if liftover_status in ("failed", "skipped_no_chain"):
            dist = None   # distance is unreliable without liftover
        classification = classify_prdm9(dist)

        records.append({
            "allele_id": rec["allele_id"],
            "chrom_hg38": chrom_hg38,
            "pos_hg38": pos_hg38,
            "chrom_hg19": chrom_hg19,
            "pos_hg19": pos_hg19,
            "liftover_status": liftover_status,
            "left_flank_bp": rec["left_len"],
            "right_flank_bp": rec["right_len"],
            "nearest_peak_center_hg19": peak_info["peak_center"],
            "distance_bp": dist,
            "peak_enrichment": peak_info["peak_enrichment"],
            "peak_pvalue": peak_info["peak_pvalue"],
            "prdm9_classification": classification,
        })

    if n_liftover_fail:
        print(f"  [WARN] {n_liftover_fail} positions failed liftover (set to NO_PRDM9).",
              file=sys.stderr)

    if n_permutations > 0:
        print(f"Running permutation test ({n_permutations} iterations) …", file=sys.stderr)
        perm_result = run_permutation_test(records, peaks_df, n_iter=n_permutations)
    else:
        perm_result = {"n_iter": 0, "p_value": None}

    return records, perm_result


# ─── Output formatters ────────────────────────────────────────────────────────

_TSV_FIELDS = [
    "allele_id",
    "chrom_hg38", "pos_hg38",
    "chrom_hg19", "pos_hg19", "liftover_status",
    "left_flank_bp", "right_flank_bp",
    "nearest_peak_center_hg19", "distance_bp",
    "peak_enrichment", "peak_pvalue",
    "prdm9_classification",
]


def format_tsv(records: list[dict]) -> str:
    lines = ["\t".join(_TSV_FIELDS)]
    for r in records:
        lines.append("\t".join(str(r.get(f, "")) for f in _TSV_FIELDS))
    return "\n".join(lines)


def format_json_output(records: list[dict], perm_result: dict) -> str:
    return json.dumps({"numt_records": records, "permutation_test": perm_result}, indent=2)


def print_console_summary(records: list[dict], perm_result: dict) -> None:
    """Print human-readable summary to stderr."""
    print("\n" + "=" * 60, file=sys.stderr)
    print("PRDM9 ChIP-seq proximity analysis — Wei et al. (2022) method", file=sys.stderr)
    print("=" * 60, file=sys.stderr)
    print(f"\nSource: Altemose et al. (2017) GSE99407 — {sum(1 for r in records if r['distance_bp'] is not None)} / {len(records)} NUMTs with peaks on same chromosome", file=sys.stderr)
    print(f"\nThresholds: SOMATIC < {PRDM9_SOMATIC_THRESHOLD_BP} bp | GERMINAL < {PRDM9_GERMINAL_THRESHOLD_BP} bp\n", file=sys.stderr)

    # Classification counts
    counts: dict[str, int] = {}
    for r in records:
        counts[r["prdm9_classification"]] = counts.get(r["prdm9_classification"], 0) + 1

    for cls in ["PRDM9_SOMATIC", "PRDM9_GERMINAL", "NO_PRDM9"]:
        n = counts.get(cls, 0)
        bar = "█" * n
        print(f"  {cls:<20} {n:>2}  {bar}", file=sys.stderr)

    print("\nPer-candidate distances (hg19 coordinates after liftover):", file=sys.stderr)
    print(f"  {'AlleleID':<12} {'Chr(hg19)':<10} {'Pos(hg19)':>12} {'NearPeak':>12} "
          f"{'Dist(bp)':>10}  {'Liftover':<10}  Classification",
          file=sys.stderr)
    print("  " + "-" * 85, file=sys.stderr)
    for r in sorted(records, key=lambda x: x.get("distance_bp") or math.inf):
        dist_str = str(r["distance_bp"]) if r["distance_bp"] is not None else "N/A"
        pos_hg19 = r.get("pos_hg19", "?")
        pos_str = f"{pos_hg19:,}" if isinstance(pos_hg19, int) else str(pos_hg19)
        print(
            f"  {r['allele_id']:<12} {r.get('chrom_hg19', '?'):<10} {pos_str:>12} "
            f"{str(r.get('nearest_peak_center_hg19') or ''):>12}  {dist_str:>8}  "
            f"{r.get('liftover_status', '?'):<10}  {r['prdm9_classification']}",
            file=sys.stderr,
        )

    # Permutation test
    if perm_result.get("p_value") is not None:
        print(f"\nPermutation test ({perm_result['n_iter']} iterations):", file=sys.stderr)
        print(f"  Observed mean distance : {perm_result['observed_mean_dist']:,.0f} bp", file=sys.stderr)
        print(f"  Permutation mean ± SD  : {perm_result['permutation_mean']:,.0f} ± {perm_result['permutation_sd']:,.0f} bp", file=sys.stderr)
        print(f"  P-value (one-tailed)   : {perm_result['p_value']:.4f}", file=sys.stderr)
        if perm_result["p_value"] < 0.05:
            print("  → NUMTs are significantly closer to PRDM9 peaks than expected by chance.", file=sys.stderr)
        else:
            print("  → No significant enrichment near PRDM9 peaks (consistent with NHEJ dominant).", file=sys.stderr)
    print(file=sys.stderr)


# ─── CLI ─────────────────────────────────────────────────────────────────────


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "PRDM9 ChIP-seq peak proximity analysis for NUMT candidates. "
            "Implements Wei et al. (2022) methodology using Altemose et al. (2017) data."
        )
    )
    parser.add_argument(
        "--peaks", required=True, type=Path,
        metavar="FILE",
        help="GSE99407 ChIP-seq peak file (.txt or .txt.gz).",
    )
    parser.add_argument(
        "--chain", type=Path, default=DEFAULT_CHAIN,
        metavar="FILE",
        help=(
            f"hg38→hg19 liftover chain file (default: {DEFAULT_CHAIN}). "
            "Required when peaks are hg19 (GSE99407). Download with --download-chain."
        ),
    )
    parser.add_argument(
        "--download-chain", type=Path, default=None,
        metavar="DEST",
        help="Download hg38ToHg19.over.chain.gz from UCSC to DEST and exit.",
    )
    parser.add_argument(
        "--flanking-dir", type=Path, default=DEFAULT_FLANKING_DIR,
        metavar="DIR",
        help=f"Directory with *_flanking.fasta files (default: {DEFAULT_FLANKING_DIR}).",
    )
    parser.add_argument(
        "--permutations", type=int, default=DEFAULT_PERMUTATIONS,
        metavar="N",
        help=f"Number of permutation iterations (default: {DEFAULT_PERMUTATIONS}; 0 = skip).",
    )
    parser.add_argument(
        "--format", choices=["tsv", "json"], default=DEFAULT_FORMAT,
        help=f"Output format (default: {DEFAULT_FORMAT}).",
    )
    parser.add_argument(
        "--no-summary", action="store_true",
        help="Suppress human-readable summary on stderr.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    # Handle --download-chain
    if args.download_chain is not None:
        download_chain(args.download_chain)
        return 0

    peaks_path = args.peaks
    if not peaks_path.exists():
        print(f"Error: peaks file not found: {peaks_path}", file=sys.stderr)
        return 1

    flanking_dir = args.flanking_dir
    if not flanking_dir.is_dir():
        print(f"Error: flanking directory not found: {flanking_dir}", file=sys.stderr)
        return 1

    chain_path = args.chain

    records, perm_result = analyse_all(
        peaks_path=peaks_path,
        flanking_dir=flanking_dir,
        chain_path=chain_path,
        n_permutations=args.permutations,
    )

    if not args.no_summary:
        print_console_summary(records, perm_result)

    if args.format == "tsv":
        print(format_tsv(records))
    else:
        print(format_json_output(records, perm_result))

    return 0


if __name__ == "__main__":
    sys.exit(main())
