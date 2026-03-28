#!/usr/bin/env python3
"""
compare_blast_last.py — Side-by-side comparison of BLAST vs LAST mito results.

Joins on allele_id (best hit per variant), emits a wide TSV with status
(BOTH / BLAST_ONLY / LAST_ONLY) and a text report summarizing discrepancies.

Usage:
  python 01_code/numt_pipeline/last_validation/compare_blast_last.py \
    --blast 02_data/processed/excel_blast_mito.tsv \
    --last  02_data/processed/excel_last_mito.tsv \
    --output-tsv    02_data/processed/blast_vs_last_comparison.tsv \
    --output-report 02_data/processed/blast_vs_last_report.txt
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from run_pipeline_excel import log


def _best_hit_per_allele(df: pd.DataFrame) -> pd.DataFrame:
    """Keep the best hit (highest combined_score) per allele_id."""
    if df.empty or "allele_id" not in df.columns:
        return pd.DataFrame()
    df = df.copy()
    df["allele_id"] = pd.to_numeric(df["allele_id"], errors="coerce")
    df = df.dropna(subset=["allele_id"])
    df["allele_id"] = df["allele_id"].astype(int)
    if "combined_score" in df.columns:
        df["combined_score"] = pd.to_numeric(df["combined_score"], errors="coerce")
        idx = df.groupby("allele_id")["combined_score"].idxmax()
        return df.loc[idx].reset_index(drop=True)
    return df.drop_duplicates(subset=["allele_id"], keep="first").reset_index(drop=True)


def compare(blast_tsv: Path, last_tsv: Path, output_tsv: Path, output_report: Path) -> pd.DataFrame:
    """Build comparison TSV and text report."""
    # Load
    df_blast = pd.read_csv(blast_tsv, sep="\t") if blast_tsv.exists() else pd.DataFrame()
    df_last = pd.read_csv(last_tsv, sep="\t") if last_tsv.exists() else pd.DataFrame()

    blast_best = _best_hit_per_allele(df_blast)
    last_best = _best_hit_per_allele(df_last)

    blast_ids = set(blast_best["allele_id"]) if not blast_best.empty else set()
    last_ids = set(last_best["allele_id"]) if not last_best.empty else set()
    all_ids = sorted(blast_ids | last_ids)

    if not all_ids:
        log("No hits in either BLAST or LAST — nothing to compare.")
        return pd.DataFrame()

    # Build wide table
    rows = []
    for aid in all_ids:
        b = blast_best[blast_best["allele_id"] == aid].iloc[0] if aid in blast_ids else None
        l = last_best[last_best["allele_id"] == aid].iloc[0] if aid in last_ids else None

        gene = (b["gene_symbol"] if b is not None else l["gene_symbol"]) if True else ""

        if b is not None and l is not None:
            status = "BOTH"
        elif b is not None:
            status = "BLAST_ONLY"
        else:
            status = "LAST_ONLY"

        blast_tier = str(b["hit_tier"]) if b is not None else ""
        last_tier = str(l["hit_tier"]) if l is not None else ""
        tier_change = ""
        if blast_tier and last_tier and blast_tier != last_tier:
            tier_change = f"{blast_tier}→{last_tier}"

        blast_region = str(b["mtdna_region"]) if b is not None else ""
        last_region = str(l["mtdna_region"]) if l is not None else ""
        region_change = blast_region != last_region and blast_region and last_region

        rows.append({
            "allele_id": aid,
            "gene_symbol": gene,
            "blast_pident": round(float(b["pident"]), 2) if b is not None else "",
            "blast_length": int(b["length"]) if b is not None else "",
            "blast_hit_tier": blast_tier,
            "blast_mtdna_region": blast_region,
            "blast_combined_score": round(float(b["combined_score"]), 4) if b is not None else "",
            "last_pident": round(float(l["pident"]), 2) if l is not None else "",
            "last_length": int(l["length"]) if l is not None else "",
            "last_hit_tier": last_tier,
            "last_mtdna_region": last_region,
            "last_combined_score": round(float(l["combined_score"]), 4) if l is not None else "",
            "last_is_split": bool(l["is_split"]) if l is not None and "is_split" in l.index else "",
            "last_mismap": round(float(l["mismap"]), 4) if l is not None and "mismap" in l.index and pd.notna(l.get("mismap")) else "",
            "status": status,
            "tier_change": tier_change,
            "region_change": region_change,
        })

    df_cmp = pd.DataFrame(rows)

    # Save TSV
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    df_cmp.to_csv(output_tsv, sep="\t", index=False)
    log(f"Comparison TSV saved: {output_tsv}")

    # Generate text report
    report_lines = []
    report_lines.append("=" * 65)
    report_lines.append("BLAST vs LAST — Comparison Report")
    report_lines.append("=" * 65)
    report_lines.append("")

    n_both = len(df_cmp[df_cmp["status"] == "BOTH"])
    n_blast_only = len(df_cmp[df_cmp["status"] == "BLAST_ONLY"])
    n_last_only = len(df_cmp[df_cmp["status"] == "LAST_ONLY"])

    report_lines.append(f"  BOTH (concordant)   : {n_both}")
    report_lines.append(f"  BLAST_ONLY          : {n_blast_only}")
    report_lines.append(f"  LAST_ONLY           : {n_last_only}")
    report_lines.append(f"  Total variants      : {len(df_cmp)}")
    report_lines.append("")

    # Tier changes
    tier_changes = df_cmp[df_cmp["tier_change"] != ""]
    if not tier_changes.empty:
        report_lines.append(f"Tier changes ({len(tier_changes)}):")
        for _, r in tier_changes.iterrows():
            report_lines.append(
                f"  AlleleID {r['allele_id']:>8}  {r['gene_symbol']:<15}  {r['tier_change']}"
            )
        report_lines.append("")

    # Top LAST_ONLY candidates
    last_only = df_cmp[df_cmp["status"] == "LAST_ONLY"].copy()
    if not last_only.empty:
        last_only["last_combined_score"] = pd.to_numeric(last_only["last_combined_score"], errors="coerce")
        last_only = last_only.sort_values("last_combined_score", ascending=False)
        report_lines.append(f"Top LAST_ONLY candidates (max 10):")
        for _, r in last_only.head(10).iterrows():
            report_lines.append(
                f"  AlleleID {r['allele_id']:>8}  {r['gene_symbol']:<15}  "
                f"pident={r['last_pident']}  region={r['last_mtdna_region']}  "
                f"score={r['last_combined_score']}"
            )
        report_lines.append("")

    # Top chimeric (split) candidates
    chimeric = df_cmp[df_cmp["last_is_split"] == True].copy()
    if not chimeric.empty:
        chimeric["last_combined_score"] = pd.to_numeric(chimeric["last_combined_score"], errors="coerce")
        chimeric = chimeric.sort_values("last_combined_score", ascending=False)
        report_lines.append(f"Chimeric candidates (is_split=True, max 10):")
        for _, r in chimeric.head(10).iterrows():
            report_lines.append(
                f"  AlleleID {r['allele_id']:>8}  {r['gene_symbol']:<15}  "
                f"pident={r['last_pident']}  region={r['last_mtdna_region']}  "
                f"score={r['last_combined_score']}  mismap={r['last_mismap']}"
            )
        report_lines.append("")

    report_text = "\n".join(report_lines)
    output_report.parent.mkdir(parents=True, exist_ok=True)
    output_report.write_text(report_text)
    log(f"Comparison report saved: {output_report}")
    print(report_text)

    return df_cmp


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare BLAST vs LAST mitochondrial alignment results.",
    )
    parser.add_argument("--blast", type=Path, required=True, help="BLAST mito TSV")
    parser.add_argument("--last", type=Path, required=True, help="LAST mito TSV")
    parser.add_argument("--output-tsv", type=Path, required=True, help="Output comparison TSV")
    parser.add_argument("--output-report", type=Path, required=True, help="Output text report")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    compare(args.blast, args.last, args.output_tsv, args.output_report)
