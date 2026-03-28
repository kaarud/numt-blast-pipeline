"""
show_alignments.py — ASCII alignment viewer for NUMT genome BLAST results.

For each of the 21 final candidates (20 NOVEL_NUMT_CANDIDATE + 1 CONFIRMED_NUMT),
reads the cached genome BLAST XML and renders the best hit alignment.

Usage:
    python 01_code/numt_pipeline/show_alignments.py
    python 01_code/numt_pipeline/show_alignments.py --width 100
    python 01_code/numt_pipeline/show_alignments.py --allele 2134240
"""

import argparse
import sys
from pathlib import Path
from xml.etree import ElementTree as ET

# ---------------------------------------------------------------------------
# Candidate metadata (from pipeline classification)
# ---------------------------------------------------------------------------

CANDIDATES = [
    # AlleleID, Gene, ClinSig, Class
    (96448,   "MSH2",    "Pathogenic",             "CONFIRMED_NUMT"),
    (2068854, "GNE",     "Pathogenic",             "NOVEL_NUMT_CANDIDATE"),
    (3780309, "KIT",     "Uncertain significance",  "NOVEL_NUMT_CANDIDATE"),
    (3853582, "CBLIF",   "Uncertain significance",  "NOVEL_NUMT_CANDIDATE"),
    (1679511, "ECM1",    "Pathogenic",             "NOVEL_NUMT_CANDIDATE"),
    (3791414, "ATM",     "Pathogenic",             "NOVEL_NUMT_CANDIDATE"),
    (3788914, "MYH6",    "Uncertain significance",  "NOVEL_NUMT_CANDIDATE"),
    (2750920, "APC",     "Pathogenic",             "NOVEL_NUMT_CANDIDATE"),
    (2853236, "COL6A2",  "Likely pathogenic",       "NOVEL_NUMT_CANDIDATE"),
    (2150417, "DCLRE1C", "Pathogenic",             "NOVEL_NUMT_CANDIDATE"),
    (2153648, "SGSH",    "Pathogenic",             "NOVEL_NUMT_CANDIDATE"),
    (2134240, "RARS2",   "Pathogenic",             "NOVEL_NUMT_CANDIDATE"),
    (2093557, "RAPSN",   "Pathogenic",             "NOVEL_NUMT_CANDIDATE"),
    (2183837, "DCLRE1C", "Pathogenic",             "NOVEL_NUMT_CANDIDATE"),
    (3030219, "SACS",    "Likely benign",           "NOVEL_NUMT_CANDIDATE"),
    (2163625, "SGSH",    "Pathogenic",             "NOVEL_NUMT_CANDIDATE"),
    (3074437, "MCOLN1",  "Likely pathogenic",       "NOVEL_NUMT_CANDIDATE"),
    (3034753, "UNC13D",  "Pathogenic",             "NOVEL_NUMT_CANDIDATE"),
    (1064280, "COL1A1",  "Pathogenic",             "NOVEL_NUMT_CANDIDATE"),
    (4837144, "CCDC88A", "Pathogenic",             "NOVEL_NUMT_CANDIDATE"),
    (1444905, "ZEB2",    "Pathogenic",             "NOVEL_NUMT_CANDIDATE"),
]

CACHE_DIR = Path(__file__).parent.parent.parent / "02_data" / "processed" / "blast_genome_cache"

CLINSIG_SYMBOL = {
    "Pathogenic":            "[P]",
    "Likely pathogenic":     "[LP]",
    "Uncertain significance": "[VUS]",
    "Likely benign":         "[LB]",
    "Benign":                "[B]",
}


def parse_hits(xml_path: Path) -> list[dict]:
    """Parse BLAST XML and return a list of hit dicts, sorted by bitscore desc."""
    tree = ET.parse(xml_path)
    root = tree.getroot()
    hits = []
    for hit in root.iter("Hit"):
        hit_def = hit.findtext("Hit_def", "")
        hit_acc = hit.findtext("Hit_accession", "")
        hit_len = int(hit.findtext("Hit_len", 0))
        for hsp in hit.iter("Hsp"):
            hits.append({
                "hit_def":   hit_def,
                "hit_acc":   hit_acc,
                "hit_len":   hit_len,
                "bitscore":  float(hsp.findtext("Hsp_bit-score", 0)),
                "evalue":    float(hsp.findtext("Hsp_evalue", 1)),
                "pident":    100 * int(hsp.findtext("Hsp_identity", 0)) / max(1, int(hsp.findtext("Hsp_align-len", 1))),
                "identity":  int(hsp.findtext("Hsp_identity", 0)),
                "gaps":      int(hsp.findtext("Hsp_gaps", 0)),
                "align_len": int(hsp.findtext("Hsp_align-len", 0)),
                "q_from":    int(hsp.findtext("Hsp_query-from", 0)),
                "q_to":      int(hsp.findtext("Hsp_query-to", 0)),
                "h_from":    int(hsp.findtext("Hsp_hit-from", 0)),
                "h_to":      int(hsp.findtext("Hsp_hit-to", 0)),
                "qseq":      hsp.findtext("Hsp_qseq", ""),
                "hseq":      hsp.findtext("Hsp_hseq", ""),
                "midline":   hsp.findtext("Hsp_midline", ""),
            })
    hits.sort(key=lambda h: -h["bitscore"])
    return hits


def is_nuclear_hit(hit_def: str) -> bool:
    """Heuristic: nuclear if the definition mentions a chromosome or genomic clone."""
    low = hit_def.lower()
    nuclear_keywords = ["chromosome", "clone", "pseudogene", "genomic"]
    mito_keywords = ["mitochondri", "mitochondrion"]
    if any(k in low for k in mito_keywords):
        return False
    if any(k in low for k in nuclear_keywords):
        return True
    return False


def format_alignment(hit: dict, width: int = 80, label: str = "hit") -> str:
    """Format a BLAST HSP as a multi-line ASCII alignment block."""
    qseq    = hit["qseq"]
    midline = hit["midline"]
    hseq    = hit["hseq"]
    align_len = len(qseq)

    q_pos = hit["q_from"]
    h_pos = hit["h_from"]
    h_strand = 1 if hit["h_from"] <= hit["h_to"] else -1

    lines = []
    chunk = width

    for i in range(0, align_len, chunk):
        qs = qseq[i:i+chunk]
        ms = midline[i:i+chunk]
        hs = hseq[i:i+chunk]

        # Count consumed (non-gap) bases
        q_consumed = sum(1 for c in qs if c != "-")
        h_consumed = sum(1 for c in hs if c != "-")

        q_end = q_pos + q_consumed - 1
        h_end = h_pos + (h_consumed - 1) * h_strand

        lines.append(f"Query  {q_pos:>8}  {qs}  {q_end}")
        lines.append(f"       {'':8}  {ms}")
        lines.append(f"{label:<6} {h_pos:>8}  {hs}  {h_end}")
        lines.append("")

        q_pos = q_end + 1
        h_pos = h_end + h_strand

    return "\n".join(lines)


def render_candidate(allele_id: int, gene: str, clinsig: str, numt_class: str,
                     width: int) -> str:
    xml_path = CACHE_DIR / f"{allele_id}.xml"
    if not xml_path.exists():
        return f"[ERROR] XML not found: {xml_path}\n"

    hits = parse_hits(xml_path)
    if not hits:
        return f"[NO HITS] AlleleID {allele_id} ({gene}) — no BLAST hits in XML\n"

    sym = CLINSIG_SYMBOL.get(clinsig, "")
    sep = "═" * (width + 30)
    out = [sep]
    out.append(f"  AlleleID : {allele_id}   Gene : {gene}   {sym} {clinsig}   Class : {numt_class}")

    if numt_class == "CONFIRMED_NUMT":
        # Show best nuclear hit (if any), then best mito hit
        nuclear_hits = [h for h in hits if is_nuclear_hit(h["hit_def"])]
        mito_hits    = [h for h in hits if not is_nuclear_hit(h["hit_def"])]

        if nuclear_hits:
            best = nuclear_hits[0]
            out.append(f"  NUCLEAR HIT → {best['hit_def'][:80]}")
            out.append(f"  Identity: {best['pident']:.1f}%  ({best['identity']}/{best['align_len']} bp)  "
                       f"Gaps: {best['gaps']}  E-value: {best['evalue']:.2e}  Bitscore: {best['bitscore']:.1f}")
            out.append(f"  Query pos: {best['q_from']}–{best['q_to']}  |  Hit pos: {best['h_from']}–{best['h_to']}")
            out.append("")
            out.append(format_alignment(best, width=width, label="Nuclear"))
        else:
            out.append("  [No nuclear hit found]")
            out.append("")

        if mito_hits:
            best = mito_hits[0]
            out.append(f"  MITO HIT → {best['hit_def'][:80]}")
            out.append(f"  Identity: {best['pident']:.1f}%  ({best['identity']}/{best['align_len']} bp)  "
                       f"Gaps: {best['gaps']}  E-value: {best['evalue']:.2e}  Bitscore: {best['bitscore']:.1f}")
            out.append(f"  Query pos: {best['q_from']}–{best['q_to']}  |  Hit pos: {best['h_from']}–{best['h_to']}")
            out.append("")
            out.append(format_alignment(best, width=width, label="Mito"))

    else:
        # NOVEL — best hit (should be mito only)
        best = hits[0]
        out.append(f"  Best hit → {best['hit_def'][:80]}")
        out.append(f"  Identity: {best['pident']:.1f}%  ({best['identity']}/{best['align_len']} bp)  "
                   f"Gaps: {best['gaps']}  E-value: {best['evalue']:.2e}  Bitscore: {best['bitscore']:.1f}")
        out.append(f"  Query pos: {best['q_from']}–{best['q_to']}  |  Hit pos: {best['h_from']}–{best['h_to']}")
        if len(hits) > 1:
            out.append(f"  (Total hits in XML: {len(hits)} — showing best only)")
        out.append("")
        out.append(format_alignment(best, width=width, label="Subj"))

    out.append(sep)
    out.append("")
    return "\n".join(out)


def main():
    parser = argparse.ArgumentParser(description="ASCII BLAST alignment viewer for NUMT candidates")
    parser.add_argument("--width", type=int, default=80, help="Alignment line width (default: 80)")
    parser.add_argument("--allele", type=int, help="Show only this AlleleID")
    args = parser.parse_args()

    candidates = CANDIDATES
    if args.allele:
        candidates = [c for c in candidates if c[0] == args.allele]
        if not candidates:
            print(f"AlleleID {args.allele} not found in candidate list.", file=sys.stderr)
            sys.exit(1)

    for allele_id, gene, clinsig, numt_class in candidates:
        print(render_candidate(allele_id, gene, clinsig, numt_class, args.width))


if __name__ == "__main__":
    main()
