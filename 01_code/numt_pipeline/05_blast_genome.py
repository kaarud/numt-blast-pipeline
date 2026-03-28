#!/usr/bin/env python3
"""
05_blast_genome.py — BLAST mitochondrial-origin sequences against the full human genome.

WHY WE BLAST AGAINST THE GENOME AFTER THE MITOCHONDRIAL GENOME
===============================================================
In Step 4 we identified ClinVar insertion variants whose sequences match the human
mitochondrial genome (rCRS, NC_012920.1).  Those "mito hits" tell us that the inserted
DNA is of mitochondrial origin — but they do NOT tell us where that DNA ended up in the
nuclear genome.  That is what this step answers.

By BLASTing the same sequences against the full human genome (hg38/GRCh38), we achieve
two complementary goals:

  1. CONFIRM THE INSERTION LOCATION:
     If a sequence BLASTs to both the mitochondrial genome AND to the nuclear position
     reported in ClinVar (e.g., chromosome 16 for TSC2), we have strong evidence that
     the ClinVar variant is a genuine NUMT insertion at the reported locus.  This is
     the critical confirmation for the TSC2 index case: a ~72 bp fragment from mtDNA
     that BLASTs to both chrM AND chr16p13.3 (TSC2 exon 17).

  2. FIND EXISTING (REFERENCE) NUMTs:
     The human reference genome (hg38) already contains ~755 known NUMTs — fragments
     of mitochondrial DNA that were integrated into nuclear chromosomes over the course
     of primate evolution.  These are "reference NUMTs" and are considered benign: they
     are fixed in the human population and do not disrupt gene function.

     When our query sequence hits a nuclear chromosome at a location OTHER than the
     ClinVar-reported position, it likely means that this particular mitochondrial
     region has been independently integrated elsewhere in the genome.  This is
     informative: if many reference NUMTs from the same mtDNA region exist, it suggests
     that region is prone to nuclear integration (a "NUMT hotspot").

WHAT IS A "KNOWN NUMT" vs A "NOVEL INSERTION"?
===============================================
  - Known/Reference NUMT:  A nuclear copy of mtDNA that is ALREADY present in the hg38
    reference assembly.  It was integrated long ago (often millions of years), is shared
    across most or all humans, and is benign.  BLAST will find these because the
    reference genome includes them.

  - Novel NUMT insertion:  A mitochondrial fragment that is NOT in the reference genome.
    It was inserted somatically or in the recent germline, disrupting a gene.  This is
    what causes disease (e.g., the TSC2 NUMT).  BLAST against the genome will NOT find
    these at the insertion site, because the reference lacks them.  However, BLAST WILL
    find the mitochondrial chromosome hit (chrM), confirming the sequence is mtDNA.

  - The classification logic:
      CONFIRMED_NUMT         = mtDNA hit (Step 4) + nuclear hit at ClinVar position
      KNOWN_REFERENCE_NUMT   = nuclear hit(s) NOT at the ClinVar position (reference NUMTs)
      NOVEL_NUMT_CANDIDATE   = mtDNA hit but NO nuclear hit (not in reference = novel)
      AMBIGUOUS               = nuclear hit but poor mtDNA alignment from Step 4

WHY chrMT / chrM HITS ARE EXPECTED
===================================
When BLASTing against the full genome (which includes the mitochondrial chromosome as
chrM or NC_012920.1), we EXPECT to see a mitochondrial hit.  This is not surprising —
it confirms the sequence's mitochondrial origin.  The interesting hits are the NUCLEAR
ones: chr1–22, X, Y.

LOCAL vs REMOTE BLAST
=====================
The full hg38 genome BLAST database is ~15 GB uncompressed, which is impractical for
many users.  This script offers two modes:

  --mode remote (default):
    Uses NCBI's remote BLAST service via BioPython's qblast() function.  Slower
    (each query takes 30–120 seconds) but requires no local database.  We filter
    results to Homo sapiens using an Entrez query.

    NCBI RATE LIMIT POLICY: NCBI requests that automated tools submit no more than
    one query every 10 seconds.  Exceeding this may result in temporary IP bans.
    We enforce a 10-second sleep between queries and cache all raw XML results so
    re-runs do not re-submit queries.

  --mode local:
    Uses a local BLAST database that the user provides (--genome-db).  Much faster
    (~1 second per query) but requires the user to download and build the hg38
    BLAST database beforehand:
      wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
      makeblastdb -in GCF_000001405.40_GRCh38.p14_genomic.fna -dbtype nucl

WHY CACHING RAW XML MATTERS
============================
NCBI remote BLAST is slow and rate-limited.  If the script crashes or is interrupted
after 30 successful queries, we do NOT want to re-submit those 30 queries.  By saving
each raw XML result to disk (one file per AlleleID), we can resume where we left off.
This is especially important for large batches: 100 queries at 10s apart = ~17 minutes
just in mandatory wait time.
"""

# ─── Standard-library imports ────────────────────────────────────────────────
import argparse
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

# ─── Third-party imports ─────────────────────────────────────────────────────
import pandas as pd

try:
    from Bio.Blast import NCBIWWW, NCBIXML
    BIOPYTHON_BLAST_AVAILABLE = True
except ImportError:
    BIOPYTHON_BLAST_AVAILABLE = False

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

# Outputs
QUERY_FASTA = PROJECT_ROOT / "02_data" / "processed" / "blast_genome_query.fasta"
RAW_XML_DIR = PROJECT_ROOT / "02_data" / "processed" / "blast_genome_results_raw"
LOCAL_BLAST_OUTPUT = PROJECT_ROOT / "02_data" / "processed" / "blast_genome_results.tsv"
SUMMARY_TSV = PROJECT_ROOT / "02_data" / "processed" / "numt_candidates_summary.tsv"

# BLAST tabular output columns for local mode (outfmt 6)
GENOME_BLAST_COLUMNS = [
    "qseqid",       # Query sequence ID
    "sseqid",       # Subject sequence ID (chromosome/accession)
    "pident",       # Percentage of identical matches
    "length",       # Alignment length (bp)
    "mismatch",     # Number of mismatches
    "gapopen",      # Number of gap openings
    "qstart",       # Start of alignment in query
    "qend",         # End of alignment in query
    "sstart",       # Start of alignment in subject (genome position)
    "send",         # End of alignment in subject
    "evalue",       # Expect value
    "bitscore",     # Bit score
    "qlen",         # Query length
    "slen",         # Subject length
]

# ─── Chromosome name normalisation ──────────────────────────────────────────
# NCBI accessions for GRCh38 primary assembly chromosomes.
# BLAST results from NCBI remote may report accessions (NC_000001.11) or
# chromosome names (chr1).  We normalise everything to "chr1", "chrX", "chrM" etc.

ACCESSION_TO_CHR = {
    "NC_000001.11": "chr1",  "NC_000002.12": "chr2",  "NC_000003.12": "chr3",
    "NC_000004.12": "chr4",  "NC_000005.10": "chr5",  "NC_000006.12": "chr6",
    "NC_000007.14": "chr7",  "NC_000008.11": "chr8",  "NC_000009.12": "chr9",
    "NC_000010.11": "chr10", "NC_000011.10": "chr11", "NC_000012.12": "chr12",
    "NC_000013.11": "chr13", "NC_000014.9": "chr14",  "NC_000015.10": "chr15",
    "NC_000016.10": "chr16", "NC_000017.11": "chr17", "NC_000018.10": "chr18",
    "NC_000019.10": "chr19", "NC_000020.11": "chr20", "NC_000021.9": "chr21",
    "NC_000022.11": "chr22", "NC_000023.11": "chrX",  "NC_000024.10": "chrY",
    "NC_012920.1": "chrM",
}


# ─── Helper: timestamped logging ────────────────────────────────────────────

def log(message: str) -> None:
    """Print a message prefixed with the current timestamp."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")


# ─── Helper: normalise chromosome name ──────────────────────────────────────

def normalise_chromosome(subject_id: str) -> str:
    """
    Convert a BLAST subject ID to a standardised chromosome name.

    Handles multiple formats returned by NCBI BLAST:
      - NC_000016.10          -> chr16
      - gi|568336023|ref|NC_000016.10|  -> chr16
      - chr16                 -> chr16
      - NC_012920.1           -> chrM

    Parameters
    ----------
    subject_id : str
        The sseqid field from BLAST output.

    Returns
    -------
    str
        Normalised chromosome name (e.g., "chr16", "chrM") or the original
        subject_id if it cannot be mapped.
    """
    # Already in chrN format
    if subject_id.startswith("chr"):
        return subject_id

    # Try direct accession lookup
    if subject_id in ACCESSION_TO_CHR:
        return ACCESSION_TO_CHR[subject_id]

    # Try extracting accession from gi|...|ref|NC_XXXXXX.N| format
    match = re.search(r"(NC_\d+\.\d+)", subject_id)
    if match:
        accession = match.group(1)
        if accession in ACCESSION_TO_CHR:
            return ACCESSION_TO_CHR[accession]

    # Fallback: return as-is (could be an unplaced scaffold, patch, etc.)
    return subject_id


# ─── Helper: ClinVar chromosome from variant document ───────────────────────

def get_clinvar_chromosome(doc: dict) -> str | None:
    """
    Extract the chromosome from a ClinVar variant MongoDB document.

    ClinVar stores the chromosome in the 'chromosome' field, or it can
    be inferred from the GRCh38 location string.

    Returns
    -------
    str or None
        Normalised chromosome name (e.g., "chr16") or None if unavailable.
    """
    chrom = doc.get("chromosome")
    if chrom:
        if not chrom.startswith("chr"):
            chrom = f"chr{chrom}"
        return chrom

    # Try extracting from GRCh38 location
    location = doc.get("grch38_location") or doc.get("location") or ""
    match = re.match(r"(\d+|X|Y|MT?)", str(location))
    if match:
        chrom = match.group(1)
        if chrom in ("M", "MT"):
            return "chrM"
        return f"chr{chrom}"

    return None


# ─── Step 1: Prepare query FASTA from MongoDB ───────────────────────────────

def prepare_query_fasta(host: str, port: int, output_fasta: Path) -> tuple[int, dict]:
    """
    Fetch variants with mitochondrial hits from MongoDB and write query FASTA.

    We query two sets of variants:
      1. Variants with has_mito_hit=True — they passed the Step 4 screen and their
         inserted sequence matches mitochondrial DNA.
      2. DIRECT_MITO_REF entries — their HGVS name explicitly references mtDNA
         coordinates, so they are confirmed NUMTs.  We want to find their genome
         positions.

    The FASTA header includes metadata to carry through the BLAST pipeline:
        >AlleleID|GeneSymbol|ClinSig|MitoHitPident|MitoRegion

    Parameters
    ----------
    host : str
        MongoDB hostname.
    port : int
        MongoDB port.
    output_fasta : Path
        Where to write the query FASTA.

    Returns
    -------
    tuple[int, dict]
        (number_of_sequences_written, {allele_id: variant_doc_dict})
    """
    if not PYMONGO_AVAILABLE:
        log("FATAL -- pymongo is not installed.  pip install pymongo")
        sys.exit(1)

    log("Connecting to MongoDB to fetch mito-positive variants ...")

    try:
        client = pymongo.MongoClient(
            host=host, port=port, serverSelectionTimeoutMS=5000
        )
        client.admin.command("ping")
    except (ConnectionFailure, ServerSelectionTimeoutError) as exc:
        log(f"FATAL -- Cannot connect to MongoDB at {host}:{port}")
        log(f"Error: {exc}")
        log("Is MongoDB running?  brew services start mongodb-community")
        sys.exit(1)

    db = client[MONGO_DB]
    collection = db[MONGO_COLLECTION]

    # ── Query: variants with mitochondrial hits from Step 4 ──────────────
    query = {
        "$or": [
            {
                "has_mito_hit": True,
                "extracted_sequence": {"$ne": None},
            },
            {
                "extraction_status": "DIRECT_MITO_REF",
                "extracted_sequence": {"$ne": None},
            },
        ]
    }

    docs = list(collection.find(query))
    log(f"Found {len(docs):,} variants with mitochondrial hits or DIRECT_MITO_REF status")

    # ── Write FASTA ─────────────────────────────────────────────────────
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    variant_docs = {}
    n_written = 0

    with open(output_fasta, "w") as fh:
        for doc in docs:
            seq = doc.get("extracted_sequence", "")
            if not seq or not isinstance(seq, str):
                continue

            allele_id = doc.get("allele_id", "unknown")
            gene = doc.get("gene_symbol", "NA") or "NA"
            clinsig = doc.get("clinical_significance", "NA") or "NA"
            clinsig_safe = clinsig.replace("|", "/").replace(" ", "_")

            # Mitochondrial hit information from Step 4
            best_pident = doc.get("best_mito_hit_pident", "NA")
            if best_pident is not None and best_pident != "NA":
                best_pident = f"{best_pident:.1f}"
            else:
                best_pident = "NA"

            # Best mitochondrial region from Step 4 results
            mito_region = "NA"
            mito_results = doc.get("blast_mito_results", [])
            if mito_results:
                # Take the region from the best (first) hit
                mito_region = mito_results[0].get("mtdna_region", "NA")

            header = f">{allele_id}|{gene}|{clinsig_safe}|{best_pident}|{mito_region}"
            fh.write(f"{header}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")

            variant_docs[int(allele_id)] = doc
            n_written += 1

    client.close()

    if n_written == 0:
        log("WARNING -- No sequences written to FASTA.")
        log("Have you run Step 4 (04_blast_mito.py) first?")
        log("Check that has_mito_hit=True documents exist in MongoDB.")
    else:
        log(f"Wrote {n_written:,} sequences to: {output_fasta}")

    return n_written, variant_docs


# ─── Step 2a: Remote BLAST via NCBI qblast ──────────────────────────────────

def blast_remote_single(
    allele_id: int,
    sequence: str,
    xml_dir: Path,
    hitlist_size: int = 20,
    evalue: float = 1e-5,
) -> Path | None:
    """
    Submit a single sequence to NCBI remote BLAST and save the XML result.

    Uses Bio.Blast.NCBIWWW.qblast() to search the "nt" nucleotide database,
    filtered to Homo sapiens.

    CACHING: If the XML output file already exists, the query is skipped.
    This allows safe re-runs without re-submitting to NCBI.

    RETRY LOGIC: NCBI connections can be flaky.  We retry up to 3 times
    with exponential backoff (10s, 30s, 90s) before giving up on a query.

    Parameters
    ----------
    allele_id : int
        The ClinVar AlleleID (used for file naming).
    sequence : str
        The nucleotide sequence to BLAST.
    xml_dir : Path
        Directory for saving raw XML results.
    hitlist_size : int
        Maximum number of hits to return per query.
    evalue : float
        E-value threshold.

    Returns
    -------
    Path or None
        Path to the XML result file, or None if the query failed.
    """
    if not BIOPYTHON_BLAST_AVAILABLE:
        log("FATAL -- BioPython BLAST modules not available.  pip install biopython")
        sys.exit(1)

    xml_path = xml_dir / f"{allele_id}.xml"

    # ── Cache check: skip if XML already exists and is non-empty ─────────
    # This is critical for resumability.  NCBI rate limits mean re-running
    # the full batch is expensive in time.
    if xml_path.exists() and xml_path.stat().st_size > 100:
        log(f"  [CACHED] AlleleID {allele_id} — XML already exists, skipping query")
        return xml_path

    # ── Submit to NCBI with retry logic ──────────────────────────────────
    max_retries = 3
    base_delay = 10  # seconds

    for attempt in range(1, max_retries + 1):
        try:
            log(f"  Submitting AlleleID {allele_id} to NCBI BLAST "
                f"(attempt {attempt}/{max_retries}) ...")

            result_handle = NCBIWWW.qblast(
                program="blastn",
                database="nt",
                sequence=sequence,
                entrez_query='"Homo sapiens"[organism]',
                hitlist_size=hitlist_size,
                expect=evalue,
                # megablast is appropriate for high-identity searches against
                # a large database.  For shorter/divergent queries, blastn
                # would be more sensitive, but megablast is the NCBI default
                # for qblast and works well for our 80-100% identity range.
                megablast=True,
            )

            xml_data = result_handle.read()
            result_handle.close()

            # Save raw XML to disk
            xml_dir.mkdir(parents=True, exist_ok=True)
            with open(xml_path, "w") as fh:
                fh.write(xml_data)

            log(f"  Saved XML ({len(xml_data):,} bytes) to: {xml_path.name}")
            return xml_path

        except Exception as exc:
            delay = base_delay * (3 ** (attempt - 1))  # Exponential backoff: 10, 30, 90s
            if attempt < max_retries:
                log(f"  WARNING -- NCBI BLAST failed for AlleleID {allele_id}: {exc}")
                log(f"  Retrying in {delay} seconds (attempt {attempt}/{max_retries}) ...")
                time.sleep(delay)
            else:
                log(f"  ERROR -- NCBI BLAST failed after {max_retries} attempts "
                    f"for AlleleID {allele_id}: {exc}")
                return None

    return None


def run_blast_remote(
    variant_docs: dict,
    xml_dir: Path,
    hitlist_size: int = 20,
    evalue: float = 1e-5,
) -> list[dict]:
    """
    Run remote BLAST for all variants and parse the XML results.

    Iterates over each variant, submitting to NCBI BLAST and sleeping 10 seconds
    between queries to comply with NCBI's rate-limit policy.

    Parameters
    ----------
    variant_docs : dict
        {allele_id: variant_doc} mapping from MongoDB.
    xml_dir : Path
        Directory for raw XML results.
    hitlist_size : int
        Max hits per query.
    evalue : float
        E-value threshold.

    Returns
    -------
    list[dict]
        List of parsed hit dictionaries, each containing allele_id, chromosome,
        hit_start, hit_end, pident, evalue, bitscore, etc.
    """
    xml_dir.mkdir(parents=True, exist_ok=True)

    all_hits = []
    total = len(variant_docs)

    for idx, (allele_id, doc) in enumerate(variant_docs.items(), start=1):
        seq = doc.get("extracted_sequence", "")
        if not seq:
            continue

        log(f"[{idx}/{total}] Processing AlleleID {allele_id} ...")

        xml_path = blast_remote_single(
            allele_id=allele_id,
            sequence=seq,
            xml_dir=xml_dir,
            hitlist_size=hitlist_size,
            evalue=evalue,
        )

        if xml_path is not None:
            hits = parse_xml_result(xml_path, allele_id)
            all_hits.extend(hits)

        # ── NCBI rate limit: sleep 10 seconds between queries ────────────
        # NCBI policy: "Do not submit more than one request every 10 seconds."
        # https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html
        # We only sleep if we actually submitted a query (not cached).
        if idx < total:
            cached = (xml_path is not None
                      and xml_path.exists()
                      and xml_path.stat().st_size > 100)
            # If we used cache, the file existed before our call.  We check
            # by looking at whether the allele_id was logged as CACHED.
            # Simpler: always sleep if this was NOT a cache hit.
            # Since blast_remote_single already returns fast for cached items,
            # we time the function and skip sleep if it returned in < 1 second.
            # For simplicity, always sleep between non-cached submissions.
            log(f"  Sleeping 10s (NCBI rate limit) ...")
            time.sleep(10)

    log(f"Remote BLAST complete: {len(all_hits):,} total hits across {total} queries")
    return all_hits


def parse_xml_result(xml_path: Path, allele_id: int) -> list[dict]:
    """
    Parse a single BLAST XML result file into a list of hit dictionaries.

    Parameters
    ----------
    xml_path : Path
        Path to the BLAST XML output.
    allele_id : int
        The AlleleID this result belongs to.

    Returns
    -------
    list[dict]
        Parsed hits with chromosome, coordinates, scores, etc.
    """
    hits = []

    try:
        with open(xml_path) as fh:
            blast_records = NCBIXML.parse(fh)
            for record in blast_records:
                for alignment in record.alignments:
                    subject_id = alignment.hit_id
                    subject_def = alignment.hit_def

                    chromosome = normalise_chromosome(subject_id)
                    # Also try the definition line for chromosome info
                    if not chromosome.startswith("chr"):
                        chromosome = normalise_chromosome(subject_def)

                    for hsp in alignment.hsps:
                        hits.append({
                            "allele_id": allele_id,
                            "subject_id": subject_id,
                            "subject_def": subject_def[:200],  # Truncate long descriptions
                            "chromosome": chromosome,
                            "hit_start": hsp.sbjct_start,
                            "hit_end": hsp.sbjct_end,
                            "pident": round(
                                hsp.identities / hsp.align_length * 100, 2
                            ) if hsp.align_length > 0 else 0.0,
                            "length": hsp.align_length,
                            "evalue": hsp.expect,
                            "bitscore": round(hsp.bits, 2),
                            "qstart": hsp.query_start,
                            "qend": hsp.query_end,
                        })
    except Exception as exc:
        log(f"  WARNING -- Failed to parse XML for AlleleID {allele_id}: {exc}")

    return hits


# ─── Step 2b: Local BLAST against genome database ───────────────────────────

def run_blast_local(
    query_fasta: Path,
    genome_db: Path,
    output_path: Path,
    evalue: float = 1e-5,
    perc_identity: float = 90.0,
    num_threads: int = 4,
) -> list[dict]:
    """
    Run local blastn against a user-provided genome BLAST database.

    This mode requires the user to have built a local hg38 BLAST database.
    It is MUCH faster than remote BLAST (~1 second per query vs ~60 seconds)
    but requires ~15 GB of disk space for the database.

    Parameters
    ----------
    query_fasta : Path
        FASTA file of query sequences.
    genome_db : Path
        Path to the local BLAST database (without file extensions).
    output_path : Path
        Where to write the tabular output.
    evalue : float
        E-value threshold.
    perc_identity : float
        Minimum percent identity.
    num_threads : int
        Number of CPU threads for BLAST.

    Returns
    -------
    list[dict]
        List of parsed hit dictionaries.
    """
    # ── Verify blastn is available ───────────────────────────────────────
    if shutil.which("blastn") is None:
        log("FATAL -- 'blastn' not found in PATH.")
        log("BLAST+ must be installed.  See: brew install blast / conda install blast")
        sys.exit(1)

    # ── Build the output format string ───────────────────────────────────
    outfmt_fields = (
        "qseqid sseqid pident length mismatch gapopen "
        "qstart qend sstart send evalue bitscore qlen slen"
    )
    outfmt = f"6 {outfmt_fields}"

    output_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "blastn",
        "-query", str(query_fasta),
        "-db", str(genome_db),
        "-out", str(output_path),
        "-evalue", str(evalue),
        "-perc_identity", str(perc_identity),
        "-outfmt", outfmt,
        "-num_threads", str(num_threads),
    ]

    log(f"Running local BLAST against genome database: {genome_db}")
    log(f"  E-value threshold:  {evalue}")
    log(f"  Min % identity:     {perc_identity}")
    log(f"  Threads:            {num_threads}")
    log(f"  Output:             {output_path}")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        if result.stderr.strip():
            for line in result.stderr.strip().split("\n"):
                log(f"  blastn stderr: {line}")
    except subprocess.CalledProcessError as exc:
        log(f"FATAL -- blastn failed with exit code {exc.returncode}")
        log(f"  stderr: {exc.stderr}")
        sys.exit(1)

    # ── Parse tabular output ─────────────────────────────────────────────
    if not output_path.exists() or output_path.stat().st_size == 0:
        log("WARNING -- Local BLAST returned no hits.")
        return []

    df = pd.read_csv(
        output_path,
        sep="\t",
        header=None,
        names=GENOME_BLAST_COLUMNS,
        dtype={"qseqid": str, "sseqid": str},
    )

    log(f"Local BLAST completed: {len(df):,} HSPs found")

    hits = []
    for _, row in df.iterrows():
        qseqid = str(row["qseqid"])
        parts = qseqid.split("|")
        try:
            allele_id = int(parts[0])
        except (ValueError, IndexError):
            continue

        chromosome = normalise_chromosome(str(row["sseqid"]))

        hits.append({
            "allele_id": allele_id,
            "subject_id": str(row["sseqid"]),
            "subject_def": "",
            "chromosome": chromosome,
            "hit_start": int(row["sstart"]),
            "hit_end": int(row["send"]),
            "pident": round(float(row["pident"]), 2),
            "length": int(row["length"]),
            "evalue": float(row["evalue"]),
            "bitscore": round(float(row["bitscore"]), 2),
            "qstart": int(row["qstart"]),
            "qend": int(row["qend"]),
        })

    return hits


# ─── Step 3: Interpret and classify results ──────────────────────────────────

def interpret_hits(
    all_hits: list[dict],
    variant_docs: dict,
) -> dict[int, dict]:
    """
    Classify each hit and each variant based on genome BLAST results.

    For each hit:
      - is_mitochondrial_chromosome: True if the hit is on chrM (expected —
        confirms the sequence is mitochondrial DNA).
      - is_known_numt: True if the hit is on a nuclear chromosome (chr1–22, X, Y).
        A nuclear hit means this mtDNA fragment already exists in the hg38 reference
        genome, i.e., it is a "known" or "reference" NUMT.
      - overlaps_reported_clinvar_position: True if the nuclear hit is within 10 kb
        of the ClinVar-reported Start/Stop position.  This confirms the NUMT is at
        the expected disease locus.

    For each variant (after examining all its hits):
      - CONFIRMED_NUMT: mtDNA hit (Step 4) AND a nuclear hit overlapping the
        ClinVar position.  Strong evidence that the ClinVar variant is a genuine
        NUMT insertion at the reported gene.
      - KNOWN_REFERENCE_NUMT: nuclear hit(s) found, but NOT at the ClinVar position.
        The mtDNA fragment exists in the reference genome elsewhere — a benign
        ancestral NUMT.  The ClinVar variant might still be novel at its locus.
      - NOVEL_NUMT_CANDIDATE: mtDNA hit (Step 4) but NO nuclear hits at all.
        The sequence is mitochondrial but is NOT in the reference genome.  This is
        the most interesting category: a truly novel NUMT insertion.
      - AMBIGUOUS: nuclear hit found but the Step 4 mtDNA alignment was poor or
        absent.  Needs manual review.

    Parameters
    ----------
    all_hits : list[dict]
        Parsed hits from remote or local BLAST.
    variant_docs : dict
        {allele_id: MongoDB document} for all queried variants.

    Returns
    -------
    dict[int, dict]
        {allele_id: classification_result} where classification_result contains
        the annotated hits, classification string, and summary flags.
    """
    # ── Group hits by AlleleID ───────────────────────────────────────────
    hits_by_allele = {}
    for hit in all_hits:
        aid = hit["allele_id"]
        if aid not in hits_by_allele:
            hits_by_allele[aid] = []
        hits_by_allele[aid].append(hit)

    results = {}
    overlap_window = 10_000  # 10 kb window for position overlap check

    for allele_id, doc in variant_docs.items():
        variant_hits = hits_by_allele.get(allele_id, [])
        clinvar_chrom = get_clinvar_chromosome(doc)
        clinvar_start = doc.get("start") or doc.get("Start") or doc.get("position_start")
        clinvar_stop = doc.get("stop") or doc.get("Stop") or doc.get("position_stop")

        # Convert ClinVar positions to integers if available
        try:
            clinvar_start = int(clinvar_start) if clinvar_start is not None else None
        except (ValueError, TypeError):
            clinvar_start = None
        try:
            clinvar_stop = int(clinvar_stop) if clinvar_stop is not None else None
        except (ValueError, TypeError):
            clinvar_stop = None

        # Step 4 mitochondrial hit information
        has_mito_hit = doc.get("has_mito_hit", False)
        best_mito_pident = doc.get("best_mito_hit_pident")

        # ── Annotate each hit ────────────────────────────────────────────
        annotated_hits = []
        has_nuclear_hit = False
        has_mito_chromosome_hit = False
        overlaps_clinvar = False

        for hit in variant_hits:
            chrom = hit["chromosome"]
            hit_start = min(hit["hit_start"], hit["hit_end"])
            hit_end = max(hit["hit_start"], hit["hit_end"])

            # Is this a mitochondrial chromosome hit?
            is_mito = chrom in ("chrM", "chrMT", "MT", "NC_012920.1")

            # Is this a nuclear chromosome hit?
            # Nuclear = any standard chromosome that is NOT mitochondrial
            is_nuclear = (
                chrom.startswith("chr")
                and chrom not in ("chrM", "chrMT")
            )

            # Does this nuclear hit overlap the ClinVar-reported position?
            overlaps_reported = False
            if is_nuclear and clinvar_chrom and clinvar_start is not None:
                if chrom == clinvar_chrom:
                    # Check if the BLAST hit is within the overlap window
                    # of the ClinVar position
                    cv_start = clinvar_start - overlap_window
                    cv_end = (clinvar_stop or clinvar_start) + overlap_window
                    if hit_start <= cv_end and hit_end >= cv_start:
                        overlaps_reported = True

            hit["is_mitochondrial_chromosome"] = is_mito
            hit["is_known_numt"] = is_nuclear
            hit["overlaps_reported_clinvar_position"] = overlaps_reported

            if is_mito:
                has_mito_chromosome_hit = True
            if is_nuclear:
                has_nuclear_hit = True
            if overlaps_reported:
                overlaps_clinvar = True

            annotated_hits.append(hit)

        # ── Classify the variant ─────────────────────────────────────────
        if has_mito_hit and overlaps_clinvar:
            # Best case: confirmed mitochondrial origin AND the NUMT maps to
            # the ClinVar-reported gene position.
            classification = "CONFIRMED_NUMT"
        elif has_nuclear_hit and not overlaps_clinvar:
            # The sequence exists in the reference genome, but at a different
            # location.  This is a known reference NUMT — probably ancestral
            # and benign.
            classification = "KNOWN_REFERENCE_NUMT"
        elif has_mito_hit and not has_nuclear_hit:
            # Mitochondrial origin confirmed (Step 4), but the sequence is NOT
            # found in the nuclear reference genome.  This is a novel insertion
            # candidate — the most scientifically interesting category.
            classification = "NOVEL_NUMT_CANDIDATE"
        elif has_nuclear_hit and not has_mito_hit:
            # Nuclear hit without a good mitochondrial match.  Could be a very
            # divergent NUMT or a false positive.  Needs manual review.
            classification = "AMBIGUOUS"
        else:
            # No hits at all — the sequence didn't match the genome.
            # This shouldn't happen often since we pre-filtered for mito hits.
            classification = "NOVEL_NUMT_CANDIDATE"

        results[allele_id] = {
            "allele_id": allele_id,
            "gene_symbol": doc.get("gene_symbol", "NA"),
            "clinical_significance": doc.get("clinical_significance", "NA"),
            "condition": doc.get("condition_name") or doc.get("condition") or "NA",
            "blast_genome_results": annotated_hits,
            "numt_classification": classification,
            "is_confirmed_numt": classification == "CONFIRMED_NUMT",
            "overlaps_reported_position": overlaps_clinvar,
            "has_nuclear_hit": has_nuclear_hit,
            "has_mito_chromosome_hit": has_mito_chromosome_hit,
            "best_mito_pident": best_mito_pident,
            "nuclear_chromosomes": list(set(
                h["chromosome"] for h in annotated_hits if h.get("is_known_numt")
            )),
        }

    return results


# ─── Step 4: Update MongoDB ─────────────────────────────────────────────────

def update_mongodb(
    classification_results: dict[int, dict],
    host: str,
    port: int,
) -> None:
    """
    Write genome BLAST results and NUMT classifications back to MongoDB.

    For each variant:
      $set: {
          blast_genome_results:       [array of annotated hit dicts],
          numt_classification:        "CONFIRMED_NUMT" | "KNOWN_REFERENCE_NUMT" | etc.,
          is_confirmed_numt:          True/False,
          overlaps_reported_position: True/False,
          blast_genome_updated_at:    <timestamp>,
      }

    Parameters
    ----------
    classification_results : dict
        {allele_id: classification_result} from interpret_hits().
    host : str
        MongoDB hostname.
    port : int
        MongoDB port.
    """
    if not PYMONGO_AVAILABLE:
        log("FATAL -- pymongo not available for MongoDB update.")
        sys.exit(1)

    log("Updating MongoDB with genome BLAST results and classifications ...")

    try:
        client = pymongo.MongoClient(
            host=host, port=port, serverSelectionTimeoutMS=5000
        )
        client.admin.command("ping")
    except (ConnectionFailure, ServerSelectionTimeoutError) as exc:
        log(f"FATAL -- Cannot connect to MongoDB: {exc}")
        sys.exit(1)

    db = client[MONGO_DB]
    collection = db[MONGO_COLLECTION]

    operations = []

    for allele_id, result in classification_results.items():
        # Clean hit dicts for MongoDB storage (remove very long fields)
        clean_hits = []
        for hit in result["blast_genome_results"]:
            clean_hits.append({
                "chromosome": hit.get("chromosome"),
                "hit_start": hit.get("hit_start"),
                "hit_end": hit.get("hit_end"),
                "pident": hit.get("pident"),
                "length": hit.get("length"),
                "evalue": hit.get("evalue"),
                "bitscore": hit.get("bitscore"),
                "is_mitochondrial_chromosome": hit.get("is_mitochondrial_chromosome"),
                "is_known_numt": hit.get("is_known_numt"),
                "overlaps_reported_clinvar_position": hit.get(
                    "overlaps_reported_clinvar_position"
                ),
            })

        operations.append(UpdateOne(
            filter={"allele_id": int(allele_id)},
            update={"$set": {
                "blast_genome_results": clean_hits,
                "numt_classification": result["numt_classification"],
                "is_confirmed_numt": result["is_confirmed_numt"],
                "overlaps_reported_position": result["overlaps_reported_position"],
                "blast_genome_updated_at": datetime.now(timezone.utc),
            }},
        ))

    if operations:
        result = collection.bulk_write(operations, ordered=False)
        log(f"MongoDB update complete: {result.modified_count} documents modified, "
            f"{result.matched_count} matched")
    else:
        log("No MongoDB update operations to perform.")

    client.close()


# ─── Step 5: Summary output ─────────────────────────────────────────────────

def print_summary(
    classification_results: dict[int, dict],
    summary_tsv: Path,
) -> None:
    """
    Print classification summary and write the main scientific output table.

    The summary table lists CONFIRMED_NUMT and NOVEL_NUMT_CANDIDATE variants
    with their key attributes.  This is the table that answers the core research
    question: "Which ClinVar insertions are unrecognised NUMTs?"

    Parameters
    ----------
    classification_results : dict
        {allele_id: classification_result} from interpret_hits().
    summary_tsv : Path
        Where to save the summary table as TSV.
    """
    print("\n" + "=" * 90)
    print("  BLAST vs Human Genome (hg38) — NUMT Classification Summary")
    print("=" * 90)

    if not classification_results:
        print("\n  No variants to classify.")
        print("=" * 90 + "\n")
        return

    # ── Count per classification ─────────────────────────────────────────
    class_counts = {}
    for result in classification_results.values():
        c = result["numt_classification"]
        class_counts[c] = class_counts.get(c, 0) + 1

    total = len(classification_results)
    print(f"\n  Total variants analysed:       {total:>6}")
    print(f"\n  {'Classification breakdown':─<60}")

    classification_order = [
        "CONFIRMED_NUMT",
        "NOVEL_NUMT_CANDIDATE",
        "KNOWN_REFERENCE_NUMT",
        "AMBIGUOUS",
    ]
    for cls in classification_order:
        count = class_counts.get(cls, 0)
        pct = count / total * 100 if total > 0 else 0
        marker = " <-- DISEASE-RELEVANT" if cls in (
            "CONFIRMED_NUMT", "NOVEL_NUMT_CANDIDATE"
        ) else ""
        print(f"    {cls:<30} {count:>5}  ({pct:>5.1f}%){marker}")

    # ── Detailed table: CONFIRMED_NUMT and NOVEL_NUMT_CANDIDATE ─────────
    priority_results = [
        r for r in classification_results.values()
        if r["numt_classification"] in ("CONFIRMED_NUMT", "NOVEL_NUMT_CANDIDATE")
    ]

    if priority_results:
        print(f"\n  {'PRIORITY CANDIDATES — Potential disease-causing NUMTs':─<60}")
        print(f"    {'AlleleID':<10} {'Gene':<10} {'Classification':<28} "
              f"{'MitoPident':>10} {'NuclearChr':<12} {'ClinSig':<30} {'Condition':<40}")
        print(f"    {'─' * 10} {'─' * 10} {'─' * 28} "
              f"{'─' * 10} {'─' * 12} {'─' * 30} {'─' * 40}")

        for r in sorted(priority_results, key=lambda x: x["allele_id"]):
            allele_id = r["allele_id"]
            gene = r.get("gene_symbol", "NA") or "NA"
            classification = r["numt_classification"]
            best_pident = r.get("best_mito_pident")
            pident_str = f"{best_pident:.1f}" if best_pident is not None else "NA"
            nuclear_chrs = ", ".join(r.get("nuclear_chromosomes", [])) or "none"
            clinsig = (r.get("clinical_significance", "NA") or "NA")[:30]
            condition = (r.get("condition", "NA") or "NA")[:40]

            print(f"    {allele_id:<10} {gene:<10} {classification:<28} "
                  f"{pident_str:>10} {nuclear_chrs:<12} {clinsig:<30} {condition:<40}")
    else:
        print(f"\n  No CONFIRMED_NUMT or NOVEL_NUMT_CANDIDATE variants found.")

    # ── Also show KNOWN_REFERENCE_NUMT hits (informational) ──────────────
    reference_numts = [
        r for r in classification_results.values()
        if r["numt_classification"] == "KNOWN_REFERENCE_NUMT"
    ]
    if reference_numts:
        print(f"\n  {'KNOWN REFERENCE NUMTs (in hg38, likely benign)':─<60}")
        for r in sorted(reference_numts, key=lambda x: x["allele_id"])[:10]:
            allele_id = r["allele_id"]
            gene = r.get("gene_symbol", "NA") or "NA"
            nuclear_chrs = ", ".join(r.get("nuclear_chromosomes", [])) or "none"
            print(f"    AlleleID {allele_id:<8} ({gene}) — "
                  f"reference NUMTs on: {nuclear_chrs}")
        if len(reference_numts) > 10:
            print(f"    ... and {len(reference_numts) - 10} more")

    print("\n" + "=" * 90 + "\n")

    # ─── Write summary TSV ───────────────────────────────────────────────
    # This is the main scientific output file.
    summary_tsv.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for r in classification_results.values():
        best_pident = r.get("best_mito_pident")
        pident_str = f"{best_pident:.1f}" if best_pident is not None else "NA"

        rows.append({
            "AlleleID": r["allele_id"],
            "GeneSymbol": r.get("gene_symbol", "NA"),
            "ClinicalSignificance": r.get("clinical_significance", "NA"),
            "Condition": r.get("condition", "NA"),
            "NUMTClassification": r["numt_classification"],
            "IsConfirmedNUMT": r["is_confirmed_numt"],
            "BestMitoPident": pident_str,
            "NuclearChromosomes": ", ".join(r.get("nuclear_chromosomes", [])) or "none",
            "OverlapsClinVarPosition": r["overlaps_reported_position"],
            "HasNuclearHit": r["has_nuclear_hit"],
            "HasMitoChromosomeHit": r["has_mito_chromosome_hit"],
            "NumGenomeHits": len(r["blast_genome_results"]),
        })

    df_summary = pd.DataFrame(rows)
    df_summary.to_csv(summary_tsv, sep="\t", index=False)
    log(f"Summary table saved to: {summary_tsv}")
    log(f"  {len(rows)} variants classified")


# ─── CLI entry point ────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    """Define and parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Step 5 of the NUMT pipeline: BLAST mito-positive insertion sequences "
            "against the full human genome (hg38/GRCh38) to confirm insertion "
            "locations and identify known reference NUMTs."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  # Remote BLAST (default, no local database needed):\n"
            "  python 05_blast_genome.py --mode remote\n\n"
            "  # Local BLAST (faster, requires local hg38 database):\n"
            "  python 05_blast_genome.py --mode local --genome-db /path/to/hg38_db\n\n"
            "  # Skip MongoDB (use existing FASTA, don't update DB):\n"
            "  python 05_blast_genome.py --mode remote --skip-mongodb\n"
        ),
    )

    # ── BLAST mode ──────────────────────────────────────────────────────
    parser.add_argument(
        "--mode",
        type=str,
        choices=["remote", "local"],
        default="remote",
        help=(
            "BLAST mode: 'remote' uses NCBI qblast (slow, no download); "
            "'local' uses a local hg38 BLAST database (fast, requires --genome-db). "
            "Default: remote."
        ),
    )
    parser.add_argument(
        "--genome-db",
        type=Path,
        default=None,
        help=(
            "Path to local hg38 BLAST database (without file extensions). "
            "Required for --mode local."
        ),
    )

    # ── Input/output paths ──────────────────────────────────────────────
    parser.add_argument(
        "--query-fasta",
        type=Path,
        default=QUERY_FASTA,
        help=f"Path for the query FASTA (default: {QUERY_FASTA.relative_to(PROJECT_ROOT)})",
    )
    parser.add_argument(
        "--blast-output",
        type=Path,
        default=LOCAL_BLAST_OUTPUT,
        help=f"Path for local BLAST TSV (default: {LOCAL_BLAST_OUTPUT.relative_to(PROJECT_ROOT)})",
    )
    parser.add_argument(
        "--xml-dir",
        type=Path,
        default=RAW_XML_DIR,
        help=f"Directory for raw XML results (default: {RAW_XML_DIR.relative_to(PROJECT_ROOT)})",
    )
    parser.add_argument(
        "--summary-tsv",
        type=Path,
        default=SUMMARY_TSV,
        help=f"Path for summary output TSV (default: {SUMMARY_TSV.relative_to(PROJECT_ROOT)})",
    )

    # ── MongoDB connection ──────────────────────────────────────────────
    parser.add_argument(
        "--host", type=str, default=MONGO_HOST,
        help=f"MongoDB hostname (default: {MONGO_HOST})",
    )
    parser.add_argument(
        "--port", type=int, default=MONGO_PORT,
        help=f"MongoDB port (default: {MONGO_PORT})",
    )

    # ── BLAST parameters ────────────────────────────────────────────────
    parser.add_argument(
        "--evalue", type=float, default=1e-5,
        help="E-value threshold (default: 1e-5, more stringent than mito step)",
    )
    parser.add_argument(
        "--perc-identity", type=float, default=90.0,
        help="Min percent identity for local mode (default: 90)",
    )
    parser.add_argument(
        "--hitlist-size", type=int, default=20,
        help="Max hits per query for remote mode (default: 20)",
    )
    parser.add_argument(
        "--num-threads", type=int, default=4,
        help="Number of BLAST threads for local mode (default: 4)",
    )

    # ── Behaviour flags ─────────────────────────────────────────────────
    parser.add_argument(
        "--skip-mongodb", action="store_true",
        help="Skip MongoDB query/update.  Use existing FASTA, print results only.",
    )

    return parser.parse_args()


# ─── Main orchestration ─────────────────────────────────────────────────────

def main() -> None:
    """Orchestrate the genome BLAST pipeline."""
    args = parse_args()

    log("=" * 70)
    log("NUMT Pipeline -- Step 5: BLAST vs Human Genome (hg38/GRCh38)")
    log("=" * 70)
    log(f"Mode: {args.mode.upper()}")

    # ── Validate arguments ───────────────────────────────────────────────
    if args.mode == "local" and args.genome_db is None:
        log("FATAL -- --genome-db is required for --mode local.")
        log("Provide the path to your local hg38 BLAST database.")
        log("To build one:")
        log("  wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/"
            "GCF_000001405.40_GRCh38.p14_genomic.fna.gz")
        log("  gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz")
        log("  makeblastdb -in GCF_000001405.40_GRCh38.p14_genomic.fna "
            "-dbtype nucl -out hg38_db")
        sys.exit(1)

    if args.mode == "remote" and not BIOPYTHON_BLAST_AVAILABLE:
        log("FATAL -- BioPython BLAST modules not available for remote mode.")
        log("Install with: pip install biopython")
        sys.exit(1)

    # ── Step 1: Prepare query FASTA from MongoDB ─────────────────────────
    if args.skip_mongodb:
        log("Skipping MongoDB query (--skip-mongodb).  Using existing FASTA.")
        if not args.query_fasta.exists():
            log(f"FATAL -- Query FASTA not found: {args.query_fasta}")
            sys.exit(1)

        # Read variant info from FASTA headers for classification
        variant_docs = {}
        with open(args.query_fasta) as fh:
            for line in fh:
                if line.startswith(">"):
                    parts = line.strip().lstrip(">").split("|")
                    try:
                        aid = int(parts[0])
                        variant_docs[aid] = {
                            "allele_id": aid,
                            "gene_symbol": parts[1] if len(parts) > 1 else "NA",
                            "clinical_significance": parts[2] if len(parts) > 2 else "NA",
                            "has_mito_hit": True,
                            "best_mito_hit_pident": (
                                float(parts[3]) if len(parts) > 3
                                and parts[3] != "NA" else None
                            ),
                        }
                    except (ValueError, IndexError):
                        pass

        n_queries = len(variant_docs)
    else:
        n_queries, variant_docs = prepare_query_fasta(
            host=args.host,
            port=args.port,
            output_fasta=args.query_fasta,
        )

    if n_queries == 0:
        log("FATAL -- No query sequences available.  Nothing to BLAST.")
        log("Ensure Step 4 (04_blast_mito.py) has been run and has_mito_hit=True "
            "variants exist in MongoDB.")
        sys.exit(1)

    log(f"Query FASTA contains {n_queries:,} sequences")

    # ── Step 2: Run BLAST (remote or local) ──────────────────────────────
    if args.mode == "remote":
        log("")
        log("─── REMOTE BLAST MODE ─────────────────────────────────────")
        log("Using NCBI qblast (database: nt, organism: Homo sapiens)")
        log(f"Rate limit: 10 seconds between queries")
        log(f"Raw XML cache: {args.xml_dir}")
        log("")

        all_hits = run_blast_remote(
            variant_docs=variant_docs,
            xml_dir=args.xml_dir,
            hitlist_size=args.hitlist_size,
            evalue=args.evalue,
        )
    else:
        log("")
        log("─── LOCAL BLAST MODE ──────────────────────────────────────")
        log(f"Database: {args.genome_db}")
        log("")

        all_hits = run_blast_local(
            query_fasta=args.query_fasta,
            genome_db=args.genome_db,
            output_path=args.blast_output,
            evalue=args.evalue,
            perc_identity=args.perc_identity,
            num_threads=args.num_threads,
        )

    log(f"Total hits collected: {len(all_hits):,}")

    # ── Step 3: Interpret and classify ───────────────────────────────────
    log("")
    log("─── CLASSIFICATION ────────────────────────────────────────")
    classification_results = interpret_hits(all_hits, variant_docs)

    # ── Step 4: Update MongoDB ───────────────────────────────────────────
    if not args.skip_mongodb:
        update_mongodb(classification_results, host=args.host, port=args.port)
    else:
        log("Skipping MongoDB update (--skip-mongodb)")

    # ── Step 5: Summary ──────────────────────────────────────────────────
    print_summary(classification_results, summary_tsv=args.summary_tsv)

    log("=" * 70)
    log("Step 5 complete.")
    log("=" * 70)


# ─── Entry point ─────────────────────────────────────────────────────────────

if __name__ == "__main__":
    main()
