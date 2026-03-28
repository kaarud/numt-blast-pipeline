#!/usr/bin/env python3
"""
04_blast_mito.py — BLAST extracted insertion sequences against the human mitochondrial genome.

BIOLOGICAL RATIONALE
====================
We are searching for unrecognised NUMTs (Nuclear Mitochondrial DNA Segments) hiding among
ClinVar pathogenic and VUS insertion variants.  The confirmed index case is a ~72 bp
fragment from the mitochondrial genome inserted into TSC2 exon 17, causing Tuberous
Sclerosis Complex.  How many similar cases remain undetected?

In the previous steps we:
  Step 1 — Downloaded ClinVar insertions > 50 bp.
  Step 2 — Extracted nucleotide sequences from HGVS notation.
  Step 3 — Stored everything in MongoDB (numt_db.clinvar_variants).

This script now BLASTs those extracted sequences against the human mitochondrial genome
to identify which insertions are of mitochondrial origin.

THE HUMAN MITOCHONDRIAL GENOME
==============================
The human mitochondrial genome is a circular, double-stranded DNA molecule of 16,569 bp
(revised Cambridge Reference Sequence, rCRS, NCBI accession NC_012920.1).  It encodes:

  - 13 protein-coding genes (components of the oxidative phosphorylation chain):
      ND1, ND2, COX1, COX2, ATP8, ATP6, COX3, ND3, ND4L, ND4, ND5, ND6, CYTB
  - 22 transfer RNAs (tRNAs)
  - 2 ribosomal RNAs (12S rRNA, 16S rRNA)
  - The D-loop (displacement loop) control region: a non-coding regulatory region
    containing the heavy-strand origin of replication and both promoters.

NUMTs can originate from ANY region of the mtDNA.  Once inserted into the nuclear genome,
they become "molecular fossils" — they accumulate mutations at the nuclear rate (slower
than mitochondrial), so older NUMTs may have diverged significantly from the current rCRS.
This is why we use permissive BLAST parameters: we expect 80-100% identity depending on
the age of the insertion event.

BLAST PARAMETER RATIONALE
=========================
For short, potentially divergent sequences (~50-500 bp) we tune parameters for sensitivity:

  -task blastn-short : Optimised for queries < 50 bp (uses smaller word size internally).
                       For queries >= 100 bp we use standard blastn.

  -evalue 1e-3      : The E-value (Expect value) is the number of alignments with an
                       equal or better score expected by chance in a database of this size.
                       1e-3 is permissive — a smaller database (16.5 kb mtDNA) means even
                       moderate scores are significant, but we keep this relaxed to catch
                       divergent NUMTs in the initial screen.

  -word_size 7      : The minimum length of an exact match ("seed") needed to initiate
                       an alignment extension.  BLAST first finds exact word matches, then
                       extends them.  Smaller word size = more sensitive but slower.
                       Default blastn word_size is 11; we use 7 for short-query sensitivity.

  -perc_identity 80 : Minimum percent identity threshold.  NUMTs that inserted millions
                       of years ago may have accumulated ~20% divergence from the current
                       rCRS.  80% catches these ancient insertions while filtering noise.

  -qcov_hsp_perc 50 : At least 50% of the query must participate in the alignment.
                       This filters out tiny partial matches that hit by chance.  For a
                       100 bp query, at least 50 bp must align.

  -outfmt 6         : Tabular output (no headers, machine-parseable).  We request 14
                       columns including query/subject coordinates and alignment stats.

  -num_threads 4    : Parallelise the search across 4 CPU cores.
"""

# ─── Standard-library imports ────────────────────────────────────────────────
import argparse
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

# ─── Third-party imports ─────────────────────────────────────────────────────
import pandas as pd

try:
    from Bio import Entrez, SeqIO
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

# Input: MongoDB collection (primary) and FASTA backup
MONGO_HOST = "localhost"
MONGO_PORT = 27017
MONGO_DB = "numt_db"
MONGO_COLLECTION = "clinvar_variants"

# Mitochondrial reference genome
MITO_FASTA = PROJECT_ROOT / "02_data" / "raw" / "NC_012920.1.fasta"
BLAST_DB_DIR = PROJECT_ROOT / "02_data" / "raw" / "blast_db" / "mito"
BLAST_DB_NAME = "mito_db"

# Outputs
QUERY_FASTA = PROJECT_ROOT / "02_data" / "processed" / "blast_query_sequences.fasta"
BLAST_OUTPUT = PROJECT_ROOT / "02_data" / "processed" / "blast_mito_results.tsv"

# BLAST tabular output column names (matching -outfmt 6 fields we request)
BLAST_COLUMNS = [
    "qseqid",       # Query sequence ID (AlleleID)
    "sseqid",       # Subject sequence ID (NC_012920.1)
    "pident",       # Percentage of identical matches
    "length",       # Alignment length (bp)
    "mismatch",     # Number of mismatches
    "gapopen",      # Number of gap openings
    "qstart",       # Start of alignment in query
    "qend",         # End of alignment in query
    "sstart",       # Start of alignment in subject (mtDNA position)
    "send",         # End of alignment in subject
    "evalue",       # Expect value — probability of hit occurring by chance
    "bitscore",     # Bit score — normalised alignment score
    "qlen",         # Query sequence length
    "slen",         # Subject (mtDNA) sequence length
]

# ─── Mitochondrial genome annotation ────────────────────────────────────────
# Coordinates from the rCRS (NC_012920.1) GenBank annotation.
# Each entry: (name, start, end, type)
# Note: mtDNA coordinates are 1-based; the genome is circular so the D-loop
# spans the origin (positions 16024-16569 and 1-576).

MTDNA_FEATURES = [
    # D-loop / control region (non-coding, contains replication origin + promoters)
    ("D-loop",         16024, 16569, "control_region"),
    ("D-loop",             1,   576, "control_region"),
    # Ribosomal RNAs
    ("12S rRNA",         648,  1601, "rRNA"),
    ("16S rRNA",        1671,  3229, "rRNA"),
    # Transfer RNAs (22 tRNAs scattered throughout the genome)
    ("tRNA-Phe",         577,   647, "tRNA"),
    ("tRNA-Val",        1602,  1670, "tRNA"),
    ("tRNA-Leu(UUR)",   3230,  3304, "tRNA"),
    ("tRNA-Ile",        4263,  4331, "tRNA"),
    ("tRNA-Gln",        4329,  4400, "tRNA"),
    ("tRNA-Met",        4402,  4469, "tRNA"),
    ("tRNA-Trp",        5512,  5579, "tRNA"),
    ("tRNA-Ala",        5587,  5655, "tRNA"),
    ("tRNA-Asn",        5657,  5729, "tRNA"),
    ("tRNA-Cys",        5761,  5826, "tRNA"),
    ("tRNA-Tyr",        5826,  5891, "tRNA"),
    ("tRNA-Ser(UCN)",   7446,  7514, "tRNA"),
    ("tRNA-Asp",        7518,  7585, "tRNA"),
    ("tRNA-Lys",        8295,  8364, "tRNA"),
    ("tRNA-Gly",        9991, 10058, "tRNA"),
    ("tRNA-Arg",       10405, 10469, "tRNA"),
    ("tRNA-His",       12138, 12206, "tRNA"),
    ("tRNA-Ser(AGY)",  12207, 12265, "tRNA"),
    ("tRNA-Leu(CUN)",  12266, 12336, "tRNA"),
    ("tRNA-Glu",       14674, 14742, "tRNA"),
    ("tRNA-Thr",       15888, 15953, "tRNA"),
    ("tRNA-Pro",       15956, 16023, "tRNA"),
    # Protein-coding genes (13 genes of the oxidative phosphorylation chain)
    ("MT-ND1",          3307,  4262, "protein_coding"),   # NADH dehydrogenase subunit 1
    ("MT-ND2",          4470,  5511, "protein_coding"),   # NADH dehydrogenase subunit 2
    ("MT-CO1",          5904,  7445, "protein_coding"),   # Cytochrome c oxidase subunit I
    ("MT-CO2",          7586,  8269, "protein_coding"),   # Cytochrome c oxidase subunit II
    ("MT-ATP8",         8366,  8572, "protein_coding"),   # ATP synthase F0 subunit 8
    ("MT-ATP6",         8527,  9207, "protein_coding"),   # ATP synthase F0 subunit 6
    ("MT-CO3",          9207,  9990, "protein_coding"),   # Cytochrome c oxidase subunit III
    ("MT-ND3",         10059, 10404, "protein_coding"),   # NADH dehydrogenase subunit 3
    ("MT-ND4L",        10470, 10766, "protein_coding"),   # NADH dehydrogenase subunit 4L
    ("MT-ND4",         10760, 12137, "protein_coding"),   # NADH dehydrogenase subunit 4
    ("MT-ND5",         12337, 14148, "protein_coding"),   # NADH dehydrogenase subunit 5
    ("MT-ND6",         14149, 14673, "protein_coding"),   # NADH dehydrogenase subunit 6 (L-strand)
    ("MT-CYB",         14747, 15887, "protein_coding"),   # Cytochrome b
]

# Genes associated with known disease — we flag hits in variants affecting these
DISEASE_GENES = {
    "TSC1", "TSC2",       # Tuberous Sclerosis Complex
    "NF1", "NF2",         # Neurofibromatosis
    "BRCA1", "BRCA2",     # Hereditary breast/ovarian cancer
    "APC",                # Familial adenomatous polyposis
    "TP53",               # Li-Fraumeni syndrome
    "RB1",                # Retinoblastoma
    "VHL",                # Von Hippel-Lindau
    "MLH1", "MSH2", "MSH6", "PMS2",  # Lynch syndrome
    "CFTR",               # Cystic fibrosis
    "DMD",                # Duchenne muscular dystrophy
    "FBN1",               # Marfan syndrome
    "PKD1", "PKD2",       # Polycystic kidney disease
    "ATM",                # Ataxia-telangiectasia
    "PTEN",               # PTEN hamartoma tumour syndrome
}


# ─── Helper: timestamped logging ────────────────────────────────────────────

def log(message: str) -> None:
    """Print a message prefixed with the current timestamp."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")


# ─── Step 1: Prepare query FASTA from MongoDB ───────────────────────────────

def prepare_query_fasta(
    host: str,
    port: int,
    output_fasta: Path,
) -> int:
    """
    Read variants from MongoDB and write a FASTA file for BLAST querying.

    We select two categories of variants:
      1. DIRECT_SEQUENCE — variants where we extracted the nucleotide sequence
         directly from the HGVS name.  These are the primary BLAST queries.
      2. DIRECT_MITO_REF — variants whose HGVS name explicitly references the
         mitochondrial genome (e.g., ins(NC_012920.1:m.1-72)).  These are already
         known NUMTs, but we include them as positive controls to validate that
         our BLAST pipeline correctly identifies mitochondrial-origin sequences.

    FASTA header format:
        >AlleleID|GeneSymbol|ClinSig|LengthBp

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
    int
        Number of sequences written.
    """
    if not PYMONGO_AVAILABLE:
        log("FATAL -- pymongo is not installed.  pip install pymongo")
        sys.exit(1)

    log("Connecting to MongoDB to fetch query sequences ...")

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

    # ── Query 1: Variants with directly extracted sequences ──────────────
    # These have extraction_status == "DIRECT_SEQUENCE" and a non-null sequence.
    query_direct = {
        "extraction_status": "DIRECT_SEQUENCE",
        "extracted_sequence": {"$ne": None},
    }
    direct_docs = list(collection.find(query_direct))
    log(f"Found {len(direct_docs):,} DIRECT_SEQUENCE variants with sequences")

    # ── Query 2: DIRECT_MITO_REF entries (positive controls) ────────────
    # These variants reference the mitochondrial genome in their HGVS name.
    # They may or may not have an extracted_sequence field (the sequence needs
    # to be fetched from the mtDNA reference using the coordinate range).
    query_mito_ref = {
        "extraction_status": "DIRECT_MITO_REF",
        "extracted_sequence": {"$ne": None},
    }
    mito_ref_docs = list(collection.find(query_mito_ref))
    log(f"Found {len(mito_ref_docs):,} DIRECT_MITO_REF variants with sequences")

    # ── Write FASTA ─────────────────────────────────────────────────────
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    n_written = 0
    with open(output_fasta, "w") as fh:
        for doc in direct_docs + mito_ref_docs:
            seq = doc.get("extracted_sequence", "")
            if not seq or not isinstance(seq, str):
                continue

            allele_id = doc.get("allele_id", "unknown")
            gene = doc.get("gene_symbol", "NA") or "NA"
            clinsig = doc.get("clinical_significance", "NA") or "NA"
            # Sanitise clinical significance for FASTA header (remove spaces/pipes)
            clinsig_safe = clinsig.replace("|", "/").replace(" ", "_")
            length_bp = len(seq)

            # FASTA header: >AlleleID|GeneSymbol|ClinSig|LengthBp
            header = f">{allele_id}|{gene}|{clinsig_safe}|{length_bp}"
            fh.write(f"{header}\n")
            # Wrap sequence at 80 characters per line (FASTA convention)
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")

            n_written += 1

    client.close()

    if n_written == 0:
        log("WARNING -- No sequences written to FASTA.  The query file is empty.")
        log("Have you run steps 01-03 first?  Check that extraction_status "
            "== 'DIRECT_SEQUENCE' documents exist in MongoDB.")
    else:
        log(f"Wrote {n_written:,} sequences to: {output_fasta}")

    return n_written


# ─── Step 2: Download rCRS and build BLAST database ─────────────────────────

def download_rcrs(mito_fasta: Path) -> None:
    """
    Download the revised Cambridge Reference Sequence (rCRS) from NCBI.

    The rCRS (NC_012920.1) is the standard reference for the human mitochondrial
    genome.  It is 16,569 bp and was established by Andrews et al. (1999).

    We use Biopython's Entrez module to fetch the sequence in FASTA format.
    The download is skipped if the file already exists on disk.

    Parameters
    ----------
    mito_fasta : Path
        Where to save the downloaded FASTA file.
    """
    if mito_fasta.exists():
        log(f"rCRS already present at: {mito_fasta}")
        log("Skipping download (delete the file to force re-download)")
        return

    if not BIOPYTHON_AVAILABLE:
        log("FATAL -- Biopython is not installed.  pip install biopython")
        log("Biopython is required to download the mitochondrial reference from NCBI.")
        sys.exit(1)

    log("Downloading rCRS (NC_012920.1) from NCBI Entrez ...")

    # NCBI requires an email for Entrez API usage — this is a courtesy identifier.
    # It is not used for authentication.
    Entrez.email = "numt_pipeline@example.com"

    try:
        # efetch retrieves a sequence record from NCBI's Nucleotide database.
        # rettype="fasta" returns the sequence in FASTA format.
        # retmode="text" returns plain text (as opposed to XML).
        handle = Entrez.efetch(
            db="nucleotide",
            id="NC_012920.1",
            rettype="fasta",
            retmode="text",
        )
        fasta_content = handle.read()
        handle.close()
    except Exception as exc:
        log(f"FATAL -- Failed to download rCRS from NCBI: {exc}")
        log("Check your network connection.  Alternatively, download manually:")
        log("  wget 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
            "db=nucleotide&id=NC_012920.1&rettype=fasta&retmode=text' "
            f"-O {mito_fasta}")
        sys.exit(1)

    if not fasta_content or len(fasta_content.strip()) < 100:
        log("FATAL -- Downloaded rCRS appears empty or truncated.")
        sys.exit(1)

    mito_fasta.parent.mkdir(parents=True, exist_ok=True)
    with open(mito_fasta, "w") as fh:
        fh.write(fasta_content)

    log(f"Saved rCRS ({len(fasta_content):,} bytes) to: {mito_fasta}")


def build_blast_db(mito_fasta: Path, db_dir: Path, db_name: str) -> Path:
    """
    Create a local BLAST nucleotide database from the mitochondrial FASTA.

    Uses the `makeblastdb` command from the BLAST+ suite.  The resulting database
    files (.nhr, .nin, .nsq) are stored in the specified directory.

    makeblastdb creates an indexed representation of the FASTA sequences that
    BLAST can search efficiently.  For a single 16.5 kb sequence this is near-
    instantaneous, but the database format is required by blastn.

    Parameters
    ----------
    mito_fasta : Path
        Path to the mitochondrial genome FASTA file.
    db_dir : Path
        Directory where BLAST database files will be stored.
    db_name : str
        Base name for the BLAST database files.

    Returns
    -------
    Path
        Full path to the BLAST database (without file extensions).
    """
    # ── Check that makeblastdb is available ──────────────────────────────
    if shutil.which("makeblastdb") is None:
        log("FATAL -- 'makeblastdb' not found in PATH.")
        log("BLAST+ must be installed.  To install:")
        log("  macOS:   brew install blast")
        log("  Ubuntu:  sudo apt-get install ncbi-blast+")
        log("  conda:   conda install -c bioconda blast")
        log("  Or download from: https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html")
        sys.exit(1)

    db_dir.mkdir(parents=True, exist_ok=True)
    db_path = db_dir / db_name

    # Check if database already exists (look for the .nsq file)
    if (db_dir / f"{db_name}.nsq").exists() or (db_dir / f"{db_name}.ndb").exists():
        log(f"BLAST database already exists at: {db_path}")
        log("Skipping makeblastdb (delete .nsq/.ndb files to force rebuild)")
        return db_path

    log(f"Building BLAST database from: {mito_fasta}")

    cmd = [
        "makeblastdb",
        "-in", str(mito_fasta),
        "-dbtype", "nucl",              # Nucleotide database
        "-out", str(db_path),
        "-title", "Human_mitochondrial_genome_rCRS",
        "-parse_seqids",                # Index sequence IDs for retrieval
    ]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
        )
        log(f"makeblastdb completed successfully")
        if result.stdout.strip():
            log(f"  stdout: {result.stdout.strip()}")
    except subprocess.CalledProcessError as exc:
        log(f"FATAL -- makeblastdb failed with exit code {exc.returncode}")
        log(f"  stderr: {exc.stderr}")
        sys.exit(1)

    return db_path


# ─── Step 3: Run BLAST ──────────────────────────────────────────────────────

def run_blast(
    query_fasta: Path,
    db_path: Path,
    output_path: Path,
    evalue: float = 1e-3,
    word_size: int = 7,
    perc_identity: float = 80.0,
    qcov_hsp_perc: float = 50.0,
    num_threads: int = 4,
    task: str | None = None,
) -> Path:
    """
    Execute blastn to align query sequences against the mitochondrial database.

    The function determines whether to use 'blastn-short' or 'blastn' task mode
    based on the query sequences (unless overridden).  blastn-short is optimised
    for queries shorter than ~50 bp, using a word size of 7 and a scoring matrix
    tuned for short alignments.  For longer queries, standard blastn is more
    appropriate.

    Parameters
    ----------
    query_fasta : Path
        FASTA file of query sequences.
    db_path : Path
        Path to the BLAST database (without extensions).
    output_path : Path
        Where to write the tabular BLAST results.
    evalue : float
        E-value threshold.
    word_size : int
        Minimum word (seed) size for initiating alignments.
    perc_identity : float
        Minimum percent identity.
    qcov_hsp_perc : float
        Minimum query coverage per HSP (%).
    num_threads : int
        Number of CPU threads.
    task : str or None
        BLAST task ('blastn', 'blastn-short', or None for auto-detect).

    Returns
    -------
    Path
        Path to the BLAST output file.
    """
    # ── Verify blastn is available ───────────────────────────────────────
    if shutil.which("blastn") is None:
        log("FATAL -- 'blastn' not found in PATH.")
        log("BLAST+ must be installed.  See: brew install blast / conda install blast")
        sys.exit(1)

    # ── Auto-detect task based on query lengths ──────────────────────────
    # Read the query FASTA to determine sequence lengths.  If the median
    # length is < 100 bp, use blastn-short for better sensitivity on short
    # queries.  Otherwise use standard blastn.
    if task is None:
        lengths = []
        with open(query_fasta) as fh:
            current_len = 0
            for line in fh:
                if line.startswith(">"):
                    if current_len > 0:
                        lengths.append(current_len)
                    current_len = 0
                else:
                    current_len += len(line.strip())
            if current_len > 0:
                lengths.append(current_len)

        if not lengths:
            log("FATAL -- Query FASTA is empty (no sequences found).")
            sys.exit(1)

        median_len = sorted(lengths)[len(lengths) // 2]
        task = "blastn-short" if median_len < 100 else "blastn"
        log(f"Auto-detected BLAST task: {task} (median query length: {median_len} bp)")

    # ── Build the BLAST output format string ─────────────────────────────
    # Format 6 = tabular, no headers.  We request specific columns:
    outfmt_fields = (
        "qseqid sseqid pident length mismatch gapopen "
        "qstart qend sstart send evalue bitscore qlen slen"
    )
    outfmt = f"6 {outfmt_fields}"

    output_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "blastn",
        "-query", str(query_fasta),
        "-db", str(db_path),
        "-out", str(output_path),
        "-task", task,
        # E-value: the expected number of chance alignments with this score
        # or better in a database of this size.  Lower = more stringent.
        "-evalue", str(evalue),
        # Word size: minimum length of an exact-match seed to initiate
        # alignment extension.  Smaller = more sensitive, slower.
        "-word_size", str(word_size),
        # Percent identity: minimum fraction of matching bases in the
        # aligned region.  80% allows for divergence accumulated since
        # the NUMT insertion event.
        "-perc_identity", str(perc_identity),
        # Query coverage per HSP: minimum fraction of the query that
        # must be included in the alignment.  Prevents tiny partial hits.
        "-qcov_hsp_perc", str(qcov_hsp_perc),
        # Tabular output with specified columns
        "-outfmt", outfmt,
        # Number of parallel threads for the search
        "-num_threads", str(num_threads),
    ]

    log(f"Running BLAST: {task}")
    log(f"  E-value threshold:  {evalue}")
    log(f"  Word size:          {word_size}")
    log(f"  Min % identity:     {perc_identity}")
    log(f"  Min query coverage: {qcov_hsp_perc}%")
    log(f"  Threads:            {num_threads}")
    log(f"  Output:             {output_path}")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
        )
        if result.stderr.strip():
            # blastn often prints warnings to stderr (e.g., about sequence IDs)
            for line in result.stderr.strip().split("\n"):
                log(f"  blastn stderr: {line}")
    except subprocess.CalledProcessError as exc:
        log(f"FATAL -- blastn failed with exit code {exc.returncode}")
        log(f"  stderr: {exc.stderr}")
        sys.exit(1)

    # ── Check if results are empty ───────────────────────────────────────
    if not output_path.exists() or output_path.stat().st_size == 0:
        log("WARNING -- BLAST returned no hits (output file is empty).")
        return output_path

    n_lines = sum(1 for _ in open(output_path))
    log(f"BLAST completed: {n_lines:,} HSPs (high-scoring pairs) found")

    return output_path


# ─── Step 4: Parse and classify BLAST results ───────────────────────────────

def annotate_mtdna_region(sstart: int, send: int) -> tuple[str, str]:
    """
    Determine which mitochondrial gene/region a BLAST hit falls in.

    Given the subject (mtDNA) start and end coordinates of an alignment, this
    function checks which annotated feature(s) overlap.  If the hit spans
    multiple features, the one with the greatest overlap is reported.

    Parameters
    ----------
    sstart : int
        Start position on the mtDNA subject (1-based).
    send : int
        End position on the mtDNA subject (1-based).

    Returns
    -------
    tuple[str, str]
        (feature_name, feature_type) e.g. ("MT-ND1", "protein_coding")
        or ("intergenic", "intergenic") if no annotated feature overlaps.
    """
    # Ensure start <= end (BLAST may report them in either order depending
    # on the strand of the alignment)
    hit_start = min(sstart, send)
    hit_end = max(sstart, send)

    best_name = "intergenic"
    best_type = "intergenic"
    best_overlap = 0

    for name, feat_start, feat_end, feat_type in MTDNA_FEATURES:
        # Calculate overlap between the hit and this feature
        overlap_start = max(hit_start, feat_start)
        overlap_end = min(hit_end, feat_end)
        overlap = max(0, overlap_end - overlap_start + 1)

        if overlap > best_overlap:
            best_overlap = overlap
            best_name = name
            best_type = feat_type

    return best_name, best_type


def classify_hit(pident: float, length: int, evalue: float) -> str:
    """
    Assign a confidence tier to a BLAST hit based on alignment statistics.

    Tiered classification:

      STRONG_HIT:   pident >= 90% AND length >= 50 bp AND evalue <= 1e-5
        → High confidence.  The sequence is very similar to the mtDNA reference
          over a substantial length.  Likely a recent NUMT or one in a conserved
          region.

      MODERATE_HIT: pident >= 80% AND length >= 40 bp AND evalue <= 1e-3
        → Medium confidence.  Some divergence is present, possibly from an older
          insertion event or one that accumulated mutations.  Worth investigating.

      WEAK_HIT:     Anything else that passed the initial BLAST filters.
        → Low confidence.  Short or divergent matches that could be real but
          need additional evidence (e.g., flanking sequence context, repeat
          element analysis).

    Parameters
    ----------
    pident : float
        Percent identity of the alignment.
    length : int
        Length of the alignment in base pairs.
    evalue : float
        E-value of the alignment.

    Returns
    -------
    str
        One of "STRONG_HIT", "MODERATE_HIT", or "WEAK_HIT".
    """
    if pident >= 90.0 and length >= 50 and evalue <= 1e-5:
        return "STRONG_HIT"
    elif pident >= 80.0 and length >= 40 and evalue <= 1e-3:
        return "MODERATE_HIT"
    else:
        return "WEAK_HIT"


def parse_blast_results(blast_output: Path) -> pd.DataFrame:
    """
    Parse BLAST tabular output into a DataFrame with annotations.

    For each HSP (High-Scoring Pair), we:
      1. Parse the 14 tab-separated columns
      2. Classify the hit into a confidence tier (STRONG/MODERATE/WEAK)
      3. Annotate which mtDNA region the hit maps to
      4. Compute a combined score: pident * (length / qlen) for ranking

    Parameters
    ----------
    blast_output : Path
        Path to the BLAST tabular output (outfmt 6).

    Returns
    -------
    pd.DataFrame
        Annotated BLAST results.  Empty DataFrame if no results.
    """
    if not blast_output.exists() or blast_output.stat().st_size == 0:
        log("No BLAST results to parse (output file empty or missing).")
        return pd.DataFrame(columns=BLAST_COLUMNS + [
            "hit_tier", "mtdna_region", "mtdna_region_type", "combined_score"
        ])

    log(f"Parsing BLAST results from: {blast_output}")

    df = pd.read_csv(
        blast_output,
        sep="\t",
        header=None,
        names=BLAST_COLUMNS,
        dtype={
            "qseqid": str,
            "sseqid": str,
        },
    )

    log(f"Parsed {len(df):,} HSPs")

    # ── Classify each hit ────────────────────────────────────────────────
    df["hit_tier"] = df.apply(
        lambda row: classify_hit(row["pident"], row["length"], row["evalue"]),
        axis=1,
    )

    # ── Annotate mtDNA region ────────────────────────────────────────────
    regions = df.apply(
        lambda row: annotate_mtdna_region(int(row["sstart"]), int(row["send"])),
        axis=1,
    )
    df["mtdna_region"] = regions.apply(lambda x: x[0])
    df["mtdna_region_type"] = regions.apply(lambda x: x[1])

    # ── Combined score for ranking ───────────────────────────────────────
    # pident * query_coverage gives a single metric that rewards both high
    # identity and good query coverage.  A 95% identity hit covering 90%
    # of the query scores 85.5, while a 100% identity hit covering only
    # 30% of the query scores 30.0.
    df["query_coverage"] = (df["length"] / df["qlen"]) * 100.0
    df["combined_score"] = df["pident"] * (df["length"] / df["qlen"]) / 100.0

    # ── Extract AlleleID from qseqid ─────────────────────────────────────
    # The qseqid format is: AlleleID|GeneSymbol|ClinSig|LengthBp
    # We extract the numeric AlleleID for MongoDB lookups.
    df["allele_id"] = df["qseqid"].apply(
        lambda x: int(x.split("|")[0]) if "|" in str(x) else None
    )
    df["gene_symbol"] = df["qseqid"].apply(
        lambda x: x.split("|")[1] if "|" in str(x) and len(x.split("|")) > 1 else None
    )

    return df


# ─── Step 5: Update MongoDB with results ────────────────────────────────────

def update_mongodb(
    df_results: pd.DataFrame,
    query_allele_ids: set,
    host: str,
    port: int,
) -> None:
    """
    Write BLAST results back into MongoDB variant documents.

    For each variant with hits:
        $set: {
            blast_mito_results: [list of hit dicts],
            has_mito_hit: True,
            best_mito_hit_pident: <highest pident across all HSPs>,
        }

    For variants queried but without hits:
        $set: {has_mito_hit: False}

    Parameters
    ----------
    df_results : pd.DataFrame
        Annotated BLAST results DataFrame.
    query_allele_ids : set
        Set of all AlleleIDs that were included in the BLAST query.
    host : str
        MongoDB hostname.
    port : int
        MongoDB port.
    """
    if not PYMONGO_AVAILABLE:
        log("FATAL -- pymongo not available for MongoDB update.")
        sys.exit(1)

    log("Updating MongoDB with BLAST results ...")

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

    # ── Group results by AlleleID ────────────────────────────────────────
    allele_ids_with_hits = set()

    if not df_results.empty:
        for allele_id, group in df_results.groupby("allele_id"):
            if allele_id is None:
                continue

            allele_ids_with_hits.add(allele_id)

            # Build a list of hit dictionaries for this variant
            hit_list = []
            for _, row in group.iterrows():
                hit_list.append({
                    "pident": round(float(row["pident"]), 2),
                    "length": int(row["length"]),
                    "mismatch": int(row["mismatch"]),
                    "gapopen": int(row["gapopen"]),
                    "qstart": int(row["qstart"]),
                    "qend": int(row["qend"]),
                    "sstart": int(row["sstart"]),
                    "send": int(row["send"]),
                    "evalue": float(row["evalue"]),
                    "bitscore": float(row["bitscore"]),
                    "qlen": int(row["qlen"]),
                    "hit_tier": row["hit_tier"],
                    "mtdna_region": row["mtdna_region"],
                    "mtdna_region_type": row["mtdna_region_type"],
                    "query_coverage": round(float(row["query_coverage"]), 2),
                    "combined_score": round(float(row["combined_score"]), 4),
                })

            best_pident = max(h["pident"] for h in hit_list)

            operations.append(UpdateOne(
                filter={"allele_id": int(allele_id)},
                update={"$set": {
                    "blast_mito_results": hit_list,
                    "has_mito_hit": True,
                    "best_mito_hit_pident": best_pident,
                    "blast_mito_updated_at": datetime.now(timezone.utc),
                }},
            ))

    # ── Mark variants without hits ───────────────────────────────────────
    allele_ids_no_hits = query_allele_ids - allele_ids_with_hits
    for allele_id in allele_ids_no_hits:
        operations.append(UpdateOne(
            filter={"allele_id": int(allele_id)},
            update={"$set": {
                "has_mito_hit": False,
                "blast_mito_results": [],
                "best_mito_hit_pident": None,
                "blast_mito_updated_at": datetime.now(timezone.utc),
            }},
        ))

    if operations:
        result = collection.bulk_write(operations, ordered=False)
        log(f"MongoDB update complete: {result.modified_count} documents modified, "
            f"{result.matched_count} matched")
    else:
        log("No MongoDB update operations to perform.")

    client.close()


# ─── Step 6: Print results summary ──────────────────────────────────────────

def print_summary(df_results: pd.DataFrame, n_queries: int) -> None:
    """
    Print a human-readable summary of the BLAST results.

    Includes:
      - Total queries and total with hits
      - Breakdown by hit tier (STRONG / MODERATE / WEAK)
      - Top 20 hits ranked by combined score
      - Disease gene flags

    Parameters
    ----------
    df_results : pd.DataFrame
        Annotated BLAST results.
    n_queries : int
        Total number of query sequences sent to BLAST.
    """
    print("\n" + "=" * 80)
    print("  BLAST vs Mitochondrial Genome — Results Summary")
    print("=" * 80)

    if df_results.empty:
        print(f"\n  Total query sequences:  {n_queries:,}")
        print("  Total with mito hits:   0")
        print("\n  No mitochondrial hits found.")
        print("=" * 80 + "\n")
        return

    n_variants_with_hits = df_results["allele_id"].nunique()

    print(f"\n  Total query sequences:          {n_queries:>8,}")
    print(f"  Variants with mito hits:        {n_variants_with_hits:>8,}")
    print(f"  Total HSPs (alignments):        {len(df_results):>8,}")
    hit_rate = n_variants_with_hits / n_queries * 100 if n_queries > 0 else 0
    print(f"  Hit rate:                       {hit_rate:>7.1f}%")

    # ── Breakdown by hit tier ────────────────────────────────────────────
    print(f"\n  {'Breakdown by hit confidence tier':─<60}")
    tier_counts = df_results["hit_tier"].value_counts()
    for tier in ["STRONG_HIT", "MODERATE_HIT", "WEAK_HIT"]:
        count = tier_counts.get(tier, 0)
        # Count unique variants per tier
        if count > 0:
            n_variants = df_results[df_results["hit_tier"] == tier]["allele_id"].nunique()
        else:
            n_variants = 0
        print(f"    {tier:<25} {count:>6} HSPs across {n_variants:>4} variants")

    # ── Breakdown by mtDNA region ────────────────────────────────────────
    print(f"\n  {'Breakdown by mtDNA region':─<60}")
    region_counts = df_results["mtdna_region"].value_counts()
    for region, count in region_counts.head(15).items():
        print(f"    {region:<30} {count:>6} HSPs")

    # ── Top 20 hits by combined score ────────────────────────────────────
    print(f"\n  {'Top 20 hits by combined score (pident x coverage)':─<60}")
    # Keep best hit per variant for the top list
    best_per_variant = df_results.sort_values("combined_score", ascending=False)
    best_per_variant = best_per_variant.drop_duplicates(subset=["allele_id"], keep="first")
    top20 = best_per_variant.head(20)

    print(f"    {'AlleleID':<12} {'Gene':<10} {'%Ident':>7} {'Len':>5} "
          f"{'E-value':>10} {'mtDNA region':<18} {'Tier':<14} {'Score':>6}")
    print(f"    {'─' * 12} {'─' * 10} {'─' * 7} {'─' * 5} "
          f"{'─' * 10} {'─' * 18} {'─' * 14} {'─' * 6}")

    for _, row in top20.iterrows():
        allele_id = row.get("allele_id", "?")
        gene = row.get("gene_symbol", "NA") or "NA"
        pident = row["pident"]
        length = row["length"]
        evalue = row["evalue"]
        region = row["mtdna_region"]
        tier = row["hit_tier"]
        score = row["combined_score"]

        # Flag disease genes with a marker
        disease_flag = ""
        if gene in DISEASE_GENES:
            disease_flag = " *** DISEASE GENE ***"

        print(f"    {str(allele_id):<12} {gene:<10} {pident:>7.1f} {length:>5} "
              f"{evalue:>10.1e} {region:<18} {tier:<14} {score:>6.2f}"
              f"{disease_flag}")

    # ── Disease gene flag summary ────────────────────────────────────────
    disease_hits = df_results[df_results["gene_symbol"].isin(DISEASE_GENES)]
    if not disease_hits.empty:
        n_disease = disease_hits["allele_id"].nunique()
        print(f"\n  *** ATTENTION: {n_disease} hit(s) found in known disease genes ***")
        for _, row in disease_hits.drop_duplicates("allele_id").iterrows():
            gene = row["gene_symbol"]
            allele = row["allele_id"]
            pident = row["pident"]
            print(f"      AlleleID {allele} in {gene} ({pident:.1f}% identity)")
        print("  These variants merit priority manual review as potential")
        print("  misclassified NUMTs causing disease through gene disruption.")

    print("\n" + "=" * 80 + "\n")


# ─── Step 7: Relaxation logic ───────────────────────────────────────────────

def should_relax(df_results: pd.DataFrame, threshold: int = 10) -> bool:
    """
    Determine whether BLAST parameters should be relaxed.

    If fewer than `threshold` unique variants have hits, the initial parameters
    may have been too stringent.  This can happen when:
      - The query sequences are very short (< 60 bp)
      - The NUMTs are ancient and have diverged significantly from rCRS
      - The insertion contains only a fragment of a mitochondrial gene

    Parameters
    ----------
    df_results : pd.DataFrame
        BLAST results from the initial run.
    threshold : int
        Minimum number of hits before relaxation is triggered.

    Returns
    -------
    bool
        True if parameters should be relaxed.
    """
    if df_results.empty:
        return True
    n_variants = df_results["allele_id"].nunique()
    return n_variants < threshold


# ─── Helper: extract AlleleIDs from FASTA ────────────────────────────────────

def get_allele_ids_from_fasta(fasta_path: Path) -> set:
    """
    Extract the set of AlleleIDs from the query FASTA headers.

    Parameters
    ----------
    fasta_path : Path
        Path to the query FASTA file.

    Returns
    -------
    set
        Set of integer AlleleIDs.
    """
    allele_ids = set()
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                header = line.strip().lstrip(">")
                parts = header.split("|")
                try:
                    allele_ids.add(int(parts[0]))
                except (ValueError, IndexError):
                    pass
    return allele_ids


# ─── CLI entry point ────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    """Define and parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Step 4 of the NUMT pipeline: BLAST extracted insertion sequences "
            "against the human mitochondrial genome (rCRS, NC_012920.1) to "
            "identify variants of mitochondrial origin."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # ── Input/output paths ───────────────────────────────────────────────
    parser.add_argument(
        "--query-fasta",
        type=Path,
        default=QUERY_FASTA,
        help=f"Path for the query FASTA (default: {QUERY_FASTA.relative_to(PROJECT_ROOT)})",
    )
    parser.add_argument(
        "--blast-output",
        type=Path,
        default=BLAST_OUTPUT,
        help=f"Path for BLAST results TSV (default: {BLAST_OUTPUT.relative_to(PROJECT_ROOT)})",
    )
    parser.add_argument(
        "--mito-fasta",
        type=Path,
        default=MITO_FASTA,
        help=f"Path to rCRS FASTA (default: {MITO_FASTA.relative_to(PROJECT_ROOT)})",
    )

    # ── MongoDB connection ───────────────────────────────────────────────
    parser.add_argument(
        "--host", type=str, default=MONGO_HOST,
        help=f"MongoDB hostname (default: {MONGO_HOST})",
    )
    parser.add_argument(
        "--port", type=int, default=MONGO_PORT,
        help=f"MongoDB port (default: {MONGO_PORT})",
    )

    # ── BLAST parameters (all configurable) ──────────────────────────────
    parser.add_argument(
        "--evalue", type=float, default=1e-3,
        help="E-value threshold (default: 1e-3, permissive for short seqs)",
    )
    parser.add_argument(
        "--word-size", type=int, default=7,
        help="Minimum word (seed) size (default: 7, sensitive for short queries)",
    )
    parser.add_argument(
        "--perc-identity", type=float, default=80.0,
        help="Minimum percent identity (default: 80, allows NUMT divergence)",
    )
    parser.add_argument(
        "--qcov-hsp-perc", type=float, default=50.0,
        help="Minimum query coverage per HSP (default: 50%%)",
    )
    parser.add_argument(
        "--num-threads", type=int, default=4,
        help="Number of BLAST threads (default: 4)",
    )
    parser.add_argument(
        "--task", type=str, default=None, choices=["blastn", "blastn-short"],
        help="BLAST task (default: auto-detect based on query lengths)",
    )

    # ── Behaviour flags ──────────────────────────────────────────────────
    parser.add_argument(
        "--relax", action="store_true",
        help=(
            "Automatically re-run BLAST with relaxed parameters if fewer than "
            "10 hits are found.  Relaxed: word_size=4, perc_identity=70."
        ),
    )
    parser.add_argument(
        "--skip-mongodb", action="store_true",
        help="Skip MongoDB query (use existing FASTA) and skip result updates.",
    )
    parser.add_argument(
        "--skip-download", action="store_true",
        help="Skip rCRS download (assume it already exists).",
    )

    return parser.parse_args()


# ─── Main orchestration ─────────────────────────────────────────────────────

def main() -> None:
    """Orchestrate the mitochondrial BLAST pipeline."""
    args = parse_args()

    log("=" * 60)
    log("NUMT Pipeline -- Step 4: BLAST vs Mitochondrial Genome")
    log("=" * 60)

    # ── Step 1: Prepare query FASTA from MongoDB ─────────────────────────
    if args.skip_mongodb:
        log("Skipping MongoDB query (--skip-mongodb).  Using existing FASTA.")
        if not args.query_fasta.exists():
            log(f"FATAL -- Query FASTA not found: {args.query_fasta}")
            log("Cannot skip MongoDB if no FASTA exists.  Run without --skip-mongodb.")
            sys.exit(1)
        n_queries = sum(1 for line in open(args.query_fasta) if line.startswith(">"))
    else:
        n_queries = prepare_query_fasta(
            host=args.host,
            port=args.port,
            output_fasta=args.query_fasta,
        )

    if n_queries == 0:
        log("FATAL -- No query sequences available.  Nothing to BLAST.")
        sys.exit(1)

    log(f"Query FASTA contains {n_queries:,} sequences")

    # ── Step 2: Download rCRS and build BLAST database ───────────────────
    if not args.skip_download:
        download_rcrs(args.mito_fasta)
    else:
        log("Skipping rCRS download (--skip-download)")
        if not args.mito_fasta.exists():
            log(f"FATAL -- rCRS FASTA not found: {args.mito_fasta}")
            sys.exit(1)

    db_path = build_blast_db(args.mito_fasta, BLAST_DB_DIR, BLAST_DB_NAME)

    # ── Step 3: Run BLAST ────────────────────────────────────────────────
    blast_output = run_blast(
        query_fasta=args.query_fasta,
        db_path=db_path,
        output_path=args.blast_output,
        evalue=args.evalue,
        word_size=args.word_size,
        perc_identity=args.perc_identity,
        qcov_hsp_perc=args.qcov_hsp_perc,
        num_threads=args.num_threads,
        task=args.task,
    )

    # ── Step 4: Parse and classify results ───────────────────────────────
    df_results = parse_blast_results(blast_output)

    # ── Step 7: Relaxation logic ─────────────────────────────────────────
    # If --relax is set and fewer than 10 variants have hits, re-run BLAST
    # with more permissive parameters to catch divergent NUMTs.
    if args.relax and should_relax(df_results):
        n_initial = df_results["allele_id"].nunique() if not df_results.empty else 0
        log("")
        log("=" * 60)
        log("RELAXATION TRIGGERED")
        log(f"Only {n_initial} variant(s) had hits in the initial run.")
        log("Re-running BLAST with relaxed parameters to catch divergent NUMTs:")
        log("  word_size:     7 -> 4  (shorter seeds = more sensitive)")
        log("  perc_identity: 80 -> 70  (allow more divergence from rCRS)")
        log("")
        log("Biological justification: NUMTs that inserted into the nuclear")
        log("genome millions of years ago may have accumulated substantial")
        log("mutations.  The nuclear mutation rate is ~2.5e-8 per site per")
        log("generation, so a NUMT from ~5 Mya could show ~15-20% divergence.")
        log("=" * 60)

        # Re-run with relaxed parameters
        relaxed_output = args.blast_output.with_stem(
            args.blast_output.stem + "_relaxed"
        )
        run_blast(
            query_fasta=args.query_fasta,
            db_path=db_path,
            output_path=relaxed_output,
            evalue=args.evalue,
            word_size=4,                   # More sensitive seed matching
            perc_identity=70.0,            # Allow greater divergence
            qcov_hsp_perc=args.qcov_hsp_perc,
            num_threads=args.num_threads,
            task=args.task,
        )

        df_relaxed = parse_blast_results(relaxed_output)

        if not df_relaxed.empty:
            n_relaxed = df_relaxed["allele_id"].nunique()
            log(f"Relaxed search found {n_relaxed} variant(s) with hits "
                f"(was {n_initial})")

            # Merge: use relaxed results if they found more
            if n_relaxed > n_initial:
                df_results = df_relaxed
                blast_output = relaxed_output
                log("Using relaxed results (more hits found).")
            else:
                log("Relaxed search did not improve results.  Keeping original.")
        else:
            log("Relaxed search also found no hits.")

    # ── Step 5: Update MongoDB ───────────────────────────────────────────
    query_allele_ids = get_allele_ids_from_fasta(args.query_fasta)

    if not args.skip_mongodb:
        update_mongodb(
            df_results=df_results,
            query_allele_ids=query_allele_ids,
            host=args.host,
            port=args.port,
        )
    else:
        log("Skipping MongoDB update (--skip-mongodb)")

    # ── Step 6: Print results summary ────────────────────────────────────
    print_summary(df_results, n_queries)

    log("Done.")


if __name__ == "__main__":
    main()
