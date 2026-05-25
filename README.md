# NUMT-ClinVar Pipeline

**Identification of unrecognised Nuclear Mitochondrial DNA Segments (NUMTs) among pathogenic ClinVar insertions**

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

---

## Scientific Background

NUMTs (Nuclear Mitochondrial DNA Segments) are fragments of mitochondrial DNA that have been transferred and integrated into the nuclear genome. When a NUMT insertion occurs within a coding exon, it can disrupt gene function and cause Mendelian disease.

**Index case:** A 138 bp mitochondrial fragment (MT-ND5, 98% identity to rCRS) inserted into exon 17 of *TSC2* causes Tuberous Sclerosis Complex (NM_000548.5(TSC2):c.1830_1831ins[138]; CHU Angers).

**Hypothesis:** Among the large insertions (≥50 bp) classified as pathogenic in ClinVar, a subset may represent unrecognised de novo NUMTs — mischaracterised by standard short-read sequencing pipelines.

---

## Pipeline Overview

```
ClinVar FTP (~350 MB)
    │
    ▼ Step 1 — 01_clinvar_fetch.py
Insertions ≥ 50 bp  →  clinvar_insertions_gt50bp.tsv
    │                   (~1,700 rows, ~875 unique AlleleIDs)
    ▼ Step 2 — 02_sequence_extract.py
Sequences from HGVS  →  extracted_sequences.fasta
    │
    ▼ Step 3 — Deduplication (one entry per AlleleID, GRCh38 preferred)
    │
    ▼ Step 4 — 04_blast_mito.py — BLAST vs rCRS (NC_012920.1)
Mitochondrial hits  →  4_BLAST_Mito sheet
    │
    ▼ Step 5 — 05_blast_genome.py — BLAST vs hg38
Nuclear hits        →  5_BLAST_Genome sheet
    │
    ▼ Step 6 — 06_mechanism.py — Classification + mechanism analysis
numt_class per variant  →  6_Classification + NOVEL_CANDIDATES sheets
```

### Classification

| Class | Definition |
|---|---|
| `NOVEL_NUMT_CANDIDATE` | mito hit + **no** nuclear hit → putative de novo pathogenic NUMT |
| `CONFIRMED_NUMT` | mito hit + nuclear hit at the ClinVar locus (ancient copy in hg38) |
| `KNOWN_REFERENCE_NUMT` | mito hit + nuclear hit elsewhere in hg38 |
| `MITO_HIT_ONLY` | mito hit, genome BLAST not performed |
| `NO_MITO_HIT` | no mitochondrial alignment |

> **Note:** `CONFIRMED_NUMT` denotes a *bioinformatic* classification — the mtDNA-like sequence is already present in the hg38 reference at the reported locus (e.g. MSH2/MTND5P11). This is distinct from *literature recognition* of NUMT origin.

---

## Quick Start

```bash
# Clone
git clone https://github.com/kaarud/numt-blast-pipeline.git
cd numt-blast-pipeline

# Set up environment (micromamba, mamba, or conda)
micromamba env create -f 01_code/numt_pipeline/environment.yml
micromamba activate NUMT

# Run the full pipeline → multi-sheet Excel output
make excel

# Or skip genome BLAST (faster, partial classification only)
make excel-skip-genome
```

Output: `02_data/processed/numt_pipeline_results.xlsx`

---

## Requirements

- Python ≥ 3.11
- BLAST+ ≥ 2.17 (`blastn`, `makeblastdb`)
- micromamba, mamba, or conda

### Python dependencies

Installed automatically via `environment.yml`:

| Package | Purpose |
|---|---|
| pandas | Data processing |
| biopython | NCBI Entrez, BLAST parsing |
| blast (BLAST+) | Sequence alignment (via bioconda) |
| openpyxl | Excel output |
| requests | ClinVar download |
| tqdm | Progress bars |

### Reference databases

The pipeline auto-builds the mitochondrial BLAST database from rCRS (NC_012920.1) on first run.

For genome BLAST (Step 5), you need either a local hg38 database or an NCBI email for remote BLAST:

```bash
# Option A — local hg38 BLAST database (~15 GB)
makeblastdb -in hg38.fa -dbtype nucl -out 02_data/raw/blast_db/hg38
make excel GENOME_DB=02_data/raw/blast_db/hg38

# Option B — remote NCBI BLAST (slower, no local db required)
python 01_code/numt_pipeline/run_pipeline_excel.py --email your@email.com
```

---

## Usage

### Makefile commands

```bash
make excel              # Full pipeline → Excel output (recommended)
make excel-skip-genome  # Skip genome BLAST (fast, mito-only classification)
make step1              # Step 1 only: ClinVar download + filter
make step2              # Step 2 only: HGVS sequence extraction
make step4              # Step 4 only: mitochondrial BLAST
make step5              # Step 5 only: genome BLAST + classification
make step6              # Step 6 only: mechanism analysis
make status             # Show status of all intermediate files
make clean              # Remove intermediate files
```

### Resume from a specific step

```bash
# Resume from step 4 (steps 1-2 already done)
python 01_code/numt_pipeline/run_pipeline_excel.py --from-step 4 --skip-genome
```

### Individual scripts

```bash
# Step 1 — Download and filter ClinVar
python 01_code/numt_pipeline/01_clinvar_fetch.py [--force] [--min-length 50]

# Step 2 — Extract sequences from HGVS notation
python 01_code/numt_pipeline/02_sequence_extract.py

# Step 4 — BLAST vs mitochondrial genome
python 01_code/numt_pipeline/04_blast_mito.py

# Step 5 — BLAST vs hg38 + classification
python 01_code/numt_pipeline/05_blast_genome.py [--genome-db /path/to/hg38] [--skip-genome]

# Step 6 — Insertion mechanism analysis
python 01_code/numt_pipeline/06_mechanism.py
```

---

## Output

The pipeline produces a multi-sheet Excel workbook at `02_data/processed/numt_pipeline_results.xlsx`:

| Sheet | Content |
|---|---|
| `1_ClinVar_Raw` | Filtered ClinVar insertions |
| `2_Sequences` | After HGVS parsing |
| `3_Deduplicated` | 1 row per AlleleID (GRCh38 preferred) |
| `4_BLAST_Mito` | Mitochondrial BLAST hits (1 row per HSP) |
| `5_BLAST_Genome` | Nuclear genome BLAST hits |
| `6_Classification` | Final colour-coded classification |
| **`NOVEL_CANDIDATES`** | mito+ / nuclear− → main targets |
| `MITO_AND_NUCLEAR` | mito+ / nuclear+ |

### Expected results (ClinVar 2025)

| Class | N |
|---|---|
| `NO_MITO_HIT` | 854 |
| `NOVEL_NUMT_CANDIDATE` | 20 |
| `CONFIRMED_NUMT` | 1 |

---

## BLAST Parameters

### Mitochondrial BLAST (Step 4)

| Parameter | Value | Rationale |
|---|---|---|
| `evalue` | 1e-3 | Permissive for small database |
| `word_size` | 7 | Sensitive for short sequences |
| `perc_identity` | 80% | Tolerates ~20% divergence (ancient NUMTs) |
| `qcov_hsp_perc` | 50% | At least 50% of query aligned |

Confidence tiers: `STRONG_HIT` (pident ≥ 90%, len ≥ 50 bp, e-value ≤ 1e-5) / `MODERATE_HIT` / `WEAK_HIT`

### Nuclear BLAST (Step 5)

Only canonical chromosomes (chr1–22, chrX, chrY) count as nuclear hits. Alternative loci (NW_/NT_ accessions) are excluded.

---

## Repository Structure

```
01_code/numt_pipeline/
    01_clinvar_fetch.py         # Step 1: ClinVar download + filter
    02_sequence_extract.py      # Step 2: HGVS sequence extraction
    04_blast_mito.py            # Step 4: mitochondrial BLAST
    05_blast_genome.py          # Step 5: genome BLAST + classification
    06_mechanism.py             # Step 6: insertion mechanism prediction
    run_pipeline_excel.py       # Main orchestrator (steps 1–6)
    environment.yml             # Conda/micromamba environment
02_data/
    raw/                        # ClinVar download, rCRS, BLAST databases (not versioned)
    processed/                  # Intermediate files + Excel output (not versioned)
Makefile                        # Pipeline commands
```

---

## Citation

If you use this pipeline in your research, please cite:

> Durand A. *et al.* Identification of unrecognised nuclear mitochondrial DNA segments (NUMTs) among pathogenic ClinVar insertions. [*in preparation*]

See also [CITATION.cff](CITATION.cff) for machine-readable citation metadata.

---

## License

MIT — see [LICENSE](LICENSE) for details.

---

## Contact

Aksel Durand — aksel.durand@chu-angers.fr
CHU Angers, Service de Génétique
