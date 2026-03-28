# 1. NUMT-ClinVar Pipeline

**Identification of unrecognised Nuclear Mitochondrial DNA Segments (NUMTs) among pathogenic ClinVar insertions**

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

---

## 1.1. Scientific Background

NUMTs (Nuclear Mitochondrial DNA Segments) are fragments of mitochondrial DNA that have been transferred and integrated into the nuclear genome over evolutionary time. When a NUMT insertion occurs within a coding exon, it can disrupt gene function and cause Mendelian disease.

**Index case:** A 138 bp mitochondrial fragment (MT-ND5, 98% identity to rCRS) inserted into exon 17 of *TSC2* causes Tuberous Sclerosis Complex (NM_000548.5(TSC2):c.1830_1831ins[138]; CHU Angers).

**Hypothesis:** Among the large insertions (≥50 bp) classified as pathogenic in ClinVar, a subset may represent unrecognised de novo NUMTs — mischaracterised by standard short-read sequencing pipelines.

**Target class — `NOVEL_NUMT_CANDIDATE`:** Variants whose inserted sequence aligns to the mitochondrial genome (rCRS) but has **no hit** in the human nuclear reference genome (hg38). These represent putative de novo mtDNA insertions.

---

## 1.2. Pipeline Overview

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
    ▼ Step 4 — BLAST vs rCRS (NC_012920.1)
Mitochondrial hits  →  4_BLAST_Mito sheet
    │
    ▼ Step 5 — BLAST vs hg38
Nuclear hits        →  5_BLAST_Genome sheet
    │
    ▼ Classification
numt_class per variant  →  6_Classification + NOVEL_CANDIDATES sheets
```

### 1.2.1. Classification

| Class | Definition |
|---|---|
| `NOVEL_NUMT_CANDIDATE` ⭐ | mito hit + **no** nuclear hit → putative de novo pathogenic NUMT |
| `CONFIRMED_NUMT` | mito hit + nuclear hit at the ClinVar locus (ancient copy in hg38) |
| `KNOWN_REFERENCE_NUMT` | mito hit + nuclear hit elsewhere in hg38 |

> **Note:** `CONFIRMED_NUMT` denotes a *bioinformatic* classification — the mtDNA-like sequence is already present in the hg38 reference at the reported locus (e.g. MSH2/MTND5P11). This is distinct from *literature recognition* of NUMT origin, tracked in the `literature_numt_mention` column — only MCOLN1 (Goldin 2004) has a published article acknowledging the mitochondrial origin of the insertion.
| `DIRECT_MITO_REF` | mitochondrial origin explicit in the HGVS name |
| `MITO_HIT_ONLY` | mito hit, genome BLAST not performed |
| `NO_MITO_HIT` | no mitochondrial alignment |

---

## 1.3. Quick Start

```bash
# Clone
git clone github.com/kaarud/numt‑blast‑pipeline
cd NUMT

# Set up environment
micromamba env create -f 01_code/numt_pipeline/environment.yml
micromamba activate NUMT

# Run the full pipeline → multi-sheet Excel output
make excel

# Fast run (skip genome BLAST)
make excel-skip-genome
```

Output: `02_data/processed/numt_pipeline_results.xlsx`

---

## 1.4. Installation

### 1.4.1. Requirements

- Python ≥ 3.11
- BLAST+ ≥ 2.17 (`blastn`, `makeblastdb`)
- micromamba / conda

### 1.4.2. Dependencies

```bash
micromamba env create -f 01_code/numt_pipeline/environment.yml
micromamba activate NUMT
```

| Package | Version tested | Purpose |
|---|---|---|
| Python | 3.11 | Runtime |
| pandas | ≥ 2.0 | Data processing |
| biopython | ≥ 1.83 | NCBI Entrez, BLAST parsing |
| blast (BLAST+) | ≥ 2.17 | Sequence alignment |
| openpyxl | ≥ 3.1 | Excel output |
| requests | — | ClinVar download |
| tqdm | — | Progress bars |

### 1.4.3. Reference databases

The pipeline auto-builds the mitochondrial BLAST database from rCRS on first run.

For genome BLAST (Step 5), a local hg38 BLAST database is required (~15 GB):
```bash
# Download hg38 FASTA, then:
makeblastdb -in hg38.fa -dbtype nucl -out 02_data/raw/blast_db/hg38

# Then run:
make excel GENOME_DB=02_data/raw/blast_db/hg38
```

Alternatively, use remote NCBI BLAST (slower, no local db required):
```bash
python 01_code/numt_pipeline/run_pipeline_excel.py --email your@email.com
```

---

## 1.5. Usage

### 1.5.1. Makefile commands

```bash
make excel                  # Full pipeline → Excel output (recommended)
make excel-skip-genome      # Skip genome BLAST (fast, partial classification)
make step1                  # Step 1 only: ClinVar download + filter
make step2                  # Step 2 only: HGVS sequence extraction
make step4                  # Step 4 only: mitochondrial BLAST
make step5                  # Step 5 only: genome BLAST + classification
make status                 # Show status of all intermediate files
make clean                  # Remove intermediate files
```

### 1.5.2. Resume from a step

```bash
# Resume from step 3 (steps 1-2 already done)
python 01_code/numt_pipeline/run_pipeline_excel.py --from-step 3 --skip-genome
```

### 1.5.3. Individual scripts

```bash
# Step 1 — Download and filter ClinVar
python 01_code/numt_pipeline/01_clinvar_fetch.py [--force] [--min-length 50]

# Step 2 — Extract sequences from HGVS notation
python 01_code/numt_pipeline/02_sequence_extract.py

# Step 4 — BLAST vs mitochondrial genome
python 01_code/numt_pipeline/04_blast_mito.py

# Step 5 — BLAST vs hg38 + classification
python 01_code/numt_pipeline/05_blast_genome.py [--genome-db /path/to/hg38] [--skip-genome]

# Step 6 — Insertion mechanism analysis (NOVEL_NUMT_CANDIDATE only)
python 01_code/numt_pipeline/06_mechanism.py
```

---

## 1.6. Output

The pipeline produces a multi-sheet Excel workbook at `02_data/processed/numt_pipeline_results.xlsx`:

| Sheet | Content |
|---|---|
| `1_ClinVar_Raw` | Filtered ClinVar insertions |
| `2_Sequences` | After HGVS parsing |
| `3_Deduplicated` | 1 row per AlleleID (GRCh38 preferred) |
| `4_BLAST_Mito` | Mitochondrial BLAST hits (1 row per HSP) |
| `5_BLAST_Genome` | Nuclear genome BLAST hits |
| `6_Classification` | Final colour-coded classification |
| **`NOVEL_CANDIDATES`** | ⭐ mito+ / nuclear− → main targets |
| `MITO_AND_NUCLEAR` | mito+ / nuclear+ |
| `DIRECT_MITO_REF` | HGVS-annotated mitochondrial origin |

### 1.6.1. Results (ClinVar 2025, genome BLAST complete)

| Class | N | Notes |
|---|---|---|
| `NO_MITO_HIT` | 854 | No mitochondrial homology — excluded |
| `NOVEL_NUMT_CANDIDATE` ⭐ | **20** | No nuclear copy in hg38 — putative de novo NUMTs |
| `CONFIRMED_NUMT` | **1** | MSH2/Lynch syndrome — ancient nuclear copy at ClinVar locus (chr5, MTND5P11) |
| `KNOWN_REFERENCE_NUMT` | 0 | — |
| `DIRECT_MITO_REF` | 0 | — |

**81% (17/21) of all mito-hit candidates carry Pathogenic or Likely pathogenic classification** (16 NOVEL + 1 CONFIRMED_NUMT — MSH2/Lynch syndrome).

See `02_data/data_interpretation.md` for full results and working interpretation.

---

## 1.7. BLAST parameters

### 1.7.1. Mitochondrial BLAST (Step 4)

| Parameter | Value | Rationale |
|---|---|---|
| `evalue` | 1e-3 | Permissive for small database |
| `word_size` | 7 | Sensitive for short sequences |
| `perc_identity` | 80% | Tolerates ~20% divergence (ancient NUMTs) |
| `qcov_hsp_perc` | 50% | At least 50% of query aligned |

Confidence tiers: `STRONG_HIT` (pident ≥ 90%, len ≥ 50 bp, e-value ≤ 1e-5) / `MODERATE_HIT` / `WEAK_HIT`

### 1.7.2. Nuclear BLAST (Step 5)

Only canonical chromosomes (chr1–22, chrX, chrY) count as nuclear hits. Alternative loci (NW_/NT_ accessions) are excluded.

---

## 1.8. Repository Structure

```
01_code/
  numt_pipeline/            # Pipeline scripts
    01_clinvar_fetch.py     # Step 1: ClinVar download + filter
    02_sequence_extract.py  # Step 2: HGVS sequence extraction
    04_blast_mito.py        # Step 4: mitochondrial BLAST
    05_blast_genome.py      # Step 5: genome BLAST + classification
    06_mechanism.py         # Step 6: insertion mechanism prediction
    run_pipeline_excel.py   # Main orchestrator
    environment.yml         # Conda environment
02_data/
  raw/                      # ClinVar download, rCRS FASTA, BLAST databases (not versioned)
  processed/                # Intermediate files + Excel output (not versioned)
03_manuscript/              # Manuscript files
Makefile                    # Pipeline commands
```

---

## 1.9. Citation

If you use this pipeline in your research, please cite:

> Durand A. *et al.* Identification of unrecognised nuclear mitochondrial DNA segments (NUMTs) among pathogenic ClinVar insertions. [*in preparation*]

See also [CITATION.cff](CITATION.cff) for machine-readable citation metadata.

---

## 1.10. License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.

---

## 1.11. Contact

Aksel Durand — aksel.durand@chu-angers.fr
CHU Angers, Service de Génétique
