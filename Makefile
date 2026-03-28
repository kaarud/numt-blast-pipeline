# ─────────────────────────────────────────────────────────────────────────────
# NUMT Pipeline — Makefile
# Usage:
#   make excel          Pipeline complète sans MongoDB → Excel multi-feuilles
#   make excel-skip-genome  Idem, sans BLAST vs génome (plus rapide)
#   make step1          Run Step 1: fetch ClinVar data
#   make step2          Run Step 2: extract sequences (HGVS regex)
#   make step3          Run Step 3: store in MongoDB  [LEGACY - non utilisé]
#   make step4          Run Step 4: BLAST vs mitochondrial genome
#   make step5          Run Step 5: BLAST vs human genome
#   make step6          Run Step 6: insertion mechanism analysis (standalone)
#   make excel          Pipeline Excel complète, inclus mécanisme (étape 7)
#   make all            Pipeline Excel complète (= make excel)
#   make status         Show current state of processed files
#   make clean          Remove all processed/intermediate files
# ─────────────────────────────────────────────────────────────────────────────

PYTHON       := python3
PIPELINE     := 01_code/numt_pipeline
PIPELINE_NEW := 01_code/numt_pipeline
EMAIL        := durand.aksel@gmail.com
DATA         := 02_data/processed

# ── Output files (used as Make targets for dependency tracking) ───────────────
CLINVAR_TSV  := $(DATA)/clinvar_insertions_gt50bp.tsv
SEQS_TSV     := $(DATA)/clinvar_insertions_with_sequences.tsv
FASTA        := $(DATA)/extracted_sequences.fasta
MITO_RESULTS := $(DATA)/blast_mito_results.tsv
NUMT_SUMMARY := $(DATA)/numt_candidates_summary.tsv


EXCEL_OUTPUT     := $(DATA)/numt_pipeline_results.xlsx
MECHANISM_JSON   := $(DATA)/mechanism_analysis.json
MECHANISM_TSV    := $(DATA)/mechanism_summary.tsv
LAST_MITO_TSV    := $(DATA)/excel_last_mito.tsv
BLAST_LAST_TSV   := $(DATA)/blast_vs_last_comparison.tsv
BLAST_LAST_RPT   := $(DATA)/blast_vs_last_report.txt
FLANKING_DIR     := $(DATA)/flanking_sequences
BLAST_GENOME_CACHE := $(DATA)/blast_genome_cache

.PHONY: all excel excel-skip-genome excel-validate-last step1 step2 step3 step4 step5 step6 prdm9 show-alignments compare-last status clean

# ── Pipeline Excel (sans MongoDB) ─────────────────────────────────────────────
all: excel

excel:
	@echo ""
	@echo "════════════════════════════════════════════════════════════"
	@echo "  Pipeline NUMT — Excel multi-feuilles (sans MongoDB)"
	@echo "════════════════════════════════════════════════════════════"
	$(PYTHON) $(PIPELINE_NEW)/run_pipeline_excel.py --email $(EMAIL)
	@echo ""
	@echo "  Résultat : $(EXCEL_OUTPUT)"
	@echo ""

excel-skip-genome:
	@echo ""
	@echo "════════════════════════════════════════════════════════════"
	@echo "  Pipeline NUMT — Excel (BLAST mito uniquement)"
	@echo "════════════════════════════════════════════════════════════"
	$(PYTHON) $(PIPELINE_NEW)/run_pipeline_excel.py --skip-genome
	@echo ""
	@echo "  Résultat : $(EXCEL_OUTPUT)"
	@echo ""

excel-validate-last:
	@echo ""
	@echo "════════════════════════════════════════════════════════════"
	@echo "  Pipeline NUMT — Excel + LAST validation (skip genome)"
	@echo "════════════════════════════════════════════════════════════"
	$(PYTHON) $(PIPELINE_NEW)/run_pipeline_excel.py --validate-last --skip-genome
	@echo ""
	@echo "  Résultat : $(EXCEL_OUTPUT)"
	@echo ""


# ── Step 1: Fetch ClinVar insertions > 50 bp ─────────────────────────────────
step1: $(CLINVAR_TSV)

$(CLINVAR_TSV):
	@echo ""
	@echo "════════════════════════════════════════════════════════════"
	@echo "  Step 1 — Fetch ClinVar insertions > 50 bp"
	@echo "════════════════════════════════════════════════════════════"
	$(PYTHON) $(PIPELINE)/01_clinvar_fetch.py
	@echo ""
	@echo "─── Vérification Step 1 ────────────────────────────────────"
	@$(PYTHON) -c "\
import pandas as pd; \
df = pd.read_csv('$(CLINVAR_TSV)', sep='\t', low_memory=False); \
print(f'  Lignes totales    : {len(df):,}'); \
print(f'  AlleleIDs uniques : {df[\"#AlleleID\"].nunique():,}'); \
print(); \
print('  Assemblies:'); \
[print(f'    {k}: {v:,}') for k,v in df['Assembly'].value_counts().items()]; \
print(); \
print('  Signification clinique (top 5):'); \
[print(f'    {k}: {v:,}') for k,v in list(df['ClinicalSignificance'].value_counts().items())[:5]]; \
"
	@echo "────────────────────────────────────────────────────────────"
	@echo "  Step 1 OK. Vérifier les comptes ci-dessus."
	@echo "  Si OK, lancer : make step2"
	@echo ""


# ── Step 2: Extract sequences from HGVS notation ─────────────────────────────
step2: $(SEQS_TSV) $(FASTA)

$(SEQS_TSV) $(FASTA): $(CLINVAR_TSV)
	@echo ""
	@echo "════════════════════════════════════════════════════════════"
	@echo "  Step 2 — Extraction séquences (regex HGVS)"
	@echo "════════════════════════════════════════════════════════════"
	$(PYTHON) $(PIPELINE)/02_sequence_extract.py
	@echo ""
	@echo "─── Vérification Step 2 ────────────────────────────────────"
	@$(PYTHON) -c "\
import pandas as pd; \
df = pd.read_csv('$(SEQS_TSV)', sep='\t', low_memory=False); \
print('  Statuts extraction:'); \
[print(f'    {k:<30} {v:>6,}') for k,v in df['extraction_status'].value_counts().items()]; \
print(); \
n_fasta = sum(1 for l in open('$(FASTA)') if l.startswith('>')); \
print(f'  Séquences dans FASTA   : {n_fasta:,}'); \
direct = df[df['extraction_status']=='DIRECT_SEQUENCE']['sequence_length_extracted'].dropna().astype(float); \
if len(direct): \
    print(f'  DIRECT_SEQUENCE longueur : min={direct.min():.0f}  max={direct.max():.0f}  médiane={direct.median():.0f}'); \
mito = df[df['extraction_status']=='DIRECT_MITO_REF']; \
if len(mito): \
    print(); \
    print(f'  DIRECT_MITO_REF ({len(mito)} variants) :'); \
    [print(f'    AlleleID {row.get(\"#AlleleID\",\"?\")}  range={row.get(\"mito_ref_range\",\"?\")}') for _,row in mito.iterrows()]; \
"
	@echo ""
	@echo "  *** Si des variants attendus sont UNPARSEABLE,"
	@echo "  *** inspecter manuellement : make inspect-unparseable"
	@echo ""
	@echo "────────────────────────────────────────────────────────────"
	@echo "  Step 2 OK. Si OK, lancer : make step4"
	@echo "  (Step 3 = MongoDB legacy, non utilisé — passer directement à step4)"
	@echo ""


# ── Cible utilitaire : inspecter les UNPARSEABLE ──────────────────────────────
.PHONY: inspect-unparseable
inspect-unparseable:
	@$(PYTHON) -c "\
import pandas as pd; \
df = pd.read_csv('$(SEQS_TSV)', sep='\t', low_memory=False); \
unp = df[df['extraction_status']=='UNPARSEABLE']; \
print(f'UNPARSEABLE : {len(unp)} variants'); \
mito_unp = unp[unp['Name'].str.contains('NC_012920|chrM|mito', case=False, na=False)]; \
print(f'  dont référence mitochondriale dans Name : {len(mito_unp)}'); \
if len(mito_unp): \
    print(mito_unp[['#AlleleID','GeneSymbol','Name']].to_string(index=False)); \
print(); \
print('30 premiers noms UNPARSEABLE:'); \
print(unp['Name'].head(30).to_string()); \
" | less


# ── Step 3: Store in MongoDB ──────────────────────────────────────────────────
step3: $(FASTA)
	@echo ""
	@echo "════════════════════════════════════════════════════════════"
	@echo "  Step 3 — Stockage MongoDB"
	@echo "════════════════════════════════════════════════════════════"
	$(PYTHON) $(PIPELINE)/03_mongodb_store.py
	@echo ""
	@echo "  Step 3 OK. Si OK, lancer : make step4"
	@echo ""


# ── Step 4: BLAST vs mitochondrial genome ─────────────────────────────────────
step4: $(MITO_RESULTS)

$(MITO_RESULTS): $(FASTA)
	@echo ""
	@echo "════════════════════════════════════════════════════════════"
	@echo "  Step 4 — BLAST vs génome mitochondrial (rCRS)"
	@echo "════════════════════════════════════════════════════════════"
	$(PYTHON) $(PIPELINE)/04_blast_mito.py
	@echo ""
	@echo "─── Vérification Step 4 ────────────────────────────────────"
	@if [ -s "$(MITO_RESULTS)" ]; then \
		$(PYTHON) -c "\
import pandas as pd; \
cols = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore','qlen','slen']; \
df = pd.read_csv('$(MITO_RESULTS)', sep='\t', header=None, names=cols); \
n = df['qseqid'].nunique(); \
print(f'  Variants avec hit mito : {n}'); \
print(f'  pident : min={df.pident.min():.1f}  max={df.pident.max():.1f}  médiane={df.pident.median():.1f}'); \
print(); \
print('  Hits (triés par pident) :'); \
best = df.sort_values('pident', ascending=False).drop_duplicates('qseqid'); \
print(best[['qseqid','pident','length','evalue']].to_string(index=False)); \
"; \
	else \
		echo "  ATTENTION : aucun hit mito trouvé."; \
		echo "  → Vérifier paramètres BLAST dans 04_blast_mito.py (evalue, perc_identity, word_size)"; \
	fi
	@echo "────────────────────────────────────────────────────────────"
	@echo "  Step 4 OK. Si OK, lancer : make step5"
	@echo ""


# ── Step 5: BLAST vs human genome ─────────────────────────────────────────────
step5: $(NUMT_SUMMARY)

$(NUMT_SUMMARY): $(MITO_RESULTS)
	@echo ""
	@echo "════════════════════════════════════════════════════════════"
	@echo "  Step 5 — BLAST vs génome humain + classification NUMT"
	@echo "════════════════════════════════════════════════════════════"
	$(PYTHON) $(PIPELINE)/05_blast_genome.py
	@echo ""
	@echo "─── Vérification Step 5 ────────────────────────────────────"
	@if [ -s "$(NUMT_SUMMARY)" ]; then \
		$(PYTHON) -c "\
import pandas as pd; \
df = pd.read_csv('$(NUMT_SUMMARY)', sep='\t'); \
print('  Classification NUMT :'); \
[print(f'    {k:<30} {v:>4}') for k,v in df['NUMTClassification'].value_counts().items()]; \
print(); \
novel = df[df['NUMTClassification']=='NOVEL_NUMT_CANDIDATE']; \
known = df[df['NUMTClassification']=='KNOWN_REFERENCE_NUMT']; \
confirmed = df[df['NUMTClassification']=='CONFIRMED_NUMT']; \
print(f'  NUMTs de novo  (NOVEL)     : {len(novel)}'); \
print(f'  NUMTs anciens  (KNOWN_REF) : {len(known)}'); \
print(f'  NUMTs confirmés (CONFIRMED): {len(confirmed)}'); \
print(); \
if len(novel): \
    print('  NOVEL_NUMT_CANDIDATE :'); \
    print(novel[['AlleleID','GeneSymbol','ClinicalSignificance','BestMitoPident']].to_string(index=False)); \
"; \
	else \
		echo "  ATTENTION : fichier summary vide."; \
	fi
	@echo "────────────────────────────────────────────────────────────"
	@echo "  Step 5 OK. Si OK, lancer : make step6"
	@echo ""


# ── Step 6: Mechanism analysis ────────────────────────────────────────────────
step6: $(NUMT_SUMMARY)
	@echo ""
	@echo "════════════════════════════════════════════════════════════"
	@echo "  Step 6 — Analyse du mécanisme d'insertion"
	@echo "════════════════════════════════════════════════════════════"
	$(PYTHON) $(PIPELINE)/06_mechanism.py
	@echo ""
	@echo "  Step 6 terminé."
	@echo ""


# ── PRDM9 scan + DSB-repair mechanism analysis ────────────────────────────────
prdm9: $(FLANKING_DIR)
	@echo ""
	@echo "════════════════════════════════════════════════════════════"
	@echo "  PRDM9 scan + analyse mécanismes DSB-repair"
	@echo "════════════════════════════════════════════════════════════"
	$(PYTHON) $(PIPELINE)/07_prdm9_scan.py --scan-all --numt-analysis --format tsv
	@echo ""
	@echo "  Résultats dans : $(FLANKING_DIR)/"
	@echo ""

$(FLANKING_DIR):
	@echo "  ATTENTION : $(FLANKING_DIR) absent."
	@echo "  Lancer d'abord 'make step6' (ou 'make excel') pour générer les séquences flanquantes."
	@exit 1


# ── Show alignments ASCII ─────────────────────────────────────────────────────
show-alignments:
	@echo ""
	@echo "════════════════════════════════════════════════════════════"
	@echo "  Visualisation alignements BLAST génome (21 candidats)"
	@echo "════════════════════════════════════════════════════════════"
	@if [ ! -d "$(BLAST_GENOME_CACHE)" ]; then \
		echo "  ATTENTION : $(BLAST_GENOME_CACHE) absent."; \
		echo "  Lancer d'abord 'make excel' (avec BLAST génome) pour générer le cache XML."; \
		exit 1; \
	fi
	$(PYTHON) $(PIPELINE)/show_alignments.py
	@echo ""


# ── Compare BLAST vs LAST results ─────────────────────────────────────────────
compare-last:
	@echo ""
	@echo "════════════════════════════════════════════════════════════"
	@echo "  Comparaison BLAST vs LAST (mito)"
	@echo "════════════════════════════════════════════════════════════"
	@if [ ! -s "$(DATA)/excel_blast_mito.tsv" ] || [ ! -s "$(LAST_MITO_TSV)" ]; then \
		echo "  ATTENTION : excel_blast_mito.tsv ou excel_last_mito.tsv manquant."; \
		echo "  Lancer 'make excel-validate-last' pour générer les deux."; \
		exit 1; \
	fi
	$(PYTHON) $(PIPELINE)/last_validation/compare_blast_last.py \
		--blast $(DATA)/excel_blast_mito.tsv \
		--last  $(LAST_MITO_TSV) \
		--output-tsv    $(BLAST_LAST_TSV) \
		--output-report $(BLAST_LAST_RPT)
	@echo ""
	@echo "  TSV    : $(BLAST_LAST_TSV)"
	@echo "  Rapport : $(BLAST_LAST_RPT)"
	@echo ""


# ── Status : état des fichiers intermédiaires ─────────────────────────────────
status:
	@echo ""
	@echo "════════ État du pipeline NUMT ════════"
	@$(PYTHON) -c "\
from pathlib import Path; \
files = { \
    'Step 1 — ClinVar TSV':         '$(CLINVAR_TSV)', \
    'Step 2 — Séquences TSV':        '$(SEQS_TSV)', \
    'Step 2 — FASTA':                '$(FASTA)', \
    'Step 4 — BLAST mito results':   '$(MITO_RESULTS)', \
    'Step 5 — NUMT summary':         '$(NUMT_SUMMARY)', \
    'Steps 1-8 — Excel output':      '$(EXCEL_OUTPUT)', \
    'Step 7 — Mechanism JSON':       '$(MECHANISM_JSON)', \
    'Step 7 — Mechanism TSV':        '$(MECHANISM_TSV)', \
    'Validation — LAST mito TSV':    '$(LAST_MITO_TSV)', \
    'Validation — BLAST vs LAST':    '$(BLAST_LAST_TSV)', \
}; \
for label, path in files.items(): \
    p = Path(path); \
    if p.exists(): \
        size = p.stat().st_size; \
        lines = sum(1 for _ in open(p, errors='ignore')) if size < 50_000_000 else '?'; \
        print(f'  [OK] {label:<35} {lines:>8} lignes  ({size/1024:.0f} Ko)'); \
    else: \
        print(f'  [--] {label:<35} (absent)'); \
"
	@echo ""


# ── Clean ─────────────────────────────────────────────────────────────────────
clean:
	@echo "Suppression des fichiers intermédiaires..."
	rm -f $(CLINVAR_TSV) $(SEQS_TSV) $(FASTA) $(MITO_RESULTS) $(NUMT_SUMMARY)
	rm -f $(EXCEL_OUTPUT)
	rm -f $(MECHANISM_JSON) $(MECHANISM_TSV)
	rm -f $(LAST_MITO_TSV) $(BLAST_LAST_TSV) $(BLAST_LAST_RPT)
	rm -f $(DATA)/blast_genome_query.fasta $(DATA)/excel_blast_mito.tsv $(DATA)/excel_blast_genome.tsv
	rm -rf $(DATA)/blast_genome_results_raw $(DATA)/last_work
	rm -f $(DATA)/.pipeline_checkpoint.json
	@echo "Nettoyage terminé."
