#!/usr/bin/env bash
# =============================================================================
#  18S rRNA V4-V5 reference database (~574 bp) – marine invertebrates and other eukaryotes – CRABS
# =============================================================================
#  Primers   : forward 18Smk1F CTGGTGCCAGCAGCCGCGGYAA (Machida & Knowlton 2012)
#              reverse 18Smk2R TCCGTCAATTYCTTTAAGTT (Machida & Knowlton 2012)
#  Sources   : SILVA SSU NR99 (all eukaryotes) + NCBI GenBank (marine invertebrates, algae) + BOLD
#  Output    : results/18S_V4V5_marine_eukaryotes.fasta (OBITools-style headers)
#              results/18S_V4V5_marine_eukaryotes.tsv   (table with column names)
#  Provenance (FAIR): metadata/ (parameters, software versions, sources, counts,
#              checksums, database README), logs/ (one log per step)
#
#  Folder    : next to crabs_coi, e.g. ~/…/CRABS/crabs_18s
#              (the NCBI taxonomy, the CRABS installation and the BOLD files of
#              crabs_coi are reused: nothing is downloaded twice)
#  Usage     : mkdir -p ~/…/CRABS/crabs_18s && cd ~/…/CRABS/crabs_18s
#              caffeinate -i bash crabs_18S.sh 2>&1 | tee -a crabs_18S.log
#
#  Resuming  : the script can be re-run at any time. Every step already done
#  (output present, inputs unchanged, same parameters) is skipped. If an input
#  or a parameter changes, that step and all downstream steps are re-run.
# =============================================================================
set -uo pipefail
# up NAME: first NAME found in this folder or in a parent folder (absolute path),
# so that the shared CRABS/taxo, CRABS/reference_database_creator and crabs_coi
# are found whatever the depth of the database folder (e.g. crabs_18S/18S_aligned)
up() { local d=$PWD; while :; do [[ -e "$d/$1" ]] && { echo "$d/$1"; return 0; }; [[ "$d" == / ]] && return 0; d=$(dirname "$d"); done; }   # nothing printed if not found
SCRIPT_VERSION="1.5.1 (2026-10-28, shared CRABS/taxo and CRABS/reference_database_creator; 1.5: from crabs_coi.sh 1.8 + taxonomy fix of 1.7.2, SILVA lineage check)"

# ================================ PARAMETERS ==================================
FWD="CTGGTGCCAGCAGCCGCGGYAA"
REV="TCCGTCAATTYCTTTAAGTT"
FWD_NAME="18Smk1F, Machida & Knowlton 2012"   # primer names / references (metadata only)
REV_NAME="18Smk2R, Machida & Knowlton 2012"
EMAIL="morgane.durand@mnhn.fr"         # required by NCBI for downloads
BASE_NAME="18S_V4V5_marine_eukaryotes"
TITLE="18S rRNA V4-V5 reference database (~574 bp) – marine invertebrates and other eukaryotes"
MARKER_DESC="18S rRNA (nuclear small subunit), V4-V5 region, ~574 bp"
# References of the primers, added to "Please cite" in metadata/README.md
PRIMER_CITES=("Machida & Knowlton (2012) PCR primers for metazoan nuclear 18S and 28S ribosomal DNA sequences. PLoS ONE 7: e46180. doi:10.1371/journal.pone.0046180")

# In silico PCR
MISMATCH=4           # mismatches allowed per primer (CRABS default)
RELAXED=false        # both primer sites sit in conserved regions of full-length 18S
# Pairwise global alignment (recovers amplicons deposited without primer sites).
# The slowest step (~15 h for COI). USE_PGA=false skips it: the in silico PCR
# amplicons are used as they are (for COI it had recovered only 0.9 % of them).
USE_PGA=true
PGA_IDENTITY=0.95
PGA_COVERAGE=95
PGA_SIZE_SELECT=10000
# Final filters (1.5). Amplicons with both primer sites, by phylum (1.4 run):
# median ~550 bp for most groups, but longer in Platyhelminthes (median 639,
# p90 707), Cnidaria (p90 683), Porifera and Tardigrada (p90 640), Apicomplexa
# (p90 657), Arthropoda (p90 836), and shorter in Ciliophora (p10 475). The
# window 500-680 bp of 1.1-1.4 removed 7,091 amplicons as too long (1,783
# arthropods, 1,195 flatworms, 709 cnidarians at 681-800 bp): widened to
# 450-900 bp. Over-long amplicons from a wrong primer site: TRIM_LONG below.
EXPECTED_LEN=554
MIN_LEN=450
MAX_LEN=900
# Over-long amplicons (1.5, as crabs_16S.sh 1.4). In a long record (genome,
# rDNA operon), cutadapt cuts at the FIRST acceptable primer site from the
# start, not the best one, so the "amplicon" can run over a wrong site and be
# removed by MAX_LEN. TRIM_LONG=true: in an amplicon longer than MAX_LEN, the
# primer sites found INSIDE it (up to MISMATCH mismatches) are used to re-cut
# it, keeping the cut with the fewest mismatches within MIN_LEN-MAX_LEN.
# Changes: logs/09_long_trimmed.tsv.
TRIM_LONG=true
MAX_N=2              # CRABS removes amplicons with MAX_N or MORE "N": 2 = at most 1 N allowed
# Short references: an amplicon between MIN_LEN_RESCUE and MIN_LEN is kept ONLY
# for species that have no amplicon of full length (MIN_LEN-MAX_LEN), e.g.
# barcodes deposited slightly shorter than the fragment. "" = off.
# Translation check (protein-coding markers): an amplicon is kept only if one of
# its three reading frames has no stop codon. Removes pseudogenes (NUMTs),
# sequencing errors and frame-shifted records, e.g. among "UNVERIFIED" ones.
# rRNA marker: not protein-coding, check off.
TRANSLATION_CHECK=false
STOP_CODONS=""
MIN_LEN_RESCUE=""
# Lowest rank a reference must be identified to: species | genus | family.
# 18S often does not separate congeneric species: genus is a sensible minimum.
MIN_ID_RANK=genus
MAX_MISSING_RANKS=1  # missing ranks tolerated ABOVE that rank (gaps in the NCBI lineage)
# Phylum names harmonised with Adl et al. (2019, J Eukaryot Microbiol 66: 4-119),
# the reference classification of eukaryotes (fixed table, no online lookup).
#  - PHYLUM_RENAME : phylum names replaced ("old=new;…"). Diatoms (Bacillariophyta
#    in NCBI) belong to Ochrophyta; dinoflagellates are Dinoflagellata.
#  - CLASS_TO_PHYLUM: phylum given from the class when NCBI has no phylum rank
#    (brown algae and other Ochrophyta, dinoflagellates: "Unassigned" otherwise).
# Animal phyla and Rhodophyta are kept as named in NCBI.
PHYLUM_RENAME="Bacillariophyta=Ochrophyta;Heterokontophyta=Ochrophyta;Myzozoa=Dinoflagellata;Dinophyta=Dinoflagellata"
CLASS_TO_PHYLUM="Phaeophyceae=Ochrophyta;Bacillariophyceae=Ochrophyta;Coscinodiscophyceae=Ochrophyta;Mediophyceae=Ochrophyta;Chrysophyceae=Ochrophyta;Synurophyceae=Ochrophyta;Raphidophyceae=Ochrophyta;Eustigmatophyceae=Ochrophyta;Phaeosacciophyceae=Ochrophyta;Pelagophyceae=Ochrophyta;Dictyochophyceae=Ochrophyta;Xanthophyceae=Ochrophyta;Bolidophyceae=Ochrophyta;Pinguiophyceae=Ochrophyta;Phaeothamniophyceae=Ochrophyta;Aurearenophyceae=Ochrophyta;Chrysomerophyceae=Ochrophyta;Schizocladiophyceae=Ochrophyta;Olisthodiscophyceae=Ochrophyta;Dinophyceae=Dinoflagellata"
DEREP_METHOD=unique_species
RANKS='superkingdom;phylum;class;order;family;genus;species'   # switched to "domain" automatically if needed

# Additional export formats (optional), e.g. (sintax blast-tax)
EXPORT_FORMATS=()

# Taxa removed from the final database (any rank of the NCBI lineage). Terrestrial
# arthropods are not targeted but some reach the database anyway (NCBI files
# insects under Pancrustacea next to Crustacea; misfiled or contaminant records).
# SILVA SSU also holds bacterial and archaeal 16S (and organelle rRNA filed under
# Bacteria): removed too.
EXCLUDE_TAXA="Bacteria;Archaea;Insecta;Collembola;Protura;Diplura;Arachnida;Chilopoda;Diplopoda;Pauropoda;Symphyla;Onychophora"
# Excluded taxa are removed as early as possible, so they do not go through the
# long steps: in the NCBI queries (NCBI_EXCLUDE_AT_DOWNLOAD=true adds
# NOT ("<taxon>"[Organism] OR …)), right after import (step 04c) and once more
# after the final filters, as a safety net
NCBI_EXCLUDE_AT_DOWNLOAD=true

# Genus homonyms. CRABS gives a sequence its NCBI lineage from its NAME when the
# record has no usable taxid (BOLD records, names unknown to NCBI). When the
# genus name also exists in another kingdom or phylum, the wrong genus may be
# picked: e.g. Fridericia (oligochaete) filed under a plant, Siphonaria
# (gastropod) under a chytrid fungus, Cephalothrix (nemertean) under a
# cyanobacterium, Isotoma (springtail) under a plant (and so missed by
# EXCLUDE_TAXA). Step 10_taxonomy_fix, BEFORE the exclusion of taxa, checks every
# sequence: (1) BOLD records: the phylum given by BOLD (original BOLD files) must
# match the phylum of the NCBI lineage; (2) if TARGET_CHECK=true, every record:
# the NCBI lineage must reach one of the target groups (TAXA_NCBI, TAXA_BOLD,
# TARGET_EXTRA). A wrong lineage is replaced by that of the genus of the same
# name that fits (BOLD phylum, else target groups), when there is exactly one;
# otherwise the sequence is removed. Every change: logs/10_taxonomy_fix.tsv.
# false = off.
TAXONOMY_FIX=true
TARGET_CHECK=false          # SILVA brings all eukaryotes on purpose: no target check
# 18S: the phylum of the SILVA lineage (SILVA headers) is checked like the BOLD
# phylum. CRABS finds the NCBI taxid of a SILVA record from its GenBank
# accession, else from its NAME (one of the homonymous genera, at random): the
# SILVA lineage tells which genus is right.
# Organelle rRNA (1.5): SILVA SSU also holds chloroplast and mitochondrion 16S,
# classified by SILVA in Bacteria (e.g. "Bacteria;Cyanobacteria;…;Chloroplast")
# but filed by NCBI under the plant or alga: 2,855 such amplicons (~370 bp, in
# Streptophyta, Rhodophyta, Chlorophyta, Euglenozoa) in the 1.4 run. Removed by
# the taxonomy fix (logged as "organelle rRNA").
# SILVA classifies each sequence from the sequence itself (phylogenetic
# placement), NCBI files it under the organism declared in the record, often the
# host. When the SILVA (or BOLD) phylum differs and no genus of that name exists
# in it, the sequence most likely comes from another organism (an epibiont,
# parasite, gut content or contaminant: e.g. ciliate 18S filed under Gammarus,
# flatworm 18S under Daphnia): such a reference would assign reads to the wrong
# organism. UNRESOLVED_MISMATCH=remove removes them (keep = kept as in 1.3).
# Exception: both phyla in the same kingdom listed in SAME_KINGDOM_OK (fungal
# phyla are delimited differently by NCBI and SILVA, e.g. Basidiobolomycota vs
# Zoopagomycota): not counted as a disagreement.
UNRESOLVED_MISMATCH=remove
SAME_KINGDOM_OK="Fungi"
# Target groups not named in the download lists: algae under their NCBI names,
# fish lice (Branchiura) and tongue worms (Pentastomida) filed by NCBI under
# Ichthyostraca, and Mystacocarida (all downloaded through "Crustacea" in NCBI)
TARGET_EXTRA="Ochrophyta;Bacillariophyta;Phaeophyceae;Ichthyostraca;Mystacocarida"
# Manual decisions, applied first (traceable, cited in metadata/parameters.tsv):
# a tab-separated file in this folder, columns taxon / action / reason, where
# taxon = a genus (applies only to its FLAGGED sequences) or a seq_id, and
# action = keep | remove. Lines starting with # are comments. Absent = none.
TAXONOMY_MANUAL="taxonomy_manual_crabs.tsv"

# Target groups for NCBI and BOLD (SILVA brings all the other eukaryotes:
# protists, fungi, vertebrates…, so that a read of a non-target organism is not
# forced onto an invertebrate). BOLD and NCBI do not always use the same names.
TAXA_NCBI=(Porifera Cnidaria Ctenophora Placozoa Xenacoelomorpha Platyhelminthes
           Nemertea Mollusca Annelida Sipuncula Brachiopoda Bryozoa Phoronida
           Entoprocta Chaetognatha Echinodermata Hemichordata Tunicata
           Cephalochordata Nematoda Gastrotricha Kinorhyncha Priapulida
           Tardigrada Crustacea Pycnogonida Xiphosura)
TAXA_BOLD=(Porifera Cnidaria Ctenophora Placozoa Xenacoelomorpha Platyhelminthes
           Nemertea Mollusca Annelida Sipuncula Brachiopoda Bryozoa Phoronida
           Entoprocta Chaetognatha Echinodermata Hemichordata
           Ascidiacea Thaliacea Appendicularia Leptocardii Nematoda Gastrotricha
           Kinorhyncha Priapulida Tardigrada
           Malacostraca Copepoda Ostracoda Thecostraca Branchiopoda Cephalocarida
           Remipedia Pycnogonida Merostomata)

# Marine algae (also amplified by these universal eukaryote primers)
INCLUDE_ALGAE=true
ALGAE_NCBI=(Rhodophyta Phaeophyceae Bacillariophyta)
ALGAE_BOLD=(Rhodophyta Ochrophyta)
if [[ "$INCLUDE_ALGAE" == true ]]; then
  TAXA_NCBI+=("${ALGAE_NCBI[@]}")
  TAXA_BOLD+=("${ALGAE_BOLD[@]}")
fi
# Nuclear 18S: complete rDNA operons are 5-10 kb
NCBI_GENE_QUERY='("18S ribosomal RNA"[All Fields] OR "18S rRNA"[All Fields] OR "small subunit ribosomal RNA"[All Fields]) NOT mitochondrion[filter] NOT chloroplast[filter] AND 300:10000[SLEN]'
NCBI_EXTRA_QUERY=""   # nuclear marker: no unannotated organelle records to add
BOLD_MARKER="18S"    # BOLD marker_code kept (regular expression); "" = all markers

# SILVA SSU Ref NR99 ("subset"; "full" = SSU Ref, much bigger). Check the latest
# release on https://www.arb-silva.de/download/arb-files/
USE_SILVA=true
SILVA_GENE=SSU
SILVA_TYPE=subset
SILVA_VERSION=138.2

USE_BOLD=true
USE_NCBI=true
THREADS=8
MAX_TRIES=4          # download attempts per group
PAUSE=60             # wait (s) before the 2nd attempt, doubled afterwards
BOLD_PAUSE=10        # wait (s) between two BOLD groups, to spare the server
# NOTE ON BOLD DOWNLOAD LIMITS
# The BOLD portal limits how much data a user (IP address) can download in a
# given period. Once the limit is reached, BOLD answers "503 Service Temporarily
# Unavailable" or "403 Forbidden" (even the home page), sometimes saved by CRABS
# as an HTML page instead of data. This is not a bug in CRABS or in this script:
# wait a few hours to a day, then simply re-run the script. Groups already
# downloaded are kept; only the missing ones are retried.
# BOLD V5 downloads contain ALL markers of a taxon: the files already downloaded
# for the COI database (BOLD_REUSE_DIR/raw/bold_<group>.tsv) are reused, and only
# the records of BOLD_MARKER are kept.
BOLD_REUSE_DIR="$(up crabs_coi)"     # COI database folder (searched in the parent folders)
# NCBI taxonomy: the one of the COI database is used when present, otherwise it
# is downloaded into TAXO
TAXO_REUSE="$(up taxo)"; [[ -s "$TAXO_REUSE/nodes.dmp" ]] || TAXO_REUSE="$(up crabs_coi)/taxo"   # shared CRABS/taxo (older layout: crabs_coi/taxo)
TAXO="taxo"

# =============================================================================
#                      NOTHING TO EDIT BELOW THIS LINE
# =============================================================================

SCRIPT_NAME=$(basename "$0")
LOG_NAME="${SCRIPT_NAME%.sh}.log"
# NCBI taxonomy shared with the COI database when it is already there
[[ -n "${TAXO_REUSE:-}" && -s "$TAXO_REUSE/nodes.dmp" && -s "$TAXO_REUSE/nucl_gb.accession2taxid" ]] && TAXO="$TAXO_REUSE"
RUN_START=$(date '+%Y-%m-%dT%H:%M:%S%z')
mkdir -p raw imported steps results logs metadata metadata/.steps "$TAXO"

log()  { echo "[$(date '+%H:%M:%S')] $*"; }
warn() { echo -e "\n>>> WARNING: $*\n" >&2; }
die()  { echo -e "\n✘ $*" >&2; exit 1; }
# CRABS exits without an error code when a step fails, so we always check that
# the expected output file was actually created
need() { [[ -s "$1" ]] || die "Step failed: $1 was not created (see the CRABS message above and logs/)"; }
fdate() { python3 -c 'import os,sys,datetime;print(datetime.datetime.fromtimestamp(os.path.getmtime(sys.argv[1])).strftime("%Y-%m-%d"))' "$1"; }
nlines() { [[ -s "$1" ]] && wc -l < "$1" | tr -d ' ' || echo 0; }

# ---------------------------------------------------------------------------
# Resume engine: a step is skipped if its output exists, is newer than all of
# its inputs, and was produced with the same "signature" (parameters + inputs).
# An existing output without a signature (made by an earlier script version)
# is adopted as is when the step is called with ADOPT=1.
# ---------------------------------------------------------------------------
sigfile() { echo "metadata/.steps/$(echo "$1" | tr '/' '_').sig"; }
# Signatures are compared without the folders of the files they name (only the
# file names), so that moving the folders does not make finished steps run again
SIGNORM='s#[^;|=]*/##g'
fresh() {   # fresh OUTPUT SIGNATURE ADOPT INPUTS...
  local out=$1 sig=$2 adopt=$3; shift 3
  [[ -s "$out" ]] || return 1
  local sf; sf=$(sigfile "$out")
  if [[ ! -e "$sf" ]]; then
    [[ "$adopt" == 1 ]] || return 1
    printf '%s\n' "$sig" > "$sf"; log "  (existing output adopted: $out)"; return 0
  fi
  [[ "$(sed "$SIGNORM" "$sf")" == "$(printf '%s\n' "$sig" | sed "$SIGNORM")" ]] || return 1
  local i; for i in "$@"; do [[ -e "$i" && "$i" -nt "$out" ]] && return 1; done
  return 0
}
step() {    # step NAME OUTPUT SIGNATURE ADOPT INPUTS... -- COMMAND...
  local name=$1 out=$2 sig=$3 adopt=$4; shift 4
  local inputs=()
  while [[ $# -gt 0 && "$1" != "--" ]]; do inputs+=("$1"); shift; done
  shift
  if fresh "$out" "$sig" "$adopt" ${inputs[@]+"${inputs[@]}"}; then
    log "✔ $name: already done → $out"; return 0
  fi
  log "▶ $name"
  rm -f "$(sigfile "$out")"
  "$@" 2>&1 | tee "logs/${name}.log"
  need "$out"
  printf '%s\n' "$sig" > "$(sigfile "$out")"
}

# -------------------------------- Tools ---------------------------------------
CRABS_FOUND=""
# CRABS installation: CRABS/reference_database_creator (shared by all markers);
# (searched in this folder and its parents; older layout: crabs_coi/reference_database_creator;
# or set CRABS_DIR=/path)
for D in "${CRABS_DIR:-}" "$(up reference_database_creator)" "$(up crabs_coi)/reference_database_creator"; do
  if [[ -n "$D" && -f "$D/crabs" ]]; then CRABS_FOUND="$D"; break; fi
done
[[ -n "$CRABS_FOUND" ]] || die "CRABS not found: put reference_database_creator in the CRABS folder (or set CRABS_DIR=/path/to/reference_database_creator)"
chmod +x "$CRABS_FOUND/crabs"
export PATH="$CRABS_FOUND:$PATH"
[[ -n "${CRABS_VENV:-}" && -d "$CRABS_VENV/bin" ]] && export PATH="$CRABS_FOUND:$CRABS_VENV/bin:$PATH"   # optional Python environment
hash -r
# CRABS declares the option --exclude twice (once for --download-taxonomy, once
# for --subset), so click prints a harmless "UserWarning: The parameter --exclude
# is used more than once" at every CRABS call. Only that warning is silenced.
export PYTHONWARNINGS="${PYTHONWARNINGS:+$PYTHONWARNINGS,}ignore:The parameter --exclude is used more than once:UserWarning"
crabs -h >/dev/null 2>&1 || { crabs -h; die "crabs -h fails"; }
for tool in cutadapt vsearch makeblastdb python3; do
  command -v "$tool" >/dev/null 2>&1 || die "Missing tool: $tool"
done
log "CRABS : $(command -v crabs)"

# Software versions (FAIR: reproducibility)
CRABS_VERSION=$(python3 -c "import re,sys;print(re.search(r'[\"\']([0-9][^\"\']*)', open(sys.argv[1]).read()).group(1))" "$CRABS_FOUND/function/__init__.py" 2>/dev/null || echo "unknown")
CRABS_COMMIT=$(git -C "$CRABS_FOUND" rev-parse --short HEAD 2>/dev/null || echo "unknown")
{
  printf 'software\tversion\n'
  printf 'script\t%s\n' "$SCRIPT_VERSION"
  printf 'CRABS\t%s (commit %s)\n' "$CRABS_VERSION" "$CRABS_COMMIT"
  printf 'cutadapt\t%s\n' "$(cutadapt --version 2>&1 | head -n 1)"
  printf 'vsearch\t%s\n' "$(vsearch --version 2>&1 | head -n 1 | cut -d, -f1)"
  printf 'makeblastdb\t%s\n' "$(makeblastdb -version 2>&1 | head -n 1)"
  printf 'python\t%s\n' "$(python3 --version 2>&1)"
  printf 'system\t%s\n' "$(uname -srm)"
} > metadata/versions.tsv

# Validity checks for downloaded files (an overloaded server sometimes returns
# an HTML error page that CRABS saves as if it were data)
bold_ok()  { [[ -s "$1" ]] && head -n 1 "$1" | awk -F'\t' '{p=s=0; for(i=1;i<=NF;i++){c=tolower($i); if(c=="processid")p=1; if(c=="nuc"||c=="nucleotides"||c=="sequence")s=1} exit !(p&&s)}'; }
fasta_ok() { [[ -s "$1" ]] && [[ "$(head -c 1 "$1")" == ">" ]]; }

: > metadata/.sources.tmp

# =============================================================================
# 1. NCBI taxonomy
# =============================================================================
if [[ -s "$TAXO/nodes.dmp" && -s "$TAXO/names.dmp" && -s "$TAXO/nucl_gb.accession2taxid" ]]; then
  log "✔ NCBI taxonomy: already downloaded in $TAXO ($(fdate "$TAXO/nodes.dmp"))"
else
  log "▶ Downloading NCBI taxonomy"
  crabs --download-taxonomy --output "$TAXO" 2>&1 | tee logs/00_taxonomy.log
  need "$TAXO/nodes.dmp"; need "$TAXO/nucl_gb.accession2taxid"
fi
printf 'NCBI Taxonomy\ttaxdump + nucl_gb.accession2taxid\t%s\t%s\t-\t-\n' "$TAXO" "$(fdate "$TAXO/nodes.dmp")" >> metadata/.sources.tmp
# Since 2025 the NCBI taxonomy names the top rank "domain" instead of
# "superkingdom". With the old name, that rank comes out empty ("NA") for every
# sequence, which also makes --rank-na remove sequences missing just one more rank.
if [[ "$RANKS" == superkingdom\;* ]] && ! grep -q -m1 $'\t|\tsuperkingdom\t|' "$TAXO/nodes.dmp" \
   && grep -q -m1 $'\t|\tdomain\t|' "$TAXO/nodes.dmp"; then
  RANKS="domain;${RANKS#superkingdom;}"
  log "  NCBI taxonomy uses 'domain' instead of 'superkingdom': ranks = $RANKS"
fi
# Filters derived from MIN_ID_RANK: CRABS --rank-na counts ALL missing ranks, so
# the ranks allowed to be empty below MIN_ID_RANK are added to the tolerance
IFS=';' read -r -a RANK_LIST <<< "$RANKS"
MIN_RANK_POS=0
for i in "${!RANK_LIST[@]}"; do [[ "${RANK_LIST[$i]}" == "$MIN_ID_RANK" ]] && MIN_RANK_POS=$((i + 1)); done
(( MIN_RANK_POS > 0 )) || die "MIN_ID_RANK=$MIN_ID_RANK is not one of the ranks: $RANKS"
BELOW=$(( ${#RANK_LIST[@]} - MIN_RANK_POS ))
RANK_NA=$(( MAX_MISSING_RANKS + BELOW + 1 ))       # CRABS removes rows with RANK_NA or more "NA"
MIN_RANK_COL=$(( 3 + MIN_RANK_POS ))                # column of that rank in the CRABS format
SPECIES_ARG=(); [[ "$MIN_ID_RANK" == species ]] && SPECIES_ARG=(--no-species-id)
log "  references identified at least to: $MIN_ID_RANK (--rank-na $RANK_NA)"
TAXARGS=(--names "$TAXO/names.dmp" --nodes "$TAXO/nodes.dmp"
         --acc2tax "$TAXO/nucl_gb.accession2taxid" --ranks "$RANKS")

IMPORTED=()

# =============================================================================
# 2a. BOLD: download per group, merge, import
# =============================================================================
# BOLD restricts the amount of data that can be downloaded over a given period
# (see NOTE ON BOLD DOWNLOAD LIMITS in the parameters). After repeated failures
# the script stops querying BOLD, tells the user to wait, and carries on with
# the other sources.
bold_limit_msg() {
  echo
  echo "=============================================================================="
  echo " BOLD DOWNLOAD LIMIT REACHED"
  echo "------------------------------------------------------------------------------"
  echo " BOLD restricts the amount of data a user can download over a given period."
  echo " The server is now refusing requests (403/503) after several attempts."
  echo " This is NOT an error in CRABS or in this script."
  echo
  echo " → Wait a few hours to a day, then re-run the same command:"
  echo "     caffeinate -i bash $SCRIPT_NAME 2>&1 | tee -a $LOG_NAME"
  echo "   Groups already downloaded are kept; only the missing ones are retried."
  echo
  echo " → To check whether BOLD accepts requests again (200 = OK, 403 = still blocked):"
  echo "     curl -s -o /dev/null -w \"%{http_code}\\n\" \"https://portal.boldsystems.org/\""
  echo
  echo " Meanwhile the database is built with the data already available."
  echo "=============================================================================="
  echo
}
if [[ "$USE_BOLD" == true ]]; then
  BOLD_FILES=(); BOLD_FAILED=(); BOLD_BLOCKED=false
  for T in "${TAXA_BOLD[@]}"; do
    F="raw/bold_${T}.tsv"
    if [[ -e "$F" ]] && ! bold_ok "$F"; then
      warn "BOLD: $F invalid (server error page?), downloading again"; rm -f "$F"
    fi
    if bold_ok "$F"; then BOLD_FILES+=("$F"); continue; fi
    # BOLD V5 gives ALL markers of a taxon: a group already downloaded for the COI
    # database is reused as is (no new download, no BOLD quota used)
    if [[ -n "${BOLD_REUSE_DIR:-}" ]] && bold_ok "$BOLD_REUSE_DIR/raw/bold_${T}.tsv"; then
      log "✔ BOLD: $T reused from $BOLD_REUSE_DIR"; BOLD_FILES+=("$BOLD_REUSE_DIR/raw/bold_${T}.tsv"); continue
    fi
    if [[ "$BOLD_BLOCKED" == true ]]; then BOLD_FAILED+=("$T"); continue; fi
    log "▶ BOLD: $T"
    TRY=1; WAIT=$PAUSE
    while ! bold_ok "$F" && (( TRY <= MAX_TRIES )); do
      (( TRY > 1 )) && { log "  BOLD $T: retry $TRY/$MAX_TRIES in ${WAIT}s"; sleep "$WAIT"; WAIT=$(( WAIT * 2 )); }
      crabs --download-bold --taxon "$T" --output "$F" 2>&1 | tee "logs/01_bold_${T}.log" || true
      bold_ok "$F" || rm -f "$F"
      if grep -q "not found in BOLD" "logs/01_bold_${T}.log"; then warn "BOLD: '$T' unknown to BOLD, skipped"; break; fi
      if grep -q "(403)" "logs/01_bold_${T}.log"; then
        warn "BOLD refuses requests (403): download limit reached, stopping BOLD downloads for this run"
        BOLD_BLOCKED=true; break
      fi
      TRY=$(( TRY + 1 )); sleep "$BOLD_PAUSE"
    done
    if bold_ok "$F"; then BOLD_FILES+=("$F"); else BOLD_FAILED+=("$T"); fi
  done
  (( ${#BOLD_FAILED[@]} > 0 )) && warn "BOLD: missing groups: ${BOLD_FAILED[*]} (re-run the script later to retry them)"
  # Blocked (403), or at least one group failed after all attempts with server
  # errors: most likely the BOLD download limit → tell the user to wait
  if [[ "$BOLD_BLOCKED" == true ]] || grep -lqE "\((403|429|503)\)|503 Service|Temporarily Unavailable" logs/01_bold_*.log 2>/dev/null \
       && (( ${#BOLD_FAILED[@]} > 0 )); then
    bold_limit_msg
  fi
  for F in ${BOLD_FILES[@]+"${BOLD_FILES[@]}"}; do
    printf 'BOLD Systems V5\t%s\t%s\t%s\t%s\ttax:%s\n' "$(basename "$F" .tsv | sed 's/^bold_//')" "$F" "$(fdate "$F")" \
      "$(grep -vc '^processid' "$F" || true)" "$(basename "$F" .tsv | sed 's/^bold_//')" >> metadata/.sources.tmp
  done

  if (( ${#BOLD_FILES[@]} > 0 )); then
    # A merged file made by an older version (original BOLD header instead of
    # the 68-column layout) is discarded
    if [[ -s raw/bold_all.tsv ]] && ! head -n 1 raw/bold_all.tsv | awk -F'\t' '{exit !(NF==68 && $67=="nuc")}'; then
      log "  Old BOLD merged file discarded"; rm -f raw/bold_all.tsv imported/bold_all.txt "$(sigfile imported/bold_all.txt)"
    fi
    SIG_BOLDMERGE="bold-merge-v6|marker=$BOLD_MARKER|$(printf '%s;' "${BOLD_FILES[@]}")"
    step 02_merge_bold raw/bold_all.tsv "$SIG_BOLDMERGE" 0 "${BOLD_FILES[@]}" -- \
      python3 - raw/bold_all.tsv "$BOLD_MARKER" "${BOLD_FILES[@]}" <<'PYEOF'
import re, sys
# Merges BOLD V5 TSV files into the layout expected by the CRABS importer, which
# reads columns by POSITION (2nd = identifier, 25th = taxon, 67th = sequence).
# BOLD downloads come in pages of 10,000 records, and a page may carry fewer or
# more (optional) columns than the header announces, WITHOUT a new header line.
# Columns are therefore located by NAME in the header, and each row is checked:
#  - same width as the header -> read as is;
#  - other width -> the sequence column is found by its content (a long DNA
#    string) and the identifier / taxon columns are shifted by the same offset;
#    the row is kept only if every value found is plausible.
#  - BOLD_MARKER (regular expression, case-insensitive) keeps only the records
#    whose marker_code matches (e.g. 18S), so the other markers are not imported.
out_path, marker, files = sys.argv[1], sys.argv[2], sys.argv[3:]
MARK = re.compile(marker, re.I) if marker else None
MARKER_COLS = ["marker_code", "markercode", "marker"]
ID_COLS  = ["record_id", "processid", "sampleid"]   # record_id: one per specimen AND per marker
SEQ_COLS = ["nuc", "nucleotides", "sequence"]
TAX_COLS = ["species", "identification", "subspecies", "genus", "subfamily", "family", "order", "class", "phylum"]
NCOL = 68   # the sequence (67th) must not be the last column
DNA = re.compile(r"^[ACGTURYKMSWBDHVN-]{50,}$")
NAME = re.compile(r"^[A-Z][a-z]+( [A-Za-z0-9.\-]+){0,4}$")
def find(h, names):
    low = [x.strip().lower() for x in h]
    for n in names:
        if n in low: return low.index(n)
    return None
def is_header(p):
    low = {x.strip().lower() for x in p}
    return "processid" in low and bool(low & set(SEQ_COLS))
def pick(parts, i_id, i_seq, i_tax, shift, i_mk=None):
    j = i_seq + shift
    if not (0 <= j < len(parts)): return None
    seq = parts[j].strip().upper()
    if not DNA.match(seq): return None
    # the extra/missing column may sit before or after a given column: try the
    # sequence's offset first, then no offset, and keep the first plausible value
    offs = [shift, 0] if shift else [0]
    def get(i, ok):
        for o in offs:
            if 0 <= i + o < len(parts) and ok(parts[i + o].strip()): return parts[i + o].strip()
        return None
    pid = get(i_id, lambda v: bool(v) and " " not in v and ("-" in v or "." in v)) if i_id is not None else None
    if not pid: return None
    sp = ""
    for i in i_tax:
        v = get(i, lambda v: bool(NAME.match(v)))
        if v: sp = v; break
    mk = get(i_mk, lambda v: bool(v) and " " not in v and len(v) < 40) if i_mk is not None else ""
    return pid, sp, seq.replace("-", "N"), mk or ""
total = 0
with open(out_path, "w", encoding="utf-8") as out:
    h = [""] * NCOL; h[1], h[24], h[66] = "processid", "species", "nuc"
    out.write("\t".join(h) + "\n")
    for path in files:
        kept = shifted = bad = noseq = other = 0; header = None
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for raw in f:
                line = raw.rstrip("\r\n"); parts = line.split("\t")
                if is_header(parts):
                    header = parts; n = len(header)
                    i_id, i_seq = find(header, ID_COLS), find(header, SEQ_COLS)
                    i_tax = [i for i in (find(header, [c]) for c in TAX_COLS) if i is not None]
                    i_mk = find(header, MARKER_COLS)
                    if MARK and i_mk is None:
                        print(f"   ⚠ {path}: no marker column, all markers kept", file=sys.stderr)
                    continue
                if header is None or not line.strip() or i_seq is None: continue
                res = pick(parts, i_id, i_seq, i_tax, 0, i_mk)
                if res is None and len(parts) != n:
                    # locate the DNA column and try the offset it implies
                    for j, v in enumerate(parts):
                        if len(v) >= 50 and DNA.match(v.strip().upper()):
                            res = pick(parts, i_id, i_seq, i_tax, j - i_seq, i_mk)
                            if res: shifted += 1
                            break
                if res is None:
                    if i_seq < len(parts) and not parts[i_seq].strip() and len(parts) == n: noseq += 1
                    else: bad += 1
                    continue
                pid, sp, seq, mk = res
                if MARK and i_mk is not None and not MARK.search(mk):
                    other += 1; continue
                row = [""] * NCOL; row[1], row[24], row[66] = pid, sp, seq
                out.write("\t".join(row) + "\n"); kept += 1
        if header is None:
            print(f"   ⚠ {path}: no BOLD header, file skipped", file=sys.stderr); continue
        total += kept
        extra = [x for x in (f"{other} other markers skipped" if other else "",
                             f"{shifted} realigned" if shifted else "", f"{noseq} without sequence" if noseq else "",
                             f"{bad} unreadable lines" if bad else "") if x]
        print(f"   {path}: {kept} sequences" + (f" ({', '.join(extra)})" if extra else ""))
print(f"   {out_path}: {total} sequences in total" + (f" (marker /{marker}/)" if marker else ""))
sys.exit(0)
PYEOF
    if (( $(nlines raw/bold_all.tsv) > 1 )); then
      step 03_import_bold imported/bold_all.txt "import-boldv5|$RANKS" 0 raw/bold_all.tsv -- \
        crabs --import --import-format boldv5 --input raw/bold_all.tsv "${TAXARGS[@]}" --output imported/bold_all.txt
      IMPORTED+=("imported/bold_all.txt")
    else
      warn "BOLD: no $BOLD_MARKER record in the BOLD groups, BOLD not used"
      rm -f imported/bold_all.txt
    fi
  fi
fi

# =============================================================================
# 2b. NCBI GenBank: download per group, merge, import
# =============================================================================
if [[ "$USE_NCBI" == true ]]; then
  NCBI_FILES=(); NCBI_FAILED=()
  NCBI_NOT=""
  if [[ "$NCBI_EXCLUDE_AT_DOWNLOAD" == true && -n "$EXCLUDE_TAXA" ]]; then
    NCBI_NOT=" NOT ($(echo "$EXCLUDE_TAXA" | sed 's/;/"[Organism] OR "/g; s/^/"/; s/$/"[Organism]/'))"
  fi
  # ncbi_get GROUP FILE QUERY LOGNAME: download unless a valid file is already there
  ncbi_get() {
    local T=$1 F=$2 Q=$3 L=$4
    if [[ -e "$F" ]] && ! fasta_ok "$F"; then
      warn "NCBI: $F invalid (server error page?), downloading again"; rm -f "$F"
    fi
    if ! fasta_ok "$F"; then
      log "▶ NCBI: $L"
      TRY=1; WAIT=$PAUSE
      while ! fasta_ok "$F" && (( TRY <= MAX_TRIES )); do
        (( TRY > 1 )) && { log "  NCBI $L: retry $TRY/$MAX_TRIES in ${WAIT}s"; sleep "$WAIT"; WAIT=$(( WAIT * 2 )); }
        crabs --download-ncbi --database nucleotide --email "$EMAIL" --query "$Q" --output "$F" 2>&1 | tee "logs/01_ncbi_${L}.log" || true
        fasta_ok "$F" || rm -f "$F"
        grep -qE "sequences downloaded: 0[[:space:]]*$" "logs/01_ncbi_${L}.log" && break   # no sequence: no point retrying
        TRY=$(( TRY + 1 ))
      done
    fi
    if fasta_ok "$F"; then
      NCBI_FILES+=("$F")
      printf 'NCBI GenBank (nucleotide)\t%s\t%s\t%s\t%s\t%s\n' "$T" "$F" "$(fdate "$F")" "$(grep -c '^>' "$F")" "$Q" >> metadata/.sources.tmp
      return 0
    fi
    return 1
  }
  for T in "${TAXA_NCBI[@]}"; do
    ncbi_get "$T" "raw/ncbi_${T}.fasta" "\"${T}\"[Organism] AND ${NCBI_GENE_QUERY}${NCBI_NOT}" "$T" || NCBI_FAILED+=("$T")
    # Organelle records WITHOUT gene annotation (often "UNVERIFIED: … mitochondrion
    # sequence"): missed by the gene query, but the in silico PCR finds the
    # fragment in them. Duplicates of the gene query are removed at the merge.
    if [[ -n "${NCBI_EXTRA_QUERY:-}" ]]; then
      ncbi_get "$T" "raw/ncbi_unannotated_${T}.fasta" "\"${T}\"[Organism] AND ${NCBI_EXTRA_QUERY}${NCBI_NOT}" "unannotated_${T}" || true
    fi
  done
  (( ${#NCBI_FAILED[@]} > 0 )) && warn "NCBI: groups with no sequence or not retrieved: ${NCBI_FAILED[*]}"

  if (( ${#NCBI_FILES[@]} > 0 )); then
    SIG_NCBICAT="concat|$(printf '%s;' "${NCBI_FILES[@]}")"
    step 02_merge_ncbi raw/ncbi_all.fasta "$SIG_NCBICAT" 1 "${NCBI_FILES[@]}" -- \
      sh -c 'cat "$@" > raw/ncbi_all.fasta; echo "   raw/ncbi_all.fasta: $(grep -c "^>" raw/ncbi_all.fasta) sequences"' _ "${NCBI_FILES[@]}"
    step 03_import_ncbi imported/ncbi_all.txt "import-ncbi|$RANKS" 1 raw/ncbi_all.fasta -- \
      crabs --import --import-format ncbi --input raw/ncbi_all.fasta "${TAXARGS[@]}" --output imported/ncbi_all.txt
    IMPORTED+=("imported/ncbi_all.txt")
  fi
fi

# =============================================================================
# 2c. SILVA (curated rRNA database, all domains of life)
# =============================================================================
if [[ "$USE_SILVA" == true ]]; then
  F="raw/silva_${SILVA_GENE}_${SILVA_TYPE}_${SILVA_VERSION}.fasta"
  [[ -e "$F" ]] && ! fasta_ok "$F" && { warn "SILVA: $F invalid, downloading again"; rm -f "$F"; }
  if fasta_ok "$F"; then
    log "✔ SILVA: already downloaded → $F"
  else
    log "▶ SILVA ${SILVA_GENE} ${SILVA_TYPE} ${SILVA_VERSION}"
    TRY=1; WAIT=$PAUSE
    while ! fasta_ok "$F" && (( TRY <= MAX_TRIES )); do
      (( TRY > 1 )) && { log "  SILVA: retry $TRY/$MAX_TRIES in ${WAIT}s"; sleep "$WAIT"; WAIT=$(( WAIT * 2 )); }
      crabs --download-silva --gene "$SILVA_GENE" --db-type "$SILVA_TYPE" --db-version "$SILVA_VERSION" --output "$F" 2>&1 | tee logs/01_silva.log || true
      fasta_ok "$F" || rm -f "$F"
      TRY=$(( TRY + 1 ))
    done
  fi
  if fasta_ok "$F"; then
    printf 'SILVA %s %s\t%s\t%s\t%s\t%s\trelease %s\n' "$SILVA_GENE" "$SILVA_TYPE" "all" "$F" "$(fdate "$F")" "$(grep -c '^>' "$F")" "$SILVA_VERSION" >> metadata/.sources.tmp
    step 03_import_silva imported/silva_all.txt "import-silva|$SILVA_VERSION|$SILVA_TYPE|$RANKS" 0 "$F" -- \
      crabs --import --import-format silva --input "$F" "${TAXARGS[@]}" --output imported/silva_all.txt
    IMPORTED+=("imported/silva_all.txt")
  else
    warn "SILVA not retrieved (check SILVA_VERSION on https://www.arb-silva.de/download/arb-files/)"
  fi
fi

set -e
(( ${#IMPORTED[@]} > 0 )) || die "No data imported."
log "Imported sources: ${IMPORTED[*]}"

# =============================================================================
# 3. Merge and pre-curation (on full sequences, to speed up the alignment)
# =============================================================================
if (( ${#IMPORTED[@]} == 1 )); then   # CRABS refuses to merge a single file
  step 04_merge steps/01_merged.txt "copy|${IMPORTED[0]}" 0 "${IMPORTED[0]}" -- cp "${IMPORTED[0]}" steps/01_merged.txt
else
  step 04_merge steps/01_merged.txt "merge-uniq|$(printf '%s;' "${IMPORTED[@]}")" 0 "${IMPORTED[@]}" -- \
    crabs --merge --uniq --input "$(IFS=';'; echo "${IMPORTED[*]}")" --output steps/01_merged.txt
fi
# Lineage rescue: a species name unknown to the NCBI taxonomy (recent, provisional
# or BOLD-only names) is imported by CRABS with an empty lineage and is then
# removed by the filters. Its lineage is rebuilt from its GENUS when the genus is
# in the NCBI taxonomy (species name kept as given by the source).
step 04b_rescue_lineage steps/01b_lineage_rescued.txt "rescue-genus-v1|$RANKS" 0 steps/01_merged.txt -- \
  python3 - steps/01_merged.txt steps/01b_lineage_rescued.txt "$TAXO/nodes.dmp" "$TAXO/names.dmp" "$RANKS" <<'PYEOF'
import sys
src, dst, f_nodes, f_names, ranks = sys.argv[1:6]
ranks = ranks.split(";")
wanted = set(ranks[:-1])                       # ranks above species
parent, rank = {}, {}
with open(f_nodes, encoding="utf-8", errors="ignore") as f:
    for line in f:
        p = line.split("\t|\t")
        try: t = int(p[0]); parent[t] = int(p[1]); rank[t] = p[2]
        except (ValueError, IndexError): continue
need = {t for t, r in rank.items() if r in wanted}
sciname, genus_ids = {}, {}
with open(f_names, encoding="utf-8", errors="ignore") as f:
    for line in f:
        p = line.rstrip("\t|\n").split("\t|\t")
        if len(p) < 4 or p[3] != "scientific name":
            continue
        try: t = int(p[0])
        except ValueError: continue
        if t in need:
            sciname[t] = p[1]
        if rank.get(t) == "genus":
            genus_ids.setdefault(p[1], []).append(t)
def lineage(t):
    out = {}
    for _ in range(100):
        if rank.get(t) in wanted and t in sciname:
            out[rank[t]] = sciname[t]
        if t == parent.get(t, t):
            break
        t = parent[t]
    return out, t
def in_eukaryota(t):
    for _ in range(100):
        if t == 2759:
            return True
        if t == parent.get(t, t):
            return False
        t = parent[t]
    return False
rescued = unresolved = 0
with open(src, encoding="utf-8", errors="ignore") as fi, open(dst, "w", encoding="utf-8") as fo:
    for line in fi:
        p = line.rstrip("\n").split("\t")
        lin = p[3:3 + len(ranks)]
        if p[2] == "NA" and all(x == "NA" for x in lin):
            unresolved += 1
            genus = p[1].split(" ")[0]
            cand = [g for g in genus_ids.get(genus, []) if in_eukaryota(g)]
            if len(cand) == 1:                 # ambiguous genus names are left alone
                l, _ = lineage(cand[0])
                p[2] = str(cand[0])
                p[3:3 + len(ranks)] = [l.get(r, "NA") for r in ranks[:-1]] + [p[1]]
                rescued += 1
        fo.write("\t".join(p) + "\n")
print(f"   {unresolved} sequences without lineage; {rescued} rescued from their genus "
      f"({unresolved - rescued} left without lineage)")
PYEOF
# Excluded taxa removed right away, before the long steps
if [[ -n "$EXCLUDE_TAXA" ]]; then
  step 04c_exclude_early steps/01c_taxa_excluded_early.txt "exclude|$EXCLUDE_TAXA" 0 steps/01b_lineage_rescued.txt -- \
    crabs --subset --input steps/01b_lineage_rescued.txt --output steps/01c_taxa_excluded_early.txt --exclude "$EXCLUDE_TAXA"
else
  step 04c_exclude_early steps/01c_taxa_excluded_early.txt "exclude|none" 0 steps/01b_lineage_rescued.txt -- \
    cp steps/01b_lineage_rescued.txt steps/01c_taxa_excluded_early.txt
fi
step 05_dereplicate_initial steps/02_dereplicated_initial.txt "derep|$DEREP_METHOD" 0 steps/01c_taxa_excluded_early.txt -- \
  crabs --dereplicate --input steps/01c_taxa_excluded_early.txt --output steps/02_dereplicated_initial.txt --dereplication-method "$DEREP_METHOD"
step 06_filter_initial steps/03_filtered_initial.txt "filter|env|min=$MIN_ID_RANK|ranka$RANK_NA" 0 steps/02_dereplicated_initial.txt -- \
  crabs --filter --input steps/02_dereplicated_initial.txt --output steps/03_filtered_initial.txt --environmental \
        ${SPECIES_ARG[@]+"${SPECIES_ARG[@]}"} --rank-na "$RANK_NA"

# =============================================================================
# 4. In silico PCR
# =============================================================================
RELAXED_ARG=(); [[ "$RELAXED" == true ]] && RELAXED_ARG=(--relaxed)
step 07_insilico_pcr steps/04_insilico_pcr.txt "pcr|$FWD|$REV|mm$MISMATCH|relaxed=$RELAXED" 0 steps/03_filtered_initial.txt -- \
  crabs --in-silico-pcr --input steps/03_filtered_initial.txt --output steps/04_insilico_pcr.txt \
        --forward "$FWD" --reverse "$REV" --mismatch "$MISMATCH" ${RELAXED_ARG[@]+"${RELAXED_ARG[@]}"} --threads "$THREADS"

# =============================================================================
# 5. Recovery of amplicons deposited without primer-binding sites
# =============================================================================
if [[ "$USE_PGA" != true ]]; then
  step 08_global_alignment steps/05_global_alignment.txt "pga|off" 0 steps/04_insilico_pcr.txt -- \
    cp steps/04_insilico_pcr.txt steps/05_global_alignment.txt
else
step 08_global_alignment steps/05_global_alignment.txt "pga|$FWD|$REV|$PGA_IDENTITY|$PGA_COVERAGE|$PGA_SIZE_SELECT" 0 \
     steps/03_filtered_initial.txt steps/04_insilico_pcr.txt -- \
  crabs --pairwise-global-alignment --input steps/03_filtered_initial.txt --amplicons steps/04_insilico_pcr.txt \
        --output steps/05_global_alignment.txt --forward "$FWD" --reverse "$REV" \
        --size-select "$PGA_SIZE_SELECT" --percent-identity "$PGA_IDENTITY" --coverage "$PGA_COVERAGE"
fi

# =============================================================================
# 5b. Over-long amplicons re-cut at the primer sites found inside them
# =============================================================================
if [[ "${TRIM_LONG:-false}" == true ]]; then
step 09_trim_long steps/05b_long_trimmed.txt "trimlong-v1|$FWD|$REV|mm$MISMATCH|$MIN_LEN|$MAX_LEN" 0 steps/05_global_alignment.txt -- \
  python3 - steps/05_global_alignment.txt steps/05b_long_trimmed.txt logs/09_long_trimmed.tsv "$FWD" "$REV" "$MISMATCH" "$MIN_LEN" "$MAX_LEN" <<'PYEOF'
import sys
src, dst, log_f, fwd, rev, mm, lo, hi = sys.argv[1:9]
mm, lo, hi = int(mm), int(lo), int(hi)
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT",
         "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
comp = str.maketrans("ACGTRYKMSWBDHVN", "TGCAYRMKSWVHDBN")
rrc = rev.upper().translate(comp)[::-1]          # reverse primer as found on the amplicon strand
fwd = fwd.upper()
def sites(p, s):
    """(start, mismatches) of every site of primer p in s with <= mm mismatches, no indels"""
    out, L = [], len(p)
    for i in range(len(s) - L + 1):
        e = 0
        for a, b in zip(p, s[i:i + L]):
            if b not in IUPAC[a]:
                e += 1
                if e > mm: break
        else:
            out.append((i, e))
    return out
n = trimmed = left = 0
with open(src) as fi, open(dst, "w") as fo, open(log_f, "w") as fr:
    fr.write("seq_id\tname\told_length\tnew_length\tinner_mismatches\taction\n")
    for line in fi:
        p = line.rstrip("\n").split("\t"); s = p[-1].upper(); n += 1
        if len(s) <= hi:
            fo.write(line); continue
        # new start = after an inner forward site (or 0), new end = before an inner reverse site (or len)
        starts = [(0, 0)] + [(i + len(fwd), e) for i, e in sites(fwd, s)]
        ends = [(len(s), 0)] + [(j, e) for j, e in sites(rrc, s)]
        best = None
        for a, ea in starts:
            for b, eb in ends:
                if lo <= b - a <= hi and (best is None or ea + eb < best[2]):
                    best = (a, b, ea + eb)
        if best:
            p[-1] = s[best[0]:best[1]]; trimmed += 1
            fo.write("\t".join(p) + "\n")
            fr.write(f"{p[0]}\t{p[1]}\t{len(s)}\t{len(p[-1])}\t{best[2]}\tre-cut at inner primer sites\n")
        else:
            left += 1; fo.write(line)     # left as is: removed by the length filter
            fr.write(f"{p[0]}\t{p[1]}\t{len(s)}\t\t\tno inner primer site giving {lo}-{hi} bp\n")
print(f"   {n} amplicons; {trimmed + left} longer than {hi} bp: {trimmed} re-cut at inner primer sites, {left} left as is (removed by the length filter)")
PYEOF
  DEREP_IN=steps/05b_long_trimmed.txt
else
  DEREP_IN=steps/05_global_alignment.txt
fi

# =============================================================================
# 6. Final curation
# =============================================================================
step 09_dereplicate_final steps/06_dereplicated_final.txt "derep|$DEREP_METHOD|in=$DEREP_IN" 0 "$DEREP_IN" -- \
  crabs --dereplicate --input "$DEREP_IN" --output steps/06_dereplicated_final.txt --dereplication-method "$DEREP_METHOD"
FILTER_MIN=$MIN_LEN
if [[ -n "${MIN_LEN_RESCUE:-}" ]] && (( MIN_LEN_RESCUE < MIN_LEN )); then FILTER_MIN=$MIN_LEN_RESCUE; fi
step 10_filter_final steps/07_filtered_final.txt "filter|$FILTER_MIN|$MAX_LEN|n$MAX_N|env|min=$MIN_ID_RANK|ranka$RANK_NA" 0 steps/06_dereplicated_final.txt -- \
  crabs --filter --input steps/06_dereplicated_final.txt --output steps/07_filtered_final.txt \
        --minimum-length "$FILTER_MIN" --maximum-length "$MAX_LEN" --maximum-n "$MAX_N" \
        --environmental ${SPECIES_ARG[@]+"${SPECIES_ARG[@]}"} --rank-na "$RANK_NA"
need steps/07_filtered_final.txt
# Short amplicons (FILTER_MIN to MIN_LEN-1 bp) kept only for taxa (name column)
# without any full-length amplicon
step 10_length_rescue steps/07b_length_rescued.txt "rescue-len|$FILTER_MIN|$MIN_LEN" 0 steps/07_filtered_final.txt -- \
  awk -F'\t' -v min="$MIN_LEN" '
    NR == FNR { if (length($NF) >= min) full[$2] = 1; next }
    length($NF) >= min { print > "steps/07b_length_rescued.txt"; k++; next }
    !($2 in full)      { print > "steps/07b_length_rescued.txt"; r++; sp[$2] = 1; next }
                       { d++ }
    END { printf "   %d full-length amplicons; %d short amplicons kept for %d taxa without full length; %d short amplicons removed\n", k, r, length(sp), d }' \
    steps/07_filtered_final.txt steps/07_filtered_final.txt
if [[ "$TRANSLATION_CHECK" == true ]]; then
  step 10_translation_check steps/07c_translation_checked.txt "translation-v1|$STOP_CODONS" 0 steps/07b_length_rescued.txt -- \
    python3 - steps/07b_length_rescued.txt steps/07c_translation_checked.txt logs/10_translation_removed.tsv "$STOP_CODONS" <<'PYEOF'
import sys
src, dst, removed_log, stops = sys.argv[1], sys.argv[2], sys.argv[3], set(sys.argv[4].upper().split())
def open_frame(seq):
    """True if at least one of the 3 forward frames has no stop codon (the
    amplicons are already oriented by the in silico PCR)."""
    for f in range(3):
        codons = (seq[i:i + 3] for i in range(f, len(seq) - 2, 3))
        if not any(c in stops for c in codons):
            return True
    return False
kept = removed = 0
with open(src) as fi, open(dst, "w") as fo, open(removed_log, "w") as fr:
    fr.write("seq_id\tname\tlength\n")
    for line in fi:
        p = line.rstrip("\n").split("\t")
        seq = p[-1].upper()
        if open_frame(seq):
            fo.write(line); kept += 1
        else:
            fr.write(f"{p[0]}\t{p[1]}\t{len(seq)}\n"); removed += 1
print(f"   {kept} amplicons kept; {removed} removed (stop codon in every reading frame, listed in {removed_log})")
PYEOF
else
  step 10_translation_check steps/07c_translation_checked.txt "translation|off" 0 steps/07b_length_rescued.txt -- \
    cp steps/07b_length_rescued.txt steps/07c_translation_checked.txt
fi
# Genus homonyms: lineage checked against BOLD and the target groups, corrected
# when one genus of the same name fits (see TAXONOMY_FIX). Done before the
# exclusion of taxa, so that a misfiled insect is excluded once corrected.
if [[ "$TAXONOMY_FIX" == true ]]; then
  FIX_TARGETS="$(IFS=';'; echo "${TAXA_NCBI[*]};${TAXA_BOLD[*]};$TARGET_EXTRA")"
  FIX_BOLD="$(IFS=';'; echo "${BOLD_FILES[*]+${BOLD_FILES[*]}}")"
  # manual decisions (TAXONOMY_MANUAL): part of the signature, so that editing
  # the file re-runs the step
  FIX_MANUAL_SUM=none
  [[ -s "$TAXONOMY_MANUAL" ]] && FIX_MANUAL_SUM=$(cksum < "$TAXONOMY_MANUAL" | tr -s ' ' '_')
  FIX_SILVA=""
  [[ "${USE_SILVA:-false}" == true ]] && FIX_SILVA="raw/silva_${SILVA_GENE}_${SILVA_TYPE}_${SILVA_VERSION}.fasta"
  step 10_taxonomy_fix steps/07d_taxonomy_fixed.txt "taxfix-v5-organelle|unres=${UNRESOLVED_MISMATCH:-keep}|sk=${SAME_KINGDOM_OK:-}|silva=$FIX_SILVA|$FIX_TARGETS|$RANKS|$CLASS_TO_PHYLUM|$FIX_BOLD|target=$TARGET_CHECK|manual=$FIX_MANUAL_SUM" 0 steps/07c_translation_checked.txt -- \
    python3 - steps/07c_translation_checked.txt steps/07d_taxonomy_fixed.txt logs/10_taxonomy_fix.tsv \
            "$TAXO/nodes.dmp" "$TAXO/names.dmp" "$FIX_TARGETS" "$RANKS" "$CLASS_TO_PHYLUM" "$FIX_BOLD" "$EXCLUDE_TAXA" \
            "$TAXONOMY_MANUAL" "$TARGET_CHECK" "$FIX_SILVA" "${UNRESOLVED_MISMATCH:-keep}" "${SAME_KINGDOM_OK:-}" <<'PYEOF'
import collections, csv, os, sys
src, dst, log_f, f_nodes, f_names, targets, ranks, c2p, bold_files, excl, f_manual, target_check, f_silva, unresolved, same_kingdom = sys.argv[1:16]
same_kingdom = {x for x in same_kingdom.split(";") if x}
target_check = target_check == "true"
# manual decisions: genus or seq_id -> (keep | remove, reason)
manual = {}
if f_manual and os.path.exists(f_manual):
    for r in csv.DictReader((l for l in open(f_manual, encoding="utf-8") if l.strip() and not l.startswith("#")), delimiter="\t"):
        k, a = (r.get("taxon") or "").strip(), (r.get("action") or "").strip().lower()
        if k and a in ("keep", "remove"): manual[k] = (a, (r.get("reason") or "").strip())
    print(f"   {len(manual)} manual decisions read from {f_manual}")
targets = {t for t in targets.split(";") if t}
excluded = {t for t in excl.split(";") if t}
ranks = ranks.split(";"); wanted = set(ranks[:-1])
C2P = dict(kv.split("=", 1) for kv in c2p.split(";") if "=" in kv)
# phylum names meaning the same thing (BOLD, NCBI, Adl et al. 2019)
SAME = {"Bacillariophyta": "Ochrophyta", "Heterokontophyta": "Ochrophyta", "Myzozoa": "Dinoflagellata",
        "Dinophyta": "Dinoflagellata", "Pyrrophycophyta": "Dinoflagellata", "Magnoliophyta": "Streptophyta",
        "Tracheophyta": "Streptophyta", "Charophyta": "Streptophyta", "Bryophyta": "Streptophyta",
        "Cryptophycophyta": "Cryptophyta", "Prasinodermophyta": "Chlorophyta"}
norm = lambda x: SAME.get(x, x) if x and x != "NA" else ""
HIGH = {"phylum", "subphylum", "superclass", "class", "subclass", "infraclass", "superorder", "order"}

parent, rank = {}, {}
with open(f_nodes, encoding="utf-8", errors="ignore") as f:
    for line in f:
        p = line.split("\t|\t")
        try: t = int(p[0]); parent[t] = int(p[1]); rank[t] = p[2].strip("\t|\n ")
        except (ValueError, IndexError): continue
name, genus_ids, high_names, phylum_names, phylum_tax = {}, collections.defaultdict(list), set(), set(), {}
with open(f_names, encoding="utf-8", errors="ignore") as f:
    for line in f:
        p = line.rstrip("\t|\n").split("\t|\t")
        if len(p) < 4 or p[3] != "scientific name": continue
        try: t = int(p[0])
        except ValueError: continue
        name[t] = p[1]
        if rank.get(t) == "genus": genus_ids[p[1]].append(t)
        if rank.get(t) in HIGH: high_names.add(p[1])
        if rank.get(t) == "phylum": phylum_names.add(p[1]); phylum_tax.setdefault(p[1], t)
memo_r = {}
def under(S, memo):
    """function t -> True if taxid t or one of its ancestors (any rank) is named in S"""
    def f(t):
        path = []
        for _ in range(200):
            if t in memo: res = memo[t]; break
            path.append(t)
            if name.get(t) in S: res = True; break
            if t == parent.get(t, t): res = False; break
            t = parent[t]
        else: res = False
        for x in path: memo[x] = res
        return res
    return f
on_target, in_excluded = under(targets, {}), under(excluded, {})
if not target_check:          # database with all eukaryotes on purpose (SILVA): BOLD check only
    on_target = lambda t: True
def ranks_of(t):
    if t in memo_r: return memo_r[t]
    out, s = {}, t
    for _ in range(200):
        if rank.get(s) in wanted: out.setdefault(rank[s], name.get(s, "NA"))
        if s == parent.get(s, s): break
        s = parent[s]
    memo_r[t] = out
    return out
def kingdom_of(ph):
    """NCBI kingdom of a phylum name (e.g. Ascomycota -> Fungi), "" if unknown"""
    t = phylum_tax.get(ph) or phylum_tax.get({v: k for k, v in SAME.items()}.get(ph, ""), None)
    for _ in range(200):
        if t is None: return ""
        if rank.get(t) == "kingdom": return name.get(t, "")
        if t == parent.get(t, t): return ""
        t = parent[t]
    return ""
def phylum_of(t):
    r = ranks_of(t)
    return norm(r.get("phylum") or C2P.get(r.get("class", ""), ""))

rows = [l.rstrip("\n").split("\t") for l in open(src, encoding="utf-8", errors="ignore")]
ids = {p[0] for p in rows}
# phylum given by BOLD, from the original BOLD files (identifier = record_id,
# or processid.marker_code)
bold = {}
for fpath in [x for x in bold_files.split(";") if x]:
    try: fi = open(fpath, encoding="utf-8", errors="ignore")
    except OSError: print(f"   BOLD file not found, skipped: {fpath}"); continue
    with fi:
        h = [x.strip().lower() for x in fi.readline().rstrip("\n").split("\t")]
        ix = lambda c: h.index(c) if c in h else None
        i_rec, i_pid, i_mk, i_ph = ix("record_id"), ix("processid"), ix("marker_code"), ix("phylum")
        if i_ph is None: continue
        for line in fi:
            p = line.rstrip("\n").split("\t")
            if i_ph >= len(p): continue
            keys = []
            if i_rec is not None and i_rec < len(p): keys.append(p[i_rec])
            if i_pid is not None and i_pid < len(p):
                keys.append(p[i_pid])
                if i_mk is not None and i_mk < len(p): keys.append(f"{p[i_pid]}.{p[i_mk]}")
            for k in keys:
                if k in ids: bold[k] = p[i_ph]; break
# phylum given by SILVA: the first name of the SILVA lineage (header of the
# SILVA file) that is a phylum in the NCBI taxonomy (e.g. Metazoa (Animalia);
# ...;Arthropoda;... -> Arthropoda). CRABS identifier = accession before the 1st dot.
n_bold = len(bold)
silva, organelle = {}, {}
if f_silva and os.path.exists(f_silva):
    with open(f_silva, encoding="utf-8", errors="ignore") as fi:
        for line in fi:
            if not line.startswith(">"): continue
            sid = line[1:].split(".")[0]
            if sid not in ids or sid in silva: continue
            path = line.rstrip("\n").split(" ", 1)[1].split(";")[:-1] if " " in line else []
            if path and path[0].strip() in ("Bacteria", "Archaea"):
                organelle[sid] = ";".join(path[:5]); continue
            for x in path:
                x = x.split(" (")[0].strip()
                if x in phylum_names: silva[sid] = x; break
    print(f"   {len(silva)} sequences of the database with a phylum in their SILVA lineage; "
          f"{len(organelle)} classified by SILVA in Bacteria/Archaea (organelle rRNA: chloroplast or mitochondrion 16S)")
for k, v in silva.items():
    bold.setdefault(k, v)                      # same check as for BOLD
src_of = lambda k: "BOLD" if k in bold and k not in silva else "SILVA"
print(f"   {n_bold} sequences found in the BOLD files ({len([x for x in bold_files.split(';') if x])} files)")

stats = collections.Counter()
with open(dst, "w", encoding="utf-8") as fo, open(log_f, "w", encoding="utf-8") as fr:
    fr.write("seq_id\tname\taction\treason\tbold_phylum\told_taxid\told_lineage\tnew_taxid\tnew_lineage\n")
    def log(p, action, reason, b, old_t, old_l):
        fr.write(f"{p[0]}\t{p[1]}\t{action}\t{reason}\t{b}\t{old_t}\t{old_l}\t"
                 f"{p[2] if action == 'corrected' else ''}\t{'|'.join(p[3:-1]) if action == 'corrected' else ''}\n")
    for p in rows:
        try: t = int(p[2])
        except ValueError: t = None
        # SILVA classifies the sequence in Bacteria/Archaea: a chloroplast or
        # mitochondrion 16S filed by NCBI under its eukaryote host, not 18S
        if p[0] in organelle:
            stats["removed (organelle rRNA)"] += 1
            log(p, "removed", "SILVA lineage " + organelle[p[0]] + ": organelle rRNA (chloroplast or mitochondrion 16S), not nuclear 18S",
                "", p[2], "|".join(p[3:-1]))
            continue
        b = norm(bold.get(p[0], ""))
        cur = phylum_of(t) if t is not None else norm(p[4] if len(p) > 4 else "")
        ok_t = t is not None and on_target(t)
        mismatch = bool(b and cur and b != cur)
        if mismatch and same_kingdom and kingdom_of(b) in same_kingdom and kingdom_of(b) == kingdom_of(cur):
            mismatch = False                       # e.g. two fungal phyla: delimitation, not an error
        # agrees with BOLD, or in the targets, or in a taxon removed by the next
        # step (EXCLUDE_TAXA, e.g. genuine insects): unchanged
        flagged = not (not mismatch and (ok_t or b or (t is not None and in_excluded(t))))
        # manual decision (TAXONOMY_MANUAL): a seq_id always; a genus only for its
        # flagged sequences (the correct sequences of that genus are untouched)
        genus0 = p[1].replace("UNVERIFIED: ", "").split(" ")[0]
        m = manual.get(p[0]) or (manual.get(genus0) if flagged else None)
        if m:
            if m[0] == "keep":
                fo.write("\t".join(p) + "\n"); stats["kept (manual)"] += 1
            else:
                stats["removed (manual)"] += 1
            log(p, "kept" if m[0] == "keep" else "removed", "manual decision: " + m[1],
                bold.get(p[0], ""), p[2], "|".join(p[3:-1]))
            continue
        if not flagged:
            fo.write("\t".join(p) + "\n"); stats["unchanged"] += 1; continue
        old_t, old_l = p[2], "|".join(p[3:-1])
        nm = p[1].replace("UNVERIFIED: ", "")
        genus = nm.split(" ")[0]
        reason = f"phylum differs from {src_of(p[0])}" if mismatch else "lineage outside the target groups"
        if " " not in nm.strip() and genus in high_names:   # a higher taxon written as a genus (e.g. Ctenophora, Oligochaeta)
            log(p, "removed", reason + "; name is a higher taxon", bold.get(p[0], ""), old_t, old_l); stats["removed"] += 1; continue
        cands = genus_ids.get(genus, [])
        fit = [g for g in cands if b and phylum_of(g) == b] if mismatch else []
        if len(fit) != 1:
            if mismatch and unresolved == "remove":   # sequence probably from another organism
                stats["removed (other organism)"] += 1
                log(p, "removed", reason + f"; no genus of that name in the {src_of(p[0])} phylum: sequence probably from another organism (epibiont, parasite, contaminant)",
                    bold.get(p[0], ""), old_t, old_l)
                continue
            if mismatch and ok_t:                  # BOLD disagrees but NCBI lineage is plausible: kept
                fo.write("\t".join(p) + "\n"); stats["kept (BOLD disagrees)"] += 1
                log(p, "kept", reason + "; no single genus of the BOLD phylum, NCBI lineage in the targets", bold.get(p[0], ""), old_t, old_l)
                continue
            fit = [g for g in cands if on_target(g)]
        if len(fit) == 1:
            l = ranks_of(fit[0])
            p[2] = str(fit[0])
            p[3:3 + len(ranks)] = [l.get(r, "NA") for r in ranks[:-1]] + [p[1]]
            fo.write("\t".join(p) + "\n"); stats["corrected"] += 1
            log(p, "corrected", reason + "; homonym", bold.get(p[0], ""), old_t, old_l)
        else:
            stats["removed"] += 1
            log(p, "removed", reason + ("; several genera of that name fit" if fit else "; no genus of that name fits"),
                bold.get(p[0], ""), old_t, old_l)
print("   " + "; ".join(f"{k}: {v}" for k, v in stats.most_common()) + f" (details: {log_f})")
PYEOF
  EXCL_IN=steps/07d_taxonomy_fixed.txt
else
  EXCL_IN=steps/07c_translation_checked.txt
fi
if [[ -n "$EXCLUDE_TAXA" ]]; then
  step 10_exclude_taxa steps/08_taxa_excluded.txt "exclude|$EXCLUDE_TAXA|in=$EXCL_IN" 0 "$EXCL_IN" -- \
    crabs --subset --input "$EXCL_IN" --output steps/08_taxa_excluded.txt --exclude "$EXCLUDE_TAXA"
else
  step 10_exclude_taxa steps/08_taxa_excluded.txt "exclude|none|in=$EXCL_IN" 0 "$EXCL_IN" -- cp "$EXCL_IN" steps/08_taxa_excluded.txt
fi
# Keep only references named at least to MIN_ID_RANK. In species mode, informal
# names written with a space ("Ballia sp. 1", "cf.", "aff.") are removed too:
# CRABS only catches them when written with an underscore ("Ballia_sp.").
FINAL_TXT="results/${BASE_NAME}.txt"
step 10_min_rank steps/09_min_rank.txt "minrank|$MIN_ID_RANK|$MIN_RANK_COL" 0 steps/08_taxa_excluded.txt -- \
  awk -F'\t' -v c="$MIN_RANK_COL" -v sp="$([[ $MIN_ID_RANK == species ]] && echo 1 || echo 0)" -v out="steps/09_min_rank.txt" '
    $c == "NA" || $c == "" { drop++; next }
    sp && $(NF-1) ~ / (sp|cf|aff)\.( |$)/ { drop++; next }
    { print > out; kept++ }
    END { printf "   %d references kept, %d removed (not identified to the minimum rank)\n", kept, drop }' steps/08_taxa_excluded.txt

# Phylum harmonisation with Adl et al. (2019): see PHYLUM_RENAME / CLASS_TO_PHYLUM.
# The GenBank flag "UNVERIFIED: " is also removed from the names (the sequences
# themselves passed the translation check).
# Every change is counted in metadata/taxonomy_harmonisation.tsv.
PHYLUM_COL=$(( 3 + $(echo "$RANKS" | tr ';' '\n' | grep -n -x phylum | cut -d: -f1) ))
CLASS_COL=$(( 3 + $(echo "$RANKS" | tr ';' '\n' | grep -n -x class | cut -d: -f1) ))
step 10_harmonise_phylum "$FINAL_TXT" "adl2019-v2|$PHYLUM_RENAME|$CLASS_TO_PHYLUM|$PHYLUM_COL|$CLASS_COL" 0 steps/09_min_rank.txt -- \
  awk -F'\t' -v OFS='\t' -v pc="$PHYLUM_COL" -v cc="$CLASS_COL" -v ren="$PHYLUM_RENAME" -v c2p="$CLASS_TO_PHYLUM" \
      -v out="$FINAL_TXT" -v log_f=metadata/taxonomy_harmonisation.tsv '
    BEGIN { n = split(ren, a, ";"); for (i = 1; i <= n; i++) { split(a[i], kv, "="); R[kv[1]] = kv[2] }
            n = split(c2p, a, ";"); for (i = 1; i <= n; i++) { split(a[i], kv, "="); C[kv[1]] = kv[2] } }
    {
      # GenBank "UNVERIFIED: Genus species" names: flag removed from the name
      if (sub(/^UNVERIFIED: /, "", $2)) nu++
      if (($pc == "NA" || $pc == "") && ($cc in C)) { k = "phylum filled from class\t" $cc "\t" C[$cc]; ch[k]++; $pc = C[$cc]; nf++ }
      else if ($pc in R)                             { k = "phylum renamed\t" $pc "\t" R[$pc];          ch[k]++; $pc = R[$pc]; nr++ }
      if ($pc == "NA" || $pc == "") na++
      print > out
    }
    END {
      print "change\tfrom\tto\tn_sequences" > log_f
      for (k in ch) print k "\t" ch[k] > log_f
      printf "   phylum: %d renamed, %d filled from the class (Adl et al. 2019); %d sequences still without phylum\n", nr, nf, na
      printf "   names: \"UNVERIFIED: \" removed from %d names\n", nu
    }' steps/09_min_rank.txt


# =============================================================================
# 7. Exports
# =============================================================================
crabs2fasta() {   # key=value headers readable by OBITools
  awk -F'\t' '{ lin=$4; for (i=5; i<NF; i++) lin=lin "|" $i
    printf ">%s taxid=%s; scientific_name=%s; lineage=%s;\n%s\n", $1, $3, $2, lin, $NF }' "$1" > "$2"
  echo "   $2: $(grep -c '^>' "$2") sequences"
}
crabs2tsv() {     # table with explicit column names
  { printf 'seq_id\tscientific_name\ttaxid\t%s\tsequence\n' "$(echo "$RANKS" | tr ';' '\t')"; cat "$1"; } > "$2"
  echo "   $2: table with header"
}
step 11_export_fasta "results/${BASE_NAME}.fasta" "fasta-obitools-v1" 0 "$FINAL_TXT" -- crabs2fasta "$FINAL_TXT" "results/${BASE_NAME}.fasta"
step 11_export_tsv   "results/${BASE_NAME}.tsv"   "tsv-v1|$RANKS" 0 "$FINAL_TXT" -- crabs2tsv "$FINAL_TXT" "results/${BASE_NAME}.tsv"
step 11_export_insilico results/insilico_pcr_raw.fasta "fasta-obitools-v1" 0 steps/04_insilico_pcr.txt -- \
  crabs2fasta steps/04_insilico_pcr.txt results/insilico_pcr_raw.fasta
for FMT in ${EXPORT_FORMATS[@]+"${EXPORT_FORMATS[@]}"}; do
  case "$FMT" in
    blast-tax|blast-notax) OUTF="results/${BASE_NAME}_${FMT}" ;;
    *-text|dada2-species)  OUTF="results/${BASE_NAME}_${FMT}.txt" ;;
    *)                     OUTF="results/${BASE_NAME}_${FMT}.fasta" ;;
  esac
  CHECK="$OUTF"; [[ "$FMT" == blast-* ]] && CHECK="${OUTF}.nsq"
  if [[ -s "$CHECK" && "$CHECK" -nt "$FINAL_TXT" ]]; then log "✔ export $FMT: already done"; continue; fi
  log "▶ export $FMT"
  crabs --export --input "$FINAL_TXT" --output "$OUTF" --export-format "$FMT" 2>&1 | tee "logs/11_export_${FMT}.log"
done

# =============================================================================
# 8. Metadata (FAIR) – regenerated at every run
# =============================================================================
log "▶ Metadata"
{
  printf 'source\tgroup\tfile\tdownload_date\tn_records\tquery\n'
  cat metadata/.sources.tmp
} > metadata/sources.tsv
rm -f metadata/.sources.tmp

{
  printf 'parameter\tvalue\n'
  printf 'marker\t%s\n' "$MARKER_DESC"
  printf 'primer_forward\t%s%s\n' "$FWD" "${FWD_NAME:+ ($FWD_NAME)}"
  printf 'primer_reverse\t%s%s\n' "$REV" "${REV_NAME:+ ($REV_NAME)}"
  printf 'expected_length\t%s\n' "$EXPECTED_LEN"
  printf 'pcr_mismatch\t%s\n' "$MISMATCH"
  printf 'pcr_relaxed\t%s\n' "$RELAXED"
  printf 'pga_used\t%s\n' "$USE_PGA"
  printf 'pga_identity\t%s\n' "$PGA_IDENTITY"
  printf 'pga_coverage\t%s\n' "$PGA_COVERAGE"
  printf 'pga_size_select\t%s\n' "$PGA_SIZE_SELECT"
  printf 'trim_long_amplicons\t%s\n' "$([[ "${TRIM_LONG:-false}" == true ]] && echo "yes (amplicons > MAX_LEN re-cut at inner primer sites; log logs/09_long_trimmed.tsv)" || echo no)"
  printf 'min_length\t%s\n' "$MIN_LEN"
  printf 'max_length\t%s\n' "$MAX_LEN"
  printf 'min_length_rescue\t%s\n' "${MIN_LEN_RESCUE:-}"
  printf 'translation_check\t%s\n' "$([[ "$TRANSLATION_CHECK" == true ]] && echo "yes (stop codons: $STOP_CODONS)" || echo no)"
  printf 'max_N\t%s\n' "$MAX_N"
  printf 'min_identification_rank\t%s\n' "$MIN_ID_RANK"
  printf 'max_missing_ranks_above\t%s\n' "$MAX_MISSING_RANKS"
  printf 'crabs_rank_na\t%s\n' "$RANK_NA"
  printf 'dereplication\t%s\n' "$DEREP_METHOD"
  printf 'filters\t--environmental %s\n' "${SPECIES_ARG[*]:-}"
  printf 'excluded_taxa\t%s (removed after import, after the final filters%s)\n' "$EXCLUDE_TAXA" "$([[ "$NCBI_EXCLUDE_AT_DOWNLOAD" == true && -n "$EXCLUDE_TAXA" ]] && echo ', and excluded from the NCBI queries')"
  printf 'taxonomy_fix\t%s\n' "$([[ "$TAXONOMY_FIX" == true ]] && echo "yes (genus homonyms: BOLD and SILVA phylum$([[ "$TARGET_CHECK" == true ]] && echo ", then target groups + $TARGET_EXTRA"); log logs/10_taxonomy_fix.tsv)" || echo no)"
  printf 'taxonomy_manual\t%s\n' "$([[ -s "$TAXONOMY_MANUAL" ]] && echo "$TAXONOMY_MANUAL ($(grep -v '^#' "$TAXONOMY_MANUAL" | awk 'NF' | tail -n +2 | wc -l | tr -d ' ') decisions): $(grep -v '^#' "$TAXONOMY_MANUAL" | awk -F'\t' 'NF && h++ {printf "%s%s=%s", (n++ ? "; " : ""), $1, $2}')" || echo none)"
  printf 'taxonomic_ranks\t%s\n' "$RANKS"
  printf 'phylum_renamed\t%s\n' "$PHYLUM_RENAME"
  printf 'phylum_from_class\t%s\n' "$CLASS_TO_PHYLUM"
  printf 'taxonomy_reference\tAdl et al. (2019) J Eukaryot Microbiol 66: 4-119\n'
  printf 'bold_groups\t%s\n' "${TAXA_BOLD[*]}"
  printf 'ncbi_groups\t%s\n' "${TAXA_NCBI[*]}"
  printf 'bold_marker_filter\t%s\n' "$BOLD_MARKER"
  printf 'silva\t%s\n' "$([[ "$USE_SILVA" == true ]] && echo "$SILVA_GENE $SILVA_TYPE $SILVA_VERSION" || echo no)"
  printf 'ncbi_query\t"<group>"[Organism] AND %s%s\n' "$NCBI_GENE_QUERY" "${NCBI_NOT:-}"
  printf 'ncbi_query_unannotated\t%s\n' "${NCBI_EXTRA_QUERY:+\"<group>\"[Organism] AND $NCBI_EXTRA_QUERY${NCBI_NOT:-}}"
} > metadata/parameters.tsv

{
  printf 'step\tfile\tn_sequences\n'
  for P in ${IMPORTED[@]+"${IMPORTED[@]}"}; do printf 'import_%s\t%s\t%s\n' "$(basename "$P" _all.txt)" "$P" "$(nlines "$P")"; done
  for P in steps/01_merged.txt steps/01b_lineage_rescued.txt steps/01c_taxa_excluded_early.txt steps/02_dereplicated_initial.txt steps/03_filtered_initial.txt \
           steps/04_insilico_pcr.txt steps/05_global_alignment.txt steps/05b_long_trimmed.txt steps/06_dereplicated_final.txt \
           steps/07_filtered_final.txt steps/07b_length_rescued.txt steps/07c_translation_checked.txt steps/07d_taxonomy_fixed.txt steps/08_taxa_excluded.txt steps/09_min_rank.txt "$FINAL_TXT"; do
    printf '%s\t%s\t%s\n' "$(basename "$P" .txt)" "$P" "$(nlines "$P")"
  done
} > metadata/counts.tsv

LEN_STATS=$(awk -F'\t' '{print length($NF)}' "$FINAL_TXT" | sort -n | awk '{a[NR]=$1} END{if(NR) printf "%d / %d / %d", a[1], a[int((NR+1)/2)], a[NR]}')
N_FINAL=$(nlines "$FINAL_TXT")
N_SPECIES=$(cut -f2 "$FINAL_TXT" | sort -u | wc -l | tr -d ' ')
RUN_END=$(date '+%Y-%m-%dT%H:%M:%S%z')

(cd results && shasum -a 256 "${BASE_NAME}".* insilico_pcr_raw.fasta 2>/dev/null || true) > metadata/checksums_sha256.txt

TITLE="$TITLE" MARKER_DESC="$MARKER_DESC" PRIMER_CITES="$(printf '%s\n' ${PRIMER_CITES[@]+"${PRIMER_CITES[@]}"})" \
python3 - "$BASE_NAME" "$RUN_START" "$RUN_END" "$N_FINAL" "$N_SPECIES" "$LEN_STATS" "$EMAIL" <<'PYEOF'
import csv, json, os, sys
base, start, end, n_final, n_sp, lens, email = sys.argv[1:8]
def tsv(p):
    with open(p, encoding="utf-8") as f: return list(csv.DictReader(f, delimiter="\t"))
params   = {r["parameter"]: r["value"] for r in tsv("metadata/parameters.tsv")}
versions = {r["software"]: r["version"] for r in tsv("metadata/versions.tsv")}
sources, counts = tsv("metadata/sources.tsv"), tsv("metadata/counts.tsv")
meta = {
  "title": f"{os.environ['TITLE']} ({base})",
  "creator": email, "run_start": start, "run_end": end,
  "marker": os.environ["MARKER_DESC"],
  "n_sequences": int(n_final), "n_taxa": int(n_sp),
  "length_min_median_max_bp": lens,
  "files": {"fasta_obitools": f"results/{base}.fasta", "table_tsv": f"results/{base}.tsv",
            "crabs_format": f"results/{base}.txt", "insilico_pcr_raw": "results/insilico_pcr_raw.fasta",
            "checksums": "metadata/checksums_sha256.txt"},
  "parameters": params, "software": versions, "sources": sources, "counts_per_step": counts,
  "cite": ["Jeunen et al. (2022) CRABS. Mol Ecol Resour. doi:10.1111/1755-0998.13741",
           "Ratnasingham & Hebert (2007) BOLD. Mol Ecol Notes. doi:10.1111/j.1471-8286.2007.01678.x",
           "Sayers et al. NCBI GenBank / Taxonomy (https://www.ncbi.nlm.nih.gov)",
           "Adl et al. (2019) Revisions to the classification, nomenclature, and diversity of eukaryotes. J Eukaryot Microbiol 66: 4-119. doi:10.1111/jeu.12691"]
}
meta["cite"] += (["Quast et al. (2013) SILVA. Nucleic Acids Res. doi:10.1093/nar/gks1219"]
                 if params.get("silva", "no") != "no" else [])
meta["cite"] += [c for c in os.environ.get("PRIMER_CITES", "").split("\n") if c.strip()]
meta["cite"] += [
           "Martin (2011) cutadapt. EMBnet.journal. doi:10.14806/ej.17.1.200",
           "Rognes et al. (2016) VSEARCH. PeerJ. doi:10.7717/peerj.2584"]
meta["terms_of_use"] = "Data retrieved from BOLD, NCBI (and SILVA): follow the terms of use and citation requirements of each source."
json.dump(meta, open("metadata/run_metadata.json", "w", encoding="utf-8"), ensure_ascii=False, indent=2)

n_bold = sum(1 for s in sources if s["source"] == "BOLD Systems V5")
n_ncbi = sum(1 for s in sources if s["source"] == "NCBI GenBank (nucleotide)")
silva = params.get("silva", "no")
lines = [f"# {meta['title']}", "",
  f"- **Created**: {end[:10]} (run {start} → {end})",
  f"- **Contact**: {email}",
  f"- **Content**: {n_final} sequences, {n_sp} taxa; length min / median / max: {lens} bp",
  f"- **Primers**: {params['primer_forward']} / {params['primer_reverse']}", "",
  "## Files", "",
  "| File | Content |", "|---|---|",
  f"| `results/{base}.fasta` | Final database, headers `taxid=…; scientific_name=…; lineage=…;` (OBITools) |",
  f"| `results/{base}.tsv` | Same content as a table with column names |",
  f"| `results/{base}.txt` | CRABS format (usable by every CRABS function) |",
  "| `results/insilico_pcr_raw.fasta` | Raw in silico PCR output (before alignment and filters) |",
  "| `metadata/run_metadata.json` | All metadata, machine-readable |",
  "| `metadata/parameters.tsv`, `versions.tsv`, `sources.tsv`, `counts.tsv` | Parameters, software, sources, counts per step |",
  "| `metadata/checksums_sha256.txt` | Checksums of the result files |",
  "| `logs/` | Log of every step |", "",
  "## Methods", "",
  f"1. Download from BOLD (V5 portal, {n_bold} groups" + (f", records with marker /{params['bold_marker_filter']}/ only" if params.get('bold_marker_filter') else "") + f") and NCBI GenBank ({n_ncbi} groups)" + (f", SILVA {silva}" if silva != "no" else "") + ", plus the NCBI taxonomy.",
  f"2. Import into CRABS format with NCBI lineage ({params['taxonomic_ranks']}).",
  f"3. Merge (duplicate accessions removed); lineage rescue (a species name unknown to the NCBI taxonomy takes the lineage of its genus); removal of the excluded taxa; dereplication ({params['dereplication']}), removal of environmental sequences and of sequences not identified at least to {params['min_identification_rank']} (or with more than {params['max_missing_ranks_above']} missing rank above it).",
  f"4. In silico PCR (cutadapt, {params['pcr_mismatch']} mismatches, relaxed={params['pcr_relaxed']}).",
  (f"5. Pairwise global alignment (VSEARCH, identity ≥ {params['pga_identity']}, coverage ≥ {params['pga_coverage']} %) to recover amplicons deposited without primer-binding sites." if params.get('pga_used') == 'true' else "5. Pairwise global alignment not used (in silico PCR amplicons only)."),
  f"6. Final dereplication and filters: {params['min_length']}–{params['max_length']} bp, fewer than {params['max_N']} N" + (f"; amplicons of {params['min_length_rescue']}-{int(params['min_length'])-1} bp kept only for taxa without a full-length amplicon" if params.get('min_length_rescue') else "") + (f"; amplicons with a stop codon in every reading frame removed ({params['translation_check']})" if params.get('translation_check', 'no') != 'no' else "") + f"; taxa removed: {params.get('excluded_taxa') or 'none'}.",
  ("Genus homonyms: " + params["taxonomy_fix"] + ("; manual decisions: " + params["taxonomy_manual"] if params.get("taxonomy_manual", "none") != "none" else "") + "." if params.get("taxonomy_fix", "no") != "no" else ""),
  "7. Phylum names harmonised with Adl et al. (2019): " + f"renamed {params.get('phylum_renamed', '')}; phylum given from the class when NCBI has none ({params.get('phylum_from_class', '')}). Details in metadata/taxonomy_harmonisation.tsv; species names and NCBI taxids unchanged.", "",
  "## Sequences per step", "", "| Step | Sequences |", "|---|---|"]
lines += [f"| {c['step']} | {c['n_sequences']} |" for c in counts]
lines += ["", "## Software", "", "| Software | Version |", "|---|---|"]
lines += [f"| {k} | {v} |" for k, v in versions.items()]
lines += ["", "## Please cite", ""] + [f"- {c}" for c in meta["cite"]]
lines += ["", "Source details (download date, number of records, exact query) are in `metadata/sources.tsv`."]
open("metadata/README.md", "w", encoding="utf-8").write("\n".join(lines) + "\n")
PYEOF

echo
log "Done."
echo "   Final database : results/${BASE_NAME}.fasta  (${N_FINAL} sequences, ${N_SPECIES} taxa)"
echo "   Lengths        : ${LEN_STATS} bp (min / median / max)"
echo "   Metadata       : metadata/README.md and metadata/run_metadata.json"
echo "   Figures        : Rscript $(up scripts)/crabs_figures.R $PWD"