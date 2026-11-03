#!/usr/bin/env bash
# =============================================================================
#  18S rRNA V4 "forward read" reference database – marine invertebrates and
#  other eukaryotes – CRABS
# =============================================================================
#  Why      : the 18S amplicon (~574 bp) is longer than what a MiSeq v2 2x250 run
#             can pair, so many reads are lost when forward and reverse reads are
#             aligned. Forward reads alone can be used instead, but compared with
#             a database of full amplicons their identity is underestimated (the
#             reference is ~3 times longer than the read). This database holds
#             only the region read by the forward read: the first FWD_LEN bp
#             right after the forward primer 18Smk1F.
#  Primer   : forward 18Smk1F CTGGTGCCAGCAGCCGCGGYAA (Machida & Knowlton 2012)
#  Input    : the curated full sequences of the 18S database built by
#             crabs_18S.sh (SOURCE_DIR/steps/03_filtered_initial.txt: SILVA +
#             NCBI + BOLD, merged, dereplicated, environmental and poorly
#             identified sequences removed). crabs_18S.sh must have run first.
#             Nothing is downloaded.
#  Output   : results/18S_V4_fwd<FWD_LEN>_marine_eukaryotes.fasta (OBITools-style headers)
#             results/18S_V4_fwd<FWD_LEN>_marine_eukaryotes.tsv   (table with column names)
#  Provenance (FAIR): metadata/<base>_* (parameters, software versions, sources,
#             counts, checksums, README), logs/<base>_* (one log per step)
#
#  Steps    : 1. forward-only in silico PCR (cutadapt, as CRABS does): forward
#                primer found → the FWD_LEN bp that follow it
#             2. rescue of sequences deposited without the primer site (e.g.
#                V4 sequences starting right after the primer): aligned with
#                VSEARCH on the regions of step 1, kept only if the WHOLE region
#                is covered (no shifted start)
#             3. dereplication, filters, removal of excluded taxa, minimum rank,
#                phylum names harmonised with Adl et al. (2019), exports
#
#  READS    : the reads must be cut the same way: forward reads with the primer
#             removed, cut to FWD_LEN bp, shorter ones discarded, e.g. with
#             OBITools4:  obigrep -l FWD_LEN reads.fasta | obiannotate --cut 1:FWD_LEN
#             (references and reads of the same length → identity comparable to
#             a paired-read assignment)
#
#  Folder   : the 18S folder itself, e.g. ~/…/CRABS/crabs_18s, next to
#             crabs_18S.sh. Every file written starts with the name of the
#             forward database (e.g. steps/18S_V4_fwd150_marine_eukaryotes_04_…,
#             metadata/18S_V4_fwd150_marine_eukaryotes_parameters.tsv,
#             logs/18S_V4_fwd150_marine_eukaryotes_…): the 18S files are untouched.
#  Usage    : cd ~/…/CRABS/crabs_18s
#             caffeinate -i bash crabs_18S_fwd.sh 2>&1 | tee -a crabs_18S_fwd.log
#
#  Resuming : the script can be re-run at any time. Every step already done
#  (output present, inputs unchanged, same parameters) is skipped. Changing
#  FWD_LEN re-runs everything from step 1 (a few minutes to hours, no download).
# =============================================================================
set -uo pipefail
# up NAME: first NAME found in this folder or in a parent folder (absolute path),
# so that the shared CRABS/taxo, CRABS/reference_database_creator and crabs_coi
# are found whatever the depth of the database folder (e.g. crabs_18S/18S_aligned)
up() { local d=$PWD; while :; do [[ -e "$d/$1" ]] && { echo "$d/$1"; return 0; }; [[ "$d" == / ]] && return 0; d=$(dirname "$d"); done; }   # nothing printed if not found
SCRIPT_VERSION="1.3.1 (2026-10-28, shared CRABS/taxo and CRABS/reference_database_creator; 1.3: from crabs_18S.sh 1.5: taxonomy fix with organelle rRNA removal; runs in the 18S folder)"

# ================================ PARAMETERS ==================================
FWD="CTGGTGCCAGCAGCCGCGGYAA"
FWD_NAME="18Smk1F, Machida & Knowlton 2012"   # primer name / reference (metadata only)
REV_NAME="none (forward read only)"
# Length of the forward reads AFTER removal of the primer, i.e. the length the
# reads are cut to. MiSeq v2 2x250: 250 bp - 22 bp primer = 228 bp at most,
# less if tags/spacers precede the primer or after quality trimming. Choose it
# from the length distribution of the forward reads, and cut the reads to it.
# DDU 18S clean forward reads (123 sequences, 11,052 reads): 34-191 bp, median
# 108 bp, none >= 200 bp. 150 bp keeps 36 sequences (56 % of the reads),
# 170 bp keeps 27 (49 %), 100 bp keeps 66 (71 %) but separates fewer taxa.
FWD_LEN=150
EMAIL="morgane.durand@mnhn.fr"
BASE_NAME="18S_V4_fwd${FWD_LEN}_marine_eukaryotes"
TITLE="18S rRNA V4 forward-read reference database (${FWD_LEN} bp after 18Smk1F) – marine invertebrates and other eukaryotes"
MARKER_DESC="18S rRNA (nuclear small subunit), V4 region, first ${FWD_LEN} bp after the forward primer 18Smk1F"
PRIMER_CITES=("Machida & Knowlton (2012) PCR primers for metazoan nuclear 18S and 28S ribosomal DNA sequences. PLoS ONE 7: e46180. doi:10.1371/journal.pone.0046180")

# 18S database built by crabs_18S.sh (its curated full sequences are the input).
# Layout: crabs_18S/18S_aligned (crabs_18S.sh) and crabs_18S/18S_forward (this
# script): SOURCE_DIR = ../18S_aligned. If this script sits in the 18S folder
# itself (older layout), SOURCE_DIR = "."; every file it writes starts with the
# name of the forward database (BASE_NAME_…), so no 18S file is overwritten.
SOURCE_DIR=""                  # "" = ../18S_aligned if it exists, else "."

# Forward-only in silico PCR
MISMATCH=4           # mismatches allowed in the primer (CRABS default)
# Rescue of sequences deposited without the primer site (VSEARCH global
# alignment on the regions found by the in silico PCR). false = primer only.
USE_RESCUE=true
RESCUE_IDENTITY=0.95
# Final length window: regions from the in silico PCR are exactly FWD_LEN bp;
# rescued ones may differ by a few bp (indels)
LEN_TOLERANCE=5
MIN_LEN=$(( FWD_LEN - LEN_TOLERANCE ))
MAX_LEN=$(( FWD_LEN + LEN_TOLERANCE ))
MAX_N=2              # CRABS removes sequences with MAX_N or MORE "N": 2 = at most 1 N allowed
# Lowest rank a reference must be identified to: species | genus | family
MIN_ID_RANK=genus
MAX_MISSING_RANKS=1  # missing ranks tolerated ABOVE that rank
# Phylum names harmonised with Adl et al. (2019, J Eukaryot Microbiol 66: 4-119)
# (same tables as crabs_18S.sh)
PHYLUM_RENAME="Bacillariophyta=Ochrophyta;Heterokontophyta=Ochrophyta;Myzozoa=Dinoflagellata;Dinophyta=Dinoflagellata"
CLASS_TO_PHYLUM="Phaeophyceae=Ochrophyta;Bacillariophyceae=Ochrophyta;Coscinodiscophyceae=Ochrophyta;Mediophyceae=Ochrophyta;Chrysophyceae=Ochrophyta;Synurophyceae=Ochrophyta;Raphidophyceae=Ochrophyta;Eustigmatophyceae=Ochrophyta;Phaeosacciophyceae=Ochrophyta;Pelagophyceae=Ochrophyta;Dictyochophyceae=Ochrophyta;Xanthophyceae=Ochrophyta;Bolidophyceae=Ochrophyta;Pinguiophyceae=Ochrophyta;Phaeothamniophyceae=Ochrophyta;Aurearenophyceae=Ochrophyta;Chrysomerophyceae=Ochrophyta;Schizocladiophyceae=Ochrophyta;Olisthodiscophyceae=Ochrophyta;Dinophyceae=Dinoflagellata"
# The short region separates fewer species than the full amplicon: several
# species may share one sequence (they are all kept, obitag then assigns their
# common ancestor)
DEREP_METHOD=unique_species
# Safety net (already removed in crabs_18S.sh)
EXCLUDE_TAXA="Bacteria;Archaea;Insecta;Collembola;Protura;Diplura;Arachnida;Chilopoda;Diplopoda;Pauropoda;Symphyla;Onychophora"
# Taxonomy fix, as in crabs_18S.sh 1.4 (step 10_taxonomy_fix): the phylum of the
# NCBI lineage is checked against BOLD and the SILVA lineage of the 18S
# database (genus homonyms corrected; sequences from another organism, e.g.
# ciliate 18S filed under Gammarus, removed; fungal phyla not compared). Applied
# to the forward regions, including those rescued here that the 18S database
# never saw. The manual decisions of the 18S database are reused. Since 1.3
# (as crabs_18S.sh 1.5), sequences that SILVA classifies in Bacteria/Archaea
# (chloroplast or mitochondrion 16S filed by NCBI under the host) are removed.
# false = off.
TAXONOMY_FIX=true
TARGET_CHECK=false
TARGET_EXTRA="Ochrophyta;Bacillariophyta;Phaeophyceae;Ichthyostraca;Mystacocarida"
UNRESOLVED_MISMATCH=remove
SAME_KINGDOM_OK="Fungi"
EXPORT_FORMATS=()
THREADS=8

# =============================================================================
#                      NOTHING TO EDIT BELOW THIS LINE
# =============================================================================

RUN_START=$(date '+%Y-%m-%dT%H:%M:%S%z')
PFX="${BASE_NAME}_"     # prefix of every file written (shared folder with the 18S database)
mkdir -p steps results logs metadata metadata/.steps

log()  { echo "[$(date '+%H:%M:%S')] $*"; }
warn() { echo -e "\n>>> WARNING: $*\n" >&2; }
die()  { echo -e "\n✘ $*" >&2; exit 1; }
need() { [[ -s "$1" ]] || die "Step failed: $1 was not created (see the message above and logs/)"; }
nlines() { [[ -s "$1" ]] && wc -l < "$1" | tr -d ' ' || echo 0; }

# Resume engine (same as crabs_18S.sh)
sigfile() { echo "metadata/.steps/$(echo "$1" | tr '/' '_').sig"; }
# Signatures are compared without the folders of the files they name (only the
# file names), so that moving the folders does not make finished steps run again
SIGNORM='s#[^;|=]*/##g'
fresh() {   # fresh OUTPUT SIGNATURE INPUTS...
  local out=$1 sig=$2; shift 2
  [[ -s "$out" ]] || return 1
  local sf; sf=$(sigfile "$out")
  [[ -e "$sf" && "$(sed "$SIGNORM" "$sf")" == "$(printf '%s\n' "$sig" | sed "$SIGNORM")" ]] || return 1
  local i; for i in "$@"; do [[ -e "$i" && "$i" -nt "$out" ]] && return 1; done
  return 0
}
step() {    # step NAME OUTPUT SIGNATURE INPUTS... -- COMMAND...
  local name=$1 out=$2 sig=$3; shift 3
  local inputs=()
  while [[ $# -gt 0 && "$1" != "--" ]]; do inputs+=("$1"); shift; done
  shift
  if fresh "$out" "$sig" ${inputs[@]+"${inputs[@]}"}; then
    log "✔ $name: already done → $out"; return 0
  fi
  log "▶ $name"
  rm -f "$(sigfile "$out")"
  "$@" 2>&1 | tee "logs/${PFX}${name}.log"
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
export PATH="$CRABS_FOUND:$PATH"
[[ -n "${CRABS_VENV:-}" && -d "$CRABS_VENV/bin" ]] && export PATH="$CRABS_FOUND:$CRABS_VENV/bin:$PATH"   # optional Python environment
hash -r
export PYTHONWARNINGS="${PYTHONWARNINGS:+$PYTHONWARNINGS,}ignore:The parameter --exclude is used more than once:UserWarning"
crabs -h >/dev/null 2>&1 || { crabs -h; die "crabs -h fails"; }
for tool in cutadapt vsearch python3; do
  command -v "$tool" >/dev/null 2>&1 || die "Missing tool: $tool"
done
log "CRABS : $(command -v crabs)"
CRABS_VERSION=$(python3 -c "import re,sys;print(re.search(r'[\"\']([0-9][^\"\']*)', open(sys.argv[1]).read()).group(1))" "$CRABS_FOUND/function/__init__.py" 2>/dev/null || echo "unknown")
CRABS_COMMIT=$(git -C "$CRABS_FOUND" rev-parse --short HEAD 2>/dev/null || echo "unknown")
{
  printf 'software\tversion\n'
  printf 'script\t%s\n' "$SCRIPT_VERSION"
  printf 'CRABS\t%s (commit %s)\n' "$CRABS_VERSION" "$CRABS_COMMIT"
  printf 'cutadapt\t%s\n' "$(cutadapt --version 2>&1 | head -n 1)"
  printf 'vsearch\t%s\n' "$(vsearch --version 2>&1 | head -n 1 | cut -d, -f1)"
  printf 'python\t%s\n' "$(python3 --version 2>&1)"
  printf 'system\t%s\n' "$(uname -srm)"
} > metadata/${PFX}versions.tsv

# ----------------------------- Input (crabs_18s) ------------------------------
if [[ -z "$SOURCE_DIR" ]]; then [[ -d ../18S_aligned ]] && SOURCE_DIR="../18S_aligned" || SOURCE_DIR="."; fi
SRC_FULL="$SOURCE_DIR/steps/03_filtered_initial.txt"
[[ -s "$SRC_FULL" ]] || die "Input not found: $SRC_FULL
   Run crabs_18S.sh first (at least up to step 06_filter_initial), or set SOURCE_DIR."
SRC_FULL=$(cd "$(dirname "$SRC_FULL")" && pwd)/$(basename "$SRC_FULL")
SOURCE_DIR=$(cd "$SOURCE_DIR" && pwd)
log "Input : $SRC_FULL ($(nlines "$SRC_FULL") curated full sequences)"

# Taxonomic ranks: those of the 18S database
RANKS=$(awk -F'\t' '$1 == "taxonomic_ranks" {print $2}' "$SOURCE_DIR/metadata/parameters.tsv" 2>/dev/null)
[[ -n "$RANKS" ]] || RANKS='domain;phylum;class;order;family;genus;species'
IFS=';' read -r -a RANK_LIST <<< "$RANKS"
MIN_RANK_POS=0
for i in "${!RANK_LIST[@]}"; do [[ "${RANK_LIST[$i]}" == "$MIN_ID_RANK" ]] && MIN_RANK_POS=$((i + 1)); done
(( MIN_RANK_POS > 0 )) || die "MIN_ID_RANK=$MIN_ID_RANK is not one of the ranks: $RANKS"
BELOW=$(( ${#RANK_LIST[@]} - MIN_RANK_POS ))
RANK_NA=$(( MAX_MISSING_RANKS + BELOW + 1 ))
MIN_RANK_COL=$(( 3 + MIN_RANK_POS ))
SPECIES_ARG=(); [[ "$MIN_ID_RANK" == species ]] && SPECIES_ARG=(--no-species-id)
log "  ranks: $RANKS | references identified at least to: $MIN_ID_RANK"

# =============================================================================
# 1. Forward-only in silico PCR: primer found → the FWD_LEN bp that follow it
#    cutadapt options as in CRABS (--no-indels, -e MISMATCH, full-primer overlap,
#    both strands); -l cuts after the primer, -m drops sequences that end earlier
# =============================================================================
step 07_insilico_pcr_fwd steps/${PFX}04_insilico_pcr.txt "pcr-fwd-v1|$FWD|mm$MISMATCH|len$FWD_LEN" "$SRC_FULL" -- \
  python3 - "$SRC_FULL" steps/${PFX}04_insilico_pcr.txt "$FWD" "$MISMATCH" "$FWD_LEN" "$THREADS" <<'PYEOF'
import os, subprocess, sys, tempfile
src, dst, fwd, mm, L, threads = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]), sys.argv[6]
rows = {}
tmp = tempfile.mkdtemp(prefix="fwdpcr_", dir="steps")
fin, fout = os.path.join(tmp, "in.fasta"), os.path.join(tmp, "out.fasta")
with open(src) as fi, open(fin, "w") as fo:
    for line in fi:
        p = line.rstrip("\n").split("\t")
        rows[p[0]] = p
        fo.write(f">{p[0]}\n{p[-1]}\n")
cmd = ["cutadapt", fin, "-g", fwd, "--no-indels", "-e", mm, "--overlap", str(len(fwd)),
       "--revcomp", "--discard-untrimmed", "-l", str(L), "-m", str(L),
       "--cores", threads, "--quiet", "-o", fout]
subprocess.run(cmd, check=True)
n = 0; rc = 0
with open(fout) as fi, open(dst, "w") as fo:
    name = None
    for line in fi:
        line = line.rstrip("\n")
        if line.startswith(">"):
            f = line[1:].split()
            name = f[0]; rc += (len(f) > 1 and f[-1] == "rc")
        else:
            p = rows[name][:]
            p[-1] = line
            fo.write("\t".join(p) + "\n"); n += 1
for f in (fin, fout): os.remove(f)
os.rmdir(tmp)
print(f"   {n} of {len(rows)} full sequences carry the forward primer and at least {L} bp after it ({rc} on the reverse strand)")
PYEOF

# =============================================================================
# 2. Rescue of sequences deposited without the primer site. Each sequence
#    without a primer hit (both orientations) is aligned (VSEARCH, identity ≥
#    RESCUE_IDENTITY) on the unique regions of step 1. Kept only if the region
#    is covered from its first to its last base: its start is then known
#    exactly (a shifted start would lower the identity with the reads).
# =============================================================================
if [[ "$USE_RESCUE" == true ]]; then
step 08_rescue_fwd steps/${PFX}05_global_alignment.txt "rescue-fwd-v1|$RESCUE_IDENTITY" "$SRC_FULL" steps/${PFX}04_insilico_pcr.txt -- \
  python3 - "$SRC_FULL" steps/${PFX}04_insilico_pcr.txt steps/${PFX}05_global_alignment.txt "$RESCUE_IDENTITY" "$THREADS" <<'PYEOF'
import os, subprocess, sys, tempfile
full, amp, dst, ident, threads = sys.argv[1:6]
comp = str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn", "TGCAYRMKSWVHDBNtgcayrmkswvhdbn")
rc = lambda s: s.translate(comp)[::-1]
found = set(); targets = {}
with open(amp) as fi, open(dst, "w") as fo:
    for line in fi:
        p = line.rstrip("\n").split("\t")
        found.add(p[0]); targets.setdefault(p[-1].upper(), len(targets))
        fo.write(line)
tmp = tempfile.mkdtemp(prefix="fwdrescue_", dir="steps")
fq, fdb, fhit = (os.path.join(tmp, x) for x in ("query.fasta", "db.fasta", "hits.tsv"))
with open(fdb, "w") as fo:
    for s, i in targets.items(): fo.write(f">t{i}\n{s}\n")
rows = {}
with open(full) as fi, open(fq, "w") as fo:
    for line in fi:
        p = line.rstrip("\n").split("\t")
        if p[0] in found: continue
        rows[p[0]] = p
        s = p[-1].upper()
        fo.write(f">{p[0]}|+\n{s}\n>{p[0]}|-\n{rc(s)}\n")   # both orientations, plus strand only
n_kept = 0
if rows:
    subprocess.run(["vsearch", "--usearch_global", fq, "--db", fdb, "--id", ident, "--strand", "plus",
                    "--userout", fhit, "--userfields", "query+target+id+qilo+qihi+tilo+tihi+tl",
                    "--maxaccepts", "1", "--threads", threads, "--quiet"], check=True)
    best = {}
    with open(fhit) as fi:
        for line in fi:
            q, t, idv, qilo, qihi, tilo, tihi, tl = line.rstrip("\n").split("\t")
            if int(tilo) != 1 or int(tihi) != int(tl): continue      # region not fully covered
            sid, strand = q.rsplit("|", 1)
            if sid not in best or float(idv) > best[sid][0]:
                best[sid] = (float(idv), strand, int(qilo), int(qihi))
    with open(dst, "a") as fo:
        for sid, (idv, strand, a, b) in best.items():
            p = rows[sid][:]
            s = p[-1].upper() if strand == "+" else rc(p[-1].upper())
            p[-1] = s[a - 1:b]
            fo.write("\t".join(p) + "\n"); n_kept += 1
for f in (fq, fdb, fhit):
    if os.path.exists(f): os.remove(f)
os.rmdir(tmp)
print(f"   {len(found)} regions from the in silico PCR; {len(rows)} sequences without primer site aligned, {n_kept} rescued (whole region covered, identity >= {ident})")
PYEOF
else
  step 08_rescue_fwd steps/${PFX}05_global_alignment.txt "rescue|off" steps/${PFX}04_insilico_pcr.txt -- \
    cp steps/${PFX}04_insilico_pcr.txt steps/${PFX}05_global_alignment.txt
fi

# =============================================================================
# 3. Final curation (as in crabs_18S.sh)
# =============================================================================
step 09_dereplicate_final steps/${PFX}06_dereplicated_final.txt "derep|$DEREP_METHOD" steps/${PFX}05_global_alignment.txt -- \
  crabs --dereplicate --input steps/${PFX}05_global_alignment.txt --output steps/${PFX}06_dereplicated_final.txt --dereplication-method "$DEREP_METHOD"
step 10_filter_final steps/${PFX}07_filtered_final.txt "filter|$MIN_LEN|$MAX_LEN|n$MAX_N|env|min=$MIN_ID_RANK|ranka$RANK_NA" steps/${PFX}06_dereplicated_final.txt -- \
  crabs --filter --input steps/${PFX}06_dereplicated_final.txt --output steps/${PFX}07_filtered_final.txt \
        --minimum-length "$MIN_LEN" --maximum-length "$MAX_LEN" --maximum-n "$MAX_N" \
        --environmental ${SPECIES_ARG[@]+"${SPECIES_ARG[@]}"} --rank-na "$RANK_NA"
# Inputs of the taxonomy fix, all taken from the 18S database
TAXO=""
for t in "$(up taxo)" "$SOURCE_DIR/taxo" "$(up crabs_coi)/taxo"; do   # shared CRABS/taxo first
  [[ -s "$t/nodes.dmp" && -s "$t/names.dmp" ]] && { TAXO="$t"; break; }
done
[[ -n "$TAXO" || "$TAXONOMY_FIX" != true ]] || die "NCBI taxonomy (nodes.dmp) not found next to $SOURCE_DIR"
TAXA_NCBI=(); TAXA_BOLD=(); BOLD_FILES=(); SILVA_FILE=""
while IFS=$'\t' read -r src grp file rest; do
  case "$src" in
    "BOLD Systems V5") f="$SOURCE_DIR/$file"; [[ -s "$f" ]] || f="$(up crabs_coi)/raw/$(basename "$file")"   # file moved with its folder
                       [[ -s "$f" ]] && BOLD_FILES+=("$f") ;;
    SILVA*)            [[ -s "$SOURCE_DIR/$file" ]] && SILVA_FILE="$SOURCE_DIR/$file" ;;
  esac
done < <(tail -n +2 "$SOURCE_DIR/metadata/sources.tsv" 2>/dev/null)
TAXONOMY_MANUAL="$SOURCE_DIR/taxonomy_manual_crabs.tsv"
log "  taxonomy fix: ${#BOLD_FILES[@]} BOLD files, SILVA ${SILVA_FILE:-none}, taxonomy $TAXO"
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
  [[ -s "${SILVA_FILE:-}" ]] && FIX_SILVA="$SILVA_FILE"
  step 10_taxonomy_fix steps/${PFX}07d_taxonomy_fixed.txt "taxfix-v5-organelle|unres=${UNRESOLVED_MISMATCH:-keep}|sk=${SAME_KINGDOM_OK:-}|silva=$FIX_SILVA|$FIX_TARGETS|$RANKS|$CLASS_TO_PHYLUM|$FIX_BOLD|target=$TARGET_CHECK|manual=$FIX_MANUAL_SUM" steps/${PFX}07_filtered_final.txt -- \
    python3 - steps/${PFX}07_filtered_final.txt steps/${PFX}07d_taxonomy_fixed.txt logs/${PFX}10_taxonomy_fix.tsv \
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
  EXCL_IN=steps/${PFX}07d_taxonomy_fixed.txt
else
  EXCL_IN=steps/${PFX}07_filtered_final.txt
fi
if [[ -n "$EXCLUDE_TAXA" ]]; then
  step 10_exclude_taxa steps/${PFX}08_taxa_excluded.txt "exclude|$EXCLUDE_TAXA|in=$EXCL_IN" "$EXCL_IN" -- \
    crabs --subset --input "$EXCL_IN" --output steps/${PFX}08_taxa_excluded.txt --exclude "$EXCLUDE_TAXA"
else
  step 10_exclude_taxa steps/${PFX}08_taxa_excluded.txt "exclude|none|in=$EXCL_IN" "$EXCL_IN" -- cp "$EXCL_IN" steps/${PFX}08_taxa_excluded.txt
fi
FINAL_TXT="results/${BASE_NAME}.txt"
step 10_min_rank steps/${PFX}09_min_rank.txt "minrank|$MIN_ID_RANK|$MIN_RANK_COL" steps/${PFX}08_taxa_excluded.txt -- \
  awk -F'\t' -v c="$MIN_RANK_COL" -v sp="$([[ $MIN_ID_RANK == species ]] && echo 1 || echo 0)" -v out="steps/${PFX}09_min_rank.txt" '
    $c == "NA" || $c == "" { drop++; next }
    sp && $(NF-1) ~ / (sp|cf|aff)\.( |$)/ { drop++; next }
    { print > out; kept++ }
    END { printf "   %d references kept, %d removed (not identified to the minimum rank)\n", kept, drop }' steps/${PFX}08_taxa_excluded.txt
PHYLUM_COL=$(( 3 + $(echo "$RANKS" | tr ';' '\n' | grep -n -x phylum | cut -d: -f1) ))
CLASS_COL=$(( 3 + $(echo "$RANKS" | tr ';' '\n' | grep -n -x class | cut -d: -f1) ))
step 10_harmonise_phylum "$FINAL_TXT" "adl2019-v2|$PHYLUM_RENAME|$CLASS_TO_PHYLUM|$PHYLUM_COL|$CLASS_COL" steps/${PFX}09_min_rank.txt -- \
  awk -F'\t' -v OFS='\t' -v pc="$PHYLUM_COL" -v cc="$CLASS_COL" -v ren="$PHYLUM_RENAME" -v c2p="$CLASS_TO_PHYLUM" \
      -v out="$FINAL_TXT" -v log_f=metadata/${PFX}taxonomy_harmonisation.tsv '
    BEGIN { n = split(ren, a, ";"); for (i = 1; i <= n; i++) { split(a[i], kv, "="); R[kv[1]] = kv[2] }
            n = split(c2p, a, ";"); for (i = 1; i <= n; i++) { split(a[i], kv, "="); C[kv[1]] = kv[2] } }
    {
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
    }' steps/${PFX}09_min_rank.txt

# =============================================================================
# 4. Exports
# =============================================================================
crabs2fasta() {
  awk -F'\t' '{ lin=$4; for (i=5; i<NF; i++) lin=lin "|" $i
    printf ">%s taxid=%s; scientific_name=%s; lineage=%s;\n%s\n", $1, $3, $2, lin, $NF }' "$1" > "$2"
  echo "   $2: $(grep -c '^>' "$2") sequences"
}
crabs2tsv() {
  { printf 'seq_id\tscientific_name\ttaxid\t%s\tsequence\n' "$(echo "$RANKS" | tr ';' '\t')"; cat "$1"; } > "$2"
  echo "   $2: table with header"
}
step 11_export_fasta "results/${BASE_NAME}.fasta" "fasta-obitools-v1" "$FINAL_TXT" -- crabs2fasta "$FINAL_TXT" "results/${BASE_NAME}.fasta"
step 11_export_tsv   "results/${BASE_NAME}.tsv"   "tsv-v1|$RANKS" "$FINAL_TXT" -- crabs2tsv "$FINAL_TXT" "results/${BASE_NAME}.tsv"
step 11_export_insilico results/${PFX}insilico_pcr_raw.fasta "fasta-obitools-v1" steps/${PFX}04_insilico_pcr.txt -- \
  crabs2fasta steps/${PFX}04_insilico_pcr.txt results/${PFX}insilico_pcr_raw.fasta
for FMT in ${EXPORT_FORMATS[@]+"${EXPORT_FORMATS[@]}"}; do
  case "$FMT" in
    blast-tax|blast-notax) OUTF="results/${BASE_NAME}_${FMT}" ;;
    *-text|dada2-species)  OUTF="results/${BASE_NAME}_${FMT}.txt" ;;
    *)                     OUTF="results/${BASE_NAME}_${FMT}.fasta" ;;
  esac
  log "▶ export $FMT"
  crabs --export --input "$FINAL_TXT" --output "$OUTF" --export-format "$FMT" 2>&1 | tee "logs/${PFX}11_export_${FMT}.log"
done

# =============================================================================
# 5. Metadata (FAIR) – regenerated at every run
# =============================================================================
log "▶ Metadata"
# Sources: those of the 18S database (the data were downloaded by crabs_18S.sh)
if [[ -s "$SOURCE_DIR/metadata/sources.tsv" ]]; then
  cp "$SOURCE_DIR/metadata/sources.tsv" metadata/${PFX}sources.tsv
else
  printf 'source\tgroup\tfile\tdownload_date\tn_records\tquery\n' > metadata/${PFX}sources.tsv
fi
printf 'CRABS 18S database\tcurated full sequences\t%s\t%s\t%s\t-\n' "$SRC_FULL" \
  "$(python3 -c 'import os,sys,datetime;print(datetime.datetime.fromtimestamp(os.path.getmtime(sys.argv[1])).strftime("%Y-%m-%d"))' "$SRC_FULL")" \
  "$(nlines "$SRC_FULL")" >> metadata/${PFX}sources.tsv
SRC_PARAM() { awk -F'\t' -v k="$1" '$1 == k {print $2}' "$SOURCE_DIR/metadata/parameters.tsv" 2>/dev/null; }

{
  printf 'parameter\tvalue\n'
  printf 'marker\t%s\n' "$MARKER_DESC"
  printf 'primer_forward\t%s (%s)\n' "$FWD" "$FWD_NAME"
  printf 'primer_reverse\t%s\n' "$REV_NAME"
  printf 'expected_length\t%s\n' "$FWD_LEN"
  printf 'forward_read_length\t%s\n' "$FWD_LEN"
  printf 'pcr_mismatch\t%s\n' "$MISMATCH"
  printf 'pcr_relaxed\tforward primer only\n'
  printf 'rescue_used\t%s\n' "$USE_RESCUE"
  printf 'rescue_identity\t%s\n' "$RESCUE_IDENTITY"
  printf 'min_length\t%s\n' "$MIN_LEN"
  printf 'max_length\t%s\n' "$MAX_LEN"
  printf 'max_N\t%s\n' "$MAX_N"
  printf 'min_identification_rank\t%s\n' "$MIN_ID_RANK"
  printf 'max_missing_ranks_above\t%s\n' "$MAX_MISSING_RANKS"
  printf 'crabs_rank_na\t%s\n' "$RANK_NA"
  printf 'dereplication\t%s\n' "$DEREP_METHOD"
  printf 'excluded_taxa\t%s\n' "$EXCLUDE_TAXA"
  printf 'taxonomy_fix\t%s\n' "$([[ "$TAXONOMY_FIX" == true ]] && echo "yes (as crabs_18S.sh 1.5: BOLD and SILVA phylum, other organisms and organelle rRNA removed, fungal phyla not compared; log logs/${PFX}10_taxonomy_fix.tsv)" || echo no)"
  printf 'taxonomic_ranks\t%s\n' "$RANKS"
  printf 'phylum_renamed\t%s\n' "$PHYLUM_RENAME"
  printf 'phylum_from_class\t%s\n' "$CLASS_TO_PHYLUM"
  printf 'taxonomy_reference\tAdl et al. (2019) J Eukaryot Microbiol 66: 4-119\n'
  printf 'source_database\t%s\n' "$SOURCE_DIR"
  printf 'bold_groups\t%s\n' "$(SRC_PARAM bold_groups)"
  printf 'ncbi_groups\t%s\n' "$(SRC_PARAM ncbi_groups)"
  printf 'bold_marker_filter\t%s\n' "$(SRC_PARAM bold_marker_filter)"
  printf 'silva\t%s\n' "$(SRC_PARAM silva)"
  printf 'ncbi_query\t%s\n' "$(SRC_PARAM ncbi_query)"
} > metadata/${PFX}parameters.tsv

{
  printf 'step\tfile\tn_sequences\n'
  for P in "$SRC_FULL" steps/${PFX}04_insilico_pcr.txt steps/${PFX}05_global_alignment.txt steps/${PFX}06_dereplicated_final.txt \
           steps/${PFX}07_filtered_final.txt steps/${PFX}07d_taxonomy_fixed.txt steps/${PFX}08_taxa_excluded.txt steps/${PFX}09_min_rank.txt "$FINAL_TXT"; do
    printf '%s\t%s\t%s\n' "$(basename "$P" .txt)" "$P" "$(nlines "$P")"
  done
} > metadata/${PFX}counts.tsv

LEN_STATS=$(awk -F'\t' '{print length($NF)}' "$FINAL_TXT" | sort -n | awk '{a[NR]=$1} END{if(NR) printf "%d / %d / %d", a[1], a[int((NR+1)/2)], a[NR]}')
N_FINAL=$(nlines "$FINAL_TXT")
N_SPECIES=$(cut -f2 "$FINAL_TXT" | sort -u | wc -l | tr -d ' ')
# Sequences shared by several taxa (the short region does not separate them)
N_SHARED=$(awk -F'\t' '{print $NF "\t" $2}' "$FINAL_TXT" | sort -u | cut -f1 | sort | uniq -d | wc -l | tr -d ' ')
RUN_END=$(date '+%Y-%m-%dT%H:%M:%S%z')
SHA=$(command -v shasum >/dev/null 2>&1 && echo "shasum -a 256" || echo sha256sum)
(cd results && $SHA "${BASE_NAME}".* "${PFX}insilico_pcr_raw.fasta" 2>/dev/null || true) > metadata/${PFX}checksums_sha256.txt

PFX="$PFX" TITLE="$TITLE" MARKER_DESC="$MARKER_DESC" PRIMER_CITES="$(printf '%s\n' ${PRIMER_CITES[@]+"${PRIMER_CITES[@]}"})" \
python3 - "$BASE_NAME" "$RUN_START" "$RUN_END" "$N_FINAL" "$N_SPECIES" "$LEN_STATS" "$EMAIL" "$N_SHARED" <<'PYEOF'
import csv, json, os, sys
base, start, end, n_final, n_sp, lens, email, n_shared = sys.argv[1:9]
X = os.environ.get("PFX", "")      # prefix of the metadata files
def tsv(p):
    with open(p, encoding="utf-8") as f: return list(csv.DictReader(f, delimiter="\t"))
params   = {r["parameter"]: r["value"] for r in tsv(f"metadata/{X}parameters.tsv")}
versions = {r["software"]: r["version"] for r in tsv(f"metadata/{X}versions.tsv")}
sources, counts = tsv(f"metadata/{X}sources.tsv"), tsv(f"metadata/{X}counts.tsv")
meta = {
  "title": f"{os.environ['TITLE']} ({base})",
  "creator": email, "run_start": start, "run_end": end,
  "marker": os.environ["MARKER_DESC"],
  "n_sequences": int(n_final), "n_taxa": int(n_sp), "n_sequences_shared_by_several_taxa": int(n_shared),
  "length_min_median_max_bp": lens,
  "files": {"fasta_obitools": f"results/{base}.fasta", "table_tsv": f"results/{base}.tsv",
            "crabs_format": f"results/{base}.txt", "insilico_pcr_raw": f"results/{X}insilico_pcr_raw.fasta",
            "checksums": f"metadata/{X}checksums_sha256.txt"},
  "parameters": params, "software": versions, "sources": sources, "counts_per_step": counts,
  "cite": ["Jeunen et al. (2022) CRABS. Mol Ecol Resour. doi:10.1111/1755-0998.13741",
           "Ratnasingham & Hebert (2007) BOLD. Mol Ecol Notes. doi:10.1111/j.1471-8286.2007.01678.x",
           "Sayers et al. NCBI GenBank / Taxonomy (https://www.ncbi.nlm.nih.gov)",
           "Quast et al. (2013) SILVA. Nucleic Acids Res. doi:10.1093/nar/gks1219",
           "Adl et al. (2019) Revisions to the classification, nomenclature, and diversity of eukaryotes. J Eukaryot Microbiol 66: 4-119. doi:10.1111/jeu.12691"]
}
meta["cite"] += [c for c in os.environ.get("PRIMER_CITES", "").split("\n") if c.strip()]
meta["cite"] += ["Martin (2011) cutadapt. EMBnet.journal. doi:10.14806/ej.17.1.200",
                 "Rognes et al. (2016) VSEARCH. PeerJ. doi:10.7717/peerj.2584"]
meta["terms_of_use"] = "Data retrieved from BOLD, NCBI and SILVA: follow the terms of use and citation requirements of each source."
json.dump(meta, open(f"metadata/{X}run_metadata.json", "w", encoding="utf-8"), ensure_ascii=False, indent=2)
L = params["forward_read_length"]
lines = [f"# {meta['title']}", "",
  f"- **Created**: {end[:10]} (run {start} → {end})",
  f"- **Contact**: {email}",
  f"- **Content**: {n_final} sequences, {n_sp} taxa; length min / median / max: {lens} bp; {n_shared} sequences shared by several taxa",
  f"- **Primer**: {params['primer_forward']}; region = the {L} bp following it (forward read)",
  f"- **Built from**: the curated full sequences of the 18S database `{params['source_database']}`", "",
  "## Use", "",
  f"For forward reads only: primer removed, reads cut to {L} bp (shorter reads discarded), e.g. "
  f"`obigrep -l {L} reads.fasta | obiannotate --cut 1:{L}`. References and reads then have the same length.", "",
  "## Files", "",
  "| File | Content |", "|---|---|",
  f"| `results/{base}.fasta` | Final database, headers `taxid=…; scientific_name=…; lineage=…;` (OBITools) |",
  f"| `results/{base}.tsv` | Same content as a table with column names |",
  f"| `results/{base}.txt` | CRABS format |",
  f"| `results/{X}insilico_pcr_raw.fasta` | Regions found by the forward-only in silico PCR (before rescue and filters) |",
  f"| `metadata/{X}run_metadata.json` | All metadata, machine-readable |",
  f"| `metadata/{X}parameters.tsv`, `…versions.tsv`, `…sources.tsv`, `…counts.tsv` | Parameters, software, sources, counts per step |",
  f"| `metadata/{X}checksums_sha256.txt` | Checksums of the result files |",
  f"| `logs/{X}*` | Log of every step |", "",
  "## Methods", "",
  f"1. Input: curated full sequences of the 18S database (SILVA {params.get('silva','')}, NCBI GenBank, BOLD; merged, lineage rescued, excluded taxa removed, dereplicated, environmental and poorly identified sequences removed; see its own README).",
  f"2. Forward-only in silico PCR (cutadapt, {params['pcr_mismatch']} mismatches, no indels, both strands): the {L} bp following the forward primer; sequences ending earlier discarded.",
  (f"3. Sequences without the primer site aligned (VSEARCH, identity ≥ {params['rescue_identity']}) on the regions of step 2; kept only when the whole region is covered, so that its start is exact." if params.get('rescue_used') == 'true' else "3. No rescue of sequences without the primer site."),
  f"4. Dereplication ({params['dereplication']}), filters: {params['min_length']}–{params['max_length']} bp, fewer than {params['max_N']} N, environmental sequences and sequences not identified at least to {params['min_identification_rank']} removed; taxa removed: {params.get('excluded_taxa') or 'none'}.",
  "5. Phylum names harmonised with Adl et al. (2019): " + f"renamed {params.get('phylum_renamed', '')}; phylum given from the class when NCBI has none. Details in metadata/{X}taxonomy_harmonisation.tsv.", "",
  "## Sequences per step", "", "| Step | Sequences |", "|---|---|"]
lines += [f"| {c['step']} | {c['n_sequences']} |" for c in counts]
lines += ["", "## Software", "", "| Software | Version |", "|---|---|"]
lines += [f"| {k} | {v} |" for k, v in versions.items()]
lines += ["", "## Please cite", ""] + [f"- {c}" for c in meta["cite"]]
open(f"metadata/{X}README.md", "w", encoding="utf-8").write("\n".join(lines) + "\n")
PYEOF

echo
log "Done."
echo "   Final database : results/${BASE_NAME}.fasta  (${N_FINAL} sequences, ${N_SPECIES} taxa)"
echo "   Lengths        : ${LEN_STATS} bp (min / median / max)"
echo "   Shared         : ${N_SHARED} sequences shared by several taxa (not separated by ${FWD_LEN} bp)"
echo "   Reads          : forward reads must be cut to ${FWD_LEN} bp too (see the header of this script)"
echo "   Figures        : Rscript $(up scripts)/crabs_figures.R $PWD ${BASE_NAME}   (figures in figures_${BASE_NAME}/)"