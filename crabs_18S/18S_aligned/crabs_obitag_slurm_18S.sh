#!/usr/bin/env bash
#SBATCH --job-name=DDU_crabs_obitag_18S_2026
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=100GB
#SBATCH --account=austral_benthic_cryptofauna
#SBATCH --mail-user=morgane.durand@mnhn.fr
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
# ============================================================================
# crabs_obitag_slurm_18S.sh - Index de référence (obirefidx) puis assignation
# (obitag) avec les deux bases 18S construites par CRABS, sur un cluster SLURM :
#   normal : 18S_V4V5_marine_eukaryotes.fasta        (amplicon entier, ~574 pb)
#   fwd    : 18S_V4_fwd150_marine_eukaryotes.fasta   (150 pb après 18Smk1F,
#            pour les lectures forward seules, coupées à 150 pb)
# Même script que pour le COI et le rbcL.
#
# Dossier CRABS_v1/, d'où le job est soumis, avec :
#   crabs_obitag_slurm_18S.sh                ce script
#   18S_V4V5_marine_eukaryotes.fasta         base 18S (crabs_18s/results/)
#   18S_V4_fwd150_marine_eukaryotes.fasta    base 18S forward (crabs_18s/results/)
#   new_taxdump.tar.gz / new_taxdump/        taxonomie NCBI (déjà là pour le COI)
#   OTU à assigner : OTU_NORMAL (nom à vérifier) et OTU_FWD
#   (DDU_18S_forward_only_renamed.fasta, coupées automatiquement à FWD_LEN pb) ;
#   absentes = index seul
#
# Usage : un job par base, en parallèle (l'index de la base normale, avec
# SILVA, est bien plus long que celui de la base forward) :
#   cd CRABS_v1
#   sbatch --job-name=DDU_crabs_obitag_18S_normal crabs_obitag_slurm_18S.sh normal
#   sbatch --job-name=DDU_crabs_obitag_18S_fwd    crabs_obitag_slurm_18S.sh fwd
#   Les deux à la suite dans un seul job : sbatch crabs_obitag_slurm_18S.sh
#   Sans SLURM : bash crabs_obitag_slurm_18S.sh normal|fwd
#
# Étapes, pour chaque base :
#   1. conversion des en-têtes CRABS au format OBITools4
#   2. taxids inconnus de la taxonomie retirés, taxids fusionnés mis à jour
#   3. obirefidx (sauté s'il est déjà fait et que la base n'a pas changé)
#   4. obitag sur chaque fichier d'OTU
#   5. obicsv : table OTU, taxon assigné, obitag_bestid
#
# Résultats dans CRABS_v1/obitools4/, nommés d'après la base et les OTU :
#   <base>.index.fasta, <base>.removed.tsv, <OTU>_taxo_crabs_<base>.fasta,
#   <OTU>_assignments_crabs_<base>.csv, journaux (*.log)
# Journaux du job : <nom du job>_<numéro>.out / .err
# ============================================================================
set -euo pipefail
SCRIPT_VERSION="2.0-18S (2026-10-23)"

# ---------------------------- PARAMÈTRES -------------------------------------
TAXDUMP_ARCHIVE="new_taxdump.tar.gz"      # taxonomie téléchargée (même époque que les bases CRABS)
TAXDUMP="new_taxdump"                     # dossier où elle est décompressée
OUT="obitools4"
BASE_NORMAL="18S_V4V5_marine_eukaryotes.fasta"
OTU_NORMAL=(DDU_18S_swarm_clust.fasta)              # OTU 18S (lectures appariées) ; nom à vérifier
BASE_FWD="18S_V4_fwd150_marine_eukaryotes.fasta"
OTU_FWD=(DDU_18S_forward_only_renamed.fasta)       # OTU 18S forward
# Les lectures forward doivent avoir la longueur de la base forward (150 pb
# après l'amorce) : les OTU forward sont coupées à FWD_LEN pb, les plus courtes
# écartées, dans obitools4/<OTU>_cut<FWD_LEN>.fasta, avant obitag. "" = pas de coupe
# (OTU déjà coupées).
FWD_LEN=150
# OBITools4 sur le cluster (module ; version la plus récente disponible)
module load obitools/4.4.0
# ----------------------------------------------------------------------------

# Sous SLURM, le script est copié ailleurs : on se place dans le dossier de soumission
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

WHICH="${1:-both}"
case "$WHICH" in normal|fwd|both) ;; *) echo "Usage : sbatch crabs_obitag_slurm_18S.sh normal|fwd   (sans argument : les deux)"; exit 1 ;; esac

# Nombre de cœurs : celui alloué par SLURM s'il est connu, sinon OBITools décide
CPU_OPT=(); [ -n "${SLURM_CPUS_PER_TASK:-}" ] && CPU_OPT=(--max-cpu "$SLURM_CPUS_PER_TASK")

for tool in obirefidx obitag obicsv; do
    command -v "$tool" >/dev/null 2>&1 || { echo "$tool introuvable : charger OBITools4 (voir PARAMÈTRES)"; exit 1; }
done
# Taxonomie : décompressée une seule fois depuis l'archive. Dans un dossier
# temporaire puis renommée, pour que deux jobs lancés en même temps ne se gênent pas.
if [ ! -s "$TAXDUMP/nodes.dmp" ] && [ -s "$TAXDUMP_ARCHIVE" ]; then
    tmpd="${TAXDUMP}.tmp.$$"
    log "Décompression de $TAXDUMP_ARCHIVE"
    mkdir -p "$tmpd" && tar -xzf "$TAXDUMP_ARCHIVE" -C "$tmpd"
    if [ ! -s "$TAXDUMP/nodes.dmp" ]; then rm -rf "$TAXDUMP"; mv "$tmpd" "$TAXDUMP"; else rm -rf "$tmpd"; fi
fi
[ -s "$TAXDUMP/nodes.dmp" ] || { echo "Taxonomie introuvable : ni $PWD/$TAXDUMP/nodes.dmp ni $PWD/$TAXDUMP_ARCHIVE"; exit 1; }
mkdir -p "$OUT"

run_base() {      # run_base BASE.fasta [OTU.fasta …] (absent OTU files are skipped: index only)
CRABS_FASTA="$1"; shift
OTU_FASTA=()
for f in "$@"; do
    if [ -s "$f" ]; then OTU_FASTA+=("$f"); else log "OTU absentes (index seul) : $f"; fi
done
[ -s "$CRABS_FASTA" ] || { log "Base CRABS introuvable : $PWD/$CRABS_FASTA"; return 1; }
base=$(basename "$CRABS_FASTA" .fasta)
conv="$OUT/${base}.obi.fasta"
taxok="$OUT/${base}.taxok.fasta"
removed="$OUT/${base}.removed.tsv"
index="$OUT/${base}.index.fasta"
log "crabs_obitag_slurm.sh $SCRIPT_VERSION | base : $CRABS_FASTA ($(grep -c '^>' "$CRABS_FASTA") séquences)"
log "Dossier : $PWD | $(obitag --version 2>&1 | head -n 1) | cœurs : ${SLURM_CPUS_PER_TASK:-auto}"

# 1-2. Conversion et contrôle des taxids (rapides, refaits si la base est plus récente)
if [ -s "$taxok" ] && [ "$taxok" -nt "$CRABS_FASTA" ]; then
    log "Référence déjà convertie : $taxok"
else
    log "Conversion OBITools4 : $CRABS_FASTA"
    awk '/^>/ {
            id = substr($1, 2); tx = ""; nm = ""
            if (match($0, /taxid=[^;]*;/))           tx = substr($0, RSTART + 6, RLENGTH - 7)
            if (match($0, /scientific_name=[^;]*;/)) nm = substr($0, RSTART + 16, RLENGTH - 17)
            gsub(/\\/, "\\\\", nm); gsub(/"/, "\\\"", nm)
            printf ">%s {\"taxid\":%s,\"scientific_name\":\"%s\"}\n", id, tx, nm
            next
         }
         { print }' "$CRABS_FASTA" > "$conv"
    log "    $(grep -c '^>' "$conv") séquences"

    log "Contrôle des taxids avec $TAXDUMP"
    awk -v removed="$removed" '
        FILENAME ~ /nodes\.dmp$/  { split($0, a, "\t"); ok[a[1]] = 1; next }
        FILENAME ~ /merged\.dmp$/ { split($0, a, "\t"); mg[a[1]] = a[3]; next }
        /^>/ {
            keep = 0; t = ""
            if (match($0, /"taxid":[0-9]+/)) t = substr($0, RSTART + 8, RLENGTH - 8)
            if (t in ok) keep = 1
            else if ((t in mg) && (mg[t] in ok)) {
                sub("\"taxid\":" t ",", "\"taxid\":" mg[t] ","); keep = 1; upd++
            }
            if (!keep) { print substr($1, 2) "\t" t > removed; rm++ }
            else n++
        }
        keep { print }
        END {
            printf "    %d gardées, %d taxids mis à jour (fusionnés), %d retirées (taxid inconnu)\n", n, upd, rm > "/dev/stderr"
        }' "$TAXDUMP/nodes.dmp" "$TAXDUMP/merged.dmp" "$conv" > "$taxok"
fi

# 3. Index (n'apparaît que complet : un job interrompu le refait au prochain lancement)
if [ -s "$index" ] && [ "$index" -nt "$taxok" ]; then
    log "Référence déjà indexée : $index"
else
    log "obirefidx : long"
    obirefidx ${CPU_OPT[@]+"${CPU_OPT[@]}"} -t "$TAXDUMP" "$taxok" > "$index.tmp" 2> "$OUT/${base}.obirefidx.log"
    mv "$index.tmp" "$index"
    log "    -> $index ($(grep -c '^>' "$index") références indexées)"
fi

# 4-5. Assignation et table, pour chaque fichier d'OTU
if [ ${#OTU_FASTA[@]} -eq 0 ]; then log "Pas de fichier d'OTU : index seul"; fi
for otu in ${OTU_FASTA[@]+"${OTU_FASTA[@]}"}; do
    ob=$(basename "$otu" .fasta)
    taxo="$OUT/${ob}_taxo_crabs_${base}.fasta"
    table="$OUT/${ob}_assignments_crabs_${base}.csv"
    log "obitag : $otu ($(grep -c '^>' "$otu") OTU)"
    obitag ${CPU_OPT[@]+"${CPU_OPT[@]}"} -t "$TAXDUMP" -R "$index" "$otu" > "$taxo" 2> "$OUT/${ob}_obitag_crabs_${base}.log"
    log "    -> $taxo"
    obicsv -i -t "$TAXDUMP" --taxon -k obitag_bestid "$taxo" > "$table" 2>> "$OUT/${ob}_obitag_crabs_${base}.log"
    log "    -> $table"
done

}

status=0
if [ "$WHICH" = normal ] || [ "$WHICH" = both ]; then
    log "===== 18S (amplicon entier)"; run_base "$BASE_NORMAL" "${OTU_NORMAL[@]}" || status=1
fi
if [ "$WHICH" = fwd ] || [ "$WHICH" = both ]; then
    log "===== 18S forward (${FWD_LEN:-?} pb)"
    OTU_FWD_RUN=()
    for f in "${OTU_FWD[@]}"; do
        if [ -s "$f" ] && [ -n "${FWD_LEN:-}" ]; then
            cut_f="$OUT/$(basename "$f" .fasta)_cut${FWD_LEN}.fasta"
            if [ ! -s "$cut_f" ] || [ "$f" -nt "$cut_f" ]; then
                # séquences sur une ligne, coupées à FWD_LEN, les plus courtes écartées
                awk -v L="$FWD_LEN" -v out="$cut_f" '
                    function flush() { if (h != "") { n++; if (length(q) >= L) { print h > out; print substr(q, 1, L) > out; k++ } } }
                    /^>/ { flush(); h = $0; q = ""; next }
                    { q = q $0 }
                    END { flush(); printf "    %d OTU sur %d gardées (>= %d pb), coupées à %d pb\n", k, n, L, L > "/dev/stderr" }' "$f"
            fi
            log "OTU forward coupées : $cut_f"
            OTU_FWD_RUN+=("$cut_f")
        else
            OTU_FWD_RUN+=("$f")
        fi
    done
    run_base "$BASE_FWD" "${OTU_FWD_RUN[@]}" || status=1
fi
log "Terminé. Résultats dans $PWD/$OUT/"
exit $status