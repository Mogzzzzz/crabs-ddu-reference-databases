#!/usr/bin/env bash
#SBATCH --job-name=DDU_crabs_obitag_COI
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=100GB
#SBATCH --account=austral_benthic_cryptofauna
#SBATCH --mail-user=morgane.durand@mnhn.fr
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output=crabs_obitag_COI_%j.out
#SBATCH --error=crabs_obitag_COI_%j.err
# ============================================================================
# crabs_obitag_slurm_coi.sh - Index de référence (obirefidx) puis assignation
# (obitag) avec la base COI construite par CRABS (version 1.7.3 : homonymes
# corrigés, amplicons trop longs recoupés ; COI_marine_invertebrates_algae.fasta),
# sur un cluster SLURM.
# Même script que pour le 16S, le 18S et le rbcL ; la base et les OTU COI sont
# les valeurs par défaut (BASE_DEFAULT, OTU_DEFAULT ci-dessous), et peuvent être
# données en arguments.
#
# Dossier CRABS_v1/, d'où le job est soumis, avec :
#   crabs_obitag_slurm.sh       ce script
#   <base>.fasta                base CRABS (dossier results/ de crabs_<marqueur>)
#   new_taxdump.tar.gz          taxonomie NCBI (décompressée au 1er lancement dans
#                               new_taxdump/ ; ou directement ce dossier)
#   <OTU>.fasta                 OTU à assigner (facultatif : sans OTU, l'index seul)
#
# Usage : cd CRABS_v1 && sbatch crabs_obitag_slurm_COI.sh
#         (autres fichiers : sbatch crabs_obitag_slurm_COI.sh <base>.fasta [<OTU>.fasta …])
#         Sans fichier d'OTU présent : seul l'index est construit.
#         Sans SLURM : bash crabs_obitag_slurm_COI.sh
#
# Étapes :
#   1. conversion des en-têtes CRABS au format OBITools4
#      >id taxid=123; scientific_name=X; ...  ->  >id {"taxid":123,"scientific_name":"X"}
#   2. taxids inconnus de la taxonomie retirés, taxids fusionnés mis à jour
#   3. obirefidx (long ; sauté s'il est déjà fait et que la base n'a pas changé :
#      un job relancé reprend à obitag)
#   4. obitag sur chaque fichier d'OTU
#   5. obicsv : table OTU, taxon assigné, obitag_bestid
# Pas d'obiuniq : la base CRABS est déjà dérépliquée (unique_species).
#
# Résultats dans CRABS_v1/obitools4/, nommés d'après la base et les OTU :
#   <base>.index.fasta                 référence indexée (obitag -R)
#   <base>.removed.tsv                 références retirées (taxid inconnu)
#   <OTU>_taxo_crabs_<base>.fasta      OTU annotées par obitag
#   <OTU>_assignments_crabs_<base>.csv table d'assignation
#   + journaux (*.log)
# Chaque base a ses propres fichiers (nommés d'après la base) : les jobs des
# autres marqueurs peuvent tourner en même temps dans le même dossier, mais
# jamais deux jobs sur la même base.
# ============================================================================
set -euo pipefail
SCRIPT_VERSION="2.1-COI (2026-10-26)"

# ---------------------------- PARAMÈTRES -------------------------------------
TAXDUMP_ARCHIVE="new_taxdump.tar.gz"      # taxonomie téléchargée (même époque que les bases CRABS)
TAXDUMP="new_taxdump"                     # dossier où elle est décompressée
OUT="obitools4"
BASE_DEFAULT="COI_marine_invertebrates_algae.fasta"      # base CRABS COI v1.7.3 (crabs_coi/results/)
OTU_DEFAULT=(DDU_COI_swarm_clust.fasta)                   # OTU COI à assigner
# OBITools4 sur le cluster (module ; version la plus récente disponible)
module load obitools/4.4.0
# ----------------------------------------------------------------------------

# Sous SLURM, le script est copié ailleurs : on se place dans le dossier de soumission
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

if [ $# -ge 1 ]; then
    CRABS_FASTA="$1"; shift; OTU_FASTA=("$@")
else
    CRABS_FASTA="$BASE_DEFAULT"; OTU_FASTA=()
    for f in "${OTU_DEFAULT[@]}"; do
        if [ -s "$f" ]; then OTU_FASTA+=("$f"); else echo "OTU absentes (index seul) : $f"; fi
    done
fi

# Nombre de cœurs : celui alloué par SLURM s'il est connu, sinon OBITools décide
CPU_OPT=(); [ -n "${SLURM_CPUS_PER_TASK:-}" ] && CPU_OPT=(--max-cpu "$SLURM_CPUS_PER_TASK")

for tool in obirefidx obitag obicsv; do
    command -v "$tool" >/dev/null 2>&1 || { echo "$tool introuvable : charger OBITools4 (voir PARAMÈTRES)"; exit 1; }
done
[ -s "$CRABS_FASTA" ] || { echo "Base CRABS introuvable : $PWD/$CRABS_FASTA"; exit 1; }
for otu in ${OTU_FASTA[@]+"${OTU_FASTA[@]}"}; do
    [ -s "$otu" ] || { echo "Fichier d'OTU introuvable : $PWD/$otu"; exit 1; }
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

log "Terminé. Résultats dans $PWD/$OUT/"