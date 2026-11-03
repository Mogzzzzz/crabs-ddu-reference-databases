#!/usr/bin/env bash
# =============================================================================
#  Trace species through the CRABS pipeline (diagnostic)
# =============================================================================
#  For each species of a list (e.g. species found by BLAST / BOLD ID engine but
#  absent from the database), counts its sequences at every step of
#  crabs_coi.sh and reports the step where it disappears.
#  Also reports which source brought insects (or any class) into the database.
#
#  Usage : cd ~/…/CRABS/crabs_coi
#          bash ../scripts/crabs_trace_species.sh [species_list.txt] [class_to_check]
#          without a list: species_gb_bold.txt (here or in a parent folder)
#          species_list.txt = one species name per line (e.g. "Mytilus edulis"),
#          or a table (TSV/CSV) whose header has a species column (Espèce, species…)
#          class_to_check   = default Insecta
#  Output: trace_species.tsv (+ printed table)
# =============================================================================
set -uo pipefail
# up NAME: NAME in this folder or in a parent folder (e.g. the shared CRABS/taxo)
up() { local d=$PWD; while :; do [[ -e "$d/$1" ]] && { echo "$d/$1"; return 0; }; [[ "$d" == / ]] && return 0; d=$(dirname "$d"); done; }
# Default list (no argument): the first one found here or in the CRABS folder
if [[ -n "${1:-}" ]]; then LIST=$1; else
  LIST=""
  for L in "$(up species_gb_bold.txt)" "$(up species_bold_gb.txt)"; do
    [[ -s "$L" ]] && { LIST=$L; break; }
  done
  [[ -n "$LIST" ]] || { echo "No species list given and species_gb_bold.txt not found (here or in a parent folder)" >&2; exit 1; }
  echo "  species list: $LIST" >&2
fi
CHECK_CLASS=${2:-Insecta}
[[ -s "$LIST" ]] || { echo "List not found or empty: $LIST" >&2; exit 1; }
BASE=$(ls results/*.txt 2>/dev/null | head -n 1)
OUT=trace_species.tsv

# Stages: label | file | how to match the species
#   bold_raw : raw/bold_all.tsv, 68-column layout, taxon name in column 25
#   ncbi_raw : raw/ncbi_all.fasta, species name inside the FASTA header
#   crabs    : CRABS format, column 2 = name given by the source, column 10 = species (NCBI)
STAGES=(
  "bold_raw|raw/bold_all.tsv|bold_raw"
  "ncbi_raw|raw/ncbi_all.fasta|ncbi_raw"
  "import_bold|imported/bold_all.txt|crabs"
  "import_ncbi|imported/ncbi_all.txt|crabs"
  "merged|steps/01_merged.txt|crabs"
  "lineage_rescued|steps/01b_lineage_rescued.txt|crabs"
  "derep_initial|steps/02_dereplicated_initial.txt|crabs"
  "filter_initial|steps/03_filtered_initial.txt|crabs"
  "insilico_pcr|steps/04_insilico_pcr.txt|crabs"
  "global_alignment|steps/05_global_alignment.txt|crabs"
  "derep_final|steps/06_dereplicated_final.txt|crabs"
  "filter_final|steps/07_filtered_final.txt|crabs"
  "length_rescued|steps/07b_length_rescued.txt|crabs"
  "translation_ok|steps/07c_translation_checked.txt|crabs"
  "taxa_excluded|steps/08_taxa_excluded.txt|crabs"
  "final|${BASE:-results/none.txt}|crabs"
)

TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
# The list can be a plain list (one species per line) or a table (TSV/CSV) with a
# header containing a species column (species, espèce, espece, scientific_name…)
tr -d '\r' < "$LIST" > "$TMP/raw_list"
FIRST=$(head -n 1 "$TMP/raw_list")
SEP=$'\t'
if [[ "$FIRST" != *$'\t'* ]]; then
  if [[ "$FIRST" == *';'* ]]; then SEP=';'
  elif [[ "$FIRST" == *','* ]]; then SEP=','
  fi
fi
COL=$(head -n 1 "$TMP/raw_list" | awk -F"$SEP" '
  { for (i = 1; i <= NF; i++) {
      c = tolower($i); gsub(/^[ "]+|[ "]+$/, "", c)
      if (c ~ /^(species|scientific_name|scientific name|taxon|species_name)$/ || c ~ /^esp/) { print i; exit }
  } }')
if [[ -n "$COL" ]]; then
  echo "  table detected: species taken from column $COL" >&2
  tail -n +2 "$TMP/raw_list" | awk -F"$SEP" -v c="$COL" '{gsub(/"/,"",$c); print $c}' > "$TMP/raw_list2"
  mv "$TMP/raw_list2" "$TMP/raw_list"
fi
sed 's/^[[:space:]]*//; s/[[:space:]]*$//' "$TMP/raw_list" | grep -v '^$' | grep -v '^#' | sort -u > "$TMP/list"
cp "$TMP/list" "$TMP/table"

for S in "${STAGES[@]}"; do
  IFS='|' read -r LABEL FILE MODE <<<"$S"
  if [[ ! -s "$FILE" ]]; then
    awk -v OFS='\t' '{print $0, "-"}' "$TMP/table" > "$TMP/t2"; mv "$TMP/t2" "$TMP/table"
    echo "  (missing: $FILE)" >&2; continue
  fi
  echo "  counting in $FILE" >&2
  case "$MODE" in
    bold_raw)
      awk -F'\t' '
        FILENAME == ARGV[1] { w[$0] = 1; next }
        FNR > 1 && ($25 in w) { c[$25]++ }
        END { for (k in c) print k "\t" c[k] }' "$TMP/list" "$FILE" ;;
    ncbi_raw)
      awk '
        FILENAME == ARGV[1] { w[$0] = 1; next }
        /^>/ {
          h = $0; sub(/^>[^ ]* /, "", h); sub(/^UNVERIFIED: /, "", h)
          split(h, a, " "); n = a[1] " " a[2]
          if (n in w) c[n]++
        }
        END { for (k in c) print k "\t" c[k] }' "$TMP/list" "$FILE" ;;
    crabs)
      awk -F'\t' '
        FILENAME == ARGV[1] { w[$0] = 1; next }
        { if ($2 in w) c[$2]++; else if ($10 in w) c[$10]++ }
        END { for (k in c) print k "\t" c[k] }' "$TMP/list" "$FILE" ;;
  esac > "$TMP/counts"
  awk -F'\t' -v OFS='\t' 'FILENAME==ARGV[1]{c[$1]=$2; next}{print $0, ($1 in c ? c[$1] : 0)}' "$TMP/counts" "$TMP/table" > "$TMP/t2"
  mv "$TMP/t2" "$TMP/table"
done

# Sequences imported WITHOUT a lineage (taxon name unknown to NCBI): CRABS keeps
# them with "NA" ranks, then --no-species-id / --rank-na remove them
cat imported/bold_all.txt imported/ncbi_all.txt 2>/dev/null | \
  awk -F'\t' '
    FILENAME == ARGV[1] { w[$0] = 1; next }
    ($2 in w) && $3 == "NA" {
      na = 1; for (i = 4; i < NF; i++) if ($i != "NA") na = 0
      if (na) c[$2]++
    }
    END { for (k in c) print k "\t" c[k] }' "$TMP/list" - > "$TMP/nolin"
awk -F'\t' -v OFS='\t' 'FILENAME==ARGV[1]{c[$1]=$2; next}{print $0, ($1 in c ? c[$1] : 0)}' "$TMP/nolin" "$TMP/table" > "$TMP/t2"
mv "$TMP/t2" "$TMP/table"

# Is the name in the NCBI taxonomy? (unresolved names get an all-NA lineage at
# import and are then removed by --no-species-id / --rank-na)
TAXO="$(up taxo)"; [[ -s "$TAXO/names.dmp" ]] || TAXO="$(up crabs_coi)/taxo"   # shared CRABS/taxo (older layout: crabs_coi/taxo)
if [[ -s "$TAXO/names.dmp" ]]; then
  echo "  checking names in $TAXO/names.dmp" >&2
  awk -F'\t[|]\t' '
    FILENAME == ARGV[1] { w[$0] = 1; next }
    ($2 in w) {
      cls = $4; sub(/\t[|]$/, "", cls)
      if (cls == "scientific name") s[$2] = 1
      else if (!($2 in s)) y[$2] = cls
    }
    END {
      for (k in s) print k "\tyes"
      for (k in y) if (!(k in s)) print k "\t" y[k]
    }' \
    "$TMP/list" "$TAXO/names.dmp" > "$TMP/tax"
else
  : > "$TMP/tax"
fi
awk -F'\t' -v OFS='\t' 'FILENAME==ARGV[1]{t[$1]=$2; next}{print $0, ($1 in t ? t[$1] : "NOT FOUND")}' "$TMP/tax" "$TMP/table" > "$TMP/t2"
mv "$TMP/t2" "$TMP/table"

# Step where each species disappears
HEADER="species"; for S in "${STAGES[@]}"; do HEADER="$HEADER	${S%%|*}"; done
HEADER="$HEADER	import_without_lineage	in_ncbi_taxonomy	lost_at"
{
  echo "$HEADER"
  awk -F'\t' -v OFS='\t' -v hdr="$HEADER" -v ns="${#STAGES[@]}" 'BEGIN{n=split(hdr,h,"\t")}
  { lost="kept"; fin=1+ns; nolin=fin+1
    # raw sources: columns 2-3, imports: 4-5, then one column per later step
    if ($fin+0 > 0) lost="kept"
    else if ($2+0==0 && $3+0==0) lost="not downloaded (no record in raw/bold_all.tsv or raw/ncbi_all.fasta)"
    else if ($4+0==0 && $5+0==0) lost="import (record not imported)"
    else if ($nolin+0 > 0 && $nolin+0 >= $4+$5) lost="import: name unknown to NCBI taxonomy (empty lineage, then removed by filters)"
    else for (i=6; i<=fin; i++) if ($i!="-" && $i+0==0) { lost=h[i]; break }
    print $0, lost }' "$TMP/table"
} > "$OUT"

# Names under which the GENUS appears in the final database: catches provisional
# GenBank names such as "Ballia sp. 1callitricha" for Ballia callitricha
if [[ -n "${BASE:-}" && -s "$BASE" ]]; then
  echo
  echo "=== Same genus in the final database (other names it may be filed under) ==="
  awk -F'\t' 'FILENAME==ARGV[1]{split($0,a," "); g[a[1]]=1; next} { split($2,b," "); if (b[1] in g) c[$2]++ }
       END{ for(k in c) print "  " k "\t" c[k] " seq" }' "$TMP/list" "$BASE" | sort
fi

echo
echo "=== Species trace (counts of sequences per step) ==="
column -t -s $'\t' "$OUT" 2>/dev/null || cat "$OUT"
echo
echo "Summary of where species are lost:"
tail -n +2 "$OUT" | awk -F'\t' '{print $NF}' | sort | uniq -c | sort -rn

# Origin of a class (default Insecta) in the final database
if [[ -n "${BASE:-}" && -s "$BASE" ]]; then
  echo
  echo "=== $CHECK_CLASS in the final database ==="
  awk -F'\t' -v k="$CHECK_CLASS" '$6==k' "$BASE" > "$TMP/cls"
  N=$(wc -l < "$TMP/cls" | tr -d ' ')
  echo "  $N sequences of class $CHECK_CLASS"
  if (( N > 0 )); then
    cut -f1 "$TMP/cls" | sort -u > "$TMP/cls_ids"
    for SRC in bold ncbi; do
      F="imported/${SRC}_all.txt"; [[ -s "$F" ]] || continue
      M=$(cut -f1 "$F" | sort -u | comm -12 - "$TMP/cls_ids" | wc -l | tr -d ' ')
      echo "  from $(echo "$SRC" | tr a-z A-Z): $M"
    done
    echo "  most frequent orders:"
    cut -f7 "$TMP/cls" | sort | uniq -c | sort -rn | head -n 5 | sed 's/^/    /'
    echo "  examples (id, name, phylum, class, order):"
    cut -f1,2,5,6,7 "$TMP/cls" | head -n 5 | sed 's/^/    /'
  fi
fi
echo
echo "Table written to $OUT"