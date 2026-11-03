#!/usr/bin/env bash
# =============================================================================
#  crabs_18S_split.sh – split the 18S folder into 18S_aligned/ and 18S_forward/
# =============================================================================
#  Before : crabs_18S/  (18S database + 18S forward database, forward files
#                        named <forward base>_…)
#  After  : crabs_18S/18S_aligned/   everything of the 18S database (crabs_18S.sh)
#           crabs_18S/18S_forward/   everything of the forward database
#                                    (crabs_18S_fwd.sh, files <forward base>_…)
#  Files are MOVED (mv), so their dates are kept and no finished step is run
#  again (needs crabs_18S.sh 1.5.1 and crabs_18S_fwd.sh 1.3.1 or later).
#  Usage : cd …/CRABS/crabs_18S && bash ../scripts/crabs_18S_split.sh
#          (add "-n" to only print what would be moved)
# =============================================================================
set -euo pipefail
DRY=false; [[ "${1:-}" == "-n" ]] && DRY=true
run() { if $DRY; then echo "  $*"; else "$@"; fi; }
[[ -d results ]] || { echo "Run this in the 18S folder (crabs_18S/), which holds results/"; exit 1; }
[[ -e 18S_aligned || -e 18S_forward ]] && { echo "18S_aligned/ or 18S_forward/ already exists: nothing done"; exit 1; }
if pgrep -f "crabs_18S(_fwd)?\.sh" >/dev/null 2>&1; then echo "A crabs_18S script is running: wait for it to finish"; exit 1; fi
# forward database: the final table named *fwd*
FWD_TSV=$(ls results/*fwd*.tsv 2>/dev/null | head -n 1 || true)
[[ -n "$FWD_TSV" ]] || { echo "No forward database (results/*fwd*.tsv) found"; exit 1; }
B=$(basename "$FWD_TSV" .tsv)
echo "Forward database: $B"
run mkdir -p 18S_forward/results 18S_forward/metadata/.steps 18S_forward/steps 18S_forward/logs 18S_forward/checks
shopt -s nullglob dotglob
# 1. forward files (named after the forward database) and its script, log, figures
for f in results/"$B"* results/"$B"_* steps/"$B"_* logs/"$B"_* metadata/"$B"_* checks/"$B"_* \
         metadata/.steps/*"$B"*; do
  [[ -e "$f" ]] && run mv "$f" "18S_forward/$f"
done
for f in crabs_18S_fwd*.sh crabs_18S_fwd*.log figures_"$B"; do
  [[ -e "$f" ]] && run mv "$f" 18S_forward/
done
# 2. everything else = the 18S database
run mkdir -p 18S_aligned
for f in *; do
  [[ "$f" == 18S_aligned || "$f" == 18S_forward ]] && continue
  run mv "$f" 18S_aligned/
done
$DRY && exit 0
echo "Done:"
echo "  18S_aligned/: $(ls 18S_aligned | tr '\n' ' ')"
echo "  18S_forward/: $(ls 18S_forward | tr '\n' ' ')"
echo "Check: cd 18S_aligned && bash crabs_18S.sh   then   cd ../18S_forward && bash crabs_18S_fwd.sh"
echo "       (both must answer \"already done\" at every step)"
