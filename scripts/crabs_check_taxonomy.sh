#!/usr/bin/env bash
# =============================================================================
#  crabs_check_taxonomy.sh – check the taxonomy of CRABS reference databases
#  for genus homonyms (read only: the databases are not modified)
# =============================================================================
#  Problem  : CRABS gives a sequence its NCBI lineage from its NAME when the
#             record has no usable NCBI taxid (BOLD records, species unknown to
#             NCBI, rescued lineages). When the genus name also exists in another
#             kingdom or phylum (e.g. Trentepohlia: green alga and crane fly;
#             Lessonia: kelp and bird; Alaria: kelp and trematode; Urospora:
#             green alga and gregarine), the wrong genus may be picked and the
#             sequence ends up in the wrong phylum.
#  Checks   : for each sequence of the final database (results/<base>.tsv)
#             1. off target: its NCBI lineage does not reach any of the target
#                groups of the database (bold_groups + ncbi_groups of
#                metadata/parameters.tsv). Not applied to databases built with
#                SILVA (all eukaryotes on purpose, e.g. 18S).
#             2. BOLD mismatch: for a BOLD record, the phylum given by BOLD
#                (original BOLD files listed in metadata/sources.tsv) differs from
#                the phylum of its NCBI lineage. Catches homonyms inside the
#                target groups too.
#             For each flagged sequence, the NCBI genera of the same name are
#             listed, and a correction is proposed when exactly one of them fits
#             (same phylum as BOLD, or inside the target groups).
#  Output   : <database>/checks/taxonomy_check_<date>_<time>/ (a new folder at
#             every run: earlier checks are never overwritten), and one line per
#             run added to <database>/checks/taxonomy_check_history_crabs.tsv
#             (date, database and its checksum, counts) to follow the changes:
#               taxonomy_check_summary_crabs.txt     what was checked, counts (trace)
#               taxonomy_check_sequences_crabs.tsv   one line per flagged sequence
#               taxonomy_check_genera_crabs.tsv      one line per genus concerned
#               taxonomy_check_crabs.log             screen output
#               taxonomy_manual_template_crabs.tsv   the genera to decide by hand, ready to
#                 fill (action = keep | remove), then saved in the database folder as
#                 taxonomy_manual_crabs.tsv (read by the database script)
#  Usage    : put this script in CRABS/scripts/ (or in CRABS/), then from CRABS:
#               bash scripts/crabs_check_taxonomy.sh                       → every database found
#               bash scripts/crabs_check_taxonomy.sh crabs_coi             → one database
#               bash scripts/crabs_check_taxonomy.sh crabs_18S/18S_forward
# =============================================================================
set -uo pipefail
SCRIPT_VERSION="1.7.1 (2026-10-28, shared CRABS/taxo, CRABS/scripts, database folders at any depth)"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# CRABS folder: parent of CRABS/scripts, else the folder of this script (or set CRABS_HOME)
if [[ -z "${CRABS_HOME:-}" ]]; then
  if [[ "$(basename "$SCRIPT_DIR")" == scripts ]]; then CRABS_HOME="$(dirname "$SCRIPT_DIR")"; else CRABS_HOME="$SCRIPT_DIR"; fi
fi
# upfrom DIR NAME: NAME in DIR or one of its parents (e.g. the shared CRABS/taxo)
upfrom() { local d=$1; while :; do [[ -e "$d/$2" ]] && { echo "$d/$2"; return 0; }; [[ "$d" == / ]] && return 0; d=$(dirname "$d"); done; }
STAMP=$(date '+%Y-%m-%d_%H%M%S')

# Databases to check
DBS=()
if [[ $# -ge 1 ]]; then
  for a in "$@"; do
    for d in "$a" "$CRABS_HOME/$a" "$SCRIPT_DIR/$a"; do
      if [[ -d "$d/results" ]]; then DBS+=("$(cd "$d" && pwd)"); break; fi
    done
  done
else
  for root in "$CRABS_HOME" "$SCRIPT_DIR"; do
    for d in "$root"/crabs_*/ "$root"/crabs_*/*/; do   # also sub-folders, e.g. crabs_18S/18S_aligned
      [[ -d "${d}results" ]] && ls "${d}results/"*.tsv >/dev/null 2>&1 && DBS+=("$(cd "$d" && pwd)")
    done
  done
fi
# unique
if (( ${#DBS[@]} )); then
  UNIQ=()
  while IFS= read -r d; do UNIQ+=("$d"); done < <(printf '%s\n' "${DBS[@]}" | awk '!seen[$0]++')
  DBS=("${UNIQ[@]}")
fi
(( ${#DBS[@]} )) || { echo "No CRABS database found (crabs_*/results/*.tsv) in $SCRIPT_DIR or $CRABS_HOME"; exit 1; }
echo "Databases: $(for d in "${DBS[@]}"; do printf '%s ' "$(basename "$d")"; done)"

for DB in "${DBS[@]}"; do
  echo
  echo "===== $(basename "$DB") ====="
  TSVS=$(ls "$DB"/results/*.tsv 2>/dev/null | grep -v insilico_pcr_raw)
  [[ -n "$TSVS" ]] || { echo "  no final .tsv in $DB/results: skipped"; continue; }
  # A folder may hold several databases (e.g. crabs_18s with the 18S forward
  # database): each one is checked; a database whose files carry its name
  # (metadata/<base>_parameters.tsv) uses its own metadata and results folder
  for TSV in $TSVS; do
  BASE=$(basename "$TSV" .tsv)
  PFX=""; [[ -s "$DB/metadata/${BASE}_parameters.tsv" ]] && PFX="${BASE}_"
  (( $(echo "$TSVS" | wc -l) > 1 )) && echo "----- $BASE"
  # NCBI taxonomy: the shared CRABS/taxo (or the one of this database), else the one of crabs_coi
  TAXO=""
  for t in "$(upfrom "$DB" taxo)" "$CRABS_HOME/taxo" "$(upfrom "$DB" crabs_coi)/taxo" "$CRABS_HOME/crabs_coi/taxo"; do
    [[ -s "$t/nodes.dmp" && -s "$t/names.dmp" ]] && { TAXO="$(cd "$t" && pwd)"; break; }
  done
  [[ -n "$TAXO" ]] || { echo "  NCBI taxonomy (taxo/nodes.dmp) not found: skipped"; continue; }
  OUT="$DB/checks/${PFX}taxonomy_check_$STAMP"
  n=1; while [[ -e "$OUT" ]]; do OUT="$DB/checks/${PFX}taxonomy_check_${STAMP}_$n"; n=$((n + 1)); done
  mkdir -p "$OUT"
  python3 - "$DB" "$TSV" "$TAXO" "$OUT" "$SCRIPT_VERSION" "$PFX" <<'PYEOF' 2>&1 | tee "$OUT/taxonomy_check_crabs.log"
import collections, csv, datetime, os, sys
csv.field_size_limit(10**9)
DB, TSV, TAXO, OUT, VERSION, PFX = sys.argv[1:7]
os.chdir(DB)
META = lambda f: f"metadata/{PFX}{f}" if os.path.exists(f"metadata/{PFX}{f}") else f"metadata/{f}"

# ---------------- parameters of the database --------------------------------
P = {}
if os.path.exists(META("parameters.tsv")):
    for r in csv.reader(open(META("parameters.tsv"), encoding="utf-8"), delimiter="\t"):
        if len(r) > 1: P[r[0]] = r[1]
# target groups not named in the download lists: algae under their NCBI names;
# Branchiura and Pentastomida (NCBI: Ichthyostraca) and Mystacocarida, downloaded
# through "Crustacea" in NCBI
EXTRA = ["Ochrophyta", "Bacillariophyta", "Phaeophyceae", "Ichthyostraca", "Mystacocarida"]
C2P = dict(kv.split("=", 1) for kv in P.get("phylum_from_class", "").split(";") if "=" in kv)
targets = set(P.get("bold_groups", "").split() + P.get("ncbi_groups", "").split() + EXTRA)
silva = P.get("silva", "no") not in ("", "no")
check_target = bool(targets - set(EXTRA)) and not silva

# phylum names that mean the same thing (BOLD vs NCBI vs Adl et al. 2019)
SAME = {"Bacillariophyta": "Ochrophyta", "Heterokontophyta": "Ochrophyta", "Ochrophyta": "Ochrophyta",
        "Myzozoa": "Dinoflagellata", "Dinophyta": "Dinoflagellata", "Pyrrophycophyta": "Dinoflagellata",
        "Magnoliophyta": "Streptophyta", "Tracheophyta": "Streptophyta", "Streptophyta": "Streptophyta",
        "Charophyta": "Streptophyta", "Bryophyta": "Streptophyta", "Cryptophycophyta": "Cryptophyta",
        "Prasinodermophyta": "Chlorophyta"}
norm = lambda x: SAME.get(x, x) if x not in ("", "NA", None) else ""

# manual decisions of the database script (taxonomy_manual_crabs.tsv, if any)
# files named in the metadata may have moved with their folder (e.g. crabs_18S
# split into 18S_aligned / 18S_forward): looked for here, in the source
# database, in crabs_coi/raw and in the sibling folders
import glob
def up(name):
    d = os.getcwd()
    while True:
        if os.path.exists(os.path.join(d, name)): return os.path.join(d, name)
        if d == os.path.dirname(d): return ""
        d = os.path.dirname(d)
def findf(f):
    b = os.path.basename(f)
    cands = [f, os.path.join(P.get("source_database", ""), f)]
    cands += [os.path.join(d, "raw", b) for d in [up("crabs_coi")] + sorted(glob.glob("../*"))] + [os.path.join(d, b) for d in sorted(glob.glob("../*"))]
    for c in cands:
        if c and os.path.exists(c): return c
    return f
manual = {}
f_manual = findf("taxonomy_manual_crabs.tsv") if PFX else "taxonomy_manual_crabs.tsv"
if os.path.exists(f_manual):
    for r in csv.DictReader((l for l in open(f_manual, encoding="utf-8") if l.strip() and not l.startswith("#")), delimiter="\t"):
        if r.get("taxon") and r.get("action"): manual[r["taxon"].strip()] = r["action"].strip().lower()

# ---------------- NCBI taxonomy ----------------------------------------------
print(f"Reading the NCBI taxonomy: {TAXO}")
parent, rank, name = {}, {}, {}
genus_ids = collections.defaultdict(list)
phylum_names, phylum_tax = set(), {}
SAME_KINGDOM_OK = {"Fungi"}    # fungal phyla are delimited differently by NCBI, SILVA and BOLD
for line in open(f"{TAXO}/nodes.dmp", encoding="utf-8", errors="ignore"):
    p = line.split("\t|\t")
    try: t = int(p[0]); parent[t] = int(p[1]); rank[t] = p[2].strip("\t|\n ")
    except (ValueError, IndexError): pass
for line in open(f"{TAXO}/names.dmp", encoding="utf-8", errors="ignore"):
    p = line.rstrip("\t|\n").split("\t|\t")
    if len(p) >= 4 and p[3] == "scientific name":
        try: t = int(p[0])
        except ValueError: continue
        name[t] = p[1]
        if rank.get(t) == "genus": genus_ids[p[1]].append(t)
        if rank.get(t) == "phylum": phylum_names.add(p[1]); phylum_tax.setdefault(p[1], t)
taxo_date = datetime.datetime.fromtimestamp(os.path.getmtime(f"{TAXO}/nodes.dmp")).strftime("%Y-%m-%d")

memo_t, memo_r = {}, {}
def on_target(t):
    path = []
    for _ in range(200):
        if t in memo_t: res = memo_t[t]; break
        path.append(t)
        if name.get(t) in targets: res = True; break
        if t == parent.get(t, t): res = False; break
        t = parent[t]
    else: res = False
    for x in path: memo_t[x] = res
    return res
def ranks_of(t):
    if t in memo_r: return memo_r[t]
    out, s = {}, t
    for _ in range(200):
        if rank.get(s) in ("phylum", "class", "order", "family"): out.setdefault(rank[s], name.get(s, ""))
        if s == parent.get(s, s): break
        s = parent[s]
    memo_r[t] = out
    return out

# ---------------- final database ---------------------------------------------
rows = list(csv.DictReader(open(TSV, encoding="utf-8"), delimiter="\t"))
ids = {r["seq_id"] for r in rows}
print(f"Database: {TSV} ({len(rows)} sequences)")

# ---------------- original BOLD lineage --------------------------------------
bold = {}
bold_files = []
if os.path.exists(META("sources.tsv")):
    for r in csv.DictReader(open(META("sources.tsv"), encoding="utf-8"), delimiter="\t"):
        if "BOLD" in r.get("source", "") and r.get("file"):
            bold_files.append(findf(r["file"]))
def col(h, *names):
    low = [x.strip().lower() for x in h]
    for n in names:
        if n in low: return low.index(n)
    return None
for f in bold_files:
    if not os.path.exists(f):
        print(f"  BOLD file not found (skipped): {f}"); continue
    with open(f, encoding="utf-8", errors="ignore") as fi:
        h = fi.readline().rstrip("\n").split("\t")
        i_rec, i_pid, i_mk = col(h, "record_id"), col(h, "processid"), col(h, "marker_code", "markercode")
        i_ph, i_cl, i_fa, i_ge = col(h, "phylum"), col(h, "class"), col(h, "family"), col(h, "genus")
        if i_ph is None: continue
        group = os.path.basename(f).replace("bold_", "").replace(".tsv", "")
        n = 0
        for line in fi:
            p = line.rstrip("\n").split("\t")
            keys = []
            if i_rec is not None and i_rec < len(p): keys.append(p[i_rec])
            if i_pid is not None and i_pid < len(p):
                keys.append(p[i_pid])
                if i_mk is not None and i_mk < len(p): keys.append(f"{p[i_pid]}.{p[i_mk]}")
            for k in keys:
                if k in ids:
                    g = lambda i: p[i] if i is not None and i < len(p) else ""
                    bold[k] = (g(i_ph), g(i_cl), g(i_fa), g(i_ge), group); n += 1
                    break
    print(f"  BOLD lineage read from {f}")
print(f"  {len(bold)} sequences of the database found in the BOLD files")
n_bold = len(bold)
# SILVA (18S): phylum = first name of the SILVA lineage (header) that is a phylum
# in the NCBI taxonomy; CRABS identifier = accession before the first dot
silva_files = []
if os.path.exists(META("sources.tsv")):
    for r in csv.DictReader(open(META("sources.tsv"), encoding="utf-8"), delimiter="\t"):
        if r.get("source", "").startswith("SILVA") and r.get("file"): silva_files.append(findf(r["file"]))
n_silva = 0
for f in silva_files:
    if not os.path.exists(f):
        print(f"  SILVA file not found (skipped): {f}"); continue
    with open(f, encoding="utf-8", errors="ignore") as fi:
        for line in fi:
            if not line.startswith(">"): continue
            sid = line[1:].split(".")[0]
            if sid not in ids or sid in bold: continue
            path = line.rstrip("\n").split(" ", 1)[1].split(";")[:-1] if " " in line else []
            for x in path:
                x = x.split(" (")[0].strip()
                if x in phylum_names:
                    bold[sid] = (x, "", "", "", "SILVA"); n_silva += 1; break
    print(f"  SILVA lineage read from {f}: {n_silva} database sequences with a phylum")

# ---------------- checks -----------------------------------------------------
flag = []
for r in rows:
    try: t = int(r["taxid"])
    except (ValueError, KeyError): t = None
    ph_db = r.get("phylum", "")
    issues = []
    if check_target and not (t is not None and on_target(t)):
        issues.append("off_target")
    b = bold.get(r["seq_id"])
    def kingdom_of(ph):
        t = phylum_tax.get(ph)
        for _ in range(200):
            if t is None: return ""
            if rank.get(t) == "kingdom": return name.get(t, "")
            if t == parent.get(t, t): return ""
            t = parent[t]
        return ""
    if b and norm(b[0]) and norm(ph_db) and norm(b[0]) != norm(ph_db) \
       and not (kingdom_of(b[0]) in SAME_KINGDOM_OK and kingdom_of(b[0]) == kingdom_of(ph_db)):
        issues.append("bold_phylum_mismatch")
    if not issues: continue
    genus = r["scientific_name"].replace("UNVERIFIED: ", "").split(" ")[0]
    cands = genus_ids.get(genus, [])
    desc = [f"{g}:{ranks_of(g).get('phylum', '?')}/{ranks_of(g).get('class', '?')}" for g in cands]
    # proposed correction: the genus of the same name that fits BOLD, else the targets
    ph = lambda g: norm(ranks_of(g).get("phylum") or C2P.get(ranks_of(g).get("class", ""), ""))
    fit = [g for g in cands if b and ph(g) == norm(b[0])] if b else []
    if not fit and check_target: fit = [g for g in cands if on_target(g)]
    m = manual.get(r["seq_id"]) or manual.get(genus)
    if m:
        verdict = f"manual decision: {m}"
        prop = ("", "", "", "", "")
    elif len(fit) == 1 and str(fit[0]) == r.get("taxid", ""):
        verdict = "kept: NCBI lineage plausible, BOLD disagrees (check)"
        prop = ("", "", "", "", "")
    elif len(fit) == 1:
        rk = ranks_of(fit[0])
        verdict = "homonym: correct lineage"
        prop = (str(fit[0]), rk.get("phylum", ""), rk.get("class", ""), rk.get("order", ""), rk.get("family", ""))
    else:
        verdict = "no fitting genus: remove or check" if not fit else "several fitting genera: check"
        prop = ("", "", "", "", "")
    flag.append({"seq_id": r["seq_id"], "source": (b[4] if b and b[4] == "SILVA" else "BOLD") if b else "NCBI/other", "name": r["scientific_name"],
                 "genus": genus, "taxid": r.get("taxid", ""), "db_phylum": ph_db, "db_class": r.get("class", ""),
                 "db_family": r.get("family", ""),
                 "bold_group": b[4] if b else "", "bold_phylum": b[0] if b else "", "bold_class": b[1] if b else "",
                 "bold_family": b[2] if b else "", "issue": "+".join(issues), "verdict": verdict,
                 "proposed_taxid": prop[0], "proposed_phylum": prop[1], "proposed_class": prop[2],
                 "proposed_order": prop[3], "proposed_family": prop[4],
                 "ncbi_genera_of_that_name": "; ".join(desc)})

# ---------------- outputs ----------------------------------------------------
cols = list(flag[0].keys()) if flag else ["seq_id", "source", "name", "genus", "taxid", "db_phylum", "issue", "verdict"]
with open(f"{OUT}/taxonomy_check_sequences_crabs.tsv", "w", encoding="utf-8") as fo:
    w = csv.DictWriter(fo, fieldnames=cols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(flag)
gen = collections.OrderedDict()
for x in sorted(flag, key=lambda x: (x["db_phylum"], x["genus"])):
    g = gen.setdefault(x["genus"], {"genus": x["genus"], "n_sequences": 0, "species": set(), "db_phylum": set(),
                                    "bold_phylum": set(), "proposed_phylum": set(), "verdict": set(),
                                    "ncbi_genera_of_that_name": x["ncbi_genera_of_that_name"]})
    g["n_sequences"] += 1; g["species"].add(x["name"]); g["db_phylum"].add(x["db_phylum"])
    if x["bold_phylum"]: g["bold_phylum"].add(x["bold_phylum"])
    if x["proposed_phylum"]: g["proposed_phylum"].add(x["proposed_phylum"])
    g["verdict"].add(x["verdict"])
with open(f"{OUT}/taxonomy_check_genera_crabs.tsv", "w", encoding="utf-8") as fo:
    fo.write("genus\tn_sequences\tn_names\tdb_phylum\tbold_phylum\tproposed_phylum\tverdict\tncbi_genera_of_that_name\n")
    for g in gen.values():
        fo.write("\t".join([g["genus"], str(g["n_sequences"]), str(len(g["species"])),
                            ", ".join(sorted(g["db_phylum"])), ", ".join(sorted(g["bold_phylum"])),
                            ", ".join(sorted(g["proposed_phylum"])), ", ".join(sorted(g["verdict"])),
                            g["ncbi_genera_of_that_name"]]) + "\n")
# template of manual decisions: genera that the automatic correction cannot
# settle (no single fitting genus), not decided yet
undecided = collections.OrderedDict()
for x in flag:
    if x["verdict"].startswith(("homonym", "manual")): continue
    u = undecided.setdefault(x["genus"], {"n": 0, "db": set(), "bold": set(), "verdict": x["verdict"], "ncbi": x["ncbi_genera_of_that_name"]})
    u["n"] += 1; u["db"].add(x["db_phylum"] or "NA")
    if x["bold_phylum"]: u["bold"].add(x["bold_phylum"])
f_tpl = f"{OUT}/taxonomy_manual_template_crabs.tsv"
with open(f_tpl, "w", encoding="utf-8") as fo:
    fo.write("# Manual taxonomy decisions (step 10_taxonomy_fix of the database script).\n"
             "# Fill 'action' with keep or remove (and adapt 'reason'), delete the lines you do not want to decide,\n"
             "# then add the lines to taxonomy_manual_crabs.tsv in the database folder (create it from this file the first time).\n"
             "# taxon = genus (applies only to its flagged sequences) or seq_id.\n")
    fo.write("taxon\taction\treason\tdecided\n")
    for g, u in undecided.items():
        ev = (f"{u['n']} seq; database phylum {', '.join(sorted(u['db']))}"
              + (f"; BOLD phylum {', '.join(sorted(u['bold']))}" if u["bold"] else "")
              + f"; {u['verdict']}" + (f"; NCBI genera of that name: {u['ncbi']}" if u["ncbi"] else ""))
        fo.write(f"{g}\t\t{ev}\t{datetime.date.today().isoformat()}\n")
cnt_issue = collections.Counter(x["issue"] for x in flag)
cnt_verdict = collections.Counter(x["verdict"] for x in flag)
cnt_ph = collections.Counter(x["db_phylum"] for x in flag)
summary = [
  f"Taxonomy check of a CRABS reference database – crabs_check_taxonomy.sh {VERSION}",
  f"Date              : {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}",
  f"Database          : {TSV} ({len(rows)} sequences)",
  f"NCBI taxonomy     : {TAXO} (nodes.dmp of {taxo_date})",
  f"BOLD files        : {len(bold_files)} ({n_bold} database sequences found in them)",
  f"SILVA lineage     : {len(silva_files)} file(s) ({n_silva} database sequences with a phylum in their SILVA lineage)",
  f"Target groups     : {' '.join(sorted(targets - set(EXTRA)))}" + ("" if check_target else "  [off-target check NOT applied: database built with SILVA]"),
  "",
  f"Flagged sequences : {len(flag)} ({100 * len(flag) / max(1, len(rows)):.2f} %), {len(gen)} genera",
  "  by issue        : " + ", ".join(f"{k} {v}" for k, v in cnt_issue.most_common()),
  "  by verdict      : " + ", ".join(f"{k} {v}" for k, v in cnt_verdict.most_common()),
  "  by phylum in the database: " + ", ".join(f"{k or 'NA'} {v}" for k, v in cnt_ph.most_common()),
  "",
  "Issues: off_target = NCBI lineage outside the target groups; bold_phylum_mismatch = phylum",
  "of the NCBI lineage differs from the phylum given by BOLD or by the SILVA lineage (equivalent names, e.g.",
  "Bacillariophyta/Ochrophyta, are not counted as differences).",
  "Verdict 'homonym: correct lineage' = exactly one NCBI genus of that name fits (BOLD phylum, or",
  "the target groups): proposed_* columns. Otherwise check by hand.",
  "",
  "Files: taxonomy_check_sequences_crabs.tsv (one line per sequence),",
  "       taxonomy_check_genera_crabs.tsv (one line per genus)."]
open(f"{OUT}/taxonomy_check_summary_crabs.txt", "w", encoding="utf-8").write("\n".join(summary) + "\n")
print()
print("\n".join(summary[6:11]))
print()
for g in list(gen.values())[:25]:
    print(f"  {g['genus']:<22} {g['n_sequences']:>5} seq  in {', '.join(sorted(g['db_phylum'])) or 'NA':<16}"
          f" BOLD/SILVA: {', '.join(sorted(g['bold_phylum'])) or '-':<14} → {', '.join(sorted(g['proposed_phylum'])) or '?'}")
if len(gen) > 25: print(f"  … {len(gen) - 25} more genera in taxonomy_check_genera_crabs.tsv")
# history: one line per run, never rewritten
import hashlib
h = hashlib.sha256()
with open(TSV, "rb") as fb:
    for chunk in iter(lambda: fb.read(1 << 20), b""): h.update(chunk)
hist = os.path.join(os.path.dirname(OUT), "taxonomy_check_history_crabs.tsv")
new = not os.path.exists(hist)
with open(hist, "a", encoding="utf-8") as fh:
    if new:
        fh.write("date\tscript_version\tdatabase\tdatabase_sha256\tdatabase_modified\tn_sequences\tn_flagged\t"
                 "n_genera\tn_correctable\tn_to_remove_or_check\tncbi_taxonomy_date\tresults_folder\n")
    fh.write("\t".join(map(str, [datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), VERSION.split()[0],
             os.path.relpath(TSV, DB), h.hexdigest()[:16],
             datetime.datetime.fromtimestamp(os.path.getmtime(TSV)).strftime("%Y-%m-%d %H:%M"),
             len(rows), len(flag), len(gen), cnt_verdict.get("homonym: correct lineage", 0),
             len(flag) - cnt_verdict.get("homonym: correct lineage", 0), taxo_date,
             os.path.basename(OUT)])) + "\n")
print(f"\nResults: {OUT}/")
if undecided:
    print(f"To decide by hand: {len(undecided)} genera → {f_tpl}")
    print("  fill 'action' (keep / remove), then save the lines in "
          f"{os.path.join(DB, 'taxonomy_manual_crabs.tsv')} and re-run the database script")
else:
    print("Nothing to decide by hand.")
print(f"History: {hist}")
PYEOF
  done
done