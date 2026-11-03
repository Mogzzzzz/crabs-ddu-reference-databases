#!/usr/bin/env bash
# =============================================================================
#  crabs_check_taxa.sh – are these taxa in the reference databases?
# =============================================================================
#  For each taxon of a list (any rank: species, genus, family, order, class,
#  phylum), counts the sequences and species in every final database of a CRABS
#  folder (e.g. crabs_18S: the 18S and the 18S forward databases). A taxon
#  matches a sequence when it is its species name or one of the names of its
#  lineage. When a taxon is missing from a database, the last step where it was
#  still present is given (imported data, then each step of that database).
#
#  Usage : bash crabs_check_taxa.sh <database folder> <taxa file>
#            taxa file: one taxon per line (first column if tab-separated;
#            lines starting with # or "Taxon" ignored), e.g.
#              Clathria
#              Mitrocomella frigida
#              Hexacorallia
#  Output: on screen, and <database folder>/checks/taxa_check_<date>_<time>.tsv
# =============================================================================
set -uo pipefail
[[ $# -eq 2 ]] || { echo "Usage: bash $0 <database folder> <taxa file>"; exit 1; }
DB="$1"; LIST="$2"
[[ -d "$DB/results" ]] || { echo "Not a CRABS database folder: $DB"; exit 1; }
DB=$(cd "$DB" && pwd)
[[ -s "$LIST" ]] || { echo "Taxa file not found: $LIST"; exit 1; }
LIST=$(cd "$(dirname "$LIST")" && pwd)/$(basename "$LIST")
mkdir -p "$DB/checks"
OUT="$DB/checks/taxa_check_$(date '+%Y-%m-%d_%H%M%S').tsv"
python3 - "$DB" "$LIST" "$OUT" <<'PYEOF'
import collections, glob, os, sys
db, lst, out = sys.argv[1:4]
os.chdir(db)
taxa = []
for line in open(lst, encoding="utf-8", errors="ignore"):
    t = line.rstrip("\n").split("\t")[0].strip()
    if t and not t.startswith("#") and t.lower() != "taxon" and t not in taxa: taxa.append(t)
# NCBI taxonomy (for ranks not kept in the database, e.g. subclass Hexacorallia)
parent, sci = {}, {}
def up(name):          # NAME in this folder or a parent folder (shared CRABS/taxo)
    d = os.getcwd()
    while True:
        if os.path.exists(os.path.join(d, name)): return os.path.join(d, name)
        if d == os.path.dirname(d): return name
        d = os.path.dirname(d)
for t in [up("taxo"), "taxo", os.path.join(up("crabs_coi"), "taxo")]:   # shared CRABS/taxo first
    if os.path.exists(f"{t}/nodes.dmp") and os.path.exists(f"{t}/names.dmp"):
        print(f"NCBI taxonomy: {t}")
        for line in open(f"{t}/nodes.dmp", encoding="utf-8", errors="ignore"):
            q = line.split("\t|\t")
            try: parent[int(q[0])] = int(q[1])
            except (ValueError, IndexError): pass
        for line in open(f"{t}/names.dmp", encoding="utf-8", errors="ignore"):
            q = line.rstrip("\t|\n").split("\t|\t")
            if len(q) >= 4 and q[3] == "scientific name":
                try: sci[int(q[0])] = q[1]
                except ValueError: pass
        break
memo = {}
def lineage_names(taxid):
    try: t = int(taxid)
    except ValueError: return set()
    if t in memo: return memo[t]
    out, s = set(), t
    for _ in range(200):
        if s in sci: out.add(sci[s])
        if s == parent.get(s, s): break
        s = parent[s]
    memo[t] = out
    return out
bases = sorted(os.path.basename(f)[:-4] for f in glob.glob("results/*.tsv") if "insilico_pcr_raw" not in f)
def hits(path, wanted):
    """{taxon: (n_sequences, set of species)} for a CRABS-format file (no header)"""
    res = {t: [0, set()] for t in wanted}
    with open(path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 5 or p[0] == "seq_id": continue
            names = set(p[3:-1]) | {p[1], p[1].replace("UNVERIFIED: ", "")} | lineage_names(p[2])
            for t in wanted & names:
                res[t][0] += 1; res[t][1].add(p[-2] if p[-2] != "NA" else p[1])
    return res
def steps_of(base):
    """files of a database in pipeline order: imports, steps, final"""
    pfx = base + "_"
    own = sorted(glob.glob(f"steps/{pfx}*.txt"))
    if own:     # database sharing the folder (files prefixed): its input is the folder's 03_filtered_initial
        return ["steps/03_filtered_initial.txt"] + own + [f"results/{base}.txt"]
    shared = sorted(f for f in glob.glob("steps/*.txt") if not any(os.path.basename(f).startswith(b + "_") for b in bases))
    return sorted(glob.glob("imported/*_all.txt")) + shared + [f"results/{base}.txt"]
rows = []
for base in bases:
    final = hits(f"results/{base}.txt", set(taxa))
    missing = {t for t in taxa if final[t][0] == 0}
    last = {}
    if missing:
        for f in steps_of(base)[:-1]:
            if not os.path.exists(f): continue
            h = hits(f, missing)
            for t in missing:
                if h[t][0]: last[t] = (f, h[t][0])
    for t in taxa:
        n, sp = final[t]
        if n: where = ""
        elif t in last: where = f"last seen in {last[t][0]} ({last[t][1]} seq)"
        else: where = "never in the downloaded data"
        rows.append((t, base, n, len(sp), "; ".join(sorted(sp)[:5]) + (" …" if len(sp) > 5 else ""), where))
with open(out, "w", encoding="utf-8") as fo:
    fo.write("taxon\tdatabase\tn_sequences\tn_species\tspecies (first 5)\tif absent\n")
    for r in rows: fo.write("\t".join(map(str, r)) + "\n")
w = max(len(t) for t in taxa)
for base in bases:
    print(f"\n== {base}")
    for t, b, n, ns, sp, where in rows:
        if b != base: continue
        print(f"  {t:<{w}}  {'OK' if n else '--'}  {n:>5} seq  {ns:>4} sp  {sp if n else where}")
print(f"\nTable: {os.path.join(db, out) if not os.path.isabs(out) else out}")
PYEOF
