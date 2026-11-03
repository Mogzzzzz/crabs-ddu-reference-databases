#!/usr/bin/env Rscript
# =============================================================================
#  Figures for a CRABS reference database (output of crabs_coi.sh, crabs_18S.sh,
#  crabs_16S.sh or crabs_rbcL.sh)
# =============================================================================
#  Purpose : describe the content of the reference database (pipeline yield,
#            taxonomic composition, data sources, amplicon lengths, reference
#            depth per species, lineage completeness, downloaded records).
#  Inputs  : <project>/results/<BASE_NAME>.tsv      final database (table)
#            <project>/metadata/counts.tsv          sequences per step
#            <project>/metadata/sources.tsv         downloaded records per group
#            <project>/metadata/parameters.tsv      filter bounds
#            <project>/steps/04_insilico_pcr.txt    amplicons before filtering
#            <project>/imported/bold_all.txt, ncbi_all.txt   source of each record
#  Outputs : <project>/figures/figNN_<name>.png and .pdf (vector)
#            <project>/figures/data/figNN_<name>.csv  data behind each figure
#            <project>/figures/README.md   captions, inputs, software
#            <project>/figures/sessionInfo.txt
#  Location: CRABS/scripts/ (or directly in CRABS/), the database folders in CRABS/
#              CRABS/
#              ├── scripts/crabs_figures.R
#              ├── taxo/                  (NCBI taxonomy, shared)
#              ├── crabs_coi/             (crabs_coi.sh output)
#              └── crabs_16s/, crabs_18s/, …
#  Usage   : Rscript CRABS/scripts/crabs_figures.R [database_folder] [base_name]
#            Without a folder (or with "all"), EVERY database found in CRABS
#            is drawn, each in its own folder: crabs_coi/figures, crabs_16s/figures,
#            crabs_18s/figures, crabs_18s/figures_<forward base>, crabs_rbcl/figures…
#            (a folder counts as a database when its results/ holds the final .tsv).
#            One database only, from CRABS:
#              Rscript scripts/crabs_figures.R crabs_16s    → CRABS/crabs_16s/figures
# 
#            A figure that fails is skipped (message), the others are still made.
#  Figure 1 gives the parameters of each filter. In the per-phylum figures, animals, algae, fungi
#  and other eukaryotes follow each other (animals first) (ALGAE_PHYLA…, see PARAMETERS);
#  taxa listed in EXCLUDE_FROM_FIGURES (18S: fungi) are left out of the figures
#  (not of the database).
#  Figures 2, 5, 9 and 11 follow the visualisation functions of CRABS
#  (Jeunen et al. 2022): diversity, amplicon length, amplification efficiency
#  (primer-binding sites) and completeness table. Figure 12 shows how many
#  species can be told apart with the fragment (identical amplicons).
# =============================================================================

SCRIPT_VERSION <- "2.7 (2026-10-28)"

# ================================ PARAMETERS ==================================
DEFAULT_DB  <- "all"         # "all" = every database folder found; or one folder, e.g. "crabs_rbcl"
CRABS_HOME  <- NA            # folder holding the crabs_* database folders; NA = the CRABS folder
                             # (parent of CRABS/scripts, where this script sits), or the folder of
                             # this script; can also be set outside R: export CRABS_HOME=/path/to/CRABS
args        <- commandArgs(trailingOnly = TRUE)
DB_ARG      <- if (length(args) >= 1) args[1] else DEFAULT_DB
BASE_NAME   <- if (length(args) >= 2) args[2] else NA   # NA = auto-detect in results/
EXPECTED_BP <- NA         # expected amplicon length (bp); NA = "expected_length" of metadata/parameters.tsv (313 for COI)
TOP_PHYLA   <- 12         # phyla shown individually; the rest = "Other animals", "Other algae"…
# Order of the phyla in the per-phylum figures (2, 3, 5, 8, 10, 12): animals
# first, then algae, fungi and other eukaryotes, so that algae are not mixed with
# animals. A phylum is an animal when its lineage contains Metazoa or when it is
# in ANIMAL_PHYLA; algae, fungi as listed; anything else = "Other eukaryotes".
ALGAE_PHYLA  <- c("Rhodophyta", "Ochrophyta", "Chlorophyta", "Haptophyta", "Cryptophyta",
                  "Glaucophyta", "Streptophyta", "Euglenozoa")   # Streptophyta: seagrasses (rbcL)
FUNGI_PHYLA  <- c("Ascomycota", "Basidiomycota", "Chytridiomycota", "Mucoromycota", "Zoopagomycota",
                  "Blastocladiomycota", "Cryptomycota", "Rozellomycota", "Microsporidia", "Glomeromycota",
                  "Neocallimastigomycota", "Olpidiomycota", "Aphelidiomycota", "Monoblepharomycota")
ANIMAL_PHYLA <- c("Porifera", "Cnidaria", "Ctenophora", "Placozoa", "Xenacoelomorpha", "Platyhelminthes",
                  "Gastrotricha", "Rotifera", "Acanthocephala", "Gnathostomulida", "Micrognathozoa",
                  "Cycliophora", "Entoprocta", "Nemertea", "Mollusca", "Annelida", "Brachiopoda",
                  "Phoronida", "Bryozoa", "Chaetognatha", "Nematoda", "Nematomorpha", "Priapulida",
                  "Kinorhyncha", "Loricifera", "Tardigrada", "Onychophora", "Arthropoda", "Echinodermata",
                  "Hemichordata", "Chordata", "Orthonectida", "Dicyemida")
SET_ORDER    <- c("Animals", "Algae", "Fungi", "Other eukaryotes")
# Taxa left out of every figure except figure 1, per marker (names matched at any
# rank of the lineage). 18S: fungi are not studied in these samples, so they are
# not drawn; they stay in the database (so that fungal reads are not assigned to
# another group). Streptophyta (seagrasses…) are kept.
EXCLUDE_FROM_FIGURES <- list(`18S` = c("Fungi", FUNGI_PHYLA))
# Focal taxa of figure 13 (sessile / pioneer colonisers), in display order; a
# sequence counts for a taxon when the name appears at any rank of its lineage.
# FOCAL_ALIASES lists synonyms: other names under which a focal group can appear.
# The current databases already file diatoms and brown algae under Ochrophyta
# (Adl et al. 2019), so the aliases change nothing for them; they keep figure 13
# right for databases built before that harmonisation (diatoms under
# Bacillariophyta, or only a class such as Bacillariophyceae or Phaeophyceae).
FOCAL_TAXA  <- c("Porifera", "Cnidaria", "Bryozoa", "Annelida", "Nematoda",
                 "Platyhelminthes", "Rhodophyta", "Ochrophyta")
FOCAL_ALIASES <- list(Ochrophyta = c("Ochrophyta", "Heterokontophyta", "Bacillariophyta", "Bacillariophyceae",
                                     "Coscinodiscophyceae", "Mediophyceae", "Phaeophyceae"))
FOCAL_SUBTITLE <- "Sessile and pioneer taxa; colour = source of the sequence"
FOCAL_NOTE     <- "Ochrophyta includes diatoms and brown algae."
# rbcL (marine algae and seagrasses): the focal taxa are the main algal groups.
# Labels on the left, names matched in the lineage (any rank) on the right; a
# sequence goes to the FIRST group it matches (diatoms and brown algae before
# the other Ochrophyta).
FOCAL_RBCL <- list(
  "Diatoms"           = c("Bacillariophyceae", "Coscinodiscophyceae", "Mediophyceae"),
  "Brown algae"       = "Phaeophyceae",
  "Other Ochrophyta"  = c("Ochrophyta", "Heterokontophyta", "Bacillariophyta"),
  "Red algae"         = "Rhodophyta",
  "Green algae"       = "Chlorophyta",
  "Haptophytes"       = c("Haptophyta", "Prymnesiophyceae", "Pavlovophyceae"),
  "Cryptophytes"      = c("Cryptophyta", "Cryptophyceae"),
  "Seagrasses"        = c("Zosteraceae", "Posidoniaceae", "Cymodoceaceae"))
FOCAL_RBCL_SUBTITLE <- "Main groups of marine algae and seagrasses; colour = source of the sequence"
FOCAL_RBCL_NOTE     <- "Other Ochrophyta = ochrophytes other than diatoms and brown algae."
LENGTH_VIEW <- NA         # x range of the length histogram, e.g. c(150, 500); NA = kept window ± 150 bp
MAX_PRIMER_SEQS <- 50000  # sequences sampled for the primer-site figures (speed)
DPI         <- 300
WIDTH_IN    <- 8          # figure width (inches)

# =============================================================================
#                      NOTHING TO EDIT BELOW THIS LINE
# =============================================================================

# ------------------------------ Packages --------------------------------------
needed  <- c("data.table", "ggplot2", "scales", "viridisLite")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  message("Installing missing packages: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(scales)
})

# ------------------------------ Paths -----------------------------------------
# Folder of this script (Rscript, RStudio Source, or else the working directory)
script_file <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(f)) return(normalizePath(f))
  of <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(of)) return(normalizePath(of))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- rstudioapi::getSourceEditorContext()$path
    if (nzchar(p)) return(normalizePath(p))
  }
  NA_character_
})
script_dir <- if (!is.na(script_file)) dirname(script_file) else getwd()
CRABS_HOME <- if (!is.na(CRABS_HOME)) path.expand(CRABS_HOME) else
              if (nzchar(Sys.getenv("CRABS_HOME"))) path.expand(Sys.getenv("CRABS_HOME")) else
              if (basename(script_dir) == "scripts") dirname(script_dir) else script_dir

# Final databases of a folder: results/<base>.tsv (several when a folder holds
# more than one database, e.g. crabs_18s with the 18S forward database)
db_bases <- function(x) {
  f <- list.files(file.path(x, "results"), pattern = "\\.tsv$")
  sub("\\.tsv$", "", f[!grepl("insilico_pcr_raw", f)])
}
# Figures folder: figures/ for the main database of a folder, figures_<base>/
# for a database whose files carry its name (metadata/<base>_parameters.tsv)
fig_dir_of <- function(x, b) if (file.exists(file.path(x, "metadata", paste0(b, "_parameters.tsv")))) paste0("figures_", b) else "figures"
run_each <- function(jobs) {     # jobs: list of c(folder, base); one fresh R process each
  if (is.na(script_file)) stop("Cannot find the path of this script: run it with Rscript, or Source it in RStudio")
  rscript <- file.path(R.home("bin"), "Rscript")
  status <- vapply(jobs, function(j) {
    message("\n===== ", basename(j[1]), " / ", j[2], " =====")
    system2(rscript, c(shQuote(script_file), shQuote(j[1]), shQuote(j[2])))
  }, integer(1))
  message("\n===== Summary =====")
  for (i in seq_along(jobs))
    message(sprintf("  %-14s %-40s %s", basename(jobs[[i]][1]), jobs[[i]][2],
                    if (status[i] == 0) paste0("OK → ", file.path(jobs[[i]][1], fig_dir_of(jobs[[i]][1], jobs[[i]][2]))) else "FAILED (see the messages above)"))
  if (interactive()) invokeRestart("abort") else quit(save = "no", status = as.integer(any(status != 0)))
}

# ---- All databases: each one is drawn by a separate run of this script, in its
#      own folder (a fresh R process each time: no mix-up between markers)
if (identical(DB_ARG, "all")) {
  roots <- unique(c(script_dir, CRABS_HOME))
  roots <- roots[dir.exists(roots)]
  # database folders: crabs_*/ with results/, or their sub-folders with results/
  # (e.g. crabs_18S/18S_aligned and crabs_18S/18S_forward)
  dbs <- unique(unlist(lapply(roots, function(r) {
    d <- list.dirs(r, recursive = FALSE)
    d <- d[grepl("^crabs_", basename(d))]
    d <- unlist(lapply(d, function(x) c(x, list.dirs(x, recursive = FALSE))))
    d[vapply(d, function(x) length(db_bases(x)) > 0, logical(1))]
  })))
  dbs <- dbs[!duplicated(normalizePath(dbs))]
  if (!length(dbs)) stop("No database found (crabs_*/results/*.tsv) in: ", paste(roots, collapse = ", "))
  jobs <- unlist(lapply(dbs, function(d) lapply(db_bases(d), function(b) c(d, b))), recursive = FALSE)
  message("Databases found: ", paste(vapply(jobs, function(j) paste0(basename(j[1]), "/", j[2]), ""), collapse = ", "))
  run_each(jobs)
}

db_path <- path.expand(DB_ARG)
# Where the database folder is looked for. In R / RStudio the working directory
# may still be the folder of a previous run (the script moves into it), so
# CRABS_HOME (see PARAMETERS) and the parent folder are tried too.
cands <- if (startsWith(db_path, "/")) db_path else
  unique(c(file.path(script_dir, DB_ARG), file.path(CRABS_HOME, DB_ARG),
           file.path(getwd(), DB_ARG), file.path(dirname(getwd()), DB_ARG)))
hit <- cands[dir.exists(cands) & dir.exists(file.path(cands, "results"))]
if (!length(hit)) stop("Database folder not found: ", DB_ARG, "\n  looked in:\n    ",
                       paste(cands, collapse = "\n    "))
PROJECT_DIR <- hit[1]
message("Database folder: ", normalizePath(PROJECT_DIR))
setwd(PROJECT_DIR)
if (is.na(BASE_NAME)) {
  cand <- db_bases(".")
  if (!length(cand)) stop("No final database (results/<base>.tsv) in ", getwd())
  if (length(cand) > 1) run_each(lapply(cand, function(b) c(getwd(), b)))   # several databases in this folder
  BASE_NAME <- cand
}
# A database whose files carry its name (e.g. 18S forward: metadata/<base>_…)
# has its own metadata, steps and figures, all named after it
PFX     <- if (file.exists(file.path("metadata", paste0(BASE_NAME, "_parameters.tsv")))) paste0(BASE_NAME, "_") else ""
FIG_DIR <- fig_dir_of(".", BASE_NAME)
message("Database: ", BASE_NAME, if (nzchar(PFX)) paste0(" (files ", PFX, "*; figures in ", FIG_DIR, "/)") else "")
dir.create(file.path(FIG_DIR, "data"), recursive = TRUE, showWarnings = FALSE)
F_FINAL    <- file.path("results", paste0(BASE_NAME, ".tsv"))
F_COUNTS   <- paste0("metadata/", PFX, "counts.tsv")
F_SOURCES  <- paste0("metadata/", PFX, "sources.tsv")
F_PARAMS   <- paste0("metadata/", PFX, "parameters.tsv")
F_INSILICO <- paste0("steps/", PFX, "04_insilico_pcr.txt")
F_RAW      <- "steps/03_filtered_initial.txt"   # full sequences, before in silico PCR
# imported records: in this folder, else (database built from another one, e.g.
# 18S forward from 18S_aligned) in the folder of the source database
IMP_DIR    <- "imported"
if (!dir.exists(IMP_DIR) && nzchar(PFX)) {
  sib <- file.path(list.dirs("..", recursive = FALSE), "imported")
  sib <- sib[dir.exists(sib)]
  if (length(sib)) IMP_DIR <- sib[1]
}
F_BOLD     <- file.path(IMP_DIR, "bold_all.txt")
F_NCBI     <- file.path(IMP_DIR, "ncbi_all.txt")
F_SILVA    <- file.path(IMP_DIR, "silva_all.txt")
if (!file.exists(F_RAW) && IMP_DIR != "imported") F_RAW <- file.path(dirname(IMP_DIR), "steps", "03_filtered_initial.txt")
# NCBI taxonomy (taxdump): the shared CRABS/taxo, else in this database folder or in crabs_coi
TAXO_DIR   <- c("../taxo", "../../taxo", "taxo", "../crabs_coi/taxo", "../../crabs_coi/taxo")   # shared CRABS/taxo first
TAXO_DIR   <- c(TAXO_DIR[file.exists(file.path(TAXO_DIR, "nodes.dmp"))], "../taxo")[1]
# Optional manual corrections of the group table: TSV with columns group, comparison_group
F_GROUP_FIX <- "metadata/group_taxonomy_manual.tsv"
# Optional list of species of interest (one name per line) for the completeness
# table/figure (CRABS --completeness-table): in this database folder or next to the script
F_SPECIES  <- c("species_of_interest.txt", file.path(script_dir, "species_of_interest.txt"))
F_SPECIES  <- F_SPECIES[file.exists(F_SPECIES)][1]
if (!file.exists(F_FINAL)) stop("Final database not found: ", F_FINAL)
# Marker-specific values read from the pipeline metadata (COI, 18S, 16S, rbcL…)
PARAMS <- if (file.exists(F_PARAMS)) fread(F_PARAMS, sep = "\t", quote = "") else data.table(parameter = character(), value = character())
par_num <- function(k) suppressWarnings(as.numeric(PARAMS[parameter == k, value][1]))
if (is.na(EXPECTED_BP)) EXPECTED_BP <- par_num("expected_length")
if (is.na(EXPECTED_BP) && grepl("COI", basename(F_FINAL), ignore.case = TRUE)) EXPECTED_BP <- 313
# Marker: from metadata/parameters.tsv ("marker"), else from the folder / base name
MARKER_TXT <- paste(PARAMS[parameter == "marker", value][1], basename(PROJECT_DIR), BASE_NAME)
MARKER <- if (grepl("rbcL", MARKER_TXT, ignore.case = TRUE)) "rbcL" else
  if (grepl("18S", MARKER_TXT, ignore.case = TRUE)) "18S" else
    if (grepl("16S", MARKER_TXT, ignore.case = TRUE)) "16S" else
      if (grepl("COI", MARKER_TXT, ignore.case = TRUE)) "COI" else "other"
message("Marker: ", MARKER)
if (MARKER == "rbcL") {
  FOCAL_TAXA <- names(FOCAL_RBCL); FOCAL_ALIASES <- FOCAL_RBCL
  FOCAL_SUBTITLE <- FOCAL_RBCL_SUBTITLE; FOCAL_NOTE <- FOCAL_RBCL_NOTE
}
if (any(is.na(LENGTH_VIEW))) {
  b <- c(par_num("min_length"), par_num("max_length"))
  LENGTH_VIEW <- if (all(!is.na(b))) c(max(0, b[1] - 150), b[2] + 150) else c(0, 1000)
}

# ------------------------------ Style -----------------------------------------
VIR_BEGIN <- 0.10; VIR_END <- 0.85
VIR_2 <- viridisLite::viridis(3, begin = 0.15, end = 0.80)   # BOLD, NCBI GenBank, SILVA
COL <- list(surface = "#fcfcfb", ink = "#0b0b0b", ink2 = "#52514e", muted = "#898781",
            grid = "#e6e5e0", s1 = viridisLite::viridis(1, begin = 0.45),
            bold = VIR_2[1], ncbi = VIR_2[2], silva = VIR_2[3], band = "#e8f4f1")
fill_ordered <- function(...) scale_fill_viridis_d(begin = VIR_BEGIN, end = VIR_END, guide = "none", ...)
theme_db <- function() {
  theme_minimal(base_size = 11) +
    theme(plot.title.position = "plot",
          plot.title    = element_text(face = "bold", colour = COL$ink, size = 13),
          plot.subtitle = element_text(colour = COL$ink2, margin = margin(b = 10)),
          plot.caption  = element_text(colour = COL$muted, hjust = 0, size = 8.5),
          axis.title    = element_text(colour = COL$ink2, size = 10),
          axis.text     = element_text(colour = COL$ink2),
          strip.text    = element_text(face = "bold", colour = COL$ink, hjust = 0),
          legend.position = "top", legend.justification = "left",
          legend.title  = element_blank(), legend.text = element_text(colour = COL$ink2),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = COL$grid, linewidth = 0.3),
          plot.background  = element_rect(fill = COL$surface, colour = NA),
          plot.margin = margin(14, 18, 10, 14))
}
CAPTION <- NULL   # no caption under the figures (base and script version are in figures/README.md)
num <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
pct <- function(x) paste0(round(100 * x, 1), "%")

FIGS <- list()   # registry for the README
# A figure that fails is skipped with a message (its half-written files removed)
save_fig <- function(p, id, name, height, data, caption) {
  stem <- sprintf("fig%02d_%s", id, name)
  ok <- tryCatch({ save_fig_files(p, stem, height, data, caption); TRUE }, error = function(e) {
    while (!is.null(dev.list())) dev.off()
    unlink(file.path(FIG_DIR, paste0(stem, c(".png", ".pdf"))))
    unlink(file.path(FIG_DIR, "data", paste0(stem, ".csv")))
    message("  ", stem, ": SKIPPED (", conditionMessage(e), ")")
    FALSE
  })
  invisible(ok)
}
save_fig_files <- function(p, stem, height, data, caption) {
  png_dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else
    if (Sys.info()[["sysname"]] == "Darwin") function(...) grDevices::png(..., type = "quartz") else "png"
  ggsave(file.path(FIG_DIR, paste0(stem, ".png")), p, width = WIDTH_IN, height = height,
         dpi = DPI, bg = COL$surface, device = png_dev)
  # PDF: the base pdf() device works everywhere (better for paper publication)
  f_pdf <- file.path(FIG_DIR, paste0(stem, ".pdf"))
  ok <- FALSE
  if (Sys.info()[["sysname"]] != "Darwin" && isTRUE(capabilities("cairo"))) {
    ok <- tryCatch({ ggsave(f_pdf, p, width = WIDTH_IN, height = height, device = cairo_pdf,
                            bg = COL$surface); TRUE }, error = function(e) FALSE)
  }
  if (!ok) ggsave(f_pdf, p, width = WIDTH_IN, height = height, device = grDevices::pdf,
                  useDingbats = FALSE, bg = COL$surface)
  fwrite(data, file.path(FIG_DIR, "data", paste0(stem, ".csv")))
  FIGS[[length(FIGS) + 1]] <<- list(stem = stem, caption = caption)
  message("  ", stem)
}

# ------------------------------ Data ------------------------------------------
NA_STR <- c("", "NA", "nan", "None", "NaN")
message("Reading ", F_FINAL)
db <- fread(F_FINAL, sep = "\t", quote = "", na.strings = NA_STR, colClasses = "character")
RANKS <- setdiff(names(db), c("seq_id", "scientific_name", "taxid", "sequence"))
db[, len := nchar(sequence)]
# Some NCBI lineages have no taxon of rank "phylum" (brown algae, dinoflagellates…).
# crabs_coi.sh fills it (Adl et al. 2019); any sequence still without a phylum
# (older database, class missing from the table) is labelled with the next rank
# down (class, then order, family), so no "Unassigned" bar appears. The detail
# is saved in figures/data/phylum_not_ranked_in_ncbi.csv.
db[, phylum_raw := phylum]
# Same Adl et al. (2019) harmonisation as crabs_coi.sh (tables read from
# metadata/parameters.tsv when present), so a database built before it matches
kv <- function(s) { s <- strsplit(strsplit(s, ";", fixed = TRUE)[[1]], "=", fixed = TRUE)
setNames(vapply(s, `[`, "", 2), vapply(s, `[`, "", 1)) }
par_chr <- function(k, default) { v <- PARAMS[parameter == k, value]; if (length(v) && nzchar(v[1])) v[1] else default }
PHYLUM_RENAME <- kv(par_chr("phylum_renamed", "Bacillariophyta=Ochrophyta;Heterokontophyta=Ochrophyta;Myzozoa=Dinoflagellata;Dinophyta=Dinoflagellata"))
CLASS_TO_PHYLUM <- kv(par_chr("phylum_from_class", paste0(
  "Phaeophyceae=Ochrophyta;Bacillariophyceae=Ochrophyta;Coscinodiscophyceae=Ochrophyta;Mediophyceae=Ochrophyta;",
  "Chrysophyceae=Ochrophyta;Synurophyceae=Ochrophyta;Raphidophyceae=Ochrophyta;Eustigmatophyceae=Ochrophyta;",
  "Phaeosacciophyceae=Ochrophyta;Pelagophyceae=Ochrophyta;Dictyochophyceae=Ochrophyta;Xanthophyceae=Ochrophyta;",
  "Bolidophyceae=Ochrophyta;Pinguiophyceae=Ochrophyta;Phaeothamniophyceae=Ochrophyta;Dinophyceae=Dinoflagellata")))
db[phylum %chin% names(PHYLUM_RENAME), phylum := PHYLUM_RENAME[phylum]]
if ("class" %in% names(db)) db[is.na(phylum) & class %chin% names(CLASS_TO_PHYLUM), phylum := CLASS_TO_PHYLUM[class]]
for (r in intersect(c("class", "order", "family"), names(db))) db[is.na(phylum), phylum := get(r)]
db[is.na(phylum), phylum := "Phylum unknown"]
fb <- db[is.na(phylum_raw), .N, by = .(label_used = phylum, class = if ("class" %in% names(db)) class else NA_character_)][order(-N)]
fwrite(fb, file.path(FIG_DIR, "data", "phylum_not_ranked_in_ncbi.csv"))
if (nrow(fb)) message("  ", num(sum(fb$N)), " sequences without a phylum rank in NCBI, labelled by class/order (figures/data/phylum_not_ranked_in_ncbi.csv)")
lin_all <- intersect(RANKS, names(db))
in_lineage <- function(names_) Reduce(`|`, lapply(lin_all, function(r) db[[r]] %chin% names_))
# Taxa left out of the figures for this marker (see EXCLUDE_FROM_FIGURES)
EXCL <- EXCLUDE_FROM_FIGURES[[MARKER]]
N_EXCL <- 0L
if (length(EXCL)) {
  drop <- in_lineage(EXCL) | db$phylum %chin% EXCL
  N_EXCL <- sum(drop)
  if (N_EXCL) db <- db[!drop]
}
# Group of each phylum (panels of the per-phylum figures)
set_of <- function(ph, metazoa = FALSE) fifelse(metazoa | ph %chin% ANIMAL_PHYLA, "Animals",
                                         fifelse(ph %chin% ALGAE_PHYLA, "Algae",
                                                 fifelse(ph %chin% FUNGI_PHYLA, "Fungi", "Other eukaryotes")))
db[, set := set_of(phylum, in_lineage("Metazoa"))]
db[, set := names(sort(table(set), decreasing = TRUE))[1], by = phylum]   # one group per phylum
SET_LOOKUP <- unique(db[, .(phylum, set)])
OTHER_LABEL <- c(Animals = "Other animals", Algae = "Other algae", Fungi = "Other fungi",
                 `Other eukaryotes` = "Other eukaryotes")
top <- db[, .N, by = phylum][order(-N)][seq_len(min(.N, TOP_PHYLA)), phylum]
db[, phylum_group := fifelse(phylum %in% top, phylum, OTHER_LABEL[set])]
is_other <- function(x) as.character(x) %chin% OTHER_LABEL
SETS <- intersect(SET_ORDER, unique(db$set))
# order: by group, then by size within a group (largest at the top), "Other …" at the bottom
go <- db[, .N, by = .(set, phylum_group)][, oth := is_other(phylum_group)]
go[, set := factor(set, levels = SETS)]
group_order <- go[order(-as.integer(set), -oth, N), as.character(phylum_group)]
db[, phylum_group := factor(phylum_group, levels = group_order)]
db[, set := factor(set, levels = SETS)]
SET_OF_GROUP <- setNames(as.character(go$set), as.character(go$phylum_group))
# No panels: the groups are only ordered (animals first, then algae, fungi,
# other eukaryotes), see group_order above
by_set <- function(p, cols = ".") if (cols == ".") p else p + facet_wrap(as.formula(paste("~", cols)), scales = "free_x")
species_col <- if ("species" %in% names(db)) "species" else "scientific_name"

# Source of each retained record (exact match on the identifiers of each import)
read_ids <- function(f) if (file.exists(f)) fread(f, sep = "\t", quote = "", header = FALSE,
                                                  select = 1, colClasses = "character")[[1]] else character(0)
ids_bold <- read_ids(F_BOLD); ids_ncbi <- read_ids(F_NCBI); ids_silva <- read_ids(F_SILVA)
HAS_SOURCE <- length(ids_bold) + length(ids_ncbi) + length(ids_silva) > 0
if (HAS_SOURCE) {
  db[, source := fifelse(seq_id %chin% ids_bold, "BOLD",
                         fifelse(seq_id %chin% ids_ncbi, "NCBI GenBank",
                                 fifelse(seq_id %chin% ids_silva, "SILVA", "Other")))]
}
rm(ids_bold, ids_ncbi, ids_silva)

message("Drawing figures:")

# =============================================================================
# 1. Pipeline yield: sequences kept at each step
# =============================================================================
if (file.exists(F_COUNTS)) {
  cn <- fread(F_COUNTS, sep = "\t", quote = "")
  if (nzchar(PFX)) cn[startsWith(step, PFX), step := substring(step, nchar(PFX) + 1)]
  cn[startsWith(file, "results/"), step := "final"]
  # Step names, with the parameters of each filter (read from metadata/parameters.tsv)
  pv <- function(k) { v <- PARAMS[parameter == k, value]; if (length(v) && nzchar(v[1]) && v[1] != "NA") sub(" .*$", "", v[1]) else NA_character_ }
  wrap <- function(name, ...) { x <- c(...); x <- x[!is.na(x) & nzchar(x)]; if (length(x)) sprintf("%s (%s)", name, paste(x, collapse = ", ")) else name }
  maybe <- function(fmt, v) if (is.na(v)) NA_character_ else sprintf(fmt, v)
  pga_id <- suppressWarnings(as.numeric(pv("pga_identity"))); if (!is.na(pga_id) && pga_id <= 1) pga_id <- 100 * pga_id
  max_n <- suppressWarnings(as.integer(pv("max_N")))
  excl <- PARAMS[parameter == "excluded_taxa", value][1]
  excl <- if (length(excl) && !is.na(excl) && nzchar(excl)) strsplit(excl, ";", fixed = TRUE)[[1]] else character(0)
  derep <- maybe("%s", pv("dereplication"))
  STEP_LABEL <- c(
    import_bold             = "BOLD import",
    import_ncbi             = "NCBI GenBank import",
    import_silva            = "SILVA import",
    "01_merged"             = "Merge (unique accessions)",
    "01c_taxa_excluded_early" = "Removal of terrestrial taxa",
    "02_dereplicated_initial" = wrap("Dereplication", "records", derep),
    # CRABS --filter --environmental --rank-na: records marked "environmental", and
    # records with too many missing ranks (all ranks counted, so only a rough test
    # of the identification level; the exact one is "Minimum rank" at the end)
    "03_filtered_initial"   = wrap("Pre-filter", "environmental samples",
                                   { k <- suppressWarnings(as.integer(pv("crabs_rank_na"))); if (is.na(k)) NA else sprintf("max. %d missing ranks", k - 1L) }),
    "04_insilico_pcr"       = wrap("In silico PCR", maybe("max. %s mismatches", pv("pcr_mismatch"))),
    "05_global_alignment"   = wrap("Global alignment", maybe("min. %s%% identity", if (is.na(pga_id)) NA else format(pga_id)),
                                   maybe("%s%% coverage", pv("pga_coverage"))),
    "05b_long_trimmed"      = "Re-cut of over-long amplicons",
    "06_dereplicated_final" = wrap("Dereplication", "amplicons", derep),
    "07_filtered_final"     = wrap("Length filter", if (all(!is.na(c(pv("min_length"), pv("max_length"))))) sprintf("%s-%s bp", pv("min_length"), pv("max_length")) else NA,
                                   if (!is.na(max_n)) sprintf("max. %d N", max_n - 1L) else NA),
    "07b_length_rescued"    = wrap("Length rescue", maybe("min. %s bp", pv("min_length_rescue"))),
    "07c_translation_checked" = "Translation check (no stop codon)",
    "07d_taxonomy_fixed"    = "Taxonomy correction (homonyms)",
    "08_taxa_excluded"      = if (length(excl)) sprintf("Taxa exclusion (%d taxa, e.g. %s)", length(excl), excl[1]) else "Taxa exclusion",
    "09_min_rank"           = wrap("Minimum rank", pv("min_identification_rank")),
    final                   = "Final database")
  if (nzchar(PFX)) STEP_LABEL["03_filtered_initial"] <- "Input (pre-filtered records of the parent database)"
  cn <- cn[step %chin% names(STEP_LABEL)]               # internal steps (e.g. 01b lineage rescue) not shown
  cn <- cn[startsWith(step, "import_") | step == "final" | file.exists(file)]   # steps not run: not shown
  # steps that changed nothing (option off) are not shown
  prev <- c(NA, head(cn$n_sequences, -1)); prev[startsWith(cn$step, "import_")] <- NA
  cn <- cn[step == "final" | is.na(prev) | n_sequences != prev]
  cn[, label := unname(STEP_LABEL[step])]
  cn[, label := factor(label, levels = rev(unique(label)))]
  merged <- cn[step == "01_merged", n_sequences]
  n_fin  <- cn[step == "final", n_sequences]
  p <- ggplot(cn, aes(n_sequences, label)) +
    geom_col(aes(fill = label), width = 0.68) +
    geom_text(aes(label = num(n_sequences)), hjust = -0.12, size = 3.3, colour = COL$ink2) +
    fill_ordered(direction = -1) +
    scale_x_continuous(labels = label_comma(), expand = expansion(mult = c(0, 0.2))) +
    coord_cartesian(clip = "off") +
    labs(title = sprintf("%s sequences in the final database", num(n_fin)),
         subtitle = if (length(merged)) sprintf("%s of the merged records reach the final database", pct(n_fin / merged)) else NULL,
         x = "Sequences", y = NULL, caption = CAPTION) +
    theme_db() + theme(panel.grid.major.y = element_blank())
  save_fig(p, 1, "pipeline_yield", 1.2 + 0.30 * nrow(cn), cn[, .(step, label, n_sequences)],
           "Number of sequences after each step of the pipeline, from import to the final database (parameters in brackets; steps that changed nothing are not shown).")
}

# =============================================================================
# 2. Taxonomic composition: species and sequences per phylum
# =============================================================================
comp <- db[, .(Sequences = .N, Species = uniqueN(get(species_col))), by = .(set, phylum_group)]
comp_long <- melt(comp, id.vars = c("set", "phylum_group"), variable.name = "measure", value.name = "n")
comp_long[, measure := factor(measure, levels = c("Species", "Sequences"))]
lead <- comp[!is_other(phylum_group)][order(-Sequences)][1]
p <- ggplot(comp_long, aes(n, phylum_group)) +
  geom_col(aes(fill = phylum_group), width = 0.68) + fill_ordered(direction = -1) +
  geom_text(aes(label = num(n)), hjust = -0.12, size = 3.1, colour = COL$ink2) +
  scale_x_continuous(labels = label_comma(), expand = expansion(mult = c(0, 0.22))) +
  labs(title = sprintf("%s holds the most reference sequences", lead$phylum_group),
       subtitle = sprintf("%s sequences and %s species across %d phyla",
                          num(nrow(db)), num(uniqueN(db[[species_col]])), uniqueN(db$phylum)),
       x = NULL, y = NULL, caption = CAPTION) +
  theme_db() + theme(panel.grid.major.y = element_blank(), panel.spacing.x = unit(18, "pt"))
p <- by_set(p, "measure")
save_fig(p, 2, "composition_phylum", 0.9 + 0.32 * length(group_order), comp,
         sprintf("Number of species and of sequences per phylum, animals first, then algae (top %d phyla; the others grouped as 'Other animals', 'Other algae'...).", TOP_PHYLA))

# =============================================================================
# 3. Data source of the retained records, per phylum
# =============================================================================
if (HAS_SOURCE) {
  src <- db[, .N, by = .(set, phylum_group, source)]
  SRC_LV <- intersect(c("BOLD", "NCBI GenBank", "SILVA", "Other"), unique(src$source))
  src[, source := factor(source, levels = rev(SRC_LV))]
  top_src <- db[, .N, by = source][order(-N)][1]
  p <- ggplot(src, aes(N, phylum_group, fill = source)) +
    geom_col(width = 0.68, colour = COL$surface, linewidth = 0.4) +
    scale_fill_manual(values = c(BOLD = COL$bold, `NCBI GenBank` = COL$ncbi, SILVA = COL$silva, Other = COL$muted),
                      breaks = SRC_LV) +
    scale_x_continuous(labels = label_comma(), expand = expansion(mult = c(0, 0.04))) +
    labs(title = sprintf("%s supplies %s of the retained sequences", top_src$source, pct(top_src$N / nrow(db))),
         subtitle = "Origin of each record kept in the final database",
         x = "Sequences", y = NULL, caption = CAPTION) +
    theme_db() + theme(panel.grid.major.y = element_blank())
  p <- by_set(p)
  save_fig(p, 3, "source_by_phylum", 1.2 + 0.32 * length(group_order), dcast(src, set + phylum_group ~ source, value.var = "N", fill = 0),
           "Number of retained sequences per phylum, by the database they come from (exact match of identifiers).")
}

if (is.na(EXPECTED_BP)) EXPECTED_BP <- as.integer(median(db$len))   # no expected length given
# =============================================================================
# 4. Amplicon lengths before the final filter, with the kept window
# =============================================================================
if (file.exists(F_INSILICO)) {
  bounds <- c(NA, NA)
  if (file.exists(F_PARAMS)) {
    pr <- fread(F_PARAMS, sep = "\t", quote = "")
    bounds <- as.numeric(pr[parameter %in% c("min_length", "max_length"), value])
  }
  ncol_ins <- length(strsplit(readLines(F_INSILICO, n = 1), "\t")[[1]])
  lens <- nchar(fread(F_INSILICO, sep = "\t", quote = "", header = FALSE, select = ncol_ins,
                      colClasses = "character")[[1]])
  hist_dt <- data.table(len = lens)[, .N, by = len][order(len)]
  inside  <- if (all(!is.na(bounds))) mean(lens >= bounds[1] & lens <= bounds[2]) else NA
  p <- ggplot(hist_dt[len >= LENGTH_VIEW[1] & len <= LENGTH_VIEW[2]], aes(len, N))
  if (!is.na(inside)) p <- p + annotate("rect", xmin = bounds[1], xmax = bounds[2], ymin = -Inf, ymax = Inf,
                                        fill = COL$band) +
    annotate("text", x = bounds[2], y = Inf, label = sprintf("kept: %d-%d bp", bounds[1], bounds[2]),
             hjust = -0.05, vjust = 1.6, size = 3.2, colour = COL$ink2)
  p <- p + geom_col(fill = COL$s1, width = 1) +
    geom_vline(xintercept = EXPECTED_BP, linetype = "22", colour = COL$ink2, linewidth = 0.4) +
    annotate("text", x = EXPECTED_BP, y = Inf, label = sprintf("expected %d bp", EXPECTED_BP),
             hjust = 1.05, vjust = 1.6, size = 3.2, colour = COL$ink2) +
    scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, 0.08))) +
    labs(title = if (!is.na(inside)) sprintf("%s of in silico amplicons fall inside the kept window", pct(inside)) else "Amplicon lengths after in silico PCR",
         subtitle = sprintf("%s amplicons, median %d bp (view limited to %d-%d bp)", num(length(lens)),
                            as.integer(median(lens)), LENGTH_VIEW[1], LENGTH_VIEW[2]),
         x = "Amplicon length without primers (bp)", y = "Amplicons", caption = CAPTION) +
    theme_db() + theme(panel.grid.major.x = element_blank())
  save_fig(p, 4, "amplicon_lengths", 4.4, hist_dt,
           "Length distribution of the amplicons extracted by in silico PCR, before the final length filter; shaded = window kept.")
  rm(lens)
}

# =============================================================================
# 5. Final amplicon lengths per phylum (small multiples)
# =============================================================================
len_ph <- db[!is_other(phylum_group), .N, by = .(phylum_group, len)]
len_ph[, phylum_group := factor(phylum_group, levels = rev(group_order))]   # animals first, largest phylum first
med_len <- as.integer(median(db$len)); near <- db[, mean(abs(len - med_len) <= 3)]
p <- ggplot(len_ph, aes(len, N)) +
  geom_col(fill = COL$s1, width = 1) +
  geom_vline(xintercept = EXPECTED_BP, linetype = "22", colour = COL$ink2, linewidth = 0.35) +
  facet_wrap(~ phylum_group, scales = "free_y", ncol = 3) +
  scale_y_continuous(labels = label_comma(), n.breaks = 3) +
  labs(title = sprintf("%s of final amplicons lie within 3 bp of the median (%d bp)", pct(near), med_len),
       subtitle = sprintf("Final database; dashed line = expected %d bp; y scales differ between panels", EXPECTED_BP),
       x = "Amplicon length (bp)", y = "Sequences", caption = CAPTION) +
  theme_db() + theme(panel.grid.major.x = element_blank())
n_panels <- uniqueN(len_ph$phylum_group)
save_fig(p, 5, "lengths_by_phylum", 1.4 + 1.25 * ceiling(n_panels / 3), len_ph,
         "Length distribution of the final amplicons for each of the most represented phyla.")

# =============================================================================
# 6. Reference depth: number of sequences per species
# =============================================================================
depth <- db[, .N, by = c(species_col)]
brks <- c(0, 1, 2, 5, 10, 50, Inf)
labs_d <- c("1", "2", "3-5", "6-10", "11-50", "> 50")
depth[, bin := factor(cut(N, brks, labels = labs_d), levels = labs_d)]
dd <- depth[, .(species = .N), by = bin][order(bin)]
single <- depth[, mean(N == 1)]
p <- ggplot(dd, aes(bin, species)) +
  geom_col(aes(fill = bin), width = 0.68) + fill_ordered(drop = FALSE) +
  geom_text(aes(label = num(species)), vjust = -0.5, size = 3.2, colour = COL$ink2) +
  scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, 0.12))) +
  labs(title = sprintf("%s of species have a single reference sequence", pct(single)),
       subtitle = sprintf("Distribution of the %s species by number of unique sequences", num(nrow(depth))),
       x = "Unique sequences per species", y = "Species", caption = CAPTION) +
  theme_db() + theme(panel.grid.major.x = element_blank())
save_fig(p, 6, "sequences_per_species", 4.2, dd,
         "Number of species by number of unique reference sequences (intraspecific variants kept by the dereplication).")

# =============================================================================
# 7. Lineage completeness: share of sequences with each rank resolved
# =============================================================================
compl <- data.table(rank = RANKS, resolved = vapply(RANKS, function(r) mean(!is.na(db[[r]])), numeric(1)))
compl[, rank := factor(rank, levels = RANKS)]
weakest <- compl[which.min(resolved)]
ttl7 <- if (weakest$resolved >= 1) "Every sequence has a complete lineage" else
  sprintf("Lineages are least complete at the %s rank (%s resolved)", weakest$rank, pct(weakest$resolved))
p <- ggplot(compl, aes(rank, resolved)) +
  geom_col(aes(fill = rank), width = 0.68) + fill_ordered() +
  geom_text(aes(label = pct(resolved)), vjust = -0.5, size = 3.2, colour = COL$ink2) +
  scale_y_continuous(labels = label_percent(), limits = c(0, 1.08), breaks = seq(0, 1, 0.25)) +
  labs(title = ttl7,
       subtitle = "Share of sequences whose NCBI lineage names each rank",
       x = NULL, y = "Sequences with the rank resolved", caption = CAPTION) +
  theme_db() + theme(panel.grid.major.x = element_blank())
save_fig(p, 7, "lineage_completeness", 4.0, compl,
         "Share of sequences in the final database with a name at each taxonomic rank of the NCBI lineage.")

# =============================================================================
# 8. Records downloaded per taxonomic group, BOLD vs NCBI on a common taxonomy
# =============================================================================
# BOLD and NCBI were queried with different group names (e.g. BOLD: Malacostraca,
# Copepoda…; NCBI: Crustacea). Each name is resolved in the NCBI taxonomy and
# assigned a COMPARISON GROUP = its highest ancestor (itself included) that is one
# of the downloaded group names; so Malacostraca (BOLD) and Crustacea (NCBI) both
# become "Crustacea", Ascidiacea (BOLD) and Tunicata (NCBI) both "Tunicata".
# The table is saved (figures/data/group_taxonomy.tsv) and can be corrected by
# hand in metadata/group_taxonomy_manual.tsv (columns: group, comparison_group).
# NCBI taxonomy, read once and shared by figures 8 and 11 (NULL if absent)
TAX <- NULL
load_taxonomy <- function(taxo_dir) {
  if (!is.null(TAX)) return(TAX)
  f_names <- file.path(taxo_dir, "names.dmp"); f_nodes <- file.path(taxo_dir, "nodes.dmp")
  if (!file.exists(f_names) || !file.exists(f_nodes)) return(NULL)
  message("  (reading the NCBI taxonomy in ", taxo_dir, ", about a minute)")
  tidy <- function(x) trimws(x, whitespace = "[\\t ]")   # .dmp fields are padded with tabs
  nm <- fread(f_names, sep = "|", quote = "", header = FALSE, select = c(1, 2, 4),
              col.names = c("taxid", "name", "class"))
  nm[, `:=`(name = tidy(name), class = tidy(class))]
  nm <- nm[class %in% c("scientific name", "synonym")]
  nd <- fread(f_nodes, sep = "|", quote = "", header = FALSE, select = 1:3,
              col.names = c("taxid", "parent", "rank"))
  nd[, rank := tidy(rank)]
  n <- max(nd$taxid)
  parent <- integer(n); parent[nd$taxid] <- nd$parent
  rnk <- character(n);  rnk[nd$taxid] <- nd$rank
  sciname <- character(n); sc <- nm[class == "scientific name"]; sciname[sc$taxid] <- sc$name
  TAX <<- list(names = nm, parent = parent, rank = rnk, sciname = sciname)
  TAX
}
tax_lineage <- function(id) { v <- integer(0); while (id > 1 && length(v) < 100) { v <- c(v, id); id <- TAX$parent[id] }; v }

resolve_groups <- function(groups, taxo_dir) {
  out <- data.table(group = groups, taxid = NA_integer_, rank = NA_character_,
                    phylum = NA_character_, lineage_ids = list(integer(0)))
  if (is.null(load_taxonomy(taxo_dir))) return(out[, comparison_group := group][, lineage_ids := NULL][])
  nm <- TAX$names[name %chin% groups]; rnk <- TAX$rank; sciname <- TAX$sciname; lineage <- tax_lineage
  for (i in seq_len(nrow(out))) {
    cand <- nm[name == out$group[i]]
    if (!nrow(cand)) next
    cand[, sci := class == "scientific name"]
    # homonyms: prefer a scientific name, then a taxon inside the Eukaryota (2759)
    cand[, euk := vapply(taxid, function(t) 2759L %in% lineage(t), logical(1))]
    best <- cand[order(-sci, -euk)][1]
    lin <- lineage(best$taxid)
    out[i, `:=`(taxid = best$taxid, rank = rnk[best$taxid],
                phylum = { p <- lin[rnk[lin] == "phylum"]; if (length(p)) sciname[p[1]] else NA_character_ },
                lineage_ids = list(lin))]
  }
  # comparison group = highest ancestor that is itself a downloaded group
  ids_of_groups <- out[!is.na(taxid), setNames(taxid, group)]
  out[, comparison_group := vapply(seq_len(.N), function(i) {
    lin <- lineage_ids[[i]]
    if (!length(lin)) return(group[i])
    hit <- lin[lin %in% ids_of_groups]
    if (!length(hit)) return(group[i])
    names(ids_of_groups)[match(tail(hit, 1), ids_of_groups)]
  }, character(1))]
  out[, lineage_ids := NULL]
  out[]
}

if (file.exists(F_SOURCES)) {
  so <- fread(F_SOURCES, sep = "\t", quote = "")[grepl("BOLD|NCBI GenBank", source)]   # per-group downloads only
  so[, n_records := suppressWarnings(as.numeric(n_records))]
  so <- so[!is.na(n_records)]
  if (nrow(so)) {
    so[, source := fifelse(grepl("BOLD", source), "BOLD", "NCBI GenBank")]
    gt <- resolve_groups(unique(so$group), TAXO_DIR)
    # BOLD and NCBI name some groups differently and the NCBI taxonomy does not
    # always join them (Crustacea is paraphyletic and no longer an ancestor of the
    # crustacean classes): fixed correspondence for the groups of these scripts
    DEFAULT_MAP <- c(Malacostraca = "Crustacea", Copepoda = "Crustacea", Ostracoda = "Crustacea",
                     Thecostraca = "Crustacea", Branchiopoda = "Crustacea", Cephalocarida = "Crustacea",
                     Remipedia = "Crustacea", Ascidiacea = "Tunicata", Thaliacea = "Tunicata",
                     Appendicularia = "Tunicata", Leptocardii = "Cephalochordata", Merostomata = "Xiphosura",
                     Bacillariophyta = "Ochrophyta", Phaeophyceae = "Ochrophyta",
                     # rbcL: NCBI downloads ochrophyte classes, BOLD the whole phylum
                     Pelagophyceae = "Ochrophyta", Dictyochophyceae = "Ochrophyta", Chrysophyceae = "Ochrophyta",
                     Raphidophyceae = "Ochrophyta", Xanthophyceae = "Ochrophyta", Eustigmatophyceae = "Ochrophyta",
                     Bolidophyceae = "Ochrophyta", Cryptophyceae = "Cryptophyta",
                     Zosteraceae = "Seagrasses", Posidoniaceae = "Seagrasses", Cymodoceaceae = "Seagrasses")
    gt[group %in% names(DEFAULT_MAP), comparison_group := DEFAULT_MAP[group]]
    if (file.exists(F_GROUP_FIX)) {           # manual corrections win
      fix <- fread(F_GROUP_FIX, sep = "\t", quote = "")
      gt[fix, on = "group", comparison_group := i.comparison_group]
      message("  (manual group corrections applied from ", F_GROUP_FIX, ")")
    }
    unresolved <- gt[is.na(taxid), group]
    if (length(unresolved)) warning("Group names not found in the NCBI taxonomy (kept as is): ",
                                    paste(unresolved, collapse = ", "), call. = FALSE)
    # phylum of each downloaded group, harmonised like the database (fixed for
    # the groups whose NCBI resolution is not a phylum)
    GROUP_PHYLUM <- c(Crustacea = "Arthropoda", Malacostraca = "Arthropoda", Copepoda = "Arthropoda",
                      Ostracoda = "Arthropoda", Thecostraca = "Arthropoda", Branchiopoda = "Arthropoda",
                      Cephalocarida = "Arthropoda", Remipedia = "Arthropoda", Pycnogonida = "Arthropoda",
                      Xiphosura = "Arthropoda", Merostomata = "Arthropoda",
                      Tunicata = "Chordata", Ascidiacea = "Chordata", Thaliacea = "Chordata",
                      Appendicularia = "Chordata", Cephalochordata = "Chordata", Leptocardii = "Chordata",
                      Ochrophyta = "Ochrophyta", Bacillariophyta = "Ochrophyta", Phaeophyceae = "Ochrophyta",
                      Rhodophyta = "Rhodophyta",
                      # rbcL groups (seagrasses are flowering plants: Streptophyta in NCBI)
                      Pelagophyceae = "Ochrophyta", Dictyochophyceae = "Ochrophyta", Chrysophyceae = "Ochrophyta",
                      Raphidophyceae = "Ochrophyta", Xanthophyceae = "Ochrophyta", Eustigmatophyceae = "Ochrophyta",
                      Bolidophyceae = "Ochrophyta", Chlorophyta = "Chlorophyta", Haptophyta = "Haptophyta",
                      Cryptophyceae = "Cryptophyta", Glaucocystophyceae = "Glaucophyta",
                      Zosteraceae = "Streptophyta", Posidoniaceae = "Streptophyta", Cymodoceaceae = "Streptophyta",
                      Seagrasses = "Streptophyta")
    gt[group %in% names(GROUP_PHYLUM), phylum := GROUP_PHYLUM[group]]
    gt[is.na(phylum), phylum := comparison_group]
    gt[phylum %chin% names(PHYLUM_RENAME), phylum := PHYLUM_RENAME[phylum]]
    so <- gt[, .(group, comparison_group, phylum)][so, on = "group"]
    fwrite(merge(unique(so[, .(source, group)]), gt, by = "group")[order(comparison_group, source, group)],
           file.path(FIG_DIR, "data", "group_taxonomy.tsv"), sep = "\t")
    if (length(EXCL)) so <- so[!(phylum %chin% EXCL | group %chin% EXCL)]   # taxa not shown
    # rows = the phyla of figure 2, in the same groups and order; the rest = "Other animals"…
    so[, set := SET_LOOKUP$set[match(phylum, SET_LOOKUP$phylum)]]
    so[is.na(set), set := set_of(phylum)]
    so[, phylum_group := fifelse(phylum %chin% group_order, phylum, OTHER_LABEL[set])]
    agg <- so[, .(n_records = sum(n_records), groups = paste(sort(unique(group)), collapse = ", ")),
              by = .(phylum_group, source)]
    # every phylum of figure 2 has a row, even without download
    lv8 <- c(group_order, setdiff(unique(agg$phylum_group), group_order))
    full <- CJ(phylum_group = lv8, source = c("BOLD", "NCBI GenBank"))
    agg <- agg[full, on = c("phylum_group", "source")][is.na(n_records), n_records := 0]
    agg[, set := factor(fifelse(phylum_group %chin% names(SET_OF_GROUP), SET_OF_GROUP[phylum_group],
                                names(OTHER_LABEL)[match(phylum_group, OTHER_LABEL)]), levels = SETS)]
    agg <- agg[!is.na(set)]
    agg[, phylum_group := factor(phylum_group, levels = lv8)]
    agg[, source := factor(source, levels = c("NCBI GenBank", "BOLD"))]   # BOLD drawn on top
    tot <- agg[, sum(n_records), by = source]
    p <- ggplot(agg, aes(n_records, phylum_group, fill = source)) +
      geom_col(position = position_dodge(width = 0.78, preserve = "single"), width = 0.72) +
      scale_fill_manual(values = c(BOLD = COL$bold, `NCBI GenBank` = COL$ncbi), breaks = c("BOLD", "NCBI GenBank")) +
      scale_x_continuous(labels = label_comma(), expand = expansion(mult = c(0, 0.05))) +
      labs(title = sprintf("%s records downloaded: %s from BOLD, %s from NCBI", num(sum(agg$n_records)),
                           num(tot[source == "BOLD", sum(V1)]), num(tot[source == "NCBI GenBank", sum(V1)])),
           subtitle = "Raw records per phylum",
           x = "Records", y = NULL, caption = CAPTION) +
      theme_db() + theme(panel.grid.major.y = element_blank())
    p <- by_set(p)
    save_fig(p, 8, "downloads_by_group", 1.2 + 0.32 * length(group_order),
             agg[order(-n_records)],
             "Raw records downloaded from BOLD and NCBI GenBank, per phylum (same phyla and order as figure 2); the groups queried in each database and their phylum are in data/group_taxonomy.tsv.")
  }
}

# =============================================================================
# 9-10. Primer-binding sites (CRABS "amplification efficiency" figure)
# =============================================================================
# For a sample of final amplicons, the bases just outside the amplicon in the
# full deposited sequence are the primer-binding sites. Sequences deposited
# without them (recovered by the global alignment) are skipped.
IUPAC <- list(A = "A", C = "C", G = "G", T = "T", R = c("A","G"), Y = c("C","T"), S = c("C","G"),
              W = c("A","T"), K = c("G","T"), M = c("A","C"), B = c("C","G","T"), D = c("A","G","T"),
              H = c("A","C","T"), V = c("A","C","G"), N = c("A","C","G","T"), I = c("A","C","G","T"))
revcomp <- function(x) vapply(strsplit(chartr("ACGTRYKMBDHVN", "TGCAYRMKVHDBN", x), ""),
                              function(v) paste(rev(v), collapse = ""), character(1))
primers <- list()
if (file.exists(F_PARAMS)) {
  pr <- fread(F_PARAMS, sep = "\t", quote = "")
  getp <- function(k) { v <- pr[parameter == k, value]; if (length(v)) toupper(sub(" .*$", "", v[1])) else "" }
  primers <- list(forward = getp("primer_forward"), reverse = getp("primer_reverse"))
}
if (file.exists(F_RAW) && length(primers) && all(nzchar(unlist(primers)))) {
  FWD <- primers$forward; REV <- primers$reverse; lf <- nchar(FWD); lr <- nchar(REV)
  smp <- db[, .(seq_id, amp = toupper(sequence), phylum_group, set)]
  if (nrow(smp) > MAX_PRIMER_SEQS) { set.seed(1); smp <- smp[sample(.N, MAX_PRIMER_SEQS)] }
  ncol_raw <- length(strsplit(readLines(F_RAW, n = 1), "\t")[[1]])
  raw <- fread(F_RAW, sep = "\t", quote = "", header = FALSE, select = c(1, ncol_raw),
               colClasses = "character", col.names = c("seq_id", "raw"))[seq_id %chin% smp$seq_id]
  smp <- raw[smp, on = "seq_id", nomatch = 0]; rm(raw)
  smp[, raw := toupper(raw)]
  if (!nrow(smp)) smp[, `:=`(start = integer(0))]
  if (nrow(smp)) smp[, start := mapply(function(a, r) regexpr(a, r, fixed = TRUE)[1], amp, raw, USE.NAMES = FALSE)]
  smp[, start := as.integer(start)]
  smp[, end := start + nchar(amp) - 1L]
  smp[, fsite := fifelse(start > lf, substr(raw, start - lf, start - 1L), NA_character_)]
  smp[, rsite := fifelse(start > 0 & end + lr <= nchar(raw), substr(raw, end + 1L, end + lr), NA_character_)]
  smp[!is.na(rsite), rsite := revcomp(rsite)]          # reverse primer read 5'->3'
  site_table <- function(sites, grp, primer, label) {
    keep <- !is.na(sites); sites <- sites[keep]; grp <- grp[keep]
    if (!length(sites)) return(NULL)
    m <- do.call(rbind, strsplit(sites, ""))
    pb <- strsplit(primer, "")[[1]]
    rbindlist(lapply(seq_along(pb), function(j) {
      b <- m[, j]; ok <- b %chin% c("A","C","G","T")
      data.table(primer = label, position = j, primer_base = pb[j], base = b, ok = ok,
                 mismatch = ok & !(b %chin% IUPAC[[pb[j]]]), phylum_group = grp,
                 set = factor(SET_OF_GROUP[as.character(grp)], levels = SETS))
    }))
  }
  st <- rbind(site_table(smp$fsite, smp$phylum_group, FWD, "Forward primer"),
              site_table(smp$rsite, smp$phylum_group, REV, "Reverse primer"))
  if (!is.null(st) && nrow(st)) {
    st[, primer := factor(primer, levels = c("Forward primer", "Reverse primer"))]
    # x key unique per primer (both primers can share a "4 A" label); shown without the prefix
    st[, pos_lab := factor(sprintf("%s%02d", substr(primer, 1, 1), position))]
    lab_map <- unique(st[, .(pos_lab = as.character(pos_lab), txt = paste0(position, "\n", primer_base))])
    lab_map <- setNames(lab_map$txt, lab_map$pos_lab)
    x_labs  <- function(k) unname(lab_map[k])
    n_f <- sum(!is.na(smp$fsite)); n_r <- sum(!is.na(smp$rsite))
    # 9. base composition per position (as in CRABS), mismatches outlined
    comp9 <- st[, .(n = .N), by = .(primer, position, pos_lab, primer_base, base = fifelse(ok, base, "other"))]
    comp9[, share := n / sum(n), by = .(primer, position)]
    comp9[, base := factor(base, levels = c("A", "C", "G", "T", "other"))]
    mm9 <- st[ok == TRUE, .(mismatch = mean(mismatch)), by = .(primer, position, pos_lab, primer_base)]
    worst <- mm9[which.max(mismatch)]
    p <- ggplot(comp9, aes(pos_lab, share, fill = base)) +
      geom_col(width = 0.82, colour = COL$surface, linewidth = 0.25) +
      geom_text(data = mm9[mismatch >= 0.05], aes(pos_lab, 1.04, label = pct(mismatch)),
                inherit.aes = FALSE, size = 2.6, colour = COL$ink2) +
      facet_wrap(~ primer, ncol = 1, scales = "free_x") + scale_x_discrete(labels = x_labs) +
      scale_fill_manual(values = c(setNames(viridisLite::viridis(4, begin = VIR_BEGIN, end = VIR_END), c("A","C","G","T")),
                                   other = COL$grid)) +
      scale_y_continuous(labels = label_percent(), breaks = c(0, 0.5, 1), expand = expansion(mult = c(0, 0.1))) +
      labs(title = sprintf("Worst primer position: %s mismatches at %s position %d", pct(worst$mismatch),
                           tolower(sub(" primer", "", worst$primer)), worst$position),
           subtitle = sprintf("Base found at each primer position (5' -> 3'; primer base under the number)\n%% above a bar = sequences mismatching the primer (if >= 5%%)\nsites read: %s fwd and %s rev",
                              num(n_f), num(n_r)),
           x = NULL, y = "Sequences", caption = CAPTION) +
      theme_db() + theme(panel.grid.major.x = element_blank(), axis.text.x = element_text(size = 7.5))
    save_fig(p, 9, "primer_sites_composition", 5.8, comp9,
             "Base composition of the primer-binding sites in the full deposited sequences (CRABS 'amplification efficiency'); percentages = share of sequences mismatching the primer at that position.")
    # 10. mismatch rate per phylum and position (heatmap)
    mm10 <- st[ok == TRUE & !is_other(phylum_group), .(mismatch = mean(mismatch), n = .N),
               by = .(primer, position, pos_lab, set, phylum_group)]
    p <- ggplot(mm10, aes(pos_lab, phylum_group, fill = mismatch)) +
      geom_tile(colour = COL$surface, linewidth = 0.5) +
      facet_wrap(~ primer, ncol = 1, scales = "free_x") + scale_x_discrete(labels = x_labs) +
      scale_fill_viridis_c(labels = label_percent(), limits = c(0, NA), name = "Mismatch") +
      labs(title = "Where each phylum mismatches the primers",
           subtitle = "Share of sequences whose base differs from the primer, by position (5'->3') and phylum",
           x = NULL, y = NULL, caption = CAPTION) +
      theme_db() + theme(panel.grid = element_blank(), axis.text.x = element_text(size = 7.5),
                         legend.position = "right", legend.title = element_text(colour = COL$ink2, size = 9))
    save_fig(p, 10, "primer_mismatch_by_phylum", 1.8 + 0.5 * uniqueN(mm10$phylum_group), mm10,
             "Mismatch rate between the primer and the binding site of each sequence, per position and phylum.")
  }
  rm(smp, st)
} else message("  (primer-site figures skipped: needs ", F_RAW, " and the primers in ", F_PARAMS, ")")

# =============================================================================
# 11. Completeness for species of interest (CRABS --completeness-table)
# =============================================================================
if (!is.na(F_SPECIES) && !is.null(load_taxonomy(TAXO_DIR))) {
  soi <- unique(trimws(readLines(F_SPECIES, warn = FALSE))); soi <- soi[nzchar(soi) & !startsWith(soi, "#")]
  anc_of_rank <- function(ids, target) {
    cur <- ids; res <- rep(NA_integer_, length(ids))
    for (k in 1:60) {
      ok <- cur > 0 & cur <= length(TAX$rank)
      hit <- is.na(res) & ok & TAX$rank[pmax(cur, 1L)] == target
      res[hit] <- cur[hit]
      nxt <- TAX$parent[pmax(cur, 1L)]; nxt[!ok | nxt == cur] <- 1L; cur <- nxt
      if (all(!is.na(res) | cur == 1L)) break
    }
    res
  }
  sp_ids <- which(TAX$rank == "species")
  informal <- grepl(" sp\\.| cf\\.| aff\\.|uncultured|environmental|unidentified", TAX$sciname[sp_ids])
  sp_tbl <- data.table(taxid = sp_ids[!informal])
  sp_tbl[, `:=`(genus = anc_of_rank(taxid, "genus"), family = anc_of_rank(taxid, "family"))]
  ncbi_g <- sp_tbl[!is.na(genus), .N, by = genus]; ncbi_f <- sp_tbl[!is.na(family), .N, by = family]
  db_sp <- unique(db[!is.na(get(species_col)), .(sp = get(species_col), genus, family)])
  ct <- rbindlist(lapply(soi, function(sp) {
    hit <- TAX$names[name == sp][order(class != "scientific name")][1]
    tid <- if (nrow(hit) && !is.na(hit$taxid)) hit$taxid else NA_integer_
    g <- if (!is.na(tid)) anc_of_rank(tid, "genus") else NA_integer_
    fa <- if (!is.na(tid)) anc_of_rank(tid, "family") else NA_integer_
    gname <- if (!is.na(g)) TAX$sciname[g] else sub(" .*", "", sp)
    fname <- if (!is.na(fa)) TAX$sciname[fa] else NA_character_
    in_g <- db_sp[genus == gname, sort(unique(sp))]; in_f <- if (!is.na(fname)) db_sp[family == fname, sort(unique(sp))] else character(0)
    ng <- if (!is.na(g)) ncbi_g[genus == g, N] else NA_integer_; nf <- if (!is.na(fa)) ncbi_f[family == fa, N] else NA_integer_
    ng <- if (length(ng)) ng else NA_integer_; nf <- if (length(nf)) nf else NA_integer_
    data.table(species = sp, n_barcodes = db[get(species_col) == sp, .N], genus = gname, family = fname,
               db_species_in_genus = length(in_g), ncbi_species_in_genus = ng,
               pct_genus = if (!is.na(ng) && ng > 0) length(in_g) / ng else NA_real_,
               db_species_in_family = length(in_f), ncbi_species_in_family = nf,
               pct_family = if (!is.na(nf) && nf > 0) length(in_f) / nf else NA_real_,
               species_in_genus = paste(in_g, collapse = ", "), species_in_family = paste(in_f, collapse = ", "),
               found_in_ncbi_taxonomy = !is.na(tid))
  }))
  fwrite(ct, file.path(FIG_DIR, "data", "completeness_table.tsv"), sep = "\t")
  long11 <- rbind(ct[, .(species, panel = "Barcodes", value = n_barcodes, lab = num(n_barcodes))],
                  ct[, .(species, panel = "Congeners (% of genus)", value = pct_genus,
                         lab = fifelse(is.na(pct_genus), "n/a", sprintf("%s (%d/%s)", pct(pct_genus), db_species_in_genus, ncbi_species_in_genus)))],
                  ct[, .(species, panel = "Confamilials (% of family)", value = pct_family,
                         lab = fifelse(is.na(pct_family), "n/a", sprintf("%s (%d/%s)", pct(pct_family), db_species_in_family, ncbi_species_in_family)))])
  long11[, panel := factor(panel, levels = unique(panel))]
  long11[, species := factor(species, levels = ct[order(n_barcodes, species), species])]
  n_missing <- ct[n_barcodes == 0, .N]
  p <- ggplot(long11, aes(value, species)) +
    geom_col(data = long11[!is.na(value)], aes(fill = panel), width = 0.66, show.legend = FALSE) +
    geom_text(aes(x = fifelse(is.na(value), 0, value), label = lab), hjust = -0.1, size = 2.9, colour = COL$ink2) +
    facet_wrap(~ panel, scales = "free_x") +
    scale_fill_manual(values = setNames(viridisLite::viridis(3, begin = 0.2, end = 0.7), levels(long11$panel))) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.6)), labels = function(x) ifelse(x <= 1 & max(x, na.rm = TRUE) <= 1, label_percent()(x), label_comma()(x))) +
    labs(title = if (n_missing) sprintf("%d of %d species of interest have no barcode in the database", n_missing, nrow(ct))
         else sprintf("All %d species of interest have at least one barcode", nrow(ct)),
         subtitle = "Barcodes per species, and share of the NCBI species of its genus and family present in the database",
         x = NULL, y = NULL, caption = CAPTION) +
    theme_db() + theme(panel.grid.major.y = element_blank(), panel.spacing.x = unit(14, "pt"),
                       strip.text = element_text(size = 9))
  save_fig(p, 11, "completeness_species_of_interest", 1.4 + 0.32 * nrow(ct), ct,
           "Completeness for the species of interest (CRABS completeness table): barcodes per species and share of the NCBI species of its genus / family present in the database. NCBI counts exclude informal names (sp., cf., aff.).")
} else message("  (completeness figure skipped: add species_of_interest.txt, one species per line, in the database folder or next to the script)")

# =============================================================================
# 12. Taxonomic resolution: can species be told apart with this fragment?
# =============================================================================
# A species is "resolved" when at least one of its sequences is shared with no
# other species; "not resolved" when every sequence it has is identical to a
# sequence of another species (reads matching it cannot be assigned to species).
res_dt <- db[!is.na(get(species_col)) & !grepl(" (sp|cf|aff)\\.( |$)", get(species_col)),
             .(sp = get(species_col), seq = toupper(sequence), phylum_group, set)]
if (nrow(res_dt)) {
  res_dt[, n_sp := uniqueN(sp), by = seq]
  per_sp <- res_dt[, .(unique_seqs = sum(n_sp == 1), shared_seqs = sum(n_sp > 1)), by = .(sp, set, phylum_group)]
  per_sp[, status := fifelse(shared_seqs == 0, "Resolved: all sequences unique",
                             fifelse(unique_seqs > 0, "Partly: some sequences shared", "Not resolved: all sequences shared"))]
  lv12 <- c("Resolved: all sequences unique", "Partly: some sequences shared", "Not resolved: all sequences shared")
  per_sp[, status := factor(status, levels = rev(lv12))]
  agg12 <- per_sp[, .N, by = .(set, phylum_group, status)]
  agg12[, share := N / sum(N), by = phylum_group]
  ok12 <- per_sp[, mean(unique_seqs > 0)]
  p <- ggplot(agg12, aes(share, phylum_group, fill = status)) +
    geom_col(width = 0.68, colour = COL$surface, linewidth = 0.4) +
    scale_fill_manual(values = setNames(rev(viridisLite::viridis(3, begin = 0.15, end = 0.8)), lv12),
                      breaks = lv12) +
    scale_x_continuous(labels = label_percent(), expand = expansion(mult = c(0, 0.02))) +
    guides(fill = guide_legend(nrow = 2)) +
    labs(title = sprintf("%s of species have at least one barcode no other species shares", pct(ok12)),
         subtitle = sprintf("Identical amplicons shared between species, per phylum (%s species; exact identity,\nso resolution at the usual 97-99%% assignment thresholds is lower)", num(nrow(per_sp))),
         x = "Species", y = NULL, caption = CAPTION) +
    theme_db() + theme(panel.grid.major.y = element_blank())
  p <- by_set(p)
  save_fig(p, 12, "species_resolution", 1.6 + 0.32 * length(group_order), per_sp[order(phylum_group, sp)],
           "Share of species whose amplicons are unique to them, partly shared or entirely shared with other species (exact identity).")
  rm(res_dt)
}

# =============================================================================
# 13. Focal taxa (sessile / pioneer fauna): species and sequences by source
# =============================================================================
lin_cols <- setdiff(intersect(RANKS, names(db)), c("species", "genus"))
db[, focal := NA_character_]
for (f in FOCAL_TAXA) {
  al <- if (!is.null(FOCAL_ALIASES[[f]])) FOCAL_ALIASES[[f]] else f
  hit <- Reduce(`|`, lapply(lin_cols, function(r) db[[r]] %chin% al))
  db[is.na(focal) & hit, focal := f]
}
fo_lv <- FOCAL_TAXA
src_col <- if (HAS_SOURCE) "source" else NULL
fo <- db[!is.na(focal), .(Species = uniqueN(get(species_col)), Sequences = .N), by = c("focal", src_col)]
if (!nrow(fo)) {
  message("  fig13_focal_taxa: skipped (none of the focal taxa in this database)")
} else {
  fo_long <- melt(fo, id.vars = c("focal", src_col), variable.name = "measure", value.name = "n")
  # sorted by number of sequences (largest at the top)
  fo_lv <- fo[, .(N = sum(Sequences)), by = focal][order(N), focal]
  fo_long[, focal := factor(focal, levels = fo_lv)]
  fo_long[, measure := factor(measure, levels = c("Species", "Sequences"))]
  n_focal_sp <- db[!is.na(focal), uniqueN(get(species_col))]
  p <- ggplot(fo_long, aes(n, focal))
  p <- if (HAS_SOURCE) {
    p + geom_col(aes(fill = source), width = 0.68, colour = COL$surface, linewidth = 0.3) +
      scale_fill_manual(values = c(BOLD = COL$bold, `NCBI GenBank` = COL$ncbi, SILVA = COL$silva, Other = COL$muted),
                        breaks = intersect(c("BOLD", "NCBI GenBank", "SILVA", "Other"), unique(fo_long$source)))
  } else p + geom_col(fill = COL$s1, width = 0.68)
  p <- p + facet_wrap(~ measure, scales = "free_x") +
    scale_x_continuous(labels = label_comma(), expand = expansion(mult = c(0, 0.05))) +
    labs(title = sprintf("%s species of the %d focal taxa in the database", num(n_focal_sp), length(fo_lv)),
         subtitle = FOCAL_SUBTITLE,
         x = NULL, y = NULL, caption = CAPTION) +
    theme_db() + theme(panel.grid.major.y = element_blank(), panel.spacing.x = unit(18, "pt"))
  save_fig(p, 13, "focal_taxa", 1.2 + 0.34 * length(fo_lv), dcast(fo_long, focal + measure ~ ..., value.var = "n", fill = 0),
           sprintf("Species and sequences of the focal taxa (%s) in the final database, by source. %s",
                   paste(FOCAL_TAXA, collapse = ", "), FOCAL_NOTE))
}

# =============================================================================
# 14. Provenance (FAIR)
# =============================================================================
inputs <- c(F_FINAL, F_COUNTS, F_SOURCES, F_PARAMS, F_INSILICO, F_RAW, F_BOLD, F_NCBI, F_SILVA, if (!is.na(F_SPECIES)) F_SPECIES)
inputs <- inputs[file.exists(inputs)]
md5 <- tools::md5sum(inputs)
pkgs <- vapply(needed, function(p) as.character(packageVersion(p)), character(1))
readme <- c(
  sprintf("# Figures – %s", BASE_NAME), "",
  sprintf("- **Created**: %s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  sprintf("- **Script**: crabs_figures.R %s", SCRIPT_VERSION),
  sprintf("- **Software**: %s; %s", R.version.string, paste(names(pkgs), pkgs, collapse = ", ")),
  "- **Formats**: PNG (300 dpi) for slides and drafts, PDF (vector) for publication; the data behind each figure is in `data/`.", "",
  "## Figures", "", "| Figure | Caption | Data |", "|---|---|---|",
  vapply(FIGS, function(f) sprintf("| `%s.png` / `.pdf` | %s | `data/%s.csv` |", f$stem, f$caption, f$stem), character(1)),
  "", "## Inputs (MD5)", "", "| File | MD5 |", "|---|---|",
  sprintf("| `%s` | `%s` |", names(md5), unname(md5)), "",
  "Full R session details: `sessionInfo.txt`.")
writeLines(readme, file.path(FIG_DIR, "README.md"))
writeLines(capture.output(sessionInfo()), file.path(FIG_DIR, "sessionInfo.txt"))
message("Done: ", length(FIGS), " figures in ", normalizePath(FIG_DIR))