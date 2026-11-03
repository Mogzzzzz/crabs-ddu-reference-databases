# CRABS reference databases for Antarctic benthic metabarcoding (DDU)

Curated, reproducible reference databases for the taxonomic assignment of
metabarcoding reads from Antarctic benthic communities (Terre Adélie, Dumont
d'Urville – ARMS, REVOLTA 2014-2017), built with
[CRABS](https://github.com/gjeunen/reference_database_creator) v1.14.0 from BOLD,
NCBI GenBank and SILVA. Five databases, one per marker:

| Database | Region | Primers | Sequences | Taxa | Length (bp) |
|---|---|---|---|---|---|
| `crabs_coi` | COI (mitochondrial) | mlCOIintF / jgHCO2198 | 364,602 | 77,086 | 250–340 (median 313) |
| `crabs_16S` | 16S rRNA (mitochondrial) | Sarri_F / 16Sbr-H | 34,880 | 13,118 | 180–480 |
| `crabs_18S/18S_aligned` | 18S rRNA V4–V5 (nuclear) | 18Smk1F / 18Smk2R | 59,079 | 36,331 | 450–900 |
| `crabs_18S/18S_forward` | 18S rRNA V4, first 150 bp after 18Smk1F | 18Smk1F | 67,746 | 41,959 | 145–155 |
| `crabs_rbcl` | rbcL (plastid) | rbcL-646F / rbcL-998R | 13,737 | 6,428 | 310–346 (median 331) |

Target groups: marine invertebrates and algae (COI, 16S, rbcL); all eukaryotes
through SILVA for the 18S. The final databases themselves are not stored in
this repository (too large for GitHub): run the scripts to rebuild them, or
download the released files (see *Data availability*).

## Repository layout

```
CRABS/
├── scripts/                     shared scripts (all markers)
│   ├── crabs_figures.R          figures of every database
│   ├── crabs_check_taxonomy.sh  genus-homonym check (BOLD / SILVA vs NCBI phylum)
│   ├── crabs_check_taxa.sh      are these taxa in the databases?
│   ├── crabs_trace_species.sh   step where a species is lost
│   ├── crabs_18S_split.sh       split crabs_18S into 18S_aligned / 18S_forward
│   └── cluster/                 SLURM scripts: OBITools4 indexing and assignment
├── crabs_coi/      crabs_coi.sh            + taxonomy_manual_crabs.tsv, metadata/, figures/
├── crabs_16S/      crabs_16S.sh
├── crabs_18S/
│   ├── 18S_aligned/  crabs_18S.sh
│   └── 18S_forward/  crabs_18S_fwd.sh     (built from 18S_aligned)
└── crabs_rbcl/     crabs_rbcL.sh
```

Not in the repository (downloaded or rebuilt by the scripts): the CRABS
installation (`reference_database_creator/`), the NCBI taxonomy (`taxo/`),
raw downloads (`raw/`, `imported/`), intermediate steps (`steps/`), logs and
the final databases (`results/`).

## Requirements

- [CRABS](https://github.com/gjeunen/reference_database_creator) v1.14.0, cloned
  into `CRABS/reference_database_creator/` (or set `CRABS_DIR=/path`)
- Python ≥ 3.9, cutadapt, VSEARCH (used by CRABS and by the scripts)
- R ≥ 4.2 with data.table, ggplot2, scales, viridisLite (figures)
- OBITools4 4.4.0 (assignment, on a SLURM cluster)

## Usage

Each database script is run from its own folder and can be stopped and
relaunched at any time: finished steps are skipped (their output, inputs and
parameters are checked), so a rerun only redoes what changed.

```bash
cd CRABS/crabs_coi && bash crabs_coi.sh 2>&1 | tee -a crabs_coi.log      # first: downloads the NCBI taxonomy and BOLD
cd ../crabs_16S    && bash crabs_16S.sh 2>&1 | tee -a crabs_16S.log      # reuse the taxonomy and BOLD files of crabs_coi
cd ../crabs_rbcl   && bash crabs_rbcL.sh 2>&1 | tee -a crabs_rbcL.log
cd ../crabs_18S/18S_aligned && bash crabs_18S.sh 2>&1 | tee -a crabs_18S.log
cd ../18S_forward  && bash crabs_18S_fwd.sh 2>&1 | tee -a crabs_18S_fwd.log
cd ../../          && Rscript scripts/crabs_figures.R && bash scripts/crabs_check_taxonomy.sh
```

All parameters (primers, length windows, mismatches, target groups, excluded
taxa, minimum identification rank) are set at the top of each script and
written to `metadata/parameters.tsv` at every run.

## Pipeline (each marker)

1. Download: BOLD (by group), NCBI GenBank (marker query), SILVA (18S);
   NCBI taxonomy.
2. Cleaning: sources merged, dereplication, removal of environmental records
   and of records with an incomplete lineage.
3. Marker region: in silico PCR (cutadapt, max. 4 mismatches), recovery by
   global alignment of records deposited without primer sites, re-cut of
   over-long amplicons (inner primer sites, else alignment), dereplication.
4. Quality filters: length window and N, translation check (COI, rbcL).
5. Taxonomy: correction of genus homonyms across kingdoms (BOLD and SILVA
   phylum vs NCBI, manual decisions in `taxonomy_manual_crabs.tsv`), removal of
   organelle rRNA in the 18S, excluded taxa (e.g. insects), minimum rank (genus).
6. Phylum names harmonised with Adl et al. (2019); export to FASTA (OBITools
   headers), TSV and CRABS format, with metadata.

## Outputs and metadata (FAIR)

Each database folder holds, after a run:

- `results/<base>.fasta|.tsv|.txt` – final database
- `metadata/README.md`, `run_metadata.json` – content, sources, dates, versions
- `metadata/parameters.tsv`, `counts.tsv` (sequences per step), `sources.tsv`
  (records downloaded per group and date), `versions.tsv`, `checksums_sha256.txt`
- `logs/` – one log per step; `checks/` – taxonomy checks (timestamped)
- `figures/` – PNG, PDF and the data behind each figure, with a README (MD5 of inputs)

## Taxonomic assignment (OBITools4, SLURM)

`scripts/cluster/crabs_obitag_slurm_<marker>.sh` converts a database to
OBITools4 headers, updates merged NCBI taxids, indexes it (`obirefidx`) and
assigns the OTUs (`obitag`, table by `obicsv`).

## Data availability

The final databases will be deposited on Zenodo (DOI to come).

## Citation

See `CITATION.cff`. Please also cite CRABS (Jeunen et al. 2023, *Molecular
Ecology Resources* 23: 725–738), the data sources (BOLD, NCBI GenBank, SILVA)
and Adl et al. (2019, *J Eukaryot Microbiol* 66: 4–119).

## License

Code: MIT (see `LICENSE`).
