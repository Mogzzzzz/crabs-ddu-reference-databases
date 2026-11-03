# COI reference database (~313 bp) – marine invertebrates and algae (COI_marine_invertebrates_algae)

- **Created**: 2026-10-28 (run 2026-10-28T17:17:29+0100 → 2026-10-28T17:18:50+0100)
- **Contact**: morgane.durand@mnhn.fr
- **Content**: 364602 sequences, 77086 taxa; length min / median / max: 250 / 313 / 340 bp
- **Primers**: GGWACWGGWTGAACWGTWTAYCCYCC (mlCOIintF) / TANACYTCNGGRTGNCCRAARAAYCA (jgHCO2198)

## Files

| File | Content |
|---|---|
| `results/COI_marine_invertebrates_algae.fasta` | Final database, headers `taxid=…; scientific_name=…; lineage=…;` (OBITools) |
| `results/COI_marine_invertebrates_algae.tsv` | Same content as a table with column names |
| `results/COI_marine_invertebrates_algae.txt` | CRABS format (usable by every CRABS function) |
| `results/insilico_pcr_raw.fasta` | Raw in silico PCR output (before alignment and filters) |
| `metadata/run_metadata.json` | All metadata, machine-readable |
| `metadata/parameters.tsv`, `versions.tsv`, `sources.tsv`, `counts.tsv` | Parameters, software, sources, counts per step |
| `metadata/checksums_sha256.txt` | Checksums of the result files |
| `logs/` | Log of every step |

## Methods

1. Download from BOLD (V5 portal, 37 groups) and NCBI GenBank (60 groups), plus the NCBI taxonomy.
2. Import into CRABS format with NCBI lineage (domain;phylum;class;order;family;genus;species).
3. Merge (duplicate accessions removed); lineage rescue (a species name unknown to the NCBI taxonomy takes the lineage of its genus); dereplication (unique_species), removal of environmental sequences and of sequences not identified at least to genus (or with more than 1 missing rank above it).
4. In silico PCR (cutadapt, 4 mismatches, relaxed=true).
5. Pairwise global alignment (VSEARCH, identity ≥ 0.95, coverage ≥ 95 %) to recover amplicons deposited without primer-binding sites.
6. Final dereplication and filters: 280–340 bp, fewer than 2 N; amplicons of 250-279 bp kept only for taxa without a full-length amplicon; amplicons with a stop codon in every reading frame removed (yes (stop codons: TAA TAG)); taxa removed: Insecta;Collembola;Protura;Diplura;Arachnida;Chilopoda;Diplopoda;Pauropoda;Symphyla;Onychophora.
Genus homonyms: yes (genus homonyms: BOLD phylum, then target groups + Ochrophyta;Bacillariophyta;Phaeophyceae;Ichthyostraca;Mystacocarida; log logs/10_taxonomy_fix.tsv); manual decisions: taxonomy_manual_crabs.tsv (5 decisions): Caledoniella=keep; Roperia=remove; Trichia=remove; Anoplostrongylus=remove; Chondrus=remove.
7. Phylum names harmonised with Adl et al. (2019): renamed Bacillariophyta=Ochrophyta;Heterokontophyta=Ochrophyta;Myzozoa=Dinoflagellata;Dinophyta=Dinoflagellata; phylum given from the class when NCBI has none (Phaeophyceae=Ochrophyta;Bacillariophyceae=Ochrophyta;Coscinodiscophyceae=Ochrophyta;Mediophyceae=Ochrophyta;Chrysophyceae=Ochrophyta;Synurophyceae=Ochrophyta;Raphidophyceae=Ochrophyta;Eustigmatophyceae=Ochrophyta;Phaeosacciophyceae=Ochrophyta;Pelagophyceae=Ochrophyta;Dictyochophyceae=Ochrophyta;Xanthophyceae=Ochrophyta;Bolidophyceae=Ochrophyta;Pinguiophyceae=Ochrophyta;Phaeothamniophyceae=Ochrophyta;Aurearenophyceae=Ochrophyta;Chrysomerophyceae=Ochrophyta;Schizocladiophyceae=Ochrophyta;Olisthodiscophyceae=Ochrophyta;Dinophyceae=Dinoflagellata). Details in metadata/taxonomy_harmonisation.tsv; species names and NCBI taxids unchanged.

## Sequences per step

| Step | Sequences |
|---|---|
| import_bold | 1103382 |
| import_ncbi | 4105620 |
| 01_merged | 5209002 |
| 01b_lineage_rescued | 5209002 |
| 02_dereplicated_initial | 3124044 |
| 03_filtered_initial | 3032845 |
| 04_insilico_pcr | 2697694 |
| 05_global_alignment | 2719512 |
| 05b_long_trimmed | 2719512 |
| 06_dereplicated_final | 2058282 |
| 07_filtered_final | 1400740 |
| 07b_length_rescued | 1178378 |
| 07c_translation_checked | 1177458 |
| 07d_taxonomy_fixed | 1177360 |
| 08_taxa_excluded | 401207 |
| 09_min_rank | 364602 |
| COI_marine_invertebrates_algae | 364602 |

## Software

| Software | Version |
|---|---|
| script | 1.7.4 (2026-10-28, shared CRABS/taxo, CRABS/reference_database_creator, CRABS/scripts) |
| CRABS | 1.14.0 (commit 97b3037) |
| cutadapt | 5.2 |
| vsearch | vsearch v2.32.0_macos_aarch64 |
| makeblastdb | makeblastdb: 2.17.0+ |
| python | Python 3.13.4 |
| system | Darwin 25.5.0 arm64 |

## Please cite

- Jeunen et al. (2022) CRABS. Mol Ecol Resour. doi:10.1111/1755-0998.13741
- Ratnasingham & Hebert (2007) BOLD. Mol Ecol Notes. doi:10.1111/j.1471-8286.2007.01678.x
- Sayers et al. NCBI GenBank / Taxonomy (https://www.ncbi.nlm.nih.gov)
- Adl et al. (2019) Revisions to the classification, nomenclature, and diversity of eukaryotes. J Eukaryot Microbiol 66: 4-119. doi:10.1111/jeu.12691
- Leray et al. (2013) Front Zool. doi:10.1186/1742-9994-10-34 (mlCOIintF)
- Geller et al. (2013) Mol Ecol Resour. doi:10.1111/1755-0998.12138 (jgHCO2198)
- Martin (2011) cutadapt. EMBnet.journal. doi:10.14806/ej.17.1.200
- Rognes et al. (2016) VSEARCH. PeerJ. doi:10.7717/peerj.2584

Source details (download date, number of records, exact query) are in `metadata/sources.tsv`.
