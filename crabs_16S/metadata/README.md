# Mitochondrial 16S rRNA reference database (180–480 bp, median 290 bp) – marine invertebrates (16S_mt_marine_invertebrates)

- **Created**: 2026-10-26 (run 2026-10-26T14:26:59+0100 → 2026-10-26T14:27:53+0100)
- **Contact**: morgane.durand@mnhn.fr
- **Content**: 34880 sequences, 13118 taxa; length min / median / max: 180 / 287 / 480 bp
- **Primers**: AYAAGACGAGAAGACCC / CCGGTCTGAACTCAGATCACGT (16Sbr, Palumbi 1996)

## Files

| File | Content |
|---|---|
| `results/16S_mt_marine_invertebrates.fasta` | Final database, headers `taxid=…; scientific_name=…; lineage=…;` (OBITools) |
| `results/16S_mt_marine_invertebrates.tsv` | Same content as a table with column names |
| `results/16S_mt_marine_invertebrates.txt` | CRABS format (usable by every CRABS function) |
| `results/insilico_pcr_raw.fasta` | Raw in silico PCR output (before alignment and filters) |
| `metadata/run_metadata.json` | All metadata, machine-readable |
| `metadata/parameters.tsv`, `versions.tsv`, `sources.tsv`, `counts.tsv` | Parameters, software, sources, counts per step |
| `metadata/checksums_sha256.txt` | Checksums of the result files |
| `logs/` | Log of every step |

## Methods

1. Download from BOLD (V5 portal, 35 groups, records with marker /16S/ only) and NCBI GenBank (53 groups), plus the NCBI taxonomy.
2. Import into CRABS format with NCBI lineage (domain;phylum;class;order;family;genus;species).
3. Merge (duplicate accessions removed); lineage rescue (a species name unknown to the NCBI taxonomy takes the lineage of its genus); removal of the excluded taxa; dereplication (unique_species), removal of environmental sequences and of sequences not identified at least to genus (or with more than 1 missing rank above it).
4. In silico PCR (cutadapt, 4 mismatches, relaxed=false).
5. Pairwise global alignment (VSEARCH, identity ≥ 0.95, coverage ≥ 95 %) to recover amplicons deposited without primer-binding sites.
6. Final dereplication and filters: 180–480 bp, fewer than 2 N; taxa removed: Insecta;Collembola;Protura;Diplura;Arachnida;Chilopoda;Diplopoda;Pauropoda;Symphyla;Onychophora (removed after import, after the final filters, and excluded from the NCBI queries).
Genus homonyms: yes (genus homonyms: BOLD phylum, then target groups + Ochrophyta;Bacillariophyta;Phaeophyceae;Ichthyostraca;Mystacocarida; log logs/10_taxonomy_fix.tsv).
7. Phylum names harmonised with Adl et al. (2019): renamed Bacillariophyta=Ochrophyta;Heterokontophyta=Ochrophyta;Myzozoa=Dinoflagellata;Dinophyta=Dinoflagellata; phylum given from the class when NCBI has none (Phaeophyceae=Ochrophyta;Bacillariophyceae=Ochrophyta;Coscinodiscophyceae=Ochrophyta;Mediophyceae=Ochrophyta;Chrysophyceae=Ochrophyta;Synurophyceae=Ochrophyta;Raphidophyceae=Ochrophyta;Eustigmatophyceae=Ochrophyta;Phaeosacciophyceae=Ochrophyta;Pelagophyceae=Ochrophyta;Dictyochophyceae=Ochrophyta;Xanthophyceae=Ochrophyta;Bolidophyceae=Ochrophyta;Pinguiophyceae=Ochrophyta;Phaeothamniophyceae=Ochrophyta;Aurearenophyceae=Ochrophyta;Chrysomerophyceae=Ochrophyta;Schizocladiophyceae=Ochrophyta;Olisthodiscophyceae=Ochrophyta;Dinophyceae=Dinoflagellata). Details in metadata/taxonomy_harmonisation.tsv; species names and NCBI taxids unchanged.

## Sequences per step

| Step | Sequences |
|---|---|
| import_bold | 18541 |
| import_ncbi | 263354 |
| 01_merged | 281895 |
| 01b_lineage_rescued | 281895 |
| 01c_taxa_excluded_early | 281895 |
| 02_dereplicated_initial | 184719 |
| 03_filtered_initial | 181591 |
| 04_insilico_pcr | 30266 |
| 05_global_alignment | 52505 |
| 05b_long_trimmed | 52505 |
| 06_dereplicated_final | 36438 |
| 07_filtered_final | 35608 |
| 07b_length_rescued | 35608 |
| 07c_translation_checked | 35608 |
| 07d_taxonomy_fixed | 35608 |
| 08_taxa_excluded | 35608 |
| 09_min_rank | 34880 |
| 16S_mt_marine_invertebrates | 34880 |

## Software

| Software | Version |
|---|---|
| script | 1.4 (2026-10-26, from crabs_coi.sh 1.8 + taxonomy fix of 1.7.2) |
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
- Palumbi (1996) Nucleic acids II: the polymerase chain reaction. In Molecular Systematics, Sinauer (16Sbr)
- Martin (2011) cutadapt. EMBnet.journal. doi:10.14806/ej.17.1.200
- Rognes et al. (2016) VSEARCH. PeerJ. doi:10.7717/peerj.2584

Source details (download date, number of records, exact query) are in `metadata/sources.tsv`.
