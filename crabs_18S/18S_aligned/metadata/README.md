# 18S rRNA V4-V5 reference database (~574 bp) – marine invertebrates and other eukaryotes (18S_V4V5_marine_eukaryotes)

- **Created**: 2026-10-28 (run 2026-10-28T17:18:56+0100 → 2026-10-28T17:19:26+0100)
- **Contact**: morgane.durand@mnhn.fr
- **Content**: 59079 sequences, 36331 taxa; length min / median / max: 450 / 553 / 900 bp
- **Primers**: CTGGTGCCAGCAGCCGCGGYAA (18Smk1F, Machida & Knowlton 2012) / TCCGTCAATTYCTTTAAGTT (18Smk2R, Machida & Knowlton 2012)

## Files

| File | Content |
|---|---|
| `results/18S_V4V5_marine_eukaryotes.fasta` | Final database, headers `taxid=…; scientific_name=…; lineage=…;` (OBITools) |
| `results/18S_V4V5_marine_eukaryotes.tsv` | Same content as a table with column names |
| `results/18S_V4V5_marine_eukaryotes.txt` | CRABS format (usable by every CRABS function) |
| `results/insilico_pcr_raw.fasta` | Raw in silico PCR output (before alignment and filters) |
| `metadata/run_metadata.json` | All metadata, machine-readable |
| `metadata/parameters.tsv`, `versions.tsv`, `sources.tsv`, `counts.tsv` | Parameters, software, sources, counts per step |
| `metadata/checksums_sha256.txt` | Checksums of the result files |
| `logs/` | Log of every step |

## Methods

1. Download from BOLD (V5 portal, 37 groups, records with marker /18S/ only) and NCBI GenBank (30 groups), SILVA SSU subset 138.2, plus the NCBI taxonomy.
2. Import into CRABS format with NCBI lineage (domain;phylum;class;order;family;genus;species).
3. Merge (duplicate accessions removed); lineage rescue (a species name unknown to the NCBI taxonomy takes the lineage of its genus); removal of the excluded taxa; dereplication (unique_species), removal of environmental sequences and of sequences not identified at least to genus (or with more than 1 missing rank above it).
4. In silico PCR (cutadapt, 4 mismatches, relaxed=false).
5. Pairwise global alignment (VSEARCH, identity ≥ 0.95, coverage ≥ 95 %) to recover amplicons deposited without primer-binding sites.
6. Final dereplication and filters: 450–900 bp, fewer than 2 N; taxa removed: Bacteria;Archaea;Insecta;Collembola;Protura;Diplura;Arachnida;Chilopoda;Diplopoda;Pauropoda;Symphyla;Onychophora (removed after import, after the final filters, and excluded from the NCBI queries).
Genus homonyms: yes (genus homonyms: BOLD and SILVA phylum; log logs/10_taxonomy_fix.tsv).
7. Phylum names harmonised with Adl et al. (2019): renamed Bacillariophyta=Ochrophyta;Heterokontophyta=Ochrophyta;Myzozoa=Dinoflagellata;Dinophyta=Dinoflagellata; phylum given from the class when NCBI has none (Phaeophyceae=Ochrophyta;Bacillariophyceae=Ochrophyta;Coscinodiscophyceae=Ochrophyta;Mediophyceae=Ochrophyta;Chrysophyceae=Ochrophyta;Synurophyceae=Ochrophyta;Raphidophyceae=Ochrophyta;Eustigmatophyceae=Ochrophyta;Phaeosacciophyceae=Ochrophyta;Pelagophyceae=Ochrophyta;Dictyochophyceae=Ochrophyta;Xanthophyceae=Ochrophyta;Bolidophyceae=Ochrophyta;Pinguiophyceae=Ochrophyta;Phaeothamniophyceae=Ochrophyta;Aurearenophyceae=Ochrophyta;Chrysomerophyceae=Ochrophyta;Schizocladiophyceae=Ochrophyta;Olisthodiscophyceae=Ochrophyta;Dinophyceae=Dinoflagellata). Details in metadata/taxonomy_harmonisation.tsv; species names and NCBI taxids unchanged.

## Sequences per step

| Step | Sequences |
|---|---|
| import_bold | 77807 |
| import_ncbi | 198536 |
| import_silva | 458470 |
| 01_merged | 722505 |
| 01b_lineage_rescued | 722505 |
| 01c_taxa_excluded_early | 346661 |
| 02_dereplicated_initial | 280247 |
| 03_filtered_initial | 209930 |
| 04_insilico_pcr | 81138 |
| 05_global_alignment | 86292 |
| 05b_long_trimmed | 86292 |
| 06_dereplicated_final | 70214 |
| 07_filtered_final | 62085 |
| 07b_length_rescued | 62085 |
| 07c_translation_checked | 62085 |
| 07d_taxonomy_fixed | 61720 |
| 08_taxa_excluded | 61668 |
| 09_min_rank | 59079 |
| 18S_V4V5_marine_eukaryotes | 59079 |

## Software

| Software | Version |
|---|---|
| script | 1.5.1 (2026-10-28, shared CRABS/taxo and CRABS/reference_database_creator; 1.5: from crabs_coi.sh 1.8 + taxonomy fix of 1.7.2, SILVA lineage check) |
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
- Quast et al. (2013) SILVA. Nucleic Acids Res. doi:10.1093/nar/gks1219
- Machida & Knowlton (2012) PCR primers for metazoan nuclear 18S and 28S ribosomal DNA sequences. PLoS ONE 7: e46180. doi:10.1371/journal.pone.0046180
- Martin (2011) cutadapt. EMBnet.journal. doi:10.14806/ej.17.1.200
- Rognes et al. (2016) VSEARCH. PeerJ. doi:10.7717/peerj.2584

Source details (download date, number of records, exact query) are in `metadata/sources.tsv`.
