# rbcL reference database (~330 bp) – marine algae and seagrasses (rbcL_marine_algae)

- **Created**: 2026-10-28 (run 2026-10-28T17:17:11+0100 → 2026-10-28T17:17:28+0100)
- **Contact**: morgane.durand@mnhn.fr
- **Content**: 13737 sequences, 6428 taxa; length min / median / max: 310 / 331 / 346 bp
- **Primers**: ATGCGTTGGAGAGARCGTTTC (rbcL-646F, Kelly et al. 2018) / GATCACCTTCTAATTTACCWACAACTG (rbcL-998R, Kelly et al. 2018)

## Files

| File | Content |
|---|---|
| `results/rbcL_marine_algae.fasta` | Final database, headers `taxid=…; scientific_name=…; lineage=…;` (OBITools) |
| `results/rbcL_marine_algae.tsv` | Same content as a table with column names |
| `results/rbcL_marine_algae.txt` | CRABS format (usable by every CRABS function) |
| `results/insilico_pcr_raw.fasta` | Raw in silico PCR output (before alignment and filters) |
| `metadata/run_metadata.json` | All metadata, machine-readable |
| `metadata/parameters.tsv`, `versions.tsv`, `sources.tsv`, `counts.tsv` | Parameters, software, sources, counts per step |
| `metadata/checksums_sha256.txt` | Checksums of the result files |
| `logs/` | Log of every step |

## Methods

1. Download from BOLD (V5 portal, 7 groups, records with marker /rbcL/ only) and NCBI GenBank (34 groups), plus the NCBI taxonomy.
2. Import into CRABS format with NCBI lineage (domain;phylum;class;order;family;genus;species).
3. Merge (duplicate accessions removed); lineage rescue (a species name unknown to the NCBI taxonomy takes the lineage of its genus); removal of the excluded taxa; dereplication (unique_species), removal of environmental sequences and of sequences not identified at least to genus (or with more than 1 missing rank above it).
4. In silico PCR (cutadapt, 4 mismatches, relaxed=false).
5. Pairwise global alignment (VSEARCH, identity ≥ 0.95, coverage ≥ 95 %) to recover amplicons deposited without primer-binding sites.
6. Final dereplication and filters: 310–360 bp, fewer than 2 N; amplicons with a stop codon in every reading frame removed (yes (stop codons: TAA TAG TGA)); taxa removed:  (removed after import, after the final filters).
Genus homonyms: yes (genus homonyms: BOLD phylum, then target groups + Ochrophyta;Bacillariophyta;Phaeophyceae; log logs/10_taxonomy_fix.tsv).
7. Phylum names harmonised with Adl et al. (2019): renamed Bacillariophyta=Ochrophyta;Heterokontophyta=Ochrophyta;Myzozoa=Dinoflagellata;Dinophyta=Dinoflagellata; phylum given from the class when NCBI has none (Phaeophyceae=Ochrophyta;Bacillariophyceae=Ochrophyta;Coscinodiscophyceae=Ochrophyta;Mediophyceae=Ochrophyta;Chrysophyceae=Ochrophyta;Synurophyceae=Ochrophyta;Raphidophyceae=Ochrophyta;Eustigmatophyceae=Ochrophyta;Phaeosacciophyceae=Ochrophyta;Pelagophyceae=Ochrophyta;Dictyochophyceae=Ochrophyta;Xanthophyceae=Ochrophyta;Bolidophyceae=Ochrophyta;Pinguiophyceae=Ochrophyta;Phaeothamniophyceae=Ochrophyta;Aurearenophyceae=Ochrophyta;Chrysomerophyceae=Ochrophyta;Schizocladiophyceae=Ochrophyta;Olisthodiscophyceae=Ochrophyta;Dinophyceae=Dinoflagellata). Details in metadata/taxonomy_harmonisation.tsv; species names and NCBI taxids unchanged.

## Sequences per step

| Step | Sequences |
|---|---|
| import_bold | 18348 |
| import_ncbi | 61996 |
| 01_merged | 80344 |
| 01b_lineage_rescued | 80344 |
| 01c_taxa_excluded_early | 80344 |
| 02_dereplicated_initial | 58804 |
| 03_filtered_initial | 54133 |
| 04_insilico_pcr | 25305 |
| 05_global_alignment | 26848 |
| 05b_long_trimmed | 26848 |
| 06_dereplicated_final | 14647 |
| 07_filtered_final | 14107 |
| 07b_length_rescued | 14107 |
| 07c_translation_checked | 14042 |
| 07d_taxonomy_fixed | 14042 |
| 08_taxa_excluded | 14042 |
| 09_min_rank | 13737 |
| rbcL_marine_algae | 13737 |

## Software

| Software | Version |
|---|---|
| script | 1.4.1 (2026-10-28, shared CRABS/taxo and CRABS/reference_database_creator; 1.4: from crabs_coi.sh 1.8 + taxonomy fix of 1.7.2 + long-amplicon re-cut of 1.7.3) |
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
- Kelly M, Boonham N, Juggins S, Kille P, Mann D, Pass D, Sapp M, Sato S, Glover R (2018) A DNA based diatom metabarcoding approach for Water Framework Directive classification of rivers. Environment Agency report SC140024/R (rbcL-646F / rbcL-998R)
- Martin (2011) cutadapt. EMBnet.journal. doi:10.14806/ej.17.1.200
- Rognes et al. (2016) VSEARCH. PeerJ. doi:10.7717/peerj.2584

Source details (download date, number of records, exact query) are in `metadata/sources.tsv`.
