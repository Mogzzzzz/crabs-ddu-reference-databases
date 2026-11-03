# 18S rRNA V4 forward-read reference database (150 bp after 18Smk1F) – marine invertebrates and other eukaryotes (18S_V4_fwd150_marine_eukaryotes)

- **Created**: 2026-10-28 (run 2026-10-28T17:19:28+0100 → 2026-10-28T17:19:34+0100)
- **Contact**: morgane.durand@mnhn.fr
- **Content**: 67746 sequences, 41959 taxa; length min / median / max: 145 / 150 / 155 bp; 6227 sequences shared by several taxa
- **Primer**: CTGGTGCCAGCAGCCGCGGYAA (18Smk1F, Machida & Knowlton 2012); region = the 150 bp following it (forward read)
- **Built from**: the curated full sequences of the 18S database `/Users/morgane/PhD/BIOINFO/DDU/NEWAGE/DNA_analyses/CRABS/crabs_18S/18S_aligned`

## Use

For forward reads only: primer removed, reads cut to 150 bp (shorter reads discarded), e.g. `obigrep -l 150 reads.fasta | obiannotate --cut 1:150`. References and reads then have the same length.

## Files

| File | Content |
|---|---|
| `results/18S_V4_fwd150_marine_eukaryotes.fasta` | Final database, headers `taxid=…; scientific_name=…; lineage=…;` (OBITools) |
| `results/18S_V4_fwd150_marine_eukaryotes.tsv` | Same content as a table with column names |
| `results/18S_V4_fwd150_marine_eukaryotes.txt` | CRABS format |
| `results/18S_V4_fwd150_marine_eukaryotes_insilico_pcr_raw.fasta` | Regions found by the forward-only in silico PCR (before rescue and filters) |
| `metadata/18S_V4_fwd150_marine_eukaryotes_run_metadata.json` | All metadata, machine-readable |
| `metadata/18S_V4_fwd150_marine_eukaryotes_parameters.tsv`, `…versions.tsv`, `…sources.tsv`, `…counts.tsv` | Parameters, software, sources, counts per step |
| `metadata/18S_V4_fwd150_marine_eukaryotes_checksums_sha256.txt` | Checksums of the result files |
| `logs/18S_V4_fwd150_marine_eukaryotes_*` | Log of every step |

## Methods

1. Input: curated full sequences of the 18S database (SILVA SSU subset 138.2, NCBI GenBank, BOLD; merged, lineage rescued, excluded taxa removed, dereplicated, environmental and poorly identified sequences removed; see its own README).
2. Forward-only in silico PCR (cutadapt, 4 mismatches, no indels, both strands): the 150 bp following the forward primer; sequences ending earlier discarded.
3. Sequences without the primer site aligned (VSEARCH, identity ≥ 0.95) on the regions of step 2; kept only when the whole region is covered, so that its start is exact.
4. Dereplication (unique_species), filters: 145–155 bp, fewer than 2 N, environmental sequences and sequences not identified at least to genus removed; taxa removed: Bacteria;Archaea;Insecta;Collembola;Protura;Diplura;Arachnida;Chilopoda;Diplopoda;Pauropoda;Symphyla;Onychophora.
5. Phylum names harmonised with Adl et al. (2019): renamed Bacillariophyta=Ochrophyta;Heterokontophyta=Ochrophyta;Myzozoa=Dinoflagellata;Dinophyta=Dinoflagellata; phylum given from the class when NCBI has none. Details in metadata/18S_V4_fwd150_marine_eukaryotes_taxonomy_harmonisation.tsv.

## Sequences per step

| Step | Sequences |
|---|---|
| 03_filtered_initial | 209930 |
| 18S_V4_fwd150_marine_eukaryotes_04_insilico_pcr | 101717 |
| 18S_V4_fwd150_marine_eukaryotes_05_global_alignment | 137803 |
| 18S_V4_fwd150_marine_eukaryotes_06_dereplicated_final | 98467 |
| 18S_V4_fwd150_marine_eukaryotes_07_filtered_final | 74659 |
| 18S_V4_fwd150_marine_eukaryotes_07d_taxonomy_fixed | 71317 |
| 18S_V4_fwd150_marine_eukaryotes_08_taxa_excluded | 71281 |
| 18S_V4_fwd150_marine_eukaryotes_09_min_rank | 67746 |
| 18S_V4_fwd150_marine_eukaryotes | 67746 |

## Software

| Software | Version |
|---|---|
| script | 1.3.1 (2026-10-28, shared CRABS/taxo and CRABS/reference_database_creator; 1.3: from crabs_18S.sh 1.5: taxonomy fix with organelle rRNA removal; runs in the 18S folder) |
| CRABS | 1.14.0 (commit 97b3037) |
| cutadapt | 5.2 |
| vsearch | vsearch v2.32.0_macos_aarch64 |
| python | Python 3.13.4 |
| system | Darwin 25.5.0 arm64 |

## Please cite

- Jeunen et al. (2022) CRABS. Mol Ecol Resour. doi:10.1111/1755-0998.13741
- Ratnasingham & Hebert (2007) BOLD. Mol Ecol Notes. doi:10.1111/j.1471-8286.2007.01678.x
- Sayers et al. NCBI GenBank / Taxonomy (https://www.ncbi.nlm.nih.gov)
- Quast et al. (2013) SILVA. Nucleic Acids Res. doi:10.1093/nar/gks1219
- Adl et al. (2019) Revisions to the classification, nomenclature, and diversity of eukaryotes. J Eukaryot Microbiol 66: 4-119. doi:10.1111/jeu.12691
- Machida & Knowlton (2012) PCR primers for metazoan nuclear 18S and 28S ribosomal DNA sequences. PLoS ONE 7: e46180. doi:10.1371/journal.pone.0046180
- Martin (2011) cutadapt. EMBnet.journal. doi:10.14806/ej.17.1.200
- Rognes et al. (2016) VSEARCH. PeerJ. doi:10.7717/peerj.2584
