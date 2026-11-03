# Figures – 18S_V4_fwd150_marine_eukaryotes

- **Created**: 2026-10-28T17:29:52+0100
- **Script**: crabs_figures.R 2.7 (2026-10-28)
- **Software**: R version 4.4.2 (2024-10-31); data.table 1.17.8, ggplot2 4.0.0, scales 1.4.0, viridisLite 0.4.2
- **Formats**: PNG (300 dpi) for slides and drafts, PDF (vector) for publication; the data behind each figure is in `data/`.

## Figures

| Figure | Caption | Data |
|---|---|---|
| `fig01_pipeline_yield.png` / `.pdf` | Number of sequences after each step of the pipeline, from import to the final database (parameters in brackets; steps that changed nothing are not shown). | `data/fig01_pipeline_yield.csv` |
| `fig02_composition_phylum.png` / `.pdf` | Number of species and of sequences per phylum, animals first, then algae (top 12 phyla; the others grouped as 'Other animals', 'Other algae'...). | `data/fig02_composition_phylum.csv` |
| `fig03_source_by_phylum.png` / `.pdf` | Number of retained sequences per phylum, by the database they come from (exact match of identifiers). | `data/fig03_source_by_phylum.csv` |
| `fig04_amplicon_lengths.png` / `.pdf` | Length distribution of the amplicons extracted by in silico PCR, before the final length filter; shaded = window kept. | `data/fig04_amplicon_lengths.csv` |
| `fig05_lengths_by_phylum.png` / `.pdf` | Length distribution of the final amplicons for each of the most represented phyla. | `data/fig05_lengths_by_phylum.csv` |
| `fig06_sequences_per_species.png` / `.pdf` | Number of species by number of unique reference sequences (intraspecific variants kept by the dereplication). | `data/fig06_sequences_per_species.csv` |
| `fig07_lineage_completeness.png` / `.pdf` | Share of sequences in the final database with a name at each taxonomic rank of the NCBI lineage. | `data/fig07_lineage_completeness.csv` |
| `fig08_downloads_by_group.png` / `.pdf` | Raw records downloaded from BOLD and NCBI GenBank, per phylum (same phyla and order as figure 2); the groups queried in each database and their phylum are in data/group_taxonomy.tsv. | `data/fig08_downloads_by_group.csv` |
| `fig09_primer_sites_composition.png` / `.pdf` | Base composition of the primer-binding sites in the full deposited sequences (CRABS 'amplification efficiency'); percentages = share of sequences mismatching the primer at that position. | `data/fig09_primer_sites_composition.csv` |
| `fig10_primer_mismatch_by_phylum.png` / `.pdf` | Mismatch rate between the primer and the binding site of each sequence, per position and phylum. | `data/fig10_primer_mismatch_by_phylum.csv` |
| `fig12_species_resolution.png` / `.pdf` | Share of species whose amplicons are unique to them, partly shared or entirely shared with other species (exact identity). | `data/fig12_species_resolution.csv` |
| `fig13_focal_taxa.png` / `.pdf` | Species and sequences of the focal taxa (Porifera, Cnidaria, Bryozoa, Annelida, Nematoda, Platyhelminthes, Rhodophyta, Ochrophyta) in the final database, by source. Ochrophyta includes diatoms and brown algae. | `data/fig13_focal_taxa.csv` |

## Inputs (MD5)

| File | MD5 |
|---|---|
| `results/18S_V4_fwd150_marine_eukaryotes.tsv` | `e21a1bbf298262634727ff0038552af6` |
| `metadata/18S_V4_fwd150_marine_eukaryotes_counts.tsv` | `bb0446ae833ea359931be79f792e28ea` |
| `metadata/18S_V4_fwd150_marine_eukaryotes_sources.tsv` | `ad4484fa819f53593b90cbeb0d727c1f` |
| `metadata/18S_V4_fwd150_marine_eukaryotes_parameters.tsv` | `d896d4b2b120ea433a3d48861961fc68` |
| `steps/18S_V4_fwd150_marine_eukaryotes_04_insilico_pcr.txt` | `ed2e3ec6cdefb3e5abfd0d7e29140a29` |
| `../18S_aligned/steps/03_filtered_initial.txt` | `7f4d6e8d96bb087228a298cf92df22ae` |
| `../18S_aligned/imported/bold_all.txt` | `a60dde4d97996ae2c1737757e66b2fbb` |
| `../18S_aligned/imported/ncbi_all.txt` | `0b0ce2aaa209a19b6126facab7612f32` |
| `../18S_aligned/imported/silva_all.txt` | `45d0ba64b62ea896d25ce9e89c24f6a2` |

Full R session details: `sessionInfo.txt`.
