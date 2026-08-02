# SIMPER analysis

## Description
Similarity Percentage (SIMPER) analysis was performed using the simper() function implemented in the vegan R package to identify ASVs contributing most to the dissimilarity between experimental groups.

## Input files:

- feature_table_RA.tsv
- taxonomy.tsv
- metadata_final.tsv

## Workflow:

Relative abundance feature table
            +
Sample metadata
            +
Taxonomic annotation
            │
            ▼
       vegan::simper()
            │
            ▼
Contribution of each ASV
            │
            ▼
Taxonomic interpretation

## R workflow

The analysis was performed using Bray-Curtis dissimilarity based on relative abundance data.

Main steps:

Select pairwise comparisons between experimental treatments.
Perform SIMPER analysis using ASV-level abundance matrix.
Calculate the contribution proportion of each ASV:
Contribution proportion=
total contribution
average contribution of ASV

Select ASVs with high contribution proportions.
Merge high-contribution ASVs with taxonomy information.
Calculate mean relative abundance of selected ASVs in each treatment.
Generate contribution-abundance plots.

## Example

library(vegan)

sim <- simper(
  otu_table,
  group
)

summary(sim)

The final script additionally:

recalculated contribution proportions;
ranked ASVs according to their contribution;
linked ASVs with taxonomic information;
exported complete SIMPER results and selected high-contribution ASVs.

## Output:

SIMPER_*_full_results.csv
SIMPER_*_top_OTUs.csv
SIMPER_*_abundances.csv
SIMPER_summary_all_comparisons.csv