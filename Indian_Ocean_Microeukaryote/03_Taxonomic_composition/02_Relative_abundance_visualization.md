# Relative abundance visualization in R
## Input
- feature_table.tsv
- taxonomy.tsv
- metadata.tsv

## The R workflow:
- merge feature table with taxonomy table
- collapse ASVs to the selected taxonomic rank (Phylum / Class)
- calculate relative abundance for sorting
- merge with sample metadata
- generate stacked bar plots

## Output
- phylum_top11_relative.pdf
- class_top15_relative.pdf

The number of displayed taxa
 (e.g. Top 11 phyla or Top 15 classes), 
 color palette, and taxonomic grouping 
  (e.g. green for Autotrophs, yellow for Mixotrophs, red for Heterotrophs, blue for parasites)
  can be adjusted depending on the specific dataset and visualization requirements.