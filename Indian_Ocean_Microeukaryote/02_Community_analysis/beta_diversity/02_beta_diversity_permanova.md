# Beta diversity analysis

## Purpose

Evaluate differences in community composition among treatments and sampling days based on Bray–Curtis dissimilarity.

## Input

- bray_curtis_distance_matrix.qza
- metadata_final.tsv

## PERMANOVA

### Treatment effect

```bash
qiime diversity beta-group-significance \
  --i-distance-matrix bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata_final.tsv \
  --m-metadata-column treatment \
  --p-method permanova \
  --p-pairwise \
  --o-visualization permanova_treatment.qzv
```

### Day effect

```bash
qiime diversity beta-group-significance \
  --i-distance-matrix bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata_final.tsv \
  --m-metadata-column day \
  --p-method permanova \
  --p-pairwise \
  --o-visualization permanova_day.qzv
```

### Treatment comparisons within each sampling day

The same command was repeated for each sampling day using day-specific metadata files.

Example (Day 20):

```bash
qiime diversity beta-group-significance \
  --i-distance-matrix bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata_d20.tsv \
  --m-metadata-column treatment \
  --p-method permanova \
  --p-pairwise \
  --o-visualization permanova_d20.qzv
```

## Output

- permanova_treatment.qzv
- permanova_day.qzv
- permanova_d*.qzv

## Notes

PERMANOVA p-values were extracted from the `.qzv` visualization (viewed in QIIME 2 View) and used to annotate significance in downstream R-generated heatmaps.