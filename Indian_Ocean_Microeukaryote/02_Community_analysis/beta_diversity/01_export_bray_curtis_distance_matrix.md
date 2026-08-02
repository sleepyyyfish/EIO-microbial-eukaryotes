# Export Bray–Curtis distance matrix

## Purpose

Export the Bray–Curtis distance matrix from QIIME2 format (`.qza`) to a tab-delimited text file for downstream analyses in R.

## Input

- bray_curtis_distance_matrix.qza

## Command

```bash
qiime tools export \
  --input-path bray_curtis_distance_matrix.qza \
  --output-path bray_dist
```

## Output

```
bray_dist/
└── distance-matrix.tsv
```

## Notes

The exported distance matrix was used for downstream analyses, including NMDS ordination and average dissimilarity calculations in R.