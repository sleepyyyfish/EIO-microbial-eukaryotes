# Core metrics calculation

## Purpose

Calculate core alpha and beta diversity metrics using the rarefied ASV table in QIIME2.

## Input

- table.qza
- metadata_final.tsv

## Command

```bash
qiime diversity core-metrics \
  --i-table table.qza \
  --p-sampling-depth 45483 \
  --m-metadata-file metadata_final.tsv \
  --output-dir core-metrics-results
```

## Output

```
core-metrics-results/
├── rarefied_table.qza
├── observed_features_vector.qza
├── shannon_vector.qza
├── evenness_vector.qza
├── jaccard_distance_matrix.qza
├── bray_curtis_distance_matrix.qza
├── jaccard_pcoa_results.qza
├── bray_curtis_pcoa_results.qza
├── jaccard_emperor.qzv
└── bray_curtis_emperor.qzv
```

## Notes

- Sampling depth: **45483 reads per sample**.
- The generated alpha diversity metrics include **Observed Features**, **Shannon**, and **Pielou's Evenness**.
- The generated beta diversity metrics include **Jaccard** and **Bray–Curtis** distance matrices, corresponding PCoA results, and Emperor visualizations.
- **For the final analysis, alpha diversity indices (including Chao1, ACE, Sobs, Simpson, Good's coverage, etc.) were recalculated using an independent R script (`alpha_diversity.R`) to match the laboratory analysis workflow.**