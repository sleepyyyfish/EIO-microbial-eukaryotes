#Workflow
Raw FASTQ
   │
   ▼
QIIME2
   │
   ├── feature-table.qza
   ├── rep-seqs.qza
   ├── taxonomy.qza
   │
   ▼
Export
   │
feature-table.tsv
metadata.tsv
taxonomy.tsv
   │
   ├─────────────────────────────────────────┐
   ▼                                         ▼
Community                                    Network
analysis                                     analysis
   ├─────────────────────┐                   ├───────────────┐───────────────┐  
   │                     │                   |               |               |
   │                     │                   |               |               | 
   ▼                     ▼                   ▼               ▼               ▼
α-diversity      relative abundance    nodes & edges       Zi-Pi        Fragmentation 
(boxplot)             (barplot)              |        degree-betweenness  Robustness
β-diversity          SIMilarity              ▼         niche breadth    Vulnerability
(NMDS figure         PERcentage            Gephi
dissimilarity heatmap)                    network

