# Step 1. Import raw sequences

## Purpose

Convert raw sequencing files into a QIIME2 artifact for downstream analyses.

## Input

- Raw paired-end FASTQ files (`*.fastq.gz`)
- Manifest file (`manifest.tsv`)

## Command

```bash
# Start the QIIME2 Docker container
docker run -it --rm \
-v /home/sleepyy/bio_thesis:/data \
quay.io/qiime2/core:2023.9 \
bash

# Import paired-end FASTQ files
qiime tools import \
...
```

## Output

- `demux.qza`
  - Demultiplexed sequence artifact in QIIME2 format.

## Notes
- Raw sequencing data were provided by the sequencing company.
- Each sample contains two compressed FASTQ files:
  - `*_1.fastq.gz` (forward reads)
  - `*_2.fastq.gz` (reverse reads)
- Copy all FASTQ files into the `00_raw_data/` directory before importing.
- (how to create new folders?
    mkdir xxx
    mkdir -p bio_thesis/00_raw_data
    cd = change directory)
- A `manifest.tsv` file is required to link sample IDs with the corresponding FASTQ files.
- Docker was used to run QIIME2 (version 2023.9).

# Step 2. Summarize sequencing quality

## Purpose

Visualize sequencing quality and determine trimming parameters for DADA2.
## Input

- `demux.qza`

## Command

```bash
qiime demux summarize \
...
```

## Output

- `demux.qzv`
  - Interactive quality report.

## Notes

- Open the `.qzv` file in **QIIME 2 view** to inspect sequencing quality. (https://view.qiime2.org)
- The quality profile is used to determine the truncation length for DADA2.

# Step 3. DADA2 denoising

## Purpose

Generate high-quality ASVs by filtering low-quality reads, merging paired-end reads, and removing chimeras.

## Input

- `demux.qza`

## Command

```bash
qiime dada2 denoise-paired \
...
```

## Output

- `table.qza`
  - Feature table (ASV abundance table)
- `rep-seqs.qza`
  - Representative ASV sequences
- `stats.qza`
  - DADA2 statistics

## Notes

- Trimming parameters were determined according to the quality profile in `demux.qzv`.
- The feature table (`table.qza`) is the main input for downstream diversity analyses.
- Representative sequences (`rep-seqs.qza`) are used for taxonomic classification.