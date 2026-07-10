# RNAseqPipeline
A modular RNA-seq differential expression analysis pipeline using DESeq2, edgeR, and limma.

## Environment

`environment.yml` pins the exact R/Bioconductor package versions this pipeline was run
with. Recreate it with:

```
conda env create -f environment.yml
conda activate dge-methods
```

## Pipeline

Run `DGE_pipeline.R` after editing `config.R`. See `config.R` for all available options
(QC, DESeq2/edgeR/limma, GO/KEGG enrichment, dose-response and Venn diagram comparisons).

## External validation

`Compare_to_Published.R` validates this study's results against the published RNA-seq
dataset from Rayl et al. 2024 (see `data/external/README.md` for citation and provenance).
It reprocesses the paper's normalized count table through this pipeline's own DESeq2
module to confirm the results are reproducible, then compares this study's dose-response
gene sets against the paper's reported significant genes. Requires `DGE_pipeline.R` to
have been run first. Small summary tables are written to `analysis/external_validation/`;
figures are written to `Output/Figures/` (not version-controlled).
