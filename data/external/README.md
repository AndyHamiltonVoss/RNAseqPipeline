# External comparison data

`Rayl2024_MolPharmacol_supp_datafile2_24h.xlsx` is Supplementary Data File 2 (24-hour
timepoint) from:

Rayl M, Nemetchek MD, Voss AH, Hughes TS. *Agonists of the nuclear receptor PPARγ can
produce biased signaling.* Molecular Pharmacology (manuscript MOLPHARM-AR-2024-000992).

Sheets used by `Compare_to_Published.R`:

- `Normalized count table` — DESeq2 size-factor-normalized counts for DMSO, GW1929,
  Rosiglitazone, and MRL24 (8 replicates each), reprocessed through this pipeline's
  DESeq2 module as a reproducibility check.
- `Rosi_vs_DMSO`, `MRL24_vs_DMSO` — the paper's own reported DESeq2 results, used both
  for the reproducibility check and for comparison against this study's dose-response
  gene sets.
