# Copy number variation heatmaps from ASCAT data

Scripts to build heatmap figures from copy number variation data from ASCAT files in [Genomic Data Commons](https://docs.gdc.cancer.gov/). 

This exercise had the following objectives: 

- To explore the ASCAT files format.
- To practice the use of  [snakemake](https://snakemake.readthedocs.io/) for workflow implementation.
- To get insight into the variability of cnv data in different types of cancer for future analysis.

## Pipeline description

The pipeline uses snakemake for a reproducible workflow to run the following steps for ovary, prostate, pancreas, bladder, skin, brain, testis, liver, esophagus, breast, lung, kidney, colorectal, uterus, and thyroid tissues.

1. Query GDC to get manifest files for ASCAT and RNASeq data for normal and tumor conditions. File: `src/queryGDC.py`
2. Download data in manifest files using the [GDC Data Transfer Tool](https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Getting_Started/), also in `bin/gdc-client`
3. Build ASCAT matrix (genes in rows and samples in columns)
4. Get CNV heatmap figure (blue for 0-1 values, white for 2, red for +2)
