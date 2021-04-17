## Snakefile for ASCAT2 files from GDC
##
## Tissue type just like in GDC, lowercase is fine
TISSUES = ["ovary", "prostate", "pancreas", "bladder", "skin", "brain", "testis", "liver", "esophagus", "breast", "lung", "kidney", "colorectal", "uterus", "thyroid"]
## Only tumor samples for the ascat pipeline. It takes normal and tumor for 
## comparison, so there are no normal samples.
TTYPE = "tumor"
## The RNA files are added so that we get cases that include
## both files
DATADIR = "/datos/ot/diana/cnvs"
FIGDIR = "figures"
biomart = "input/Biomart_Ensembl102_GRCh38_p13.txt"
heatmaps = [] 
for t in TISSUES:
  heatmaps.append(FIGDIR+ "/" + t + "/" + t + "-tumor-ascat-heatmap.png")

rule all:
  input:
    heatmaps

rule get_heatmap:
  input: 
    DATADIR+"/{tissue}/{tissue}-{type}-ascat-matrix.tsv"
  output:
    FIGDIR+"/{tissue}/{tissue}-{type}-ascat-heatmap.png"
  shell:
    """
    mkdir -p {FIGDIR}/{wildcards.tissue}
    Rscript src/getHeatmap.R {wildcards.tissue} {biomart} {input} {FIGDIR}/{wildcards.tissue}
    """

## We need to run these two together because the output of the download_files
## tasks depends on the manifest and there is no easy way to specify this on 
## snakemake
rule download_files_and_get_ascat_matrix:
  input:
    ## Manifest file
    DATADIR+"/{tissue}/manifests/{tissue}-{type}-ascat.txt",
  output: 
    DATADIR+"/{tissue}/{tissue}-{type}-ascat-files.tsv",
    DATADIR+"/{tissue}/{tissue}-{type}-ascat-matrix.tsv"
  shell:
    """
    mkdir -p {DATADIR}/{wildcards.tissue}/raw/{wildcards.tissue}-{wildcards.type}-ascat
    ./bin/gdc-client download -d {DATADIR}/{wildcards.tissue}/raw/{wildcards.tissue}-{wildcards.type}-ascat -m {input} --retry-amount 3
    Rscript src/getMatrix.R {wildcards.tissue} {wildcards.type} {DATADIR}
    """

rule get_manifest:
  output:
    ## Example: data/breast/manifests/breast-tumor-ascat.txt"
    DATADIR+"/{tissue}/manifests/{tissue}-{type}-files.tsv",
    DATADIR+"/{tissue}/manifests/{tissue}-{type}-ascat.txt",
    DATADIR+"/{tissue}/manifests/{tissue}-{type}-rna_counts.txt"
  shell:
    """
    mkdir -p {DATADIR}/{wildcards.tissue}/manifests
    python src/queryGDC.py {wildcards.tissue} {wildcards.type} {DATADIR}/{wildcards.tissue}/manifests 
    """
