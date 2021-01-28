## Snakefile for ASCAT2 files from GDC
##
## Tissue type just like in GDC, lowercase is fine
TISSUE = "lung"
## Only tumor samples for the ascat pipeline. It takes normal and tumor for 
## comparison, so there are no normal samples.
TTYPE = "tumor"
## The RNA files are added so that we get cases that include
## both files
DATADIR = "/datos/ot/diana/cnvs/"+TISSUE
FIGDIR = "figures/"+TISSUE
MDIR = DATADIR+"/manifests"
RAWDIR = DATADIR+"/raw"
biomart = "data/Biomart_Ensembl102_GRCh38_p13.txt"

rule get_heatmap:
  input: 
    DATADIR+"/"+TISSUE+"-"+TTYPE+"-ascat-matrix.tsv"
  output:
    FIGDIR+"/"+TISSUE+"-ascat-heatmap.png"
  shell:
    """
    mkdir -p {FIGDIR}
    Rscript src/getHeatmap.R {TISSUE} {biomart} {input} {FIGDIR}
    """

## We need to run these two together because the output of the download_files
## tasks depends on the manifest and there is no easy way to specify this on 
## snakemake
rule download_files_and_get_ascat_matrix:
  input:
    ## Manifest file
    MDIR+"/"+TISSUE+"-"+TTYPE+"-ascat.txt",
  output: 
    expand(DATADIR+"/"+TISSUE+"-"+TTYPE+"-ascat-{ft}.tsv", ft=["files", "matrix"])
  shell:
    """
    mkdir -p {RAWDIR}/{TISSUE}-{TTYPE}-ascat
    ./bin/gdc-client download -d {RAWDIR}/{TISSUE}-{TTYPE}-ascat -m {input} --retry-amount 3
    Rscript src/getMatrix.R {TISSUE} {TTYPE} {DATADIR} {RAWDIR}
    """

rule get_manifest:
  output:
    ## Example: data/manifests/breast-tumor-ascat.txt"
    MDIR+"/"+TISSUE+"-"+TTYPE+"-files.tsv",
    expand(MDIR+"/"+TISSUE+"-"+TTYPE+"-{fts}.txt", fts=["ascat", "rna_counts"] )
  shell:
    "mkdir -p {MDIR} ; python src/queryGDC.py {TISSUE} {TTYPE} {MDIR} "

rule clean:
  shell:
    "rm -rf {DATADIR} {FIGDIR}"
