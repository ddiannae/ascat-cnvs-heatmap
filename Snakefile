import pandas as pd
 
# Tissue type just like in GDC, lowercase is fine
TISSUE = "breast"
TTYPES = ["normal", "tumor"]
FTYPES = ["ascat", "rna_counts"]
MDIR = "data/manifests/"
#  manifests_files.append("data/" + "-".join([tissue, tt, "files"]) + ".tsv")

#def get_raw_filenames(manifest):
#  mani = pd.read_csv(manifest, sep="\t")
#  rawfolder = manifest.split("/")[2].replace(".txt", "")
#  file_dirs = mani.apply(lambda x: "data/raw/" + rawfolder + "/" + x["id"] + "/" + x["filename"],
#        axis=1) 
#  return(file_dirs.tolist())

rule all:
  input:
    expand("data/{tissue}-{tt}-ascat-{ft}.tsv", tissue=TISSUE, tt=TTYPES[0], ft=["files", "matrix"])

rule get_ascat_matrix:
  input: directory("data/raw/"+TISSUE+"-"+TTYPES[0]+"-ascat") 
  output:
    expand("data/{tissue}-{tt}-ascat-{ft}.tsv", tissue=TISSUE, tt=TTYPES[0], ft=["files", "matrix"])
  shell:
    "Rscript src/getMatrix.R {TISSUE} {TTYPES[0]}"

rule download_files:
  input: expand(MDIR + "{tissue}-{tt}-ascat.txt", tissue=TISSUE, tt=TTYPES[0])
  output: directory("data/raw/"+TISSUE+"-"+TTYPES[0]+"-ascat") 
  shell:
    "mkdir -p data/raw/{TISSUE}-{TTYPES[0]}-ascat ; ./bin/gdc-client download -d data/raw/{TISSUE}-{TTYPES[0]}-ascat -m {input} --retry-amount 3"

rule get_manifest:
  output:
    expand(MDIR+"{tissue}-{tts}-{fts}.txt", tissue=TISSUE, tts=TTYPES, fts=FTYPES)
  run:
    for tt in TTYPES:
      shell("mkdir -p {MDIR} ; python src/queryGDC.py {TISSUE} {tt}")

rule clean:
  shell:
    "rm -rf data"
