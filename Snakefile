import pandas as pd

tissue = "ovary"
ttype = "cancer"
manifest = "data/" + tissue + "_" + ttype + "_manifest.txt"

mani = pd.read_csv(manifest, sep="\t")
file_dirs = mani.apply(lambda x: "data/raw/" + x["id"] + "/" + x["filename"],
        axis=1) 

print(file_dirs[0])
rule download_files:
    input:
        manifest
    output:
        file_dirs
    shell:
        "./bin/gdc-client download -d data/raw -m {manifest} --retry-amount 3"

rule clean:
    shell:
        "rm -rf data/raw/*"
