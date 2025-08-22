import pandas as pd
import warnings
from pathlib import Path 
import sys
warnings.filterwarnings("ignore")

gwas_filename = sys.argv[1]
gwas = pd.read_csv(gwas_filename)
num_phenotypes = gwas['phenotype'].nunique()

gwas = gwas[gwas["#chr"] != "intercept"]

sig_threshold = 0.05/(len(gwas)/num_phenotypes) # bonferonni correction

gwas = gwas[gwas["pvalue"] < sig_threshold]

gff_filename = sys.argv[2]
gff_window_size = int(sys.argv[3])
gff = pd.read_csv(gff_filename, sep="\t", header=None, comment='#')

gff.columns = ["#seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]

significant_annotations = []
for idx, row in gwas.iterrows():
    chr = row["#chr"]
    pos = row["pos"]

    gff_temp = gff[gff["#seqid"] == chr]
    gff_temp = gff[gff["source"] != "RefSeq"]
    mask1 = gff_temp["start"] <= (pos + gff_window_size)
    mask2 = gff_temp["end"] >= (pos - gff_window_size)
    mask = mask1 & mask2
    gff_temp = gff_temp[mask]
    gff_temp["phenotype"] = row["phenotype"]

    significant_annotations.append(gff_temp)

output_gff = pd.concat(significant_annotations, ignore_index=True)

output_path = Path(gff_filename)
output_path = output_path.with_name(output_path.stem + "_filtered" + output_path.suffix)
output_gff.to_csv(output_path, sep="\t", index=False, quoting=3)

print("File created: " + output_path.name, end="")
