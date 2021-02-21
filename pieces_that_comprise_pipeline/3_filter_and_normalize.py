import json
import numpy as np
from util import convert_gene
from config import study_dir, results_dir


with open(f"{results_dir}/all_grouped_filtered.sf", "w") as outfile:
    with open(f"{study_dir}/all_grouped.sf") as infile:
        lines = infile.readlines()
        dropped_count = 0
        for line in lines:
            line = line.strip()
            fields = line.split("\t")
            ensembl = fields[0]
            counts = fields[1:]
            min10 = 0
            keep_gene = False
            for count in counts:
                if float(count) >= 10.0:
                    min10 += 1
                    if min10 >= 3:
                        keep_gene = True
                        break
            if keep_gene:
                outfile.write(f"{ensembl}")
                for count in counts:
                    outfile.write(f"\t{count}")
                outfile.write("\n")
            else:
                dropped_count += 1
print(f"Done filtering! Dropped {dropped_count} transcripts")


print("Normalizing values")
samples = {}
for i in range(20):
    samples[f"{i}"] = []
# print(f"Empty samples dict: {samples}")

with open(f"{results_dir}/all_grouped_filtered_normalized.sf", "w") as outfile:
    with open(f"{results_dir}/all_grouped_filtered.sf", "r") as infile:
        lines = infile.readlines()
        for line in lines:
            line = line.strip()
            vals = line.split("\t")[1:]
            for i in range(20):
                samples[f"{i}"].append(float(vals[i]))

        for k, v in samples.items():
            print(f"sample {int(k)+1} has {len(v)} gene counts")

        percentiles = []
        for k, v in samples.items():
            percentile = np.percentile(v, 75)
            print(f"Calculated 75th percentile for sample {int(k)+1}: {percentile}")
            percentiles.append(percentile)

        for line in lines:
            line = line.strip()
            fields = line.split("\t")
            ensembl = fields[0]
            counts = [float(count) for count in fields[1:]]
            counts = np.array(counts)
            normalized_counts = counts / percentiles
            outfile.write(f"{ensembl}")
            for norm_count in normalized_counts:
                outfile.write(f"\t{norm_count}")
            outfile.write("\n")

print("QC step. Recalculating 75th percentile of normalized counts.")
with open(f"{results_dir}/all_grouped_filtered_normalized.sf", "r") as infile:
    samples = {}
    for i in range(20):
        samples[f"{i}"] = []
    normalized_lines = infile.readlines()
    for line in normalized_lines:
        line = line.strip()
        vals = line.split("\t")[1:]
        for i in range(20):
            samples[f"{i}"].append(float(vals[i]))

    for k, v in samples.items():
        print(f"sample {int(k)+1} has {len(v)} gene counts")

    percentiles = []
    for k, v in samples.items():
        percentile = np.percentile(v, 75)
        print(f"Recalculated 75th percentile for sample {int(k)+1}: {percentile}")
        percentiles.append(percentile)


print("Done normalizing!")
