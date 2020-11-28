import json
import numpy as np
from util import convert_gene


with open("all_grouped_filtered.sf", "w") as outfile:
    with open("all_grouped.sf") as infile:
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
values = {}
for i in range(20):
    values[f"{i}"] = []
# print(f"Empty values dict: {values}")

with open("all_grouped_filtered_normalized.sf", "w") as outfile:
    with open("all_grouped_filtered.sf", "r") as infile:
        lines = infile.readlines()
        for line in lines:
            line = line.strip()
            vals = line.split("\t")[1:]
            for i in range(20):
                values[f"{i}"].append(float(vals[i]))

        for k, v in values.items():
            print(f"sample {int(k)+1} has {len(v)} gene counts")

        percentiles = []
        for k, v in values.items():
            percentile = np.percentile(v, 75)
            print(f"Calculated 75th percentile for sample {int(k)+1}: {percentile}")
            percentiles.append(percentile)

        for line in lines:
            line = line.strip()
            fields = line.split("\t")
            ensembl = fields[0]
            counts = [float(count) for count in fields[1:]]
            normalized_counts = []
            for i in range(20):
                normalized_counts.append(counts[i] / percentiles[i])
            outfile.write(f"{ensembl}")
            for norm_count in normalized_counts:
                outfile.write(f"\t{norm_count}")
            outfile.write("\n")
print("Done normalizing!")
