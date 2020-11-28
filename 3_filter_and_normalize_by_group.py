import json
import numpy as np
from util import convert_gene

grouped_files = ["F_KO_grouped", "F_WT_grouped", "M_KO_grouped", "M_WT_grouped"]
comparisons = {"F_WT_KO": [1, 0], "KO_FM": [0, 2], "M_WT_KO": [3, 2], "WT_FM": [1, 3]}


for group, comp_list in comparisons.items():
    print(f"Working on group {group}")
    file1_idx = comparisons[group][0]
    file2_idx = comparisons[group][1]
    with open(f"{group}_filtered_keep2.sf", "w") as outfile:
        with open(f"{grouped_files[file1_idx]}.sf") as sample1_fp:
            with open(f"{grouped_files[file2_idx]}.sf") as sample2_fp:
                samp1_lines = sample1_fp.readlines()
                samp2_lines = sample2_fp.readlines()
                if len(samp1_lines) != len(samp2_lines):
                    outfile.write(f"Files have different number of lines: {group}\n")
                    continue
                dropped_count = 0
                for idx in range(len(samp1_lines)):
                    samp1_line = samp1_lines[idx].strip()
                    samp2_line = samp2_lines[idx].strip()
                    samp1_fields = samp1_line.split("\t")
                    samp2_fields = samp2_line.split("\t")
                    ensembl1 = samp1_fields[0]
                    ensembl2 = samp2_fields[0]
                    if ensembl1 != ensembl2:
                        break_str = f"Lines between files are mismatched for group {group}, ensembl1: {ensembl1}, ensembl2: {ensembl2}"
                        outfile.write(break_str)
                        print(break_str)
                        break
                    combined_samps = samp1_fields[1:] + samp2_fields[1:]
                    min10 = 0
                    keep_gene = False
                    for samp in combined_samps:
                        if float(samp) >= 10.0:
                            min10 += 1
                            if min10 >= 2:
                                keep_gene = True
                                break
                    if keep_gene:
                        outfile.write(f"{ensembl1}")
                        for samp in combined_samps:
                            outfile.write(f"\t{samp}")
                        outfile.write("\n")
                    else:
                        dropped_count += 1
                print(f"Dropped {dropped_count} transcripts from group {group}")
print(f"Done filtering!")

for group, comp_list in comparisons.items():
    print(f"Working on filtered group {group}")
    values = {}
    for i in range(10):
        values[f"{i}"] = []
    # print(f"Empty values dict: {values}")

    with open(f"{group}_filtered_normalized.sf", "w") as outfile:
        with open(f"{group}_filtered_keep2.sf", "r") as infile:
            lines = infile.readlines()
            for line in lines:
                line = line.strip()
                vals = line.split("\t")[1:]
                for i in range(10):
                    values[str(i)].append(float(vals[i]))

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
                for i in range(10):
                    normalized_counts.append(counts[i] / percentiles[i])
                outfile.write(f"{ensembl}")
                for norm_count in normalized_counts:
                    outfile.write(f"\t{norm_count}")
                outfile.write("\n")
print("Done normalizing!")
