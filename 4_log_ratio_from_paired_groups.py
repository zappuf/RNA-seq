import json
import math
import numpy as np
import statistics as stat
from util import convert_gene


F_KO = "F_KO_grouped.sf"
F_WT = "F_WT_grouped.sf"
M_KO = "M_KO_grouped.sf"
M_WT = "M_WT_grouped.sf"

F_WT_KO = "F_WT_KO_filtered_normalized.sf"
KO_FM = "KO_FM_filtered_normalized.sf"
M_WT_KO = "M_WT_KO_filtered_normalized.sf"
WT_FM = "WT_FM_filtered_normalized.sf"


def get_log2_list(list1):
    out_list = []
    for line in list1:
        line = line.strip()
        ensembl = line.split("\t")[0]
        gene = convert_gene(ensembl)

        group1_fields = line.split("\t")[1:6]
        group2_fields = line.split("\t")[6:]

        group1_floats = [float(i) for i in group1_fields]
        group2_floats = [float(i) for i in group2_fields]

        group1_avg = sum(group1_floats) / len(group1_floats)
        group2_avg = sum(group2_floats) / len(group2_floats)

        log2ratio = math.log2(((group2_avg + 1) / (group1_avg + 1)))

        out_str = f"{gene}\t{round(log2ratio, 2)}\n"
        out_list.append(out_str)
    out_list.sort()
    return out_list


def combine_same_gene(log2_list):
    dup_mapping = {}
    last = False
    last_idx = len(log2_list) - 1
    for idx in range(len(log2_list)):
        # print(f"{idx+1} of {last_idx+1}")
        if idx == last_idx:
            last = True
            break
        current_line = log2_list[idx].strip()
        next_line = log2_list[idx + 1].strip()
        current_fields = current_line.split("\t")
        if len(current_fields) < 2:
            continue
        next_fields = next_line.split("\t")

        try:
            current_gene = current_fields[0]
            if current_gene not in dup_mapping.keys() and current_gene != "":
                dup_mapping[current_gene] = [float(current_fields[1])]
            next_gene = next_fields[0]
        except Exception as e:
            print(f"Error, can't assign gene: {e}")
            continue

        if current_gene == next_gene and current_gene != "":
            dup_mapping[current_gene].append(float(next_fields[1]))
    with open("duplicates.json", "w") as fp:
        fp.write(json.dumps(dup_mapping))
    gene_sums = {}
    for k, v in dup_mapping.items():
        gene_sums[k] = sum(v)
    with open("summed.json", "w") as fp:
        fp.write(json.dumps(gene_sums))
    uniq_list = []
    for k, v in gene_sums.items():
        this_line = f"{k}\t{round(v, 2)}"
        uniq_list.append(this_line)
    return uniq_list


##############################################################
# log2 ratio for F_WT & F_KO
with open("F_WT_KO_log2.txt", "w") as outfile:
    with open(F_WT_KO, "r") as infile:
        lines_list = list(infile.readlines())
        out_list = get_log2_list(lines_list)
        uniq_list = combine_same_gene(out_list)
        for elem in uniq_list:
            outfile.write(f"{elem}\n")


# log2 ratio for M_WT & M_KO
with open("M_WT_KO_log2.txt", "w") as outfile:
    with open(M_WT_KO, "r") as infile:
        lines_list = list(infile.readlines())
        out_list = get_log2_list(lines_list)
        uniq_list = combine_same_gene(out_list)
        for elem in uniq_list:
            outfile.write(f"{elem}\n")


# log2 ratio for M_WT & F_WT
with open("WT_FM_log2.txt", "w") as outfile:
    with open(WT_FM, "r") as infile:
        lines_list = list(infile.readlines())
        out_list = get_log2_list(lines_list)
        uniq_list = combine_same_gene(out_list)
        for elem in uniq_list:
            outfile.write(f"{elem}\n")


# log2 ratio for M_KO & F_KO
with open("KO_FM_log2.txt", "w") as outfile:
    with open(KO_FM, "r") as infile:
        lines_list = list(infile.readlines())
        out_list = get_log2_list(lines_list)
        uniq_list = combine_same_gene(out_list)
        for elem in uniq_list:
            outfile.write(f"{elem}\n")
