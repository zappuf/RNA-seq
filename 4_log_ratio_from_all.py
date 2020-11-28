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
        ensembl = line[0]
        gene = convert_gene(ensembl)

        group1_fields = line[1:6]
        group2_fields = line[6:]

        group1_floats = [float(i) for i in group1_fields]
        group2_floats = [float(i) for i in group2_fields]

        group1_avg = sum(group1_floats) / len(group1_floats)
        group2_avg = sum(group2_floats) / len(group2_floats)

        log2ratio = math.log2(((group2_avg + 1) / (group1_avg + 1)))

        out_str = f"{gene}\t{round(log2ratio, 2)}\t{group1_avg}\t{group2_avg}\n"
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
        if len(current_fields) < 4:
            continue
        next_fields = next_line.split("\t")

        try:
            current_gene = current_fields[0]
            if current_gene not in dup_mapping.keys() and current_gene != "":
                dup_mapping[current_gene] = {}
                dup_mapping[current_gene]["log2"] = [float(current_fields[1])]
                dup_mapping[current_gene]["denominator"] = [float(current_fields[2])]
                dup_mapping[current_gene]["numerator"] = [float(current_fields[3])]
            next_gene = next_fields[0]
        except Exception as e:
            print(f"Error, can't assign gene: {e}")
            continue

        if current_gene == next_gene and current_gene != "":
            dup_mapping[current_gene]["log2"].append(float(next_fields[1]))
            dup_mapping[current_gene]["denominator"].append(float(next_fields[2]))
            dup_mapping[current_gene]["numerator"].append(float(next_fields[3]))
    with open("duplicates.json", "w") as fp:
        fp.write(json.dumps(dup_mapping))
    gene_sums = {}
    uniq_list = []
    for k, v in dup_mapping.items():
        gene_sums[k] = {}
        for gene_k, gene_v in v.items():
            # gene_sums[k] = {
            #     "log2": sum(gene_sums[k]["log2"]),
            #     "denominator": sum(gene_sums[k]["denominator"]),
            #     "numerator": sum(gene_sums[k]["numerator"]),
            # }
            print(f"k: {k}, gene_k: {gene_k}, gene_v: {gene_v}")
            gene_sums[k][gene_k] = sum(gene_v)
    with open("summed.json", "w") as fp:
        fp.write(json.dumps(gene_sums))
    for k, v in gene_sums.items():
        gene = k
        log2 = gene_sums[k]["log2"]
        denominator = gene_sums[k]["denominator"]
        numerator = gene_sums[k]["numerator"]
        this_line = f"{gene}\t{round(log2, 2)}\t{denominator}\t{numerator}"
        # this_line = f"{k}\t{round(v, 2)}"
        uniq_list.append(this_line)
    return uniq_list


##############################################################
with open("all_grouped_filtered_normalized.sf", "r") as infile:
    lines = infile.readlines()
    lines_list = list(lines)
    F_WT_KO_list = []
    M_WT_KO_list = []
    KO_FM_list = []
    WT_FM_list = []
    for line in lines:
        line = line.strip()
        fields = line.split("\t")
        ensembl = [fields[0]]
        F_KO = fields[1:6]
        F_WT = fields[6:11]
        M_KO = fields[11:16]
        M_WT = fields[16:]

        F_WT_KO_line = ensembl + F_WT + F_KO
        F_WT_KO_list.append(F_WT_KO_line)

        M_WT_KO_line = ensembl + M_WT + M_KO
        M_WT_KO_list.append(M_WT_KO_line)

        KO_FM_line = ensembl + F_KO + M_KO
        KO_FM_list.append(KO_FM_line)

        WT_FM_line = ensembl + F_WT + M_WT
        WT_FM_list.append(WT_FM_line)

    # log2 ratio for F_WT & F_KO
    F_WT_KO_outlist = get_log2_list(F_WT_KO_list)
    F_WT_KO_uniq_list = combine_same_gene(F_WT_KO_outlist)
    with open("log2_F_WT_KO_from_all.txt", "w") as outfile:
        outfile.write("Gene\tLog2\tDenominator\tNumerator\n")
        for elem in F_WT_KO_uniq_list:
            outfile.write(f"{elem}\n")

    # log2 ratio for M_WT & M_KO
    M_WT_KO_outlist = get_log2_list(M_WT_KO_list)
    M_WT_KO_uniq_list = combine_same_gene(M_WT_KO_outlist)
    with open("log2_M_WT_KO_from_all.txt", "w") as outfile:
        outfile.write("Gene\tLog2\tDenominator/tNumerator\n")
        for elem in M_WT_KO_uniq_list:
            outfile.write(f"{elem}\n")

    # log2 ratio for M_KO & F_KO
    KO_FM_outlist = get_log2_list(KO_FM_list)
    KO_FM_uniq_list = combine_same_gene(KO_FM_outlist)
    with open("log2_KO_FM_from_all.txt", "w") as outfile:
        outfile.write("Gene\tLog2\tDenominator\tNumerator\n")
        for elem in KO_FM_uniq_list:
            outfile.write(f"{elem}\n")

    # log2 ratio for M_WT & F_WT
    WT_FM_outlist = get_log2_list(WT_FM_list)
    WT_FM_uniq_list = combine_same_gene(WT_FM_outlist)
    with open("log2_WT_FM_from_all.txt", "w") as outfile:
        outfile.write("Gene\tLog2\tDenominator\tNumerator\n")
        for elem in WT_FM_uniq_list:
            outfile.write(f"{elem}\n")
