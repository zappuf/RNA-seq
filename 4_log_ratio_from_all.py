import json
import math
import numpy as np
import statistics as stat
from util import convert_gene
from config import study_dir, results_dir

pseudo_count = 0.01

F_KO = "F_KO_grouped.sf"
F_WT = "F_WT_grouped.sf"
M_KO = "M_KO_grouped.sf"
M_WT = "M_WT_grouped.sf"

F_WT_KO = "F_WT_KO_filtered_normalized.sf"
KO_FM = "KO_FM_filtered_normalized.sf"
M_WT_KO = "M_WT_KO_filtered_normalized.sf"
WT_FM = "WT_FM_filtered_normalized.sf"


def get_min(list1):
    min = list1[0] if list1[0] > 0.0 else float("inf")
    for elem in list1:
        if elem < min and elem > 0.0:
            min = elem
    return min


def get_log2_list(list1, outname):
    """
    list1 is a list where each element is a list
    """
    out_dict = {}
    for idx in range(len(list1)):
        ensembl = list1[idx][0]
        gene = convert_gene(ensembl)

        group1_fields = list1[idx][1:6]
        group2_fields = list1[idx][6:]

        group1_floats = [float(i) for i in group1_fields]
        group2_floats = [float(i) for i in group2_fields]

        group_floats = group1_floats + group2_floats

        min = get_min(group_floats)
        pseudo_count = min * 0.1 if min else 0.01
        print(f"Min for {ensembl} / {gene}: {min}, pseudo_count: {pseudo_count}")

        group1_avg = sum(group1_floats) / len(group1_floats)
        group2_avg = sum(group2_floats) / len(group2_floats)

        if gene in out_dict.keys():
            out_dict[gene]["denominator"].append(group1_avg)
            out_dict[gene]["numerator"].append(group2_avg)
            out_dict[gene]["log2"] = math.log2(
                (sum(out_dict[gene]["numerator"]) + pseudo_count)
                / (sum(out_dict[gene]["denominator"]) + pseudo_count)
            )
        else:
            log2ratio = math.log2(
                ((group2_avg + pseudo_count) / (group1_avg + pseudo_count))
            )
            out_dict[gene] = {
                "log2": log2ratio,
                "numerator": [group2_avg],
                "denominator": [group1_avg],
            }
    with open(f"{results_dir}/{outname}_counts.json", "w") as fp:
        fp.write(json.dumps(out_dict))

    print(f"out_dict has {len(out_dict)} entries")
    return out_dict


##############################################################
with open(f"{results_dir}/all_grouped_filtered_normalized.sf", "r") as infile:
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
        M_WT = fields[16:21]

        F_WT_KO_line = ensembl + F_WT + F_KO
        F_WT_KO_list.append(F_WT_KO_line)

        M_WT_KO_line = ensembl + M_WT + M_KO
        M_WT_KO_list.append(M_WT_KO_line)

        KO_FM_line = ensembl + F_KO + M_KO
        KO_FM_list.append(KO_FM_line)

        WT_FM_line = ensembl + F_WT + M_WT
        WT_FM_list.append(WT_FM_line)

    # log2 ratio for F_WT & F_KO
    F_WT_KO_outdict = get_log2_list(F_WT_KO_list, "F_WT_KO")
    F_WT_KO_keys = list(F_WT_KO_outdict.keys())
    F_WT_KO_keys.sort()
    with open(f"{results_dir}/log2_F_WT_KO_from_all.txt", "w") as outfile:
        outfile.write("Gene\tLog2\tDenominator\tNumerator\n")
        for key in F_WT_KO_keys:
            gene = key
            log2 = F_WT_KO_outdict[key]["log2"]
            denominator = F_WT_KO_outdict[key]["denominator"]
            numerator = F_WT_KO_outdict[key]["numerator"]
            outfile.write(
                f"{gene}\t{round(log2, 2)}\t{round(sum(denominator), 3)}\t{round(sum(numerator), 3)}\n"
            )

    # log2 ratio for M_WT & M_KO
    M_WT_KO_outdict = get_log2_list(M_WT_KO_list, "M_WT_KO")
    M_WT_KO_keys = list(M_WT_KO_outdict.keys())
    M_WT_KO_keys.sort()
    with open(f"{results_dir}/log2_M_WT_KO_from_all.txt", "w") as outfile:
        outfile.write("Gene\tLog2\tDenominator\tNumerator\n")
        for key in M_WT_KO_keys:
            gene = key
            log2 = M_WT_KO_outdict[key]["log2"]
            denominator = M_WT_KO_outdict[key]["denominator"]
            numerator = M_WT_KO_outdict[key]["numerator"]
            outfile.write(
                f"{gene}\t{round(log2, 2)}\t{round(sum(denominator), 3)}\t{round(sum(numerator), 3)}\n"
            )

    # log2 ratio for M_KO & F_KO
    KO_FM_outdict = get_log2_list(KO_FM_list, "KO_FM")
    KO_FM_keys = list(KO_FM_outdict.keys())
    KO_FM_keys.sort()
    with open(f"{results_dir}/log2_KO_FM_from_all.txt", "w") as outfile:
        outfile.write("Gene\tLog2\tDenominator\tNumerator\n")
        # outfile.write("Gene\tLog2\tDenominator\tNumerator\n")
        for key in KO_FM_keys:
            gene = key
            log2 = KO_FM_outdict[key]["log2"]
            denominator = KO_FM_outdict[key]["denominator"]
            numerator = KO_FM_outdict[key]["numerator"]
            outfile.write(
                f"{gene}\t{round(log2, 2)}\t{round(sum(denominator), 3)}\t{round(sum(numerator), 3)}\n"
            )

    # log2 ratio for M_WT & F_WT
    WT_FM_outdict = get_log2_list(WT_FM_list, "WT_FM")
    WT_FM_keys = list(WT_FM_outdict.keys())
    WT_FM_keys.sort()
    with open(f"{results_dir}/log2_WT_FM_from_all.txt", "w") as outfile:
        outfile.write("Gene\tLog2\tDenominator\tNumerator\n")
        # outfile.write("Gene\tLog2\tDenominator\tNumerator\n")
        for key in WT_FM_keys:
            gene = key
            log2 = WT_FM_outdict[key]["log2"]
            denominator = WT_FM_outdict[key]["denominator"]
            numerator = WT_FM_outdict[key]["numerator"]
            outfile.write(
                f"{gene}\t{round(log2, 2)}\t{round(sum(denominator), 3)}\t{round(sum(numerator), 3)}\n"
            )
