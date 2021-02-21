import json
import math
import numpy as np
from util import convert_gene
from config import study_dir, results_dir, samples, pseudo_count, dynamic_pseudo_count


for sample_idx in range(len(samples)):
    print(f"Working on sample {samples[sample_idx]}")
    with open(f"salmon_mapped32/{samples[sample_idx]}/quant.sf", "r") as infile:
        in_lines = infile.readlines()[1:]
        for line_idx in range(len(in_lines)):
            in_lines[line_idx] = in_lines[line_idx].strip()
            fields = in_lines[line_idx].split("\t")
            gene = fields[0]
            count = fields[-1]
            in_lines[line_idx] = (
                f"{gene}\t{count}\n" if not sample_idx else f"{count}\n"
            )
    if not sample_idx:
        with open("all_grouped.sf", "w+") as outfile:
            for in_line in in_lines:
                in_line = in_line.strip()
                outfile.write(f"{in_line}\n")
    else:
        with open("all_grouped.sf", "r") as outfile:
            out_lines = outfile.readlines()
            for out_line_idx in range(len(out_lines)):
                out_lines[out_line_idx] = out_lines[out_line_idx].strip()
                fields = out_lines[out_line_idx].split("\t")
                fields.append(in_lines[out_line_idx])
                out_lines[out_line_idx] = "\t".join(fields)
        with open("all_grouped.sf", "w") as outfile:
            outfile.writelines(out_lines)

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
            # Keep only genes where at least 3 samples have count at least 10
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


F_KO = "F_KO_grouped.sf"
F_WT = "F_WT_grouped.sf"
M_KO = "M_KO_grouped.sf"
M_WT = "M_WT_grouped.sf"

F_WT_KO = "F_WT_KO_filtered_normalized.sf"
KO_FM = "KO_FM_filtered_normalized.sf"
M_WT_KO = "M_WT_KO_filtered_normalized.sf"
WT_FM = "WT_FM_filtered_normalized.sf"


def get_min(list1):
    try:
        min = list1[0] if list1[0] > 0.0 else float("inf")
    except Exception as e:
        print(f"Error assigning min: {e}\nlist1: {list1}")
    if type(list1) == list:
        for elem in list1:
            if elem < min and elem > 0.0:
                min = elem
    else:
        min = 0
    return min


def get_log2_list(list1, outname):
    """
    list1 is a list where each element is a list.
    out_dict structure
{
   "gene":{
      "ensembl1":{
         "denominator":[
            1,
            2,
            3,
            20
         ],
         "numerator":[
            1,
            2,
            3,
            20
         ]
      },
      "ensembl2":{
         "denominator":[
            1,
            2,
            3,
            20
         ],
         "numerator":[
            1,
            2,
            3,
            20
         ]
      }
   }
}
    """

    out_dict = {}
    for idx in range(len(list1)):
        ensembl = list1[idx][0]
        gene = convert_gene(ensembl)

        group1_fields = list1[idx][1:6]
        group2_fields = list1[idx][6:]

        assert len(group1_fields) == len(
            group2_fields
        ), "input list contains incorrect number of columns"

        group1_floats = [float(i) for i in group1_fields]
        group2_floats = [float(i) for i in group2_fields]

        if gene in out_dict.keys():
            out_dict[gene][ensembl] = {
                "denominator": group1_floats,
                "numerator": group2_floats,
            }
        else:
            out_dict[gene] = {
                ensembl: {"numerator": group2_floats, "denominator": group1_floats,}
            }

    for k1, v1 in out_dict.items():
        denominator = []
        numerator = []
        for k2 in v1.keys():
            # print(
            #     f"denominator: {denominator}, numerator: {numerator}, v1[{k2}]['denominator']: {v1[k2]['denominator']}, v1[{k2}]['numerator']: {v1[k2]['numerator']}"
            # )
            denominator += v1[k2]["denominator"]
            numerator += v1[k2]["numerator"]
        out_dict[k1]["denominator_tc"] = sum(denominator)
        out_dict[k1]["numerator_tc"] = sum(numerator)
        if dynamic_pseudo_count:
            gene_counts = denominator + numerator
            min = get_min(gene_counts)
            pseudo_count = min * 0.1
            # print(f"min and pseudo_count for gene {k1}: {min}, {pseudo_count}")
        group1_avg = sum(denominator) / len(denominator)
        group2_avg = sum(numerator) / len(numerator)
        log2ratio = math.log2((group2_avg + pseudo_count) / (group1_avg + pseudo_count))
        v1["log2"] = log2ratio

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
    print("Writing log2 file for F_WT_KO")
    with open(f"{results_dir}/log2_F_WT_KO_from_all.txt", "w") as outfile:
        outfile.write("Gene\tLog2\tControl\tTest\n")
        for key in F_WT_KO_keys:
            gene = key
            log2 = F_WT_KO_outdict[key]["log2"]
            denominator = F_WT_KO_outdict[key]["denominator_tc"]
            numerator = F_WT_KO_outdict[key]["numerator_tc"]

            outfile.write(
                f"{gene}\t{round(log2, 2)}\t{round(denominator, 3)}\t{round(numerator, 3)}\n"
            )

    # log2 ratio for M_WT & M_KO
    M_WT_KO_outdict = get_log2_list(M_WT_KO_list, "M_WT_KO")
    M_WT_KO_keys = list(M_WT_KO_outdict.keys())
    M_WT_KO_keys.sort()
    print("Writing log2 file for M_WT_KO")
    with open(f"{results_dir}/log2_M_WT_KO_from_all.txt", "w") as outfile:
        outfile.write("Gene\tLog2\tControl\tTest\n")
        for key in M_WT_KO_keys:
            gene = key
            log2 = M_WT_KO_outdict[key]["log2"]
            denominator = M_WT_KO_outdict[key]["denominator_tc"]
            numerator = M_WT_KO_outdict[key]["numerator_tc"]
            outfile.write(
                f"{gene}\t{round(log2, 2)}\t{round(denominator, 3)}\t{round(numerator, 3)}\n"
            )

    # log2 ratio for M_KO & F_KO
    KO_FM_outdict = get_log2_list(KO_FM_list, "KO_FM")
    KO_FM_keys = list(KO_FM_outdict.keys())
    KO_FM_keys.sort()
    print("Writing log2 file for KO_FM")
    with open(f"{results_dir}/log2_KO_FM_from_all.txt", "w") as outfile:
        outfile.write("Gene\tLog2\tControl\tTest\n")
        # outfile.write("Gene\tLog2\tDenominator\tNumerator\n")
        for key in KO_FM_keys:
            gene = key
            log2 = KO_FM_outdict[key]["log2"]
            denominator = KO_FM_outdict[key]["denominator_tc"]
            numerator = KO_FM_outdict[key]["numerator_tc"]
            outfile.write(
                f"{gene}\t{round(log2, 2)}\t{round(denominator, 3)}\t{round(numerator, 3)}\n"
            )

    # log2 ratio for M_WT & F_WT
    WT_FM_outdict = get_log2_list(WT_FM_list, "WT_FM")
    WT_FM_keys = list(WT_FM_outdict.keys())
    WT_FM_keys.sort()
    print("Writing log2 file for WT_FM")
    with open(f"{results_dir}/log2_WT_FM_from_all.txt", "w") as outfile:
        outfile.write("Gene\tLog2\tControl\tTest\n")
        # outfile.write("Gene\tLog2\tDenominator\tNumerator\n")
        for key in WT_FM_keys:
            gene = key
            log2 = WT_FM_outdict[key]["log2"]
            denominator = WT_FM_outdict[key]["denominator_tc"]
            numerator = WT_FM_outdict[key]["numerator_tc"]
            outfile.write(
                f"{gene}\t{round(log2, 2)}\t{round(denominator, 3)}\t{round(numerator, 3)}\n"
            )
