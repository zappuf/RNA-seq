samples = {
    "F_KO": [
        "F-KO-32R_S16",
        "F-KO-41L_S17",
        "F-KO-41N_S18",
        "F-KO-44L_S19",
        "F-KO-53L_S20",
    ],
    "F_WT": [
        "F-WT-59B_S11",
        "F-WT-63B_S12",
        "F-WT-63R2_S13",
        "F-WT-65B_S14",
        "F-WT-71B_S15",
    ],
    "M_KO": [
        "M-KO-33N_S6",
        "M-KO-43R_S7",
        "M-KO-45N_S8",
        "M-KO-51L_S9",
        "M-KO-54N_S10",
    ],
    "M_WT": ["M-WT-40N_S2", "M-WT-40R_S1", "M-WT-45R_S3", "M-WT-56R_S4", "M-WT-58B_S5"],
}

for condition, samples_list in samples.items():
    print(f"working on condition {condition}")
    for idx in range(len(samples_list)):
        print(f"Working on sample {samples[condition][idx]}")
        with open(f"salmon_mapped32/{samples[condition][idx]}/quant.sf", "r") as infile:
            in_lines = infile.readlines()[1:]
            for line_idx in range(len(in_lines)):
                in_lines[line_idx] = in_lines[line_idx].strip()
                fields = in_lines[line_idx].split("\t")
                gene = fields[0]  # if idx = 0 else None
                count = fields[-1]
                in_lines[line_idx] = f"{gene}\t{count}\n" if not idx else f"{count}\n"
        if not idx:
            with open(f"{condition}_grouped.sf", "w+") as outfile:
                for in_line in in_lines:
                    in_line = in_line.strip()
                    outfile.write(f"{in_line}\n")
        else:
            with open(f"{condition}_grouped.sf", "r") as outfile:
                out_lines = outfile.readlines()
                for out_line_idx in range(len(out_lines)):
                    out_lines[out_line_idx] = out_lines[out_line_idx].strip()
                    fields = out_lines[out_line_idx].split("\t")
                    fields.append(in_lines[out_line_idx])
                    out_lines[out_line_idx] = "\t".join(fields)
            with open(f"{condition}_grouped.sf", "w") as outfile:
                outfile.writelines(out_lines)
