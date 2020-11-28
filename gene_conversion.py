import json

conversion_dict = {
    "Gene_stable_ID":[],
    "Gene_stable_ID_version":[],
    "Transcript_stable_ID":[],
    "Transcript_stable_ID_version":[],
    "Gene_symbol": [],
}

with open("mart_export.txt", "r") as mart_file:
    with open("all_gene_symbol.txt", "r") as gene_file:
        gene_lines = gene_file.readlines()
        mart_lines = mart_file.readlines()[1:]

        for idx in range(len(mart_lines)):
            mart_line = mart_lines[idx].strip()
            mart_fields = mart_line.split("\t")
            conversion_dict["Gene_stable_ID"].append(mart_fields[0])
            conversion_dict["Gene_stable_ID_version"].append(mart_fields[1])
            conversion_dict["Transcript_stable_ID"].append(mart_fields[2])
            conversion_dict["Transcript_stable_ID_version"].append(mart_fields[3])
            try:
                gene_symbol = gene_lines[idx].strip()
            except Exception as e:
                print(f"Error with assigning gene_symbol: {e}")
                gene_symbol = "N/A"
            conversion_dict["Gene_symbol"].append(gene_symbol)
            # print(fields)
            # break
        with open("gene_conversion.json", "w") as fp:
            fp.write(json.dumps(conversion_dict))
