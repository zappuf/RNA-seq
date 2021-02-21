with open("gene_symbol.txt", "w") as outfile:
    with open("F_KO_grouped.sf", "r") as infile:
        lines = infile.readlines()
        for line in lines:
            line = line.strip()
            ensembl = line.split("\t")[0]
            outfile.write(f"{ensembl}\n")
