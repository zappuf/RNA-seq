import json

conversion_list = [
    "Gene_stable_ID",
    "Gene_stable_ID_version",
    "Transcript_stable_ID",
    "Transcript_stable_ID_version",
    "Gene_symbol",
]


def get_gene_data():
    with open("gene_conversion.json", "r") as infile:
        gene_data = json.load(infile)
    return gene_data


gene_data = get_gene_data()


def convert_gene(
    id,
    id_from="Transcript_stable_ID_version",
    id_to="Gene_symbol",
    gene_data=gene_data,
):
    if id_from not in conversion_list or id_to not in conversion_list:
        return None
    try:
        index = gene_data[id_from].index(id)
    except:
        index = None
    return gene_data[id_to][index] if index else id
