import json
from config import study_dir

conversion_list = [
    "Gene_stable_ID",
    "Gene_stable_ID_version",
    "Transcript_stable_ID",
    "Transcript_stable_ID_version",
    "Gene_symbol",
]


def get_gene_data():
    with open(f"{study_dir}/gene_conversion.json", "r") as infile:
        gene_data = json.load(infile)
    return gene_data


# gene_data = get_gene_data()


def convert_gene(
    id,
    id_from="Transcript_stable_ID_version",
    id_to="Gene_symbol",
    gene_data=get_gene_data(),
):
    """
    Will return the gene symbol when given a Transcript_stable_ID_version if
    able to resolve Transcript_stable_ID_version to gene symbol.
    Otherwise returns the Transcript_stable_ID_version.
    """
    if id_from not in conversion_list or id_to not in conversion_list:
        return None
    try:
        index = gene_data[id_from].index(id)
    except:
        index = None
    return gene_data[id_to][index] if index else id
