"""Update database relation `GENE_NAMES` to the database."""

import psycopg2
import pandas as pd
from db_inserts import DatabaseInserts

DBInsert = DatabaseInserts()


def run_update_pipeline():
    """Run the pipline that inserts IDs and names into the database."""
    genes = get_gene_ids()
    df = assign_names(genes)
    DBInsert.insert_into_gene_names(df)


def get_gene_ids() -> list:
    """Retrieve all gene IDs from the database.
    :return: list of gene IDs
    """
    genes = {}
    try:
        conn = psycopg2.connect(user="lfxnboorsrhevw",
                                password="4322a747d9e7b86cb62c2ef1e44b338a9fe059ce99cf6d662278fd21ed06e388",
                                host="ec2-54-216-159-235.eu-west-1.compute.amazonaws.com",
                                port="5432",
                                database="dft1uk8fl9qnpb")
        cursor = conn.cursor()
        select_Query = "SELECT distinct human_gene from enrichment_results"
        cursor.execute(select_Query)
        gene_ids = cursor.fetchall()
        cursor.close()
        conn.close()
        genes = [id[0] for id in gene_ids]
        return genes

    except Exception as e:
        print(e)
        print("Error in `get_gene_ids`")
        return genes


def assign_names(genes: list) -> pd.DataFrame:
    """Based on the `Homo_sapiens.GRCh38.104.chr.gtf` assign names to the IDs.
    For downloading the `Homo_sapiens.GRCh38.104.chr.gtf` go to Enseble db.
    :param genes: a list of IDs
    :return: a dataframe associating IDs with names
    """
    gene_ids = {id for id in genes}
    names = []
    ids = []
    with open("/home/ivana/DARTpaths-vis/Phenotype_Enrichment/ontology_data/Homo_sapiens.GRCh38.104.chr.gtf",
              "r") as file:
        for line in file:
            if "#" not in line:
                line_list = line.split(";")
                ens_id = line_list[0].split('"')[1].strip()
                if ens_id in gene_ids and ens_id not in ids:
                    name = line_list[2].split('"')[1].strip()
                    ids.append(ens_id)
                    names.append(name)

    df = pd.DataFrame()
    df["Gene_IDs"] = ids
    df["Names"] = names
    return df


if __name__ == '__main__':
    run_update_pipeline()