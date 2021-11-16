import psycopg2
import pandas as pd


def select_from_enrichment_results(pathway_name: str) -> pd.DataFrame:
    """Select rows corresponding to the records for `pathway_name`.

    :param pathway_name: the name of the pathway to be passed in query
    :return: df corresponding to the records for `pathway_name`
    """
    human_genes = []
    ortholog_genes = []
    organism = []
    enriched_phenotypes = []
    result = pd.DataFrame()
    try:
        conn = psycopg2.connect(user="scfhbchnxiyzrp",
                                password="146d74fa3f155b51373badc4fdf76a315d36f0db7b4a0cc61bb69bfae4fcba0c",
                                host="ec2-54-74-95-84.eu-west-1.compute.amazonaws.com",
                                port="5432",
                                database="dcfbco2ll2unns")
        cursor = conn.cursor()
        postgreSQL_select_Query = "SELECT * FROM ENRICHMENT_RESULTS WHERE LOW_LEVEL_PATHWAY = %s"

        cursor.execute(postgreSQL_select_Query, (pathway_name,))
        enrichment_results = cursor.fetchall()

        print("Print each row and it's columns values")
        for row in enrichment_results:
            human_genes.append(row[3])
            ortholog_genes.append(row[4])
            organism.append(row[5])
            enriched_phenotypes.append(row[6])

        cursor.close()
        conn.close()
        print("PostgreSQL connection is closed")

        result["Orthlog_Genes"] = ortholog_genes
        result["Human_Gene"] = human_genes
        result["Organism"] = organism
        result["Enriched_Phenotypes"] = enriched_phenotypes
        return result

    except (Exception, psycopg2.Error) as error:
        print("Error while fetching data from PostgreSQL", error)
        return result


def find_top_level_pathway(pathway_name: str) -> str:
    """Find the name of the higher level pathway - if the pathway is the highest level None is returned.

    :param pathway_name: the name of the pathway to be passed in query
    :return: str corresponding to the higher level pathway name
    """
    top_level_pathway = ""
    try:
        conn = psycopg2.connect(user="scfhbchnxiyzrp",
                                password="146d74fa3f155b51373badc4fdf76a315d36f0db7b4a0cc61bb69bfae4fcba0c",
                                host="ec2-54-74-95-84.eu-west-1.compute.amazonaws.com",
                                port="5432",
                                database="dcfbco2ll2unns")
        cursor = conn.cursor()
        postgreSQL_select_Query = "SELECT TOP_LEVEL_PATHWAY FROM ENRICHMENT_RESULTS WHERE LOW_LEVEL_PATHWAY = %s " \
                                  "LIMIT 1 "
        cursor.execute(postgreSQL_select_Query, (pathway_name,))
        top_level_pathway = cursor.fetchall()
        # json_metadata is a list object with one record of type tuple which contains the pathway name and metadata dict
        cursor.close()
        conn.close()
        print("PostgreSQL connection is closed")
        return top_level_pathway

    except (Exception, psycopg2.Error) as error:
        print("Error while fetching data from PostgreSQL", error)
        return top_level_pathway


def select_from_metadata(pathway_name: str) -> dict:
    """Select rows corresponding to the records for `pathway_name`.

    :param pathway_name: the name of the pathway to be passed in query
    :return: dictionary of metadata corresponding `pathway_name`
    """
    dict_metadata = {}
    try:
        conn = psycopg2.connect(user="scfhbchnxiyzrp",
                                password="146d74fa3f155b51373badc4fdf76a315d36f0db7b4a0cc61bb69bfae4fcba0c",
                                host="ec2-54-74-95-84.eu-west-1.compute.amazonaws.com",
                                port="5432",
                                database="dcfbco2ll2unns")
        cursor = conn.cursor()
        postgreSQL_select_Query = "SELECT * FROM PHENOTYPE_METADATA WHERE PATHWAY = %s"
        cursor.execute(postgreSQL_select_Query, (pathway_name,))
        json_metadata = cursor.fetchall()
        # json_metadata is a list object with one record of type tuple which contains the pathway name and metadata dict
        tuple = json_metadata[0]
        dict_metadata = tuple[1]

        cursor.close()
        conn.close()
        print("PostgreSQL connection is closed")
        return dict_metadata

    except (Exception, psycopg2.Error) as error:
        print("Error while fetching data from PostgreSQL", error)
        return dict_metadata
