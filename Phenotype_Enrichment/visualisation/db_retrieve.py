from typing import List

import psycopg2
import pandas as pd


def select_from_pathway_genes(pathway_name: str) -> List[str]:
    genes = []
    try:
        conn = psycopg2.connect(user="wahiiuuseanslh",
                                password="3370a0e2c90b0d8eb13192a4de38b57556f0e98bc083e50d9157dd82b4d12619",
                                host="ec2-54-171-25-232.eu-west-1.compute.amazonaws.com",
                                port="5432",
                                database="d8re051vcjq8v1")
        cursor = conn.cursor()
        postgreSQL_select_Query = "SELECT GENES FROM PATHWAY_GENES WHERE PATHWAY_NAME = %s"
        cursor.execute(postgreSQL_select_Query, (pathway_name,))
        genes = cursor.fetchall()
        genes = genes[0][0]
        genes = [gene for gene in genes if str(gene) != 'nan']
        return genes
    except Exception as e:
        print(e)
        return genes


def select_from_enrichment_results(pathway_name: str) -> pd.DataFrame:
    """Select rows corresponding to the records for `pathway_name`.

    :param pathway_name: the name of the pathway to be passed in query
    :return: df corresponding to the records for `pathway_name`
    """
    genes = select_from_pathway_genes(pathway_name=pathway_name)
    human_genes = []
    ortholog_genes = []
    organism = []
    enriched_phenotypes = []
    result = pd.DataFrame()
    try:
        conn = psycopg2.connect(user="wahiiuuseanslh",
                                     password="3370a0e2c90b0d8eb13192a4de38b57556f0e98bc083e50d9157dd82b4d12619",
                                     host="ec2-54-171-25-232.eu-west-1.compute.amazonaws.com",
                                     port="5432",
                                     database="d8re051vcjq8v1")
        cursor = conn.cursor()
        postgreSQL_select_Query = "SELECT * FROM ENRICHMENT_RESULTS WHERE HUMAN_GENE = %s"
        for gene in genes:
            cursor.execute(postgreSQL_select_Query, (gene,))
            enrichment_results = cursor.fetchall()

            for row in enrichment_results:
                human_genes.append(gene)
                ortholog_genes.append(row[2])
                organism.append(row[3])
                enriched_phenotypes.append(row[4])

        cursor.close()
        conn.close()

        result["Ortholog_Genes"] = ortholog_genes
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
        conn = psycopg2.connect(user="wahiiuuseanslh",
                                     password="3370a0e2c90b0d8eb13192a4de38b57556f0e98bc083e50d9157dd82b4d12619",
                                     host="ec2-54-171-25-232.eu-west-1.compute.amazonaws.com",
                                     port="5432",
                                     database="d8re051vcjq8v1")
        cursor = conn.cursor()
        postgreSQL_select_Query = "SELECT TOP_LEVEL_PATHWAY FROM PATHWAY_HIERARCHY WHERE %s = ANY(LOW_LEVEL_PATHWAY)"
        cursor.execute(postgreSQL_select_Query, (pathway_name,))
        top_level_pathway = cursor.fetchall()
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
        conn = psycopg2.connect(user="wahiiuuseanslh",
                                     password="3370a0e2c90b0d8eb13192a4de38b57556f0e98bc083e50d9157dd82b4d12619",
                                     host="ec2-54-171-25-232.eu-west-1.compute.amazonaws.com",
                                     port="5432",
                                     database="d8re051vcjq8v1")
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
