from typing import List

import psycopg2
import pandas as pd


def select_from_pathway_genes(pathway_name: str) -> List[str]:
    """Retrieve a list of genes belonging to a pathway.
    :param pathway_name: the name of the pathway
    :return: a list of gene IDs belonging to a pathway
    """
    genes = []
    try:
        conn = psycopg2.connect(user="lfxnboorsrhevw",
                                password="4322a747d9e7b86cb62c2ef1e44b338a9fe059ce99cf6d662278fd21ed06e388",
                                host="ec2-54-216-159-235.eu-west-1.compute.amazonaws.com",
                                port="5432",
                                database="dft1uk8fl9qnpb")
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


def select_from_gene_names(human_genes: list) -> list:
    """Select names for the list of ENSG IDs.
    :param human_genes: a list of ENSG IDs
    :return: a list of names
    """
    names = []
    try:
        conn = psycopg2.connect(user="lfxnboorsrhevw",
                                password="4322a747d9e7b86cb62c2ef1e44b338a9fe059ce99cf6d662278fd21ed06e388",
                                host="ec2-54-216-159-235.eu-west-1.compute.amazonaws.com",
                                port="5432",
                                database="dft1uk8fl9qnpb")
        cursor = conn.cursor()
        for id in human_genes:
            postgreSQL_select_Query = "SELECT NAME FROM GENE_NAMES WHERE GENE_ID = %s"
            cursor.execute(postgreSQL_select_Query, (id,))
            name = cursor.fetchall()
            name = name[0][0]
            names.append(name)
        return names
    except Exception as e:
        print(e)
        return names


def select_from_enrichment_results(pathway_name: str) -> pd.DataFrame:
    """Select rows corresponding to the records for `pathway_name`.
    :param pathway_name: the name of the pathway to be passed in query
    :return: df corresponding to the records for `pathway_name`
    """
    metadata = select_from_metadata(pathway_name)
    genes = select_from_pathway_genes(pathway_name=pathway_name)
    human_genes = []
    ortholog_genes = []
    organism = []
    enriched_phenotypes = []
    result = pd.DataFrame()
    try:
        conn = psycopg2.connect(user="lfxnboorsrhevw",
                                password="4322a747d9e7b86cb62c2ef1e44b338a9fe059ce99cf6d662278fd21ed06e388",
                                host="ec2-54-216-159-235.eu-west-1.compute.amazonaws.com",
                                port="5432",
                                database="dft1uk8fl9qnpb")
        cursor = conn.cursor()
        postgreSQL_select_Query = "SELECT * FROM ENRICHMENT_RESULTS WHERE HUMAN_GENE = %s"
        for gene in genes:
            cursor.execute(postgreSQL_select_Query, (gene,))
            enrichment_results = cursor.fetchall()

            for row in enrichment_results:
                new_phen_list = []
                human_genes.append(gene)
                ortholog_genes.append(row[2])
                organism.append(row[3])
                phenotype_list = row[4]
                for phen in phenotype_list:
                    if metadata.get(phen) is not None and metadata.get(phen)[3] <= 10:
                        new_phen_list.append(phen)

                enriched_phenotypes.append(new_phen_list)

        cursor.close()
        conn.close()

        result["Ortholog_Genes"] = ortholog_genes
        result["Human_ID"] = human_genes  # select_from_gene_names(human_genes)
        result["Organism"] = organism
        result["Enriched_Phenotypes"] = enriched_phenotypes
        result["Human_Gene"] = select_from_gene_names(human_genes)
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
        conn = psycopg2.connect(user="lfxnboorsrhevw",
                                password="4322a747d9e7b86cb62c2ef1e44b338a9fe059ce99cf6d662278fd21ed06e388",
                                host="ec2-54-216-159-235.eu-west-1.compute.amazonaws.com",
                                port="5432",
                                database="dft1uk8fl9qnpb")
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
        conn = psycopg2.connect(user="lfxnboorsrhevw",
                                password="4322a747d9e7b86cb62c2ef1e44b338a9fe059ce99cf6d662278fd21ed06e388",
                                host="ec2-54-216-159-235.eu-west-1.compute.amazonaws.com",
                                port="5432",
                                database="dft1uk8fl9qnpb")
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