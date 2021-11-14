import json
import pandas as pd
import psycopg2

####################################################################################
# Define functions for inserting data into database ################################
####################################################################################


def insert_into_enrichment_results(df: pd.DataFrame) -> None:
    """Insert a dataframe of enrichment results for a given pathway.

    :param df: dataframe of the enrichment results with cols corresponding to the cols in the DB table
    """
    global cursor, conn
    try:
        conn = psycopg2.connect(user="scfhbchnxiyzrp",
                                password="146d74fa3f155b51373badc4fdf76a315d36f0db7b4a0cc61bb69bfae4fcba0c",
                                host="ec2-54-74-95-84.eu-west-1.compute.amazonaws.com",
                                port="5432",
                                database="dcfbco2ll2unns")
        cursor = conn.cursor()
        for row in range(0, len(df)):
            postgres_insert_query = """INSERT INTO ENRICHMENT_RESULTS (LOW_LEVEL_PATHWAY, TOP_LEVEL_PATHWAY, 
            HUMAN_GENE, ORTHOLOG_GENE, ORGANISM, ENRICHED_PHENOTYPES) VALUES (%s,%s,%s,%s,%s,%s)"""
            record_to_insert = (df.iloc[row, 4], df.iloc[row, 5],
                                df.iloc[row, 1], df.iloc[row, 0],
                                df.iloc[row, 2], df.iloc[row, 3])
            cursor.execute(postgres_insert_query, record_to_insert)

            conn.commit()
            count = cursor.rowcount
            print(count, "Record inserted successfully into mobile table")

        # Close the DB connection
        cursor.close()
        conn.close()
        print("PostgreSQL connection is closed")

    except (Exception, psycopg2.Error) as error:
        print("Failed to insert record into mobile table", error)


def insert_into_phenotype_metadata(pathway_name: str, metadata: json) -> None:
    """Insert a record corresponding to the json obj of metadata for app phenotypes of a given pathway.

    :param pathway_name: the name of the pathway to which will be the phenotypes associated
    :param metadata: a json obj of metadata for all phenotypes of the pathway
    """
    global conn, cursor
    try:
        conn = psycopg2.connect(user="scfhbchnxiyzrp",
                                password="146d74fa3f155b51373badc4fdf76a315d36f0db7b4a0cc61bb69bfae4fcba0c",
                                host="ec2-54-74-95-84.eu-west-1.compute.amazonaws.com",
                                port="5432",
                                database="dcfbco2ll2unns")
        cursor = conn.cursor()

        postgres_insert_query = """INSERT INTO PHENOTYPE_METADATA (PATHWAY, METADATA) VALUES (%s,%s)"""
        record_to_insert = (pathway_name, metadata)
        cursor.execute(postgres_insert_query, record_to_insert)

        conn.commit()

        # Close the DB connection
        cursor.close()
        conn.close()
        print("PostgreSQL connection is closed")

    except (Exception, psycopg2.Error) as error:
        print("Failed to insert record into mobile table", error)
