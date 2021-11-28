import json
import pandas as pd
import psycopg2


####################################################################################
# Define functions for inserting data into database ################################
####################################################################################

class DatabaseInserts:

    def __init__(self):
        self.conn = psycopg2.connect(user="wahiiuuseanslh",
                                     password="3370a0e2c90b0d8eb13192a4de38b57556f0e98bc083e50d9157dd82b4d12619",
                                     host="ec2-54-171-25-232.eu-west-1.compute.amazonaws.com",
                                     port="5432",
                                     database="d8re051vcjq8v1")
        self.cursor = self.conn.cursor()

    def close_connection(self):
        # Close the DB connection
        self.cursor.close()
        self.conn.close()
        print("PostgreSQL connection is closed")

    def insert_into_enrichment_results(self, df: pd.DataFrame) -> None:
        """Insert a dataframe of enrichment results for a given pathway.

        :param df: dataframe of the enrichment results with cols corresponding to the cols in the DB table
        """
        for row in range(0, len(df)):
            try:
                postgres_insert_query = """INSERT INTO ENRICHMENT_RESULTS (HUMAN_GENE, ORTHOLOG_GENE, 
                ORGANISM, ENRICHED_PHENOTYPES) VALUES (%s,%s,%s,%s)"""
                record_to_insert = (df.iloc[row, 1], df.iloc[row, 0],
                                    df.iloc[row, 2], df.iloc[row, 3])
                self.cursor.execute(postgres_insert_query, record_to_insert)
                self.conn.commit()
            except:
                print("Record already exists.")
                self.conn.rollback()

    def insert_into_phenotype_metadata(self, pathway_name: str, metadata: json) -> None:
        """Insert a record corresponding to the json obj of metadata for app phenotypes of a given pathway.

        :param pathway_name: the name of the pathway to which will be the phenotypes associated
        :param metadata: a json obj of metadata for all phenotypes of the pathway
        """
        postgres_insert_query = """INSERT INTO PHENOTYPE_METADATA (PATHWAY, METADATA) VALUES (%s,%s)"""
        record_to_insert = (pathway_name, metadata)
        self.cursor.execute(postgres_insert_query, record_to_insert)

        self.conn.commit()

    def insert_into_pathway_hierarchy(self, top_level_pathway: str, low_level_patway: str) -> None:
        """Insert into pathway hierarchy table.

        :param top_level_pathway: the higher level pathway
        :param low_level_patway: the lower level pathway"""
        self.conn.rollback()
        try:
            postgreSQL_select_Query = "SELECT * FROM PATHWAY_HIERARCHY WHERE TOP_LEVEL_PATHWAY = %s"
            self.cursor.execute(postgreSQL_select_Query, (top_level_pathway,))
            low_level_pathways_list = self.cursor.fetchall()
            if len(low_level_pathways_list) != 0:
                postgres_update_query = """ UPDATE PATHWAY_HIERARCHY 
                                                            SET LOW_LEVEL_PATHWAY = %s  
                                                            WHERE TOP_LEVEL_PATHWAY = %s"""
                record_to_insert = (low_level_pathways_list.append(low_level_patway), top_level_pathway)
                self.cursor.execute(postgres_update_query, record_to_insert)
                self.conn.commit()
            else:
                postgres_insert_query = """INSERT INTO PATHWAY_HIERARCHY (TOP_LEVEL_PATHWAY, LOW_LEVEL_PATHWAY) VALUES (%s,
                %s) """
                record_to_insert = (top_level_pathway, [low_level_patway])
                self.cursor.execute(postgres_insert_query, record_to_insert)
                self.conn.commit()
        except Exception as e:
            print(e)
            print("Failed")

    def insert_into_pathway_genes(self, pathway_name: str, genes_names: list) -> None:
        """Insert into pathway genes table.
        :param pathway_name: the name of the pathway
        :param genes_names: a list of unique human gene names that make up the pathway"""
        try:
            postgres_insert_query = """INSERT INTO PATHWAY_GENES (PATHWAY_NAME, GENES) VALUES (%s,%s)"""
            record_to_insert = (pathway_name, genes_names)
            self.cursor.execute(postgres_insert_query, record_to_insert)
            self.conn.commit()
        except:
            print("already exists")

    def close_conn(self):
        """Close connection after all inserts are finished."""
        self.conn.close()
