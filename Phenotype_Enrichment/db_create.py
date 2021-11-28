import psycopg2


def create_tables():
    """ create tables in the PostgreSQL database"""
    commands = (
        """
        CREATE TABLE PATHWAY_GENES (
            PATHWAY_NAME VARCHAR(255) NOT NULL PRIMARY KEY,
            GENES VARCHAR(255)[])
            """,

        """
        CREATE TABLE PATHWAY_HIERARCHY (
            TOP_LEVEL_PATHWAY VARCHAR(255) NOT NULL PRIMARY KEY,
            LOW_LEVEL_PATHWAY VARCHAR(255)[])
        """,

        """
        CREATE TABLE ENRICHMENT_RESULTS ( 
            ID SERIAL PRIMARY KEY,
            HUMAN_GENE VARCHAR(255) NOT NULL,
            ORTHOLOG_GENE VARCHAR(255),
            ORGANISM VARCHAR(255),
            ENRICHED_PHENOTYPES VARCHAR(255)[],
            UNIQUE (ID))
        """,

        """
        CREATE TABLE PHENOTYPE_METADATA (
            PATHWAY VARCHAR(255) NOT NULL PRIMARY KEY,
            METADATA JSON)
        """)
    try:
        conn = psycopg2.connect(user="wahiiuuseanslh",
                                password="3370a0e2c90b0d8eb13192a4de38b57556f0e98bc083e50d9157dd82b4d12619",
                                host="ec2-54-171-25-232.eu-west-1.compute.amazonaws.com",
                                port="5432",
                                database="d8re051vcjq8v1")
        cursor = conn.cursor()
        for command in commands:
            cursor.execute(command)
        cursor.close()
        conn.commit()
        conn.close()
    except (Exception, psycopg2.Error) as error:
        print("Failed to insert record into mobile table", error)



if __name__ == '__main__':
   create_tables()