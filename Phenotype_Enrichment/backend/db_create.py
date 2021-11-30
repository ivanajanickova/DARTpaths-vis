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
        conn = psycopg2.connect(user="lfxnboorsrhevw",
                                password="4322a747d9e7b86cb62c2ef1e44b338a9fe059ce99cf6d662278fd21ed06e388",
                                host="ec2-54-216-159-235.eu-west-1.compute.amazonaws.com",
                                port="5432",
                                database="dft1uk8fl9qnpb")
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
