import psycopg2


def create_tables():
    """ create tables in the PostgreSQL database"""
    commands = (
        """
        CREATE TABLE ENRICHMENT_RESULTS (
            ID SERIAL PRIMARY KEY,
            LOW_LEVEL_PATHWAY VARCHAR(255) NOT NULL,
            TOP_LEVEL_PATHWAY VARCHAR(255),
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
        conn = psycopg2.connect(user="scfhbchnxiyzrp",
                                password="146d74fa3f155b51373badc4fdf76a315d36f0db7b4a0cc61bb69bfae4fcba0c",
                                host="ec2-54-74-95-84.eu-west-1.compute.amazonaws.com",
                                port="5432",
                                database="dcfbco2ll2unns")
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
