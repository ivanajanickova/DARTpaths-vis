import pickle
from typing import Dict

import pandas as pd
import os


class DataExtraction():
    """Handle extraction and saving of data from `Python_APRIL_2021_phenotype_enrichment.py`.
    A new instance of the object is created upon running the `Python_APRIL_2021_phenotype_enrichment.py `
    for a given pathway. It saves the data. The saved file is name contains the `pathway_name`.
    This class is called from the `Python_APRIL_2021_phenotype_enrichment.py` script and all of the files are saved
    once the script finishes.
    """

    def __init__(self):
        self.orthologs_phenotype_df = pd.DataFrame()
        self.genes_orthologs_df = pd.DataFrame()
        self.pathway_name = ""

    def add_ortholog_vs_phenotype_data(self, uniquegenes: pd.DataFrame, organism: str) -> None:
        """Add df containing gene IDs of orthologs and the enriched phenotypes.

        :param uniquegenes: df containing gene IDs of orthologs and the enriched phenotypes
        :param organims: the name of the organism
        """
        if uniquegenes is not None and organism is not None:
            uniquegenes["Organism"] = [organism for i in range(0, len(uniquegenes))]
            data_frames = [uniquegenes, self.orthologs_phenotype_df]
            self.orthologs_phenotype_df = pd.concat(data_frames)

    def add_genes_vs_orthologs_data(self, df_ortholog: pd.DataFrame, organism: str) -> None:
        """Add df containing gene IDs of orthologs and the enriched phenotypes.

        :param df_ortholog: df containing human Ensembl gene IDs and the gene IDs of orthologs from different organisms
        :param organims: the name of the organism
        """
        if df_ortholog is not None and organism is not None:
            df_ortholog["Organism"] = [organism for i in range(0, len(df_ortholog))]
            data_frames = [df_ortholog, self.genes_orthologs_df]
            self.genes_orthologs_df = pd.concat(data_frames)

    def save_dfs_to_files(self, pathway_name: str) -> None:
        """Save extracted data from two DataFrames as csv files.
        The file will be saved in the pathway enrichment folder eg `AHR_R-HSA-8937144_Enrichment_Results`.

        :param pathway_name: the name of the given pathway. Will be included in the file names
        """
        path_phenotypes_file = os.getcwd() + "/" + pathway_name + "_orthologs_to_phenotype_data.csv"
        path_orthologs_file = os.getcwd() + "/" + pathway_name + "_genes_orthologs_data.csv"

        self.orthologs_phenotype_df.to_csv(pathway_name + "_orthologs_to_phenotype_data.csv", index=False)
        self.genes_orthologs_df.to_csv(pathway_name + "_genes_orthologs_data.csv", index=False)

