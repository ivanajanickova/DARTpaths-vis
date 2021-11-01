import pickle
from typing import Dict, List

import pandas as pd
import pickle
from pathlib import Path


def extract_phenotypes_info(ontology_files_list: List[str]) -> None:
    """Extract and save related phenotypes into a pickled dictionary.
    Extract and save the names of phenotypes into a pickled dictionary.
    """
    phenotypes_dict = {}
    names_dict = {}

    for ontology_file in ontology_files_list:

        with open(ontology_file, "r") as file:
            for line in file:

                if "id:" in line:
                    phen_id = line.split(" ")[1].strip()
                    phenotypes_dict[phen_id] = []
                    names_dict[phen_id] = ""  # add id of the main phenotype

                if "name:" in line:  # add name of the main phenotype
                    names_dict[phen_id] = line.split(":")[1].strip()

                if "is_a:" in line:
                    phenotypes_dict[phen_id].append(line.split(" ")[1])
                    names_dict[line.split(" ")[1]] = line.split("!")[1].strip()  # add id and name of related phenotype

    # pickle the dictionaries
    with open("ontology_data/related_phenotypes.pkl", "wb") as file:
        pickle.dump(phenotypes_dict, file)

    with open("ontology_data/phenotype_names.pkl", "wb") as file:
        pickle.dump(names_dict, file)


class DataExtraction:
    """Handle extraction and saving of data from `Python_APRIL_2021_phenotype_enrichment.py`.
    A new instance of the object is created upon running the `Python_APRIL_2021_phenotype_enrichment.py `
    for a given pathway. It saves the data. The saved file is name contains the `pathway_name`.
    This class is called from the `Python_APRIL_2021_phenotype_enrichment.py` script and all of the files are saved
    once the script finishes.
    """

    def __init__(self):
        self.orthologs_phenotype_dict = {}
        self.genes_orthologs_df = pd.DataFrame()
        self.go_annot_df = pd.DataFrame()
        self.enriched_phenotypes_set = set()

    def add_enrichment_phenotypes_set(self, enrichment_df: pd.DataFrame) -> None:
        for i in range(0, len(enrichment_df)):
            self.enriched_phenotypes_set.add(enrichment_df.iloc[i, 0])

        for gene_id in self.orthologs_phenotype_dict.keys():
            enriched_phenotypes_list = []
            for phen_id in self.orthologs_phenotype_dict.get(gene_id):
                if phen_id in self.enriched_phenotypes_set:
                    enriched_phenotypes_list.append(phen_id)

            self.orthologs_phenotype_dict[gene_id] = enriched_phenotypes_list
            print(f"gene id eriched list: {gene_id}; enriched list: {enriched_phenotypes_list}")
            print(self.orthologs_phenotype_dict)

    def add_ortholog_vs_phenotype_data(self, genes_phen_df: pd.DataFrame) -> None:
        """Add df containing gene IDs of orthologs and the enriched phenotypes.

        :param genes_phen_df: df containing gene IDs of orthologs and the enriched phenotypes
        :param organims: the name of the organism
        """
        dict = {}
        if genes_phen_df is not None:
            for i in range(0, len(genes_phen_df)):
                gene_id = genes_phen_df.iloc[i, 0]
                phenotype = genes_phen_df.iloc[i, 1]
                if gene_id not in dict.keys():
                    dict[gene_id] = [phenotype]
                else:
                    dict[gene_id].append(phenotype)

            self.orthologs_phenotype_dict = dict
            print(f"gene_ids unfiltered: {self.orthologs_phenotype_dict.keys()}")

    def add_genes_vs_orthologs_data(self, df_ortholog: pd.DataFrame, organism: str) -> None:
        """Add df containing gene IDs of orthologs and the enriched phenotypes.

        :param df_ortholog: df containing human Ensembl gene IDs and the gene IDs of orthologs from different organisms
        :param organims: the name of the organism
        """
        if df_ortholog is not None and organism is not None:
            df_ortholog["Organism"] = [organism for i in range(0, len(df_ortholog))]
            data_frames = [df_ortholog, self.genes_orthologs_df]
            self.genes_orthologs_df = pd.concat(data_frames)

    def save_data_to_files(self) -> None:
        """Save extracted data from two DataFrames as csv files.
        The file will be saved in the pathway enrichment folder e.g.: `AHR_R-HSA-8937144_Enrichment_Results`.
        """
        self.genes_orthologs_df.to_csv("genes_orthologs_data.csv", index=False)

        with open("orthologs_to_phenotype_data.pkl", "wb") as file:
            pickle.dump(self.orthologs_phenotype_dict, file)
            print(f"final dict: {self.orthologs_phenotype_dict}")

    def add_related_phenotypes(self) -> None:  # maybe this when df is created
        pass

    def make_phenotypes_dict(self) -> None:
        pass


def get_combined_df(path_to_pathway_enrichment: str) -> pd.DataFrame:
    genes_orthologs_df = pd.read_csv(path_to_pathway_enrichment + "/genes_orthologs_data.csv")
    with open(path_to_pathway_enrichment + "/orthologs_to_phenotype_data.pkl", "rb") as file:
        orthologs_phenotype_dict = pickle.load(file)

    associated_phenotypes = []
    for i in range(0, len(genes_orthologs_df)):
        associated_phenotypes.append(
            orthologs_phenotype_dict.get(genes_orthologs_df.iloc[i, 0]))

    genes_orthologs_df["associated_phenotype"] = associated_phenotypes
    return genes_orthologs_df


def get_related_phenotypes(phenotype_id: str) -> List[str]:
    with open("ontology_data/related_phenotypes.pkl", "rb") as file:
        related_phenotypes = pickle.load(file)

    return related_phenotypes.get(phenotype_id)


def get_phenotype_name(phenotype_id: str) -> str:
    with open("ontology_data/phenotype_names.pkl", "rb") as file:
        phenotype_name = pickle.load(file)

    return phenotype_name.get(phenotype_id)


print(get_combined_df(
    "/home/ivana/DARTpaths-vis/Python_Phenotype_Enrichment-version_April_2021/AHR_R-HSA-8937144_Enrichment_Results"))

print(get_phenotype_name("WBPhenotype:0000625"))
print(get_related_phenotypes("WBPhenotype:0000625"))
print([get_phenotype_name(phenotype) for phenotype in get_related_phenotypes("WBPhenotype:0000625")])

# TODO
""""
Check if the combining all phenotypes would not ovrercrowd the graph;
precompute the phenotypes whole network

Multiple ways to to go down the hierarchy -> how to go to the lower when there are multiple options

Phenotype ontology is hierarchical -> connections of phenotypes (extract from phenotype ontology)

.obo file -> hierarchy of phenotypes (draw a links between the terms)
"""
