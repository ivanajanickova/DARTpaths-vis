import preprocessing
import os

# this file is just for demonstration you can delete it once you understand the output
# FIRSTLY RUN THE ENRICHMENT SCRIPT FOR THE PATHWAY !!!!

# get dataframe for a given pathway
dataframe = preprocessing.get_combined_df(str(os.getcwd()) + "/AHR_R-HSA-8937144_Enrichment_Results")
print(dataframe)

# get the name of the phenotype id
name_1 = preprocessing.get_phenotype_name("MP:0005377")
print(name_1)
# get the list of related phenotypes 
list_of_related_phenotypes_1 = preprocessing.get_related_phenotypes("MP:0005377")
print(list_of_related_phenotypes_1)
print([preprocessing.get_phenotype_name(phen) for phen in preprocessing.get_related_phenotypes("MP:0005377")])
# get the name of the phenotype id
name_2 = preprocessing.get_phenotype_name("WBPhenotype:0000625")
print(name_2)
# get the list of related phenotypes
list_of_related_phenotypes_2 = preprocessing.get_related_phenotypes("WBPhenotype:0000625")
print(list_of_related_phenotypes_2)
print([preprocessing.get_phenotype_name(phen) for phen in preprocessing.get_related_phenotypes("WBPhenotype:0000625")])

