from visualisation import db_retrieve

#######################################################################################################################
# Interface for retrieving the data from DB ###########################################################################
#######################################################################################################################

# THIS SCRIPT IS JUST FOR DEMONSTRATION

##########################
# Retrieve from DB #######
##########################
# Retrieve a dataframe of genes orthologs and enriched phenotypes for a given pathway
# The list of pathways we have so far:
# AmineOxidase
# Phase2ConjugationOfCompounds
# Phase1CompoundFunctionalization
# AminoAcidConjugation
# EthanolOxidation
# AHR

df = db_retrieve.select_from_enrichment_results("AHR")
print(df)

# Retrieve phenotype metadata
metadata = db_retrieve.select_from_metadata("AHR")
print(metadata)

# Retrieve the name of the higher level pathway
name = db_retrieve.find_top_level_pathway("AHR")
name = name[0][0]
print(name)  # print pure string

# You can then analogously retrieve info for the higher level data
df_2 = db_retrieve.select_from_enrichment_results(name)
print(df_2)

############################
# Get Phenotype Metadata ###
############################
# retrieve data for a phenotype - If the phenotype data are 'nice' - such as for c. elegans
# 'Nice' means that the database in the preprocessing contained name of phenotype & list of related phenotypes

random_phenotype_list = df.iloc[40, 3]
for random_phenotype in random_phenotype_list:
    print(random_phenotype)
    print(metadata.get(random_phenotype))

# retrieve data retrieve data for a phenotype - If the phenotype data are 'nice' - such as for d. meganogaster
random_phenotype_list = df.iloc[3, 3]
for random_phenotype in random_phenotype_list:
    print(random_phenotype)
    print(metadata.get(random_phenotype))
