### This tool extracts D.Melanogaster orthologs from Human genes
### by searching each Notung file, then placing the orthologs in a dataframe and finally by appending that that to a TSV file.
### If a TSV file already exists, it will be eliminated when running this script.

#Library installation
import pandas as pd
import os
from pathlib import Path

#States the path where the script is stored
path = os.path.abspath(os.path.dirname("__file__"))
print(path)
#Function for crawling the folders, loading the Notung file data to a dataframe, filter its content to only store the intended information (Drosophila and Human),
#tranform the data to the intended form for output, and finally send the dataframe to the TSV file. This process occurs for each Notung file and the outputs are appended
#to the same TSV file.

my_file = Path("Notung_DMelanogaster_orthologs.txt")
if my_file.exists():
    os.remove("Notung_DMelanogaster_orthologs.txt")
for root, dirs, files in os.walk(path):
    for name in files:
        #print(files)
        if name.endswith("ReadyForNotung.nw.gtpruned.rooting.0.rearrange.0.homologs.txt"):
            filename = os.path.abspath(os.path.join(root, name))
            print(filename)
            lines = open(filename).readlines()
            strToFind = "DrosophilaMelanosgaster"
            count = len(open(filename).readlines(  ))
            if count > 10:
                ortho_data = pd.read_csv(filename, delimiter='\t', sep="\t", header = 12, index_col=0)

                Dmelanogaster_cols = [col for col in ortho_data.columns if 'DrosophilaMelanogaster' in col]
                ortho_data = ortho_data[Dmelanogaster_cols]
                ortho_data = ortho_data[(ortho_data == 'O')]
                ortho_data = ortho_data.dropna(how="all")
                ortho_data = ortho_data.replace("O", pd.Series(ortho_data.columns, ortho_data.columns))
                ortho_data[Dmelanogaster_cols] = ortho_data[Dmelanogaster_cols].replace({'7227.':''}, regex=True)
                ortho_data[Dmelanogaster_cols] = ortho_data[Dmelanogaster_cols].replace({'__DrosophilaMelanogaster':''}, regex=True)
                ortho_data.columns = range(ortho_data.shape[1])
                ortho_data = ortho_data.transpose()

                hSapiens_cols = [col for col in ortho_data.columns if 'HomoS' in col]
                ortho_data = ortho_data[hSapiens_cols]
                ortho_data.columns = ortho_data.columns.str.replace(r"9606.", "")
                ortho_data.columns = ortho_data.columns.str.replace(r"__HomoSapiens", "")
                ortho_data = ortho_data.dropna(how="all")

                ortho_data = ortho_data.transpose()
                ortho_data
                print(ortho_data)
                ortho_data.to_csv("Notung_DMelanogaster_orthologs.txt", header=False, mode='a', sep = "\t")
