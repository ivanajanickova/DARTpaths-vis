import random

import pandas as pd

import json
from visualisation import load_dataframe
# # Opening JSON file
# f = open('pathway_data/metadataAHR.json')
#
# # returns JSON object as
# # a dictionary
# data = json.load(f)


def load_metadata(file_path: str):
    """" Function to load metadata from json file.
    Returns a dictionary

    :param file_path: the path to the json file we wish to load
    """
    file = open(file_path)
    metadata = json.load(file)
    file.close()
    return metadata


metadata = load_metadata('pathway_data/metadataAHR.json')

df_path = 'pathway_data/AHR.csv'
genedf, n = load_dataframe(df_path)


for i in metadata.keys():
    stri = str(i)  # change to string
    for index, row in genedf.iterrows():
        phenotype = str(row['phenotype'])
        if phenotype == stri:
            print("something")
        else:
            pass









# def metadata_extraction(pathway_metadata: dict):
#     for key, value in pathway_metadata.items():
#         phenotype_id = key
#         phenotype_info = value
#         # phenotype_name = str(value[0])
#         # related_phenotypes = value[1]
#         # pval = value[2]
#         # rank = value[3]
#         # qval = value[4]
#     return phenotype_id, phenotype_info   # phenotype_name, related_phenotypes, pval, rank, qval
#
#
# ids, values = metadata_extraction(data)
# print(values[0])
# # Iterating through the json
# list
# for i in data:
#     print(i)















# df = pd.read_csv('data.csv', delimiter=',')
# df_1 = df.assign(phenotype=df['associated_phenotype'].astype(str).str.split(',')).explode('phenotype')
# df_full = df_1[['1', '2', 'Organism', 'phenotype']]
# df_full['phenotype'] = df_full['phenotype'].map(lambda x: str(x)[:-1])
# df_full['phenotype'] = df_full['phenotype'].map(lambda x: str(x)[2:])
# #gene_df['associated_phenotype'] = gene_df.associated_phenotype.apply(lambda x: x[1:-1].split('"'))
# print(df_full)
# for index, row in df_full.iterrows():
#     phenotype = row['phenotype']
    #print(phenotype)
#print(gene_df.applymap(type))
#gene_df.explode('associated_phenotype')

# def check_x_position(position_list, x1, x2, x3, x4):
#     """ Creates x coordinate (position) for a node using randomisation.
#     Checks if created x coordinate if present in provided position list.
#     If it is, function is recursively called back.
#     Finally, it appends x coordinate to the list
#
#     :param position_list: list storing all x coordinates (positions) of nodes
#      """
#     if bool(random.getrandbits(1)) is True:
#         pos_x_ort = random.uniform(x1, x2)
#     else:
#         pos_x_ort = random.uniform(x3, x4)
#
#     if pos_x_ort in position_list:
#         check_x_position(position_list)
#     else:
#         position_list.append(pos_x_ort)
#         return pos_x_ort
#
#
# def check_y_position(position_list, position_float):
#     """ Checks if node y coordinate (position) is already present in list of positions
#     If it is not true, the position is returned and added to the list of positions
#     otherwise, the function is called upon itself.
#
#     :param position_list: list where the non-overlapping positions are stored
#     :param position_float: the y coordinate (position) of related gene node
#     """
#     # since during the first iterration through df, position_float=1 and we don't wnt all to be too close
#     if position_float == 1:
#         position = random.uniform(0, 100)
#     else:
#         position = position_float + random.uniform(0, 20)
#
#     if position not in position_list:
#         position_list.append(position)
#     else:
#         check_y_position(position_list, position_float)
#     return position

coordinates_ortologs = []


# def check_coordinates(coordinates, coordinate_float,x1,x2,x3,x4):
#     """ Checks if node y coordinate (position) is already present in list of positions
#     If it is not true, the position is returned and added to the list of positions
#     otherwise, the function is called upon itself.
#
#     :param coordinates: list where the non-overlapping positions are stored
#     :param coordinate_float: the y coordinate (position) of related gene node
#     """
#
#     # x coordinate - left or right side of gene nodes:
#     if bool(random.getrandbits(1)) is True:
#         x_coordinate = random.uniform(x1, x2)
#     else:
#         x_coordinate = random.uniform(x3, x4)
#
#     # since during the first iterration through df, position_float=1 and we don't wnt all to be too close
#     # y coordinate
#     if coordinate_float == 1:
#         y_coordinate = random.uniform(0, 100)
#     else:
#         y_coordinate = coordinate_float + random.uniform(0, 20)
#
#     coordinate = (x_coordinate, y_coordinate)
#
#     # chek if given coordinates exist already
#     if coordinate not in coordinates:
#         coordinates.append(coordinate)
#     else:
#         check_coordinates(coordinates, coordinate_float, x1,x2,x3,x4)  # recursioon
#
#     return x_coordinate, y_coordinate

# df = pd.read_csv('data.csv', delimiter=',')
# df_1 = df.assign(phenotype=df['associated_phenotype'].astype(str).str.split(',')).explode('phenotype')
# gene_df = df_1[['1', '2', 'Organism', 'associated_phenotype']]
# print(gene_df)



# df = pd.read_csv('data.csv', delimiter=',')
# df_1 = df.assign(phenotype=df['associated_phenotype'].astype(str).str.split(',')).explode('phenotype')
# df_full = df_1[['1', '2', 'Organism', 'phenotype']]
# df_full.to_csv('output.csv')