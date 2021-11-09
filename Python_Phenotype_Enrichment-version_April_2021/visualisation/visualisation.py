"""
;===================================================================================================
; Title:   Visualization of enriched phenotypes
; Authors: Alexandra, Ivana, Irene, Thao
;===================================================================================================

"""
# TODO get rid of Nan nodes
# TODO import phenotypes as a different type of nodes, add it to nodes set+ cy nodes
# TODO specify allowed locations of the nodes (eg genes in the middle of page)
# TODO play with the layouts (I don't like grids)

# Import modules
import math

import numpy
import pandas as pd
import dash
import dash_cytoscape as cyto
from dash import dcc
from dash import html
import dash_bootstrap_components as dbc
import random

# Load data
df = pd.read_csv('data.csv', delimiter=',')
df_1 = df.assign(phenotype=df['associated_phenotype'].astype(str).str.split(',')).explode('phenotype')
gene_df = df_1[['1', '2', 'Organism', 'associated_phenotype']]

# make nodes and edges
nodes = []
edges = []
nodes_set = set()  # to avoid duplication

N = len(df.index)  # total number of genes TODO maybe do this as set

# for node positions
count_n = 1  # count for x_pos of gene nodes
x_all_orthologs = []  # all ortholog x positions
y_all_orthologs = []  # all ortholog y positions
x_all_phenotypes = []
y_all_phenotypes = []

# Load data into nodes and edges
for index, row in gene_df.iterrows():
    # source node is human gene
    gene = str(row['2'])
    # target node is ortholog
    # TODO get rid of Nans
    ortholog = str(row['1'])

    phenotype = str(row['associated_phenotype'])

    organism = row['Organism']

    # add data and class info to the nodes
    # gene node positions
    pos_x = 50  # always the same

    def gene_y_position(y, total_genes):
        """ Creates y coordinate of a gene node
         Y coordinate must be placed lower on the y axis as a previous one.
         As long as the count is lower than the total number of genes, the y position is produced.

         :param y: y coordinate/position of a previous node
         :param total_genes: the total number of genes to be visualised as nodes
         """
        if count_n <= total_genes:
            y = 1000 - math.log(count_n / total_genes) * 25  # normalise count_n/N, otherwise logarithmical
        else:
            print("something went wrong")
        return y

    pos_y = gene_y_position(1, N)

    cy_gene = {'data': {'id': gene, 'label': gene, 'size': 4, 'fontsize': '1.5px'}, 'classes': 'blue', 'position': {'x': pos_x, 'y': pos_y}}

    # ortholog node
    # randomize left or right side position of node

    def check_x_position(position_list, x1, x2, x3, x4):
        """ Creates x coordinate (position) for a node using randomisation.
        Checks if created x coordinate if present in provided position list.
        If it is, function is recursively called back.
        Finally, it appends x coordinate to the list

        :param position_list: list storing all x coordinates (positions) of nodes
         """
        if bool(random.getrandbits(1)) is True:
            pos_x_ort = random.uniform(x1, x2)
        else:
            pos_x_ort = random.uniform(x3, x4)

        if pos_x_ort in position_list:
            check_x_position(position_list)
        else:
            position_list.append(pos_x_ort)
            return pos_x_ort


    def check_y_position(position_list, position_float):
        """ Checks if node y coordinate (position) is already present in list of positions
        If it is not true, the position is returned and added to the list of positions
        otherwise, the function is called upon itself.

        :param position_list: list where the non-overlapping positions are stored
        :param position_float: the y coordinate (position) of related gene node
        """
        position = position_float + random.uniform(0, 2)
        if position not in position_list:
            position_list.append(position)
        else:
            check_y_position(position_list, position_float)
        return position


    pos_x_ortholog = check_x_position(x_all_orthologs, 20, 45, 55, 80)
    pos_y_ortholog = check_y_position(y_all_orthologs, pos_y)  # call the function

    cy_ortholog = {'data': {'id': ortholog, 'label': ortholog, 'size': 2, 'fontsize': '1px', 'organism': organism}, 'classes': 'red',
                   'position': {'x': pos_x_ortholog, 'y': pos_y_ortholog}}

    pos_x_phenotype = check_x_position(x_all_phenotypes, 0, 20, 80, 100)
    pos_y_phenotype = check_y_position(y_all_phenotypes, pos_y)  # call the function
    # TODO add label to these nodes
    cy_phenotype = {'data': {'id': phenotype, 'size': 1, 'fontsize': '0.5px'}, 'classes': 'purple',
                    'position': {'x': pos_x_phenotype, 'y': pos_y_phenotype}}

    cy_edge = {'data': {'id': gene + ortholog, 'source': gene, 'target': ortholog, 'width': '0.25'}}
    cy_edge_2 = {'data': {'id': ortholog + phenotype, 'source': ortholog, 'target': phenotype, 'width': '0.10'}}

    # add genes and orthologs to set nodes and cy nodes
    if gene not in nodes_set:
        nodes_set.add(gene)
        nodes.append(cy_gene)
        count_n += 1

    if ortholog not in nodes_set:
        # TODO not entirely sure this will work, since multiple genes can have single ortholog
        # it should, since this does not affect the edges
        nodes_set.add(ortholog)
        nodes.append(cy_ortholog)

    if phenotype not in nodes_set:
        nodes_set.add(phenotype)
        nodes.append(cy_phenotype)

    edges.append(cy_edge)
    edges.append(cy_edge_2)
    # TODO figure out why not all edges show
#########################
#        Graph          #
#########################

# style for the graph: node colours, names
graph_stylesheet = [
    {
        'selector': 'nodes',
        'style': {
            'content': 'data(label)',
            'width': 'data(size)',
            'height': 'data(size)',
            'font-size': 'data(fontsize)',
            # 'colour': 'data(organism)'
        }
    },
    {
        'selector': 'edges',
        'style': {
            'width': 'data(width)',
            'line-color': 'grey'
        }
    },
    {
        'selector': '.red',
        'style': {
            'background-color': 'red',
            'line-color': 'red'
        }
    },
    {
        'selector': '.purple',
        'style': {
            'background-color': 'purple',
            'line-color': 'purple'
        }
    },
    # TODO work on this
    #  {
    #     'selector': '["organism" *= "dmelanogaster"]',
    #     'background-color': '#FF4136',
    #     'line-colour': '#FF4136',
    # },
    {
        'selector': '.blue',
        'style': {
            'background-color': 'blue',
            'line-color': 'blue'
        }
    }

]

# graph itself
node_graph = dbc.Row([
    dbc.Col(
        html.Div(children=[
            cyto.Cytoscape(
                id='cytoscape-phenotype',
                layout={'name': 'preset'},
                elements=edges + nodes,
                stylesheet=graph_stylesheet,
                style={'width': '90%', 'height': '95vh'},

                responsive=True  # Changes size cystocape graph is browser window changes size
            )
        ])
    )
]
)

# specify window
# TODO add windows for metadata and all the funky stuff
window = html.Div(
    id='window',
    children=[
        node_graph
    ]
)

# run app
app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])
app.layout = html.Div([dcc.Location(id="url"), window])
if __name__ == '__main__':
    app.run_server(debug=True)
