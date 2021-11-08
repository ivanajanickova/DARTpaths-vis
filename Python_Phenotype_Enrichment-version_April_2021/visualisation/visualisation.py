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
gene_df = pd.read_csv('data.csv', delimiter=',')

# make nodes and edges
nodes = []
edges = []
nodes_set = set()  # to avoid duplication

N = len(gene_df.index)  # total number of genes TODO maybe do this as set

# for node positions
count_n = 1  # count for x_pos of gene nodes
x_all = []
y_all = []

# Load data into nodes and edges
for index, row in gene_df.iterrows():
    # source node is human gene
    gene = str(row['2'])
    # target node is ortholog
    # TODO get rid of Nans
    ortholog_all = row['1']
    if not ortholog_all == ["FBgn0003513"]:
        ortholog = str(ortholog_all)
    # TODO add phenotype nodes

    phenotype = row['associated_phenotype']
    organism = row['Organism']

    # add data and class info to the nodes
    # gene node
    # position
    pos_x = 50
    pos_y = 1  # intialise outside if/else block
    if count_n <= N:
        pos_y = 1000 - numpy.log(count_n / N) * 50  # normalise count_n/N, otherwise logarithmical
    else:
        print("something went wrong")

    cy_gene = {'data': {'id': gene, 'label': gene}, 'classes': 'blue', 'position': {'x': pos_x, 'y': pos_y}}

    # ortholog node
    # randomize left or right side position of node

    def check_x_position(position_list):
        """ Creates x coordinate (position) for a node using randomisation.
        Checks if created x coordinate if present in provided position list.
        If it is, function is recursively called back.
        Finally, it appends x coordinate to the list
        :param position_list: list storing all x coordintaes (positions) of nodes
         """
        if bool(random.getrandbits(1)) is True:
            pos_x_ort = random.randint(0, 45)
        else:
            pos_x_ort = random.randint(55, 100)
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
        position = position_float + random.randint(0, 2)
        if position not in position_list:
            position_list.append(position)
        else:
            check_y_position(position_list, position_float)
        return position

    pos_x_ortholog = check_x_position(x_all)
    pos_y_ortholog = check_y_position(y_all, pos_y)  # call the function

    cy_ortholog = {'data': {'id': ortholog, 'label': ortholog, 'organism': organism}, 'classes': 'red',
                   'position': {'x': pos_x_ortholog, 'y': pos_y_ortholog}}
    cy_edge = {'data': {'id': gene + ortholog, 'source': gene, 'target': ortholog}}

    # add genes and orthologs to set nodes and cy nodes
    if gene not in nodes_set:
        nodes_set.add(gene)
        nodes.append(cy_gene)

    if ortholog not in nodes_set:
        # TODO not entirely sure this will work, since multiple genes can have single ortholog
        nodes_set.add(ortholog)
        nodes.append(cy_ortholog)

    edges.append(cy_edge)

    count_n += 1


#########################
#        Graph          #
#########################

# style for the graph: node colours, names
graph_stylesheet = [
    {
        'selector': 'nodes',
        'style': {
            'content': 'data(label)',
            'width': '4',
            'height': '4',
            'font-size': '2px',
            # 'colour': 'data(organism)'
        }
    },
    {
        'selector': 'edges',
        'style': {
            'width': '0.25',
            'line-color': 'grey'
        }
    },
    # {
    #     'selector': '.red',
    #     'style': {
    #         'background-color': 'red',
    #         'line-color': 'red'
    #     }
    # },
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
