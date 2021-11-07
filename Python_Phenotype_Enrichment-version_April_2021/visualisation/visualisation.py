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

import pandas as pd
import dash
import dash_cytoscape as cyto
from dash import dcc
from dash import html
import dash_bootstrap_components as dbc

# Load data
gene_df = pd.read_csv('data.csv', delimiter=',')

# make nodes and edges
nodes = []
edges = []
nodes_set = set()  # to avoid duplication

# Load data into nodes and edges
for index, row in gene_df.iterrows():
    # source node is human gene
    gene = str(row['2'])
    # target node is ortholog
    ortholog = str(row['1'])
    # TODO add phenotype nodes

    phenotype = row['associated_phenotype']
    organism = row['Organism']

    cy_gene = {'data': {'id': gene, 'label': gene}, 'classes': 'blue'}
    cy_ortholog = {'data': {'id': ortholog, 'label': ortholog, 'organism': organism}, 'classes': 'red'}
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

#########################
#        Graph          #
#########################

graph_stylesheet = [
    {
        'selector': 'nodes',
        'style': {
            'content': 'data(label)'
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
        'selector': '.blue',
        'style': {
            'background-color': 'blue',
            'line-color': 'blue'
        }
    }

]

node_graph = dbc.Row([
    dbc.Col(
        html.Div(children=[
            cyto.Cytoscape(
                id='cytoscape-phenotype',
                elements=edges + nodes,
                stylesheet=graph_stylesheet,
                style={
                    'width': '100%',  # Take up 100% of the width of the space it has been assigned
                    'height': '87vh'  # 87vh = 87% of the total screen height
                },
                responsive=True  # Changes size cystocape graph is browser window changes size
            )
        ])
    )
]
)

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
