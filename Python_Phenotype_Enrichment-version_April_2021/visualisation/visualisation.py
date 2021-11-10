"""
;===================================================================================================
; Title:   Visualization of enriched phenotypes
; Authors: Alexandra, Ivana, Irene, Thao
;===================================================================================================

"""
# TODO get rid of Nan nodes
# TODO specify allowed locations of the nodes (eg genes in the middle of page)

# Import modules
import pandas as pd
import dash
import dash_cytoscape as cyto
from dash import dcc
from dash import html
import dash_bootstrap_components as dbc
from coordinates import check_coordinates, gene_y_coordinate


# Load data
df = pd.read_csv('data.csv', delimiter=',')
df_1 = df.assign(phenotype=df['associated_phenotype'].astype(str).str.split(',')).explode('phenotype')
gene_df = df_1[['1', '2', 'Organism', 'phenotype']]


# make nodes and edges
nodes = []
edges = []
nodes_set = set()  # to avoid duplication

# get total number of genes
N_genes = set(df['2'].tolist())
N = len(N_genes)  # total number of genes TODO maybe do this as set

# for node positions
count_n = 1  # count for x_pos of gene nodes

coordinates_ort = []
coordinates_phenotype = []

# Load data into nodes and edges
for index, row in gene_df.iterrows():
    # source node is human gene
    gene = str(row['2'])
    # target node is ortholog
    # TODO get rid of Nans
    ortholog = str(row['1'])
    phenotype = str(row['phenotype'])
    organism = str(row['Organism'])

    # add data and class info to the nodes
    # gene node positions
    pos_x = 50  # always the same
    pos_y = gene_y_coordinate(1, N, count_n)

    cy_gene = {'data': {'id': gene, 'label': gene, 'size': 4, 'fontsize': '1.5px'}, 'classes': 'blue', 'position': {'x': pos_x, 'y': pos_y}}

    # ortholog node
    pos_x_ortholog, pos_y_ortholog = check_coordinates(coordinates_ort, pos_y, 20, 45, 55, 80)

    cy_ortholog = {'data': {'id': ortholog, 'label': ortholog, 'size': 2, 'fontsize': '1px', 'organism': organism},
                   'position': {'x': pos_x_ortholog, 'y': pos_y_ortholog}}

    # phenotype node
    pos_x_phenotype, pos_y_phenotype = check_coordinates(coordinates_phenotype, pos_y, 0, 20, 80, 100)
    # TODO add label to these nodes
    cy_phenotype = {'data': {'id': phenotype, 'size': 1, 'fontsize': '0.5px'}, 'classes': 'purple',
                    'position': {'x': pos_x_phenotype, 'y': pos_y_phenotype}}

    cy_edge = {'data': {'id': gene + ortholog, 'source': gene, 'target': ortholog, 'width': '0.25', 'color':'#696969'}}
    cy_edge_2 = {'data': {'id': ortholog + phenotype, 'source': ortholog, 'target': phenotype, 'width': '0.10', 'color':'#B8B8B8'}}

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
            'line-color': 'data(color)',
        }
    },
    {
        'selector': '.purple',
        'style': {
            'background-color': 'purple',
            'line-color': 'purple'
        }
    },
    # colours based on organism
    {
        'selector': '[organism *= "dmelanogaster"]',
        'style': {
            'background-color': 'yellow',
            'line-colour': 'yellow',
        }
    },
    {
        'selector': '[organism *= "mouse"]',
        'style': {
            'background-color': 'green',
            'line-colour': 'green',
        }
    },
    {
        'selector': '[organism *= "celegans"]',
        'style': {
            'background-color': 'pink',
            'line-colour': 'pink',
        }
    },
    {
        'selector': '[organism *= "zebrafish"]',
        'style': {
            'background-color': 'orange',
            'line-colour': 'orange',
        }
    },
    {
        'selector': '[organism *= "slimemould"]',
        'style': {
            'background-color': 'brown',
            'line-colour': 'brown',
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

# graph itself
node_graph = dbc.Row([
    dbc.Col(
        html.Div(children=[
            cyto.Cytoscape(
                id='cytoscape-phenotype',
                layout={'name': 'preset'},
                elements=edges + nodes,
                stylesheet=graph_stylesheet,
                style={'width': '100%', 'height': '95vh'},

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
