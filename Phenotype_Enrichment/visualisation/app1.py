
import pandas as pd
import numpy as np
import plotly.express as px
import dash
import dash_cytoscape as cyto
import dash_bootstrap_components as dbc
from dash import html
from dash import dcc
from dash.dependencies import Input, Output
import json
from time import perf_counter
import coordinates
import db_retrieve


#--------------------------------------------
def expand_dataframe(df: pd.DataFrame):
    """Function to expand dataframe, so each phenotype is in a separate row.
    Returns new dataframe and the total number of genes

    :param df: dataframe to expand
    """
    df_1 = df.assign(phenotype=df['Enriched_Phenotypes'].astype(str).str.split(',')).explode('Enriched_Phenotypes')
    df_2 = df_1[['Orthlog_Genes', 'Human_Gene', 'Organism', 'Enriched_Phenotypes']]
    # get total number of genes
    n_genes = set(df['Human_Gene'].tolist())
    n = len(n_genes)  # total number of genes

    return df_2, n
#--------------------------------------------
def load_info_to_graph(dataframe: pd.DataFrame, metadata: pd.DataFrame, n_genes: int):
    """ Function used to extract information from dataframe and load it into
    nodes and edges list of dictionaries. For phenotype nodes, metadata is added as well

    :param dataframe: dataframe with information about nodes and edges
    :param metadata: dataframe with metadata of given pathway
    :param n_genes: number of unique genes, necessary for y position determination
    """

    # make nodes and edges
    nodes_list = []
    edges_list = []
    nodes_set = set()  # to avoid duplication

    # for node positions
    count_n = 1  # count for x_pos of gene nodes
    coordinates_ort = []
    coordinates_phenotype = []

    # Load data into nodes and edges
    for index, row in dataframe.iterrows():
        # source node is human gene
        gene = str(row['Human_Gene'])
        # target node is ortholog
        ortholog = str(row['Orthlog_Genes'])
        phenotype = str(row['Enriched_Phenotypes'])
        organism = str(row['Organism'])
        # add data and class info to the nodes
        # gene node positions
        # pos_x = 50  # always the same; set to 500 in case of large network
        pos_x, pos_y = coordinates.gene_coordinate(1, n_genes, count_n)
        cy_gene = {'data': {'id': gene,
                            'indent': 'gene',   #add
                            'label': gene,
                            'size': 4,
                            'fontsize': '1.5px'},
                   'classes': 'blue',
                   'selectable': True, #Allow selected
                   'position': {'x': pos_x, 'y': pos_y}}

        # ortholog node
        pos_x_ortholog, pos_y_ortholog = coordinates.check_coordinates(coordinates_ort, pos_y, n_genes, 20, 45, 55, 80)
        if ortholog != "NaN":
            cy_ortholog = {'data': {'id': ortholog,
                                    'label': ortholog,
                                    'indent': 'ortholog', #add
                                    'size': 2,
                                    'fontsize': '1px',
                                    'organism': organism},
                           'classes': 'ortholog', #add
                           'selectable': True,    #Allow selected
                           'position': {'x': pos_x_ortholog, 'y': pos_y_ortholog}}
            cy_edge = {
                'data': {'id': gene + ortholog,
                         'source': gene,
                         'target': ortholog,
                         'width': '0.25',
                         'color': '#696969'}}

        # add genes and orthologs to set nodes and cy nodes
        if gene not in nodes_set:
            nodes_set.add(gene)
            nodes_list.append(cy_gene)
            count_n += 1

        if ortholog not in nodes_set:
            nodes_set.add(ortholog)
            nodes_list.append(cy_ortholog)
        edges_list.append(cy_edge)

        if phenotype != "None": # get rid of None node
            metadata_phenotype = metadata.get(phenotype)  # load metadata
            # phenotype node
            pos_x_phenotype, pos_y_phenotype = coordinates.check_coordinates(coordinates_phenotype, pos_y, n_genes, 0, 20, 80, 100)

            cy_phenotype = {'data': {'id': phenotype,
                                     'label': phenotype,
                                     'indent': 'phenotype',
                                     'metadata': metadata_phenotype,
                                     'size': 1,
                                     'fontsize': '0.5px'},
                            'classes': 'purple',
                            'selectable': True,  #Allow selected
                            'grabbable': False,   #Not allow grab for now
                            'position': {'x': pos_x_phenotype, 'y': pos_y_phenotype}}
            cy_edge_2 = {'data': {'id': ortholog + phenotype,
                                  'source': ortholog,
                                  'target': phenotype,
                                  'width': '0.10',
                                  'color': '#B8B8B8'}}

            if phenotype not in nodes_set:
                nodes_set.add(phenotype)
                nodes_list.append(cy_phenotype)

            edges_list.append(cy_edge_2)
    return nodes_list, edges_list
#--------------------------------------------
# load data -------------------------------
gene_df1 = db_retrieve.select_from_enrichment_results("AHR")               #gene_df1
gene_df_1, N_genes_1 = expand_dataframe(gene_df1)

gene_df2 = db_retrieve.select_from_enrichment_results("AmineOxidase")     #gene_df2
gene_df_2, N_genes_2 = expand_dataframe(gene_df2)
metadata2 = db_retrieve.select_from_metadata("AmineOxidase")              #metadata2
nodes, edges = load_info_to_graph(gene_df_2, metadata2, N_genes_2)
#element = edges + nodes

#organism = ['celegans', 'zebrafish', 'mouse', 'dmelanogaster']
organism = ['dmelanogaster','mouse','celegans','zebrafish']
# https://plotly.com/python/discrete-color/
col_swatch = [px.colors.qualitative.Safe[8], px.colors.qualitative.Set1[2],
              px.colors.qualitative.Set1[7],px.colors.qualitative.Set1[4]]

# calculate the right phenotype number
no_celegans = (gene_df_2.groupby(['Organism']).count())['Enriched_Phenotypes'][0]
no_dmelanogaster = (gene_df_2.groupby(['Organism']).count())['Enriched_Phenotypes'][1]
no_mouse = (gene_df_2.groupby(['Organism']).count())['Enriched_Phenotypes'][2]
no_zebrafish = (gene_df_2.groupby(['Organism']).count())['Enriched_Phenotypes'][3]
totalPhenotype = [no_celegans,no_dmelanogaster,no_mouse,no_zebrafish]
org_totalPhenotype = dict(zip(organism,totalPhenotype))
slider = ["humangene", "organism_ortholog", "organism_phenotype"]

# --------------------------------------------------------------------------------

graph_stylesheet = [
    # {
    #     'selector': organism[i],       # colours based on organism
    #     'style': {"background-color": col_swatch[i],
    #               "line-color": col_swatch[i]
    #               },
    # }
    # for i in range(len(organism))

# colours based on organism
    {
        'selector': '[organism *= "dmelanogaster"]',
        'style': {
            'background-color': col_swatch[0],
            'line-colour': col_swatch[0],
        }
    },
    {
        'selector': '[organism *= "mouse"]',
        'style': {
            'background-color': col_swatch[1],
            'line-colour': col_swatch[1],
        }
    },
    {
        'selector': '[organism *= "celegans"]',
        'style': {
            'background-color': col_swatch[2],
            'line-colour': col_swatch[2],
        }
    },
    {
        'selector': '[organism *= "zebrafish"]',
        'style': {
            'background-color': col_swatch[3],
            'line-colour': col_swatch[3],
        }
    }
]

graph_stylesheet += [
    {
        'selector': 'nodes',
        'style': {
            'content': 'data(label)',
            'width': 'data(size)',
            'height': 'data(size)',
            "curve-style": "bezier",
            'font-size': 'data(fontsize)',
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
            'background-color': '#C918F3',
            'line-color': '#C918F3'
        }
    },

    {
        'selector': '.blue',
        'style': {'background-color': '#18B6F3',
                  'line-color': '#18B6F3'
                  }
    }

]

# Define app
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
# server = app.server


# --------------------------------------------------------------------------------
navbar = dbc.NavbarSimple(
    children=[
        dbc.NavItem(
            dbc.NavLink(
                "Article",
                href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5829367/",
            )
        ),
        dbc.NavItem(
            dbc.NavLink(
                "Source Code",
                href="https://github.com/ivanajanickova/DARTpaths-vis",
            )
        ),
    ],
    brand="Dart - Translation of Genotype to Phenotype by a Hierarchy of Model Organism - simple demo",
    brand_href="#",
    color="dark",
    dark=True,
)

topics_html = list()
for topic_html in [
    html.Span([str(i+1) + ": " + organism[i]], style={"color": col_swatch[i]})
    for i in range(len(organism))
]:
    topics_html.append(topic_html)
    topics_html.append(html.Br())

body_layout = dbc.Container(
    [
        dbc.Row(
            [
                dbc.Col(
                    [
                        dcc.Markdown(
                            f"""
                -----
                ##### Description:
                -----
                For this demonstration, phenotypes from 4 model organism were used to visualize
                the relationship between phenotypes from different organism related to a particular human
                pathway. The pathways are preference from [REACTOME](https://reactome.org/PathwayBrowser/)
                Phenotype from each organism is shown in different color on the map, as shown on the right.
                """
                        )
                    ],
                    sm=12,
                    md=8,
                ),
                dbc.Col(
                    [
                        dcc.Markdown(
                            """
                -----
                ##### Model oganisms:
                -----
                """
                        ),
                        html.Div(
                            topics_html,
                            style={
                                "fontSize": 11,
                                "height": "100px",
                                "overflow": "auto",
                            },
                        ),
                    ],
                    sm=12,
                    md=4,
                )
            ]
        ),
        dbc.Row(
            [
                dcc.Markdown(
                    """
            -----
            ##### Node representation
            There are three kinds of node in the graph. The biggest nodes in the center indicate human genes in a chosen 
            pathway. These nodes are linked directly with the second type of node which indicate orthologs from a model 
            organism. Node color indicates which organism the orthologs belong to. 
            
            The third type of node in two sides is the phenotype enrichment of organism model
            
            The edges show the relationship of 1. phenotypes and its related ortholog gene in organism 2. ortholog genes in 
            organism and its human gene in the same molecular pathway in Reactome.

            Use the filters to choose the phenotype with:
            * certain numbers of pathway/ pathway level 
            * certain numbers of model organism you want to investigate
            
            The slider is used to explore 3 level of visualisation: human gene only, human gene and orthologs or the whole graph. Toggle to explore different visualisation options.
            
            -----
            """
                )
            ]
        ),
        dbc.Row(
            [
                dbc.Col(
                    [
                        dbc.Row(
                            [
                                cyto.Cytoscape(
                                    id='cytoscape-phenotype',
                                    layout={"name": "preset"},
                                    style={"width": "100%", "height": "800px"},
                                    elements= edges + nodes,   #change
                                    stylesheet=graph_stylesheet,
                                    #style={'width': '100%', 'height': '95vh'}
                                    minZoom = 0.6
                                )
                            ]
                        ),
                        dbc.Row(
                            [
                                dbc.Alert(
                                    id='phenotype-data',
                                    children='Click on a node to see its details here',
                                    color='secondary'
                                )
                            ]
                        )
                    ],
                    sm=12,
                    md=8
                ),
                dbc.Col(
                    [
                        dbc.Label('Search phenotype for a gene: ', html_for='gene_search'),
                        dbc.Input(type='password', id='gene_search', placeholder="Enter a gene ID here"),
                        dbc.FormText('The gene is not available in the pathway', color='secondary'),
                        html.Div(style={'padding': 20}), #Space

                        html.H4(dbc.Badge(
                            "Level of pathway(s):",
                            color="success",
                            className="ms-1",

                        )),
                        dcc.Dropdown(
                            id="pathway_level",
                            #options=[{"label": k, "value": k} for k in range(1,3)],
                            options=[{"label": 'lower level', "value": 1},{"label": 'upper level', "value": 2}],
                            clearable=False,
                            value= 2,
                            style={"width": "150px"},
                        ),
                        html.H4(dbc.Badge(
                            "Model organism(s):", color="success", className="mr-1"
                        )),
                        dcc.Dropdown(
                            id="organism_dropdown",
                            options=[{"label": i + " (" + str(v)+ " phenotype(s))", "value": i} for i, v in org_totalPhenotype.items()],
                            value= [{"label": i + " (" + str(v)+ " phenotype(s))", "value": i} for i, v in org_totalPhenotype.items()] ,    #check??
                            multi=True,
                            style={"width": "100%"},
                        ),
                        html.Div(style={'padding': 7}),
                        dbc.Button("Search", id="search-button", className="me-2", n_clicks=0, size="sm",color="secondary"),
                        html.Span(id="search-status", style={"verticalAlign": "middle"}),
                        html.Div(style={'padding': 20}),  #Space

                        html.H4(dbc.Badge(
                            "Choose level of visualisation:", color="success", className="mr-1"
                        )),
                        html.Div(style={'padding': 10}), #Space
                        dcc.Slider(
                            id='detail_slider',
                            min=1,
                            max=3,
                            value=3,
                            marks={1: "human_gene",2:"organism_ortholog", 3:"organism_phenotype"},
                            step=1
                        )
                    ]
                )
                ]
        )

    ]
)


app.layout = html.Div([navbar, body_layout])


# @app.callback(
#     Output('phenotype-data', 'children'),
#     [Input('cytoscape-phenotype', 'tapNodeData')]
# )
# def display_phenotype(data):
#     contents = 'Click on a node to see its details here'
#     if data:
#         contents = []
#         #contents.append(html.H5('Phenotype: ' + str(data.get('metadata')[0])))
#         contents.append(html.H5('Phenotype: ' + data.get('id')))
#         contents.append(html.P('Hierchically related phenotype(s): ' + str(data.get('metadata')[1])))
#         contents.append(html.P('Phenotype enrichment Q-value: ' + str(data.get('metadata')[3])))
#         contents.append(html.P('Phenotype enrichment P-value: ' + str(data.get('metadata')[4])))
#         return contents
#     else:
#         return contents

# call back Node information -----------------
@app.callback(
    Output('phenotype-data', 'children'),
    [Input('cytoscape-phenotype', 'tapNodeData')]
)
def display_nodeData(data):
    contents = 'Click on a node to see its details here'
    if data:
        contents = []
        if data.get('indent') == 'phenotype':
            contents.append(html.H5('Phenotype: ' + data.get('label')))
            contents.append(html.P('Hierchically related phenotype(s): ' + str(data.get('metadata')[1])))
            contents.append(html.P('Phenotype enrichment Q-value: ' + str(data.get('metadata')[3])))
            contents.append(html.P('Phenotype enrichment P-value: ' + str(data.get('metadata')[4])))
            return contents
        elif data.get('indent') == 'gene':
            contents.append(html.H5('Gene ID: ' + data.get('label')))
            return contents
        elif data.get('indent') == 'ortholog':
            contents.append(html.H5('Ortholog ID: ' + data.get('label')))
            contents.append(html.P('Organism: ' + data.get('organism')))
            return contents
    else:
        return contents
#-----------------------------------------------

# call back Slider + Dropdown -----------------
# @app.callback(
#     Output('cytoscape-phenotype', "elements"),
#     [
#         Input('pathway_level', 'value'),
#         Input('organism_dropdown', 'value'),
#         Input('detail_slider', 'value'),
#     ]
# )
# def filter_nodes(pathway_level, organism_list, detail_slider):
#     ctx = dash.callback_context
#     if not ctx.triggered:
#         return element
#
#     else:
#         trigger_id = ctx.triggered[0]['prop_id'].split(".")[0]
#         val = ctx.triggered[0]['value'].split(".")[0]
#
#         if trigger_id == 'pathway_level':
#             pathway = val
#             df = gene_df1 if val == 1 else gene_df2
#         elif trigger_id == 'organism_dropdown':

@app.callback(
    Output('cytoscape-phenotype', "elements"),
    [Input('pathway_level', 'value')]
)
def filter_nodes(value):
    if value == 1:
        metadata1 = db_retrieve.select_from_metadata("AHR")
        new_nodes, new_edges = load_info_to_graph(gene_df_1, metadata1, N_genes_1)
        new_element = new_edges + new_nodes
        return new_element

    return element
#----------------------------------------------
# @app.callback(
#     Output('cytoscape-phenotype', 'elements'),
#     [Input('pathway_level', 'value'),
#      Input('organism_dropdown', 'value')]
# )
# def filter_nodes(pathway_level):
#     if pathway_level == 1:
#         metadata1 = db_retrieve.select_from_metadata("AHR")
#         new_nodes, new_edges = load_info_to_graph(gene_df_1, metadata1, N_genes_1)
#         new_element = new_edges + new_nodes
#         return new_element
#
#     return element
#---------------------------------------------
# @app.callback(
#     Output('cytoscape-phenotype', 'elements'),
#     [Input('detail_slider', 'value')]
# )
# def level_visualisation(slider_value):
#     if slider_value == int(1):
#
#         metadata1 = db_retrieve.select_from_metadata("AHR")
#         new_nodes, new_edges = load_info_to_graph(gene_df_1, metadata1, N_genes_1)
#         new_nodes, new_edges = load_info_to_graph(gene_df_1, metadata1, N_genes_1)
#         new_element = new_edges + new_nodes
#         return new_element
#
#     return element
#----------------------------------------------
if __name__ == '__main__':
    app.run_server(debug=False)  # set false so you can load bigger networks
