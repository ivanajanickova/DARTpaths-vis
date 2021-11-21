
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
                                    'size': 2,
                                    'fontsize': '1px',
                                    'organism': organism},
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
gene_df = db_retrieve.select_from_enrichment_results("AHR")               #gene_df
gene_df, N_genes = expand_dataframe(gene_df)
gene_df2 = db_retrieve.select_from_enrichment_results("AmineOxidase")     #gene_df2
gene_df_2, N_genes_2 = expand_dataframe(gene_df2)
metadata2 = db_retrieve.select_from_metadata("AmineOxidase")              #metadata2
nodes, edges = load_info_to_graph(gene_df_2, metadata2, N_genes_2)


organism = ['celegans', 'zebrafish', 'mouse', 'dmelanogaster']
# calculate the right phenotype number
no_celegans = (gene_df_2.groupby(['Organism']).count())['Enriched_Phenotypes'][0]
no_dmelanogaster = (gene_df_2.groupby(['Organism']).count())['Enriched_Phenotypes'][1]
no_mouse = (gene_df_2.groupby(['Organism']).count())['Enriched_Phenotypes'][2]
no_zebrafish = (gene_df_2.groupby(['Organism']).count())['Enriched_Phenotypes'][3]
totalPhenotype = [no_celegans,no_dmelanogaster,no_mouse,no_zebrafish]
org_totalPhenotype = dict(zip(organism,totalPhenotype))
slider = ["humangene", "organism_ortholog", "organism_phenotype"]
#---------------------------------------------------------------------------------
col_swatch = px.colors.qualitative.Dark24
# --------------------------------------------------------------------------------

graph_stylesheet = [# colours based on organism
    # {
    #     "selector": "." + organism[i],
    #     "style": {"background-color": col_swatch[i], "line-color": col_swatch[i]
    #               }
    # }
    # for i in range(len(organism))
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
        'style': {'background-color': '#18B6F3','line-color': '#18B6F3'
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
    html.Span([str(i) + ": " + organism[i]], style={"color": col_swatch[i]})
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
                )  # remove ,
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
                )  # remove ,
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
                                    elements=edges + nodes,
                                    stylesheet=graph_stylesheet,
                                    #style={'width': '100%', 'height': '95vh'}
                                    minZoom = 0.6
                                )
                            ]
                        ),
                        dbc.Row(
                            [
                                dbc.Alert(
                                    id="node-data",
                                    children="Click on a node to see its details here",
                                    color="secondary"
                                )
                            ]
                        )
                    ],
                    sm=12,
                    md=8
                ),
                dbc.Col(
                    [
                        dbc.Badge(
                            "Level of pathway(s):",
                            color="success",
                            className="ms-1",

                        ),
                        dcc.Dropdown(
                            id="pathway_level",
                            options=[{"label": k, "value": k} for k in range(1,3)],
                            clearable=False,
                            value=1,
                            style={"width": "50px"},
                        ),
                        dbc.Badge(
                            "Model organism(s):", color="success", className="mr-1"
                        ),
                        dcc.Dropdown(
                            id="organism_dropdown",
                            options=[{"label": i + " (" + str(v)+ " phenotype(s))", "value": i} for i, v in org_totalPhenotype.items()],
                            #value= options[1],
                            multi=True,
                            style={"width": "100%"},
                        ),
                        dbc.Badge(
                            "Choose level of visualisation:", color="success", className="mr-1"
                        ),
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


@app.callback(
    Output("node-data", "children"), [Input("cytoscape-phenotype", "selectedNodeData")]
)
def display_nodedata(selectedNodeData):
    contents = "Click on a node to see its details here"
    if selectedNodeData is not None:
        if len(selectedNodeData) > 0:
            data = selectedNodeData[-1]
            contents = []
            contents.append(html.H5("Phenotype: " + data["id"].title()))
            contents.append(
                html.P(
                    "Lable: "
                    + data["label"].title()
                    + ", Size: "
                    + data["size"]
                )
            )
            #contents.append(
            #    html.P(
            #        "Color: "
            #        + str(data["classes"])
            #        + ", Position: "
            #        + str(data["position"])
            #    )
            #)

    return contents

@app.callback(
    Output("cytoscape-phenotype", "elements"),
    [
        Input("pathway_level", "value"),
        Input("organism_dropdown", "value"),
        Input("detail_slider", "value"),
    ],
)
def filter_nodes(pathway_level, organism_list, detail_slider):
    if (detail_slider == 3
        and organism_list == organism
        and pathway_lever == 1
    ):
        logger.info("Using the default element list")
        return nodes+edges

    else:
        # Generate node list
        if org_list is not None and org_list != []:
            new_df = df[df["Organism"].isin(org_list)]

        new_gene_df, new_N_genes = load_dataframe(new_df)
        new_nodes, new_edges = load_info_to_graph(new_gene_df, new_N_genes)

        #new_element = new_nodes + new_edges

    return new_nodes + new_edges



# end_time = perf_counter()
# print("\nElapsed time: " + str(end_time - sdk_start_time) + "s")  # get run time
if __name__ == '__main__':
    app.run_server(debug=False)  # set false so you can load bigger networks
