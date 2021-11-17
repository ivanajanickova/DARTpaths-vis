
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
from coordinates import check_coordinates, gene_y_coordinate


df_path = 'pathway_data/AmineOxidase.csv'
df = pd.read_csv(df_path, delimiter=',')
#def load_dataframe(pathway: str):
def load_dataframe(df):
    """ Function to load information from input csv file to pandas dataframe.
    It returns dataframe with loaded information and the total number of unique human genes
    in the dataframe as integer.
    :param pathway: relative pathway to input .csv file
    """
    #df = pd.read_csv(pathway, delimiter=',')
    df_1 = df.assign(phenotype=df['associated_phenotype'].astype(str).str.split(',')).explode('phenotype')
    df_2 = df_1[['1', '2', 'Organism', 'phenotype']]
    df_3 = df_2[:-1]

    # get total number of genes
    n_genes = set(df['2'].tolist())
    n = len(n_genes)  # total number of genes

    return df_3, n

gene_df, N_genes = load_dataframe(df)

#--------------------------------------------
def load_info_to_graph(dataframe: pd.DataFrame, n_genes: int):
    """ Function used to extract information from dataframe and load it into
    nodes and edges list of dictionaries

    :param dataframe: dataframe with information about nodes and edges
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
        gene = str(row['2'])
        # target node is ortholog
        ortholog = str(row['1'])
        phenotype = str(row['phenotype'])

        organism = str(row['Organism'])

        # add data and class info to the nodes
        # gene node positions
        pos_x = 50  # always the same; set to 500 in case of large network
        pos_y = gene_y_coordinate(1, n_genes, count_n)

        cy_gene = {'data': {'id': gene,
                            'label': gene,
                            'size': 4,
                            'fontsize': '1.5px'},
                   'classes': 'blue',
                   'position': {'x': pos_x, 'y': pos_y},
                   "selectable": True
                   }

        # ortholog node
        pos_x_ortholog, pos_y_ortholog = check_coordinates(coordinates_ort, pos_y, 20, 45, 55, 80)

        cy_ortholog = {'data': {'id': ortholog,
                                'label': ortholog,
                                'size': 2,
                                'fontsize': '1px',
                                'organism': organism},
                       'classes': organism,
                       'position': {'x': pos_x_ortholog, 'y': pos_y_ortholog},
                       "selectable": True
                       }
        cy_edge = {'data': {'id': gene + ortholog,
                            'source': gene,
                            'target': ortholog,
                            'width': '0.25',
                            'color': '#696969'}
                   }

        # add genes and orthologs to set nodes and cy nodes
        if gene not in nodes_set:
            nodes_set.add(gene)
            nodes_list.append(cy_gene)
            count_n += 1

        if ortholog not in nodes_set:
            # TODO not entirely sure this will work, since multiple genes can have single ortholog
            # it should, since this does not affect the edges
            nodes_set.add(ortholog)
            nodes_list.append(cy_ortholog)

        edges_list.append(cy_edge)

        # phenotype node
        pos_x_phenotype, pos_y_phenotype = check_coordinates(coordinates_phenotype, pos_y, 0, 20, 80, 100)

        if phenotype != "nan":
            # remove "nan" phenotype
            cy_phenotype = {'data': {'id': phenotype,
                                     'label': phenotype,
                                     'size': 1,
                                     'fontsize': '0.5px'},
                            'classes': 'purple',
                            'position': {'x': pos_x_phenotype, 'y': pos_y_phenotype},
                            "selectable": True,
                            "grabbable": False
                            }
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

nodes, edges = load_info_to_graph(gene_df, N_genes)
# Define app
organism = ["celegans", "zebrafish", "mouse", "slimemould", "dmelanogaster"]
#TODO calculate the right phenotype number
totalPhenotype = [12,34,56,78]
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
                For this demonstration, ... phenotypes from 5 model organism were used to visualize
                the relationship between phenotypes from different organism related to a particular
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
            Node size indicates number of mutual organism (the node bigger if many organisms shrare the same phenotype ), 
            and color indicates its organism.

            The edges show the relationship of phenotypes related in the same molecular 
            pathway in Reactome.

            Use the filters to choose the phenotype with:
            * certain numbers of pathway/ pathway level (to be decided)
            * certain numbers of model organism you want to investigatet
            Try showing or hiding ..... connections with the toggle button, and explore different visualisation options.
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
                            marks={1: "humangene",2:"organism_ortholog", 3:"organism_phenotype"},
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
