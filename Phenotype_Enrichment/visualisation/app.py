#!/usr/bin/env python
# coding: utf-8
"""
;===================================================================================================
; Title:   Visualization of enriched phenotypes
; Authors: Alexandra, Ivana, Irene, Thao
;===================================================================================================

"""
# # Visualisation & Application
# ## Needed packages:

import pandas as pd
import numpy as np
import plotly.express as px
import dash
import dash_cytoscape as cyto
import dash_bootstrap_components as dbc
from dash import html
from dash import dcc
from dash.dependencies import Input, Output, State
import json
import coordinates
import db_retrieve


def expand_dataframe(df: pd.DataFrame):
    """Function to expand dataframe, so each phenotype is in a separate row.
    Returns new dataframe and the total number of genes

    :param df: dataframe to expand
    """
    df_1 = df.assign(phenotype=df['Enriched_Phenotypes'].astype(str).str.split(',')).explode('Enriched_Phenotypes')
    df_2 = df_1[['Ortholog_Genes', 'Human_ID', 'Organism', 'Enriched_Phenotypes', 'Human_Gene']]
    # get total number of genes
    n_genes = set(df['Human_ID'].tolist())
    n = len(n_genes)  # total number of genes

    return df_2, n


# ### Load the dataframe into a graph (nodes and edges)

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
        gene = str(row['Human_ID'])
        gene_name = str(row['Human_Gene'])
        # target node is ortholog
        ortholog = str(row['Ortholog_Genes'])
        phenotype = str(row['Enriched_Phenotypes'])
        organism = str(row['Organism'])
        # add data and class info to the nodes
        # gene node positions
        pos_x, pos_y = coordinates.gene_coordinate(1, n_genes, count_n)
        cy_gene = {'data': {'id': gene,
                            'indent': 'gene',
                            'label': gene_name,
                            'size': 4,
                            'fontsize': '1.5px'},
                   'classes': 'blue',
                   'selectable': True,  # Allow selected
                   'position': {'x': pos_x, 'y': pos_y}}

        # ortholog node
        pos_x_ortholog, pos_y_ortholog = coordinates.check_coordinates(coordinates_ort, pos_y, n_genes, 20, 45, 55, 80)
        if ortholog != "NaN":
            cy_ortholog = {'data': {'id': ortholog,
                                    'label': ortholog,
                                    'indent': 'ortholog',
                                    'size': 2,
                                    'fontsize': '1px',
                                    'organism': organism},
                           'classes': 'ortholog',
                           'selectable': True,  # Allow selected
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

        if phenotype != "nan":  # get rid of None node
            metadata_phenotype = metadata.get(phenotype)  # load metadata
            # phenotype node
            pos_x_phenotype, pos_y_phenotype = coordinates.check_coordinates(coordinates_phenotype, pos_y, n_genes, 0,
                                                                             20, 80, 100)

            cy_phenotype = {'data': {'id': phenotype,
                                     'label': phenotype,
                                     'indent': 'phenotype',
                                     'metadata': metadata_phenotype,
                                     'size': 1,
                                     'fontsize': '0.5px'},
                            'classes': 'purple',
                            'selectable': True,  # Allow selected
                            'grabbable': False,  # Not allow grab for now
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


# ## Load Data: first get the data from the database (db_retrieve);
# then create a df with the function expand_dataframe
# then load the metadata
# * The list of pathways we have so far:
#     * Phase2ConjugationOfCompounds (highest)
#     * Phase1CompundFunctionalization (1 level down)
#     * AmineOxidase (lowest)
#     * AminoAcidConjugation (lowest)
#     * EthanolOxidation (lowest)
#     * AHR (lowest)

# upper level
upper_df1 = db_retrieve.select_from_enrichment_results("Phase1CompoundFunctionalization")
upper_df_1, upper_N_genes_1 = expand_dataframe(upper_df1)
metadata = db_retrieve.select_from_metadata("Phase1CompoundFunctionalization")

# AHR
gene_df1 = db_retrieve.select_from_enrichment_results("AHR")  # gene_df1 (type = pandas.core.frame.DataFrame)
gene_df_1, N_genes_1 = expand_dataframe(gene_df1)
metadata1 = db_retrieve.select_from_metadata("AHR")  # metadata1

# AmineOxidase
gene_df2 = db_retrieve.select_from_enrichment_results("AmineOxidase")  # gene_df2
gene_df_2, N_genes_2 = expand_dataframe(gene_df2)
metadata2 = db_retrieve.select_from_metadata("AmineOxidase")  # metadata2

# EthanolOxidation
gene_df3 = db_retrieve.select_from_enrichment_results("EthanolOxidation")  # gene_df4
gene_df_3, N_genes_3 = expand_dataframe(gene_df3)
metadata3 = db_retrieve.select_from_metadata("EthanolOxidation")  # metadata4

# list of the possible pathways:
# add pathways to this list, when you add more pathways
path_list = ['AHR', 'Amino Oxidase', 'Ethanol Oxidation']

# nodes & edges: for each pathway
# for each dataframe and metadata, create the nodes and edges with the function load_info_to_graph
nodes, edges = load_info_to_graph(upper_df_1, metadata, upper_N_genes_1)  # upper level
nodes1, edges1 = load_info_to_graph(gene_df_1, metadata1, N_genes_1)
nodes2, edges2 = load_info_to_graph(gene_df_2, metadata2, N_genes_2)
nodes3, edges3 = load_info_to_graph(gene_df_3, metadata3, N_genes_3)

# create an element for each pathway, contains the nodes and elements:
# this is 1 element for each possible graph, with both edges and nodes
element = nodes + edges
new_element1 = nodes1 + edges1
new_element2 = nodes2 + edges2
new_element3 = nodes3 + edges3


# list of all the elements:
new_elements = [element, new_element1, new_element2, new_element3]


def find_human_genes(elements):
    """
    Find all human genes -> in list

    :param elements: edges+nodes
    """
    human_genes_list = []
    for i in list(range(0, len(elements))):
        x = str(elements[i].get('data').get('id'))
        z = str(elements[i].get('data').get('label'))
        if x[0:3] == 'ENS':
            human_genes_list.append(x)
            human_genes_list.append(z)


    return human_genes_list


def edgesNnodes_HG(els, human_gene):
    """
    Function to find the nodes and edges for a human gene -> needed for the search function

    :param els: elements = edges + nodes
    :param human gene: input form the search bar -> string with 'ENSG'-number
    """
    set_edges = []
    set_nodes = []
    ortholog_names = []
    phenotype_names = []
    # find edges and nodes the the human gene
    for i in list(range(0, len(els))):
        x = els[i]

        # when human_gene is gene_name:
        if els[i].get('data').get('label') == str(human_gene):
            ens_gene = str(els[i].get('data').get('id'))
            #nodes
            set_nodes.append(x)
            #edges
            for j in list(range(0, len(els))):
                t = els[j]
                if els[j].get('data').get('source') == ens_gene:
                    set_edges.append(t)


        # when human_gene = ensg number
        # nodes
        if els[i].get('data').get('id') == str(human_gene):
            set_nodes.append(x)

        # edges
        if els[i].get('data').get('source') == str(human_gene):
            set_edges.append(x)


    # find the orthologs names from the human genes
    for j in list(range(0, len(set_edges))):
        y = set_edges[j].get('data').get('target')
        if y not in ortholog_names:
            ortholog_names.append(y)

    # find the edges and nodes from the orthologs
    for i in list(range(0, len(els))):
        z = els[i]
        for j in list(range(0, len(ortholog_names))):
            # edges
            if els[i].get('data').get('source') == ortholog_names[j]:
                set_edges.append(z)
            # nodes
            if els[i].get('data').get('id') == ortholog_names[j]:
                set_nodes.append(z)

    # find the phenotyes names
    for j in list(range(0, len(set_edges))):
        y = set_edges[j].get('data').get('target')
        if y not in phenotype_names:
            phenotype_names.append(y)

    # find the nodes from the phenotypes
    for i in list(range(0, len(els))):
        z = els[i]
        for j in list(range(0, len(phenotype_names))):
            # nodes
            if els[i].get('data').get('id') == phenotype_names[j]:
                set_nodes.append(z)

    return set_edges, set_nodes


def edgesNnodes_Ortholog(els, ortholog):
    """
    Function to find the nodes and edges for ortholog -> needed for the ortholog selection

    :param els: element = edges + nodes
    :param ortholog: possible ortholog (mouse, zebrafish, celegans, dmelanogaster)
    """
    ortholog_nodes = []
    ortholog_edges = []

    ortholog_names = []
    hg_names = []
    phenotype_names = []

    # find the nodes of all orthologs
    for i in list(range(0, len(els))):
        x = els[i]
        # nodes
        if els[i].get('data').get('organism') == ortholog:
            # get names of orthologs
            y = els[i].get('data').get('label')
            ortholog_names.append(y)
            ortholog_nodes.append(x)  # only nodes here

    # find the human genes names from the otholog names and edges between human gene and ortholog & also for phenotype
    for i in list(range(0, len(els))):
        z = els[i]  # data in element (nodes and edges)
        for j in list(range(0, len(ortholog_names))):
            if els[i].get('data').get('target') == ortholog_names[j]:  # source is HG, target is Ortholog
                hg = els[i].get('data').get('source')  # human genes
                if hg not in hg_names:  # unique values
                    hg_names.append(hg)
                ortholog_edges.append(z)  # if the name of the ortholog is the target, then add the edge into the list
            elif els[i].get('data').get('source') == ortholog_names[j]:  # source is ortholog, target is phenotype
                pt = els[i].get('data').get('target')  # phenotypes
                if pt not in phenotype_names:  # unique values
                    phenotype_names.append(pt)
                ortholog_edges.append(z)

    for i in list(range(0, len(els))):
        a = els[i]
        for j in list(range(0, len(hg_names))):
            if els[i].get('data').get('id') == hg_names[j]:
                ortholog_nodes.append(a)  # nodes of hg
        for h in list(range(0, len(phenotype_names))):
            if els[i].get('data').get('label') == phenotype_names[h]:
                ortholog_nodes.append(a)  # nodes of phenotypes

    return ortholog_edges, ortholog_nodes


# ## make list of the nodes for each organism -> less time consuming than calculate when app is running
# dmelanogaster
ED, ND = edgesNnodes_Ortholog(new_elements[0], 'dmelanogaster')
ED1, ND1 = edgesNnodes_Ortholog(new_elements[1], 'dmelanogaster')
ED2, ND2 = edgesNnodes_Ortholog(new_elements[2], 'dmelanogaster')
ED3, ND3 = edgesNnodes_Ortholog(new_elements[3], 'dmelanogaster')

# mouse
EM, NM = edgesNnodes_Ortholog(new_elements[0], 'mouse')
EM1, NM1 = edgesNnodes_Ortholog(new_elements[1], 'mouse')
EM2, NM2 = edgesNnodes_Ortholog(new_elements[2], 'mouse')
EM3, NM3 = edgesNnodes_Ortholog(new_elements[3], 'mouse')

# celegans
EC, NC = edgesNnodes_Ortholog(new_elements[0], 'celegans')
EC1, NC1 = edgesNnodes_Ortholog(new_elements[1], 'celegans')
EC2, NC2 = edgesNnodes_Ortholog(new_elements[2], 'celegans')
EC3, NC3 = edgesNnodes_Ortholog(new_elements[3], 'celegans')

# zebrafish
EZ, NZ = edgesNnodes_Ortholog(new_elements[0], 'zebrafish')
EZ1, NZ1 = edgesNnodes_Ortholog(new_elements[1], 'zebrafish')
EZ2, NZ2 = edgesNnodes_Ortholog(new_elements[2], 'zebrafish')
EZ3, NZ3 = edgesNnodes_Ortholog(new_elements[3], 'zebrafish')

# ## Define organism and colors
# organisms list:
organism = ['dmelanogaster', 'mouse', 'celegans', 'zebrafish']
# https://plotly.com/python/discrete-color/
col_swatch = [px.colors.qualitative.Safe[8], px.colors.qualitative.Set1[2],
              px.colors.qualitative.Set1[7], px.colors.qualitative.Set1[4]]

# ## Define app:
graph_stylesheet = [

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
    brand="Dart - Translation of Genotype to Phenotype by a Hierarchy of Model Organism",  # title
    brand_href="#",
    color="dark",
    dark=True,
)

topics_html = list()
for topic_html in [
    html.Span([str(i + 1) + ": " + organism[i]], style={"color": col_swatch[i]})
    for i in range(len(organism))
]:
    topics_html.append(topic_html)
    topics_html.append(html.Br())

# ## Body layout
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
                For this demonstration, orthologs of four model organisms were used to visualize the
                 relationship between phenotypes of different organisms related to a particular human gene. The input pathways can be found on [REACTOME](https://reactome.org/PathwayBrowser/).

                The gene orthologs of each model organism are shown in different colours , as shown in the legend on the right.
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
                ##### Model organisms:
                -----
                """
                        ),
                        html.Div(
                            topics_html,
                            style={
                                "fontSize": 15,
                                "height": "100px",
                                "overflow": "auto",
                            },
                        ),
                    ],
                    sm=15,
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
            There are three types of nodes in the graph. The blue nodes in the center indicate human genes in a chosen pathway. These nodes are directly linked to the second type of node indicating the orthologs of a model organism. The node colour indicates which organism the ortholog belongs to. The third type of node (purple) on the two sides is the phenotype enrichment of the model organism. The edges show the relationship of:
            1) human genes involved in Reactome pathway and their orthologs in selected model organisms, and
            2) The orthologs in the organism and the related enriched phenotypes.

            Initially, a specific gene can be searched for using its 'ENSG' ID or gene name. This will display the specific graph for that gene.
            The pathway is divided into upper and lower levels. This can be chosen at Level of pathway(s). The lower pathway level consists of several sub-pathways, which can also be specified.
            Furthermore, filtering can be done on the basis of the model organism. Thus it is possible to only look at the graph of specific organism.

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
                                    elements=element,
                                    stylesheet=graph_stylesheet,
                                    minZoom=0.6
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
                        dbc.Label('Search for a gene: ', html_for='gene_search'),
                        html.Div(dbc.Input(id='gene_search',
                                           type='text',
                                           value='',
                                           placeholder="Enter a gene ID here")),

                        html.Button("Search", id="search-button", n_clicks=0),

                        html.Div(id='output-container-button',
                                 children='Fill in the gene name and press Search'),

                        html.Span(id="search-status", style={"verticalAlign": "middle"}),

                        html.Div(style={'padding': 20}),  # Space

                        dbc.Row(
                            [
                                dbc.Col(
                                    [
                                        html.H4(dbc.Badge(
                                            "Level of pathway(s):",
                                            color="success",
                                            className="ms-1",

                                        )),

                                        dcc.Dropdown(
                                            id="pathway_level",
                                            options=[{"label": 'lower level', "value": 1},
                                                     {"label": 'upper level', "value": 2}],
                                            clearable=False,
                                            value=2,
                                            style={"width": "150px"},
                                        )
                                    ]),
                                dbc.Col(
                                    [
                                        dbc.Badge(  # dropdown to choose pathway of lowest level
                                            "Available sub-pathway:    ",
                                            color="white",
                                            text_color="muted",
                                            class_name="border me-1",
                                        ),
                                        html.Div(style={'padding': 7}),  # Space
                                        dcc.Dropdown(
                                            id="Lowest_Level",
                                            options=[{"label": name, "value": name}
                                                     for name in path_list
                                                     ],
                                            clearable=False,
                                            value=path_list[0],  # 'AHR'
                                            style={"width": "187px"},
                                        )])]),

                        html.H4(dbc.Badge(
                            "Model organism(s):", color="success", className="mr-1"
                        )),

                        dcc.Dropdown(
                            id="organism_select",
                            options=[{"label": i, "value": i}
                                     for i in ['all organisms', 'dmelanogaster', 'mouse', 'celegans', 'zebrafish']],
                            clearable=False,
                            value='all organisms',
                            style={"width": "100%"},
                        ),

                        html.Div(style={'padding': 20}),  # Space

                    ]
                )
            ]
        )

    ]
)

app.layout = html.Div([navbar, body_layout])


# ## Dash Callbacks

# ### Click on a node to see its details here
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
            contents.append(html.P('Phenotype enrichment Q-value: ' + str(data.get('metadata')[2])))
            contents.append(html.P('Phenotype enrichment P-value: ' + str(data.get('metadata')[4])))
            return contents
        elif data.get('indent') == 'gene':
            contents.append(html.H5('Gene name: ' + data.get('label')))
            contents.append(html.P('Gene ENSG ID: ' + data.get('id')))
            return contents
        elif data.get('indent') == 'ortholog':
            contents.append(html.H5('Ortholog ID: ' + data.get('label')))
            contents.append(html.P('Organism: ' + data.get('organism')))
            return contents
    else:
        return contents


# ### Search human gene and see the linked orthologs and phenotypes &
# Level of pathway's (upper vs lower) & Lower pathways (options) & Model organisms:

@app.callback(
    [Output('cytoscape-phenotype', "elements"),
     Output('output-container-button', 'children')],
    [Input('organism_select', 'value'),
     Input('pathway_level', 'value'),
     Input('Lowest_Level', 'value'),
     Input('search-button', 'n_clicks'),
     Input('gene_search', 'value')]
)
def select_organism(value_orth, UpLow, LowLevel, n, hg_value):
    # change the click on the search button:
    changed_id = [p['prop_id'] for p in dash.callback_context.triggered][0]  # To determine if n_clicks is changed.

#overview:
    # first start by looking if there you look for 1 organism or more
    # then see if you are in the upper or the lower level
    # UPPER -> you cannot choose a lower level (in sub level) -> only use the graph of the upper leve
    # if there is no search -> give whole graph
    # if there is a search -> give the graph by using the function edgesNnodes_HG (give the edges and the nodes for that human gene)
    # if the gene is no found -> give: "The gene you were looking for was not found"

    # LOWER -> first check in which pathway you are -> go through pathway_list
    # then the options are the same as in UPPER LEVEL:
    # 1) no search = whole graph, 2) search = edgesNnodes_HG; 3) not found = text

    if value_orth == 'all organisms':
        if UpLow == 1:
            if LowLevel == path_list[0]:
                elements = edges1 + nodes1
                if 'search-button' not in changed_id:
                    text = "Enter a value and press Search"
                    return elements, text
                else:
                    hg = "{}".format(hg_value)
                    if str(hg) in find_human_genes(elements):
                        edgesHG, nodesHG = edgesNnodes_HG(elements, hg)
                        el = edgesHG + nodesHG
                        text = "you searched for: " + hg
                        return el, text
                    else:
                        text = "The gene you were looking for was not found"
                        return elements, text
            elif LowLevel == path_list[1]:
                elements = edges2 + nodes2
                if 'search-button' not in changed_id:
                    text = "Enter a value and press Search"
                    return elements, text
                else:
                    hg = "{}".format(hg_value)
                    if str(hg) in find_human_genes(elements):
                        edgesHG, nodesHG = edgesNnodes_HG(elements, hg)
                        el = edgesHG + nodesHG
                        text = "you searched for: " + hg
                        return el, text
                    else:
                        text = "The gene you were looking for was not found"
                        return elements, text
            elif LowLevel == path_list[2]:
                elements = edges3 + nodes3
                if 'search-button' not in changed_id:
                    text = "Enter a value and press Search"
                    return elements, text
                else:
                    hg = "{}".format(hg_value)
                    if str(hg) in find_human_genes(elements):
                        edgesHG, nodesHG = edgesNnodes_HG(elements, hg)
                        el = edgesHG + nodesHG
                        text = "you searched for: " + hg
                        return el, text
                    else:
                        text = "The gene you were looking for was not found"
                        return elements, text
            elements = edges1 + nodes1
            if 'search-button' not in changed_id:
                text = "Enter a value and press Search"
                return elements, text
            else:
                hg = "{}".format(hg_value)
                if str(hg) in find_human_genes(elements):
                    edgesHG, nodesHG = edgesNnodes_HG(elements, hg)
                    el = edgesHG + nodesHG
                    text = "you searched for: " + hg
                    return el, text
                else:
                    text = "The gene you were looking for was not found"
                    return elements, text
        else:
            elements = edges + nodes
            if 'search-button' not in changed_id:
                text = "Enter a value and press Search"
                return elements, text
            else:
                hg = "{}".format(hg_value)
                if str(hg) in find_human_genes(elements):
                    edgesHG, nodesHG = edgesNnodes_HG(elements, hg)
                    el = edgesHG + nodesHG
                    text = "you searched for: " + hg
                    return el, text
                else:
                    text = "The gene you were looking for was not found"
                    return elements, text
    elif value_orth == 'dmelanogaster':
        if UpLow == 1:
            if LowLevel == path_list[0]:
                elements = ED1 + ND1
                if 'search-button' not in changed_id:
                    text = "Enter a value and press Search"
                    return elements, text
                else:
                    hg = "{}".format(hg_value)
                    if str(hg) in find_human_genes(elements):
                        edgesHG, nodesHG = edgesNnodes_HG(elements, hg)
                        el = edgesHG + nodesHG
                        text = "you searched for: " + hg
                        return el, text
                    else:
                        text = "The gene you were looking for was not found"
                        return elements, text
            elif LowLevel == path_list[1]:
                elements = ED2 + ND2
                if 'search-button' not in changed_id:
                    text = "Enter a value and press Search"
                    return elements, text
                else:
                    hg = "{}".format(hg_value)
                    if str(hg) in find_human_genes(elements):
                        edgesHG, nodesHG = edgesNnodes_HG(elements, hg)
                        el = edgesHG + nodesHG
                        text = "you searched for: " + hg
                        return el, text
                    else:
                        text = "The gene you were looking for was not found"
                        return elements, text
            elif LowLevel == path_list[2]:
                elements = ED3 + ND3
                if 'search-button' not in changed_id:
                    text = "Enter a value and press Search"
                    return elements, text
                else:
                    hg = "{}".format(hg_value)
                    if str(hg) in find_human_genes(elements):
                        edgesHG, nodesHG = edgesNnodes_HG(elements, hg)
                        el = edgesHG + nodesHG
                        text = "you searched for: " + hg
                        return el, text
                    else:
                        text = "The gene you were looking for was not found"
                        return elements, text
            elements = ED1 + ND1
            if 'search-button' not in changed_id:
                text = "Enter a value and press Search"
                return elements, text
            else:
                hg = "{}".format(hg_value)
                if str(hg) in find_human_genes(elements):
                    edgesHG, nodesHG = edgesNnodes_HG(elements, hg)
                    el = edgesHG + nodesHG
                    text = "you searched for: " + hg
                    return el, text
                else:
                    text = "The gene you were looking for was not found"
                    return elements, text
        else:
            elements = ED + ND
            if 'search-button' not in changed_id:
                text = "Enter a value and press Search"
                return elements, text
            else:
                hg = "{}".format(hg_value)
                if str(hg) in find_human_genes(elements):
                    edgesHG, nodesHG = edgesNnodes_HG(elements, hg)
                    el = edgesHG + nodesHG
                    text = "you searched for: " + hg
                    return el, text
                else:
                    text = "The gene you were looking for was not found"
                    return elements, text
    elif value_orth == 'mouse':
        if UpLow == 1:
            if LowLevel == path_list[0]:
                elements = EM1 + NM1
                if 'search-button' not in changed_id:
                    text = "Enter a value and press Search"
                    return elements, text
                else:
                    hg = "{}".format(hg_value)
                    if str(hg) in find_human_genes(elements):
                        edgesHG, nodesHG = edgesNnodes_HG(elements, hg)
                        el = edgesHG + nodesHG
                        text = "you searched for: " + hg
                        return el, text
                    else:
                        text = "The gene you were looking for was not found"
                        return elements, text
            elif LowLevel == path_list[1]:
                elements = EM2 + NM2
                if 'search-button' not in changed_id:
                    text = "Enter a value and press Search"
                    return elements, text
                else:
                    hg = "{}".format(hg_value)
                    if str(hg) in find_human_genes(elements):
                        edgesHG, nodesHG = edgesNnodes_HG(elements, hg)
                        el = edgesHG + nodesHG
                        text = "you searched for: " + hg
                        return el, text
                    else:
                        text = "The gene you were looking for was not found"
                        return elements, text
            elif LowLevel == path_list[2]:
                elements = EM3 + NM3
                if 'search-button' not in changed_id:
                    text = "Enter a value and press Search"
                    return elements, text
                else:
                    hg = "{}".format(hg_value)
                    if str(hg) in find_human_genes(elements):
                        edgesHG, nodesHG = edgesNnodes_HG(elements, hg)
                        el = edgesHG + nodesHG
                        text = "you searched for: " + hg
                        return el, text
                    else:
                        text = "The gene you were looking for was not found"
                        return elements, text
            elements = EM1 + NM1
            if 'search-button' not in changed_id:
                text = "Enter a value and press Search"
                return elements, text
            else:
                hg = "{}".format(hg_value)
                if str(hg) in find_human_genes(elements):
                    edgesHG, nodesHG = edgesNnodes_HG(elements, hg)
                    el = edgesHG + nodesHG
                    text = "you searched for: " + hg
                    return el, text
                else:
                    text = "The gene you were looking for was not found"
                    return elements, text
        else:
            elements = EM + NM
            if 'search-button' not in changed_id:
                text = "Enter a value and press Search"
                return elements, text
            else:
                hg = "{}".format(hg_value)
                if str(hg) in find_human_genes(elements):
                    edgesHG, nodesHG = edgesNnodes_HG(elements, hg)
                    el = edgesHG + nodesHG
                    text = "you searched for: " + hg
                    return el, text
                else:
                    text = "The gene you were looking for was not found"
                    return elements, text
    elif value_orth == 'celegans':
        if UpLow == 1:
            if LowLevel == path_list[0]:
                elements = EC1 + NC1
                if 'search-button' not in changed_id:
                    text = "Enter a value and press Search"
                    return elements, text
                else:
                    hg = "{}".format(hg_value)
                    if str(hg) in find_human_genes(elements):
                        edgesHG, nodesHG = edgesNnodes_HG(elements, hg)
                        el = edgesHG + nodesHG
                        text = "you searched for: " + hg
                        return el, text
                    else:
                        text = "The gene you were looking for was not found"
                        return elements, text
            elif LowLevel == path_list[1]:
                elements = EC2 + NC2
                if 'search-button' not in changed_id:
                    text = "Enter a value and press Search"
                    return elements, text
                else:
                    hg = "{}".format(hg_value)
                    if str(hg) in find_human_genes(elements):
                        edgesHG, nodesHG = edgesNnodes_HG(elements, hg)
                        el = edgesHG + nodesHG
                        text = "you searched for: " + hg
                        return el, text
                    else:
                        text = "The gene you were looking for was not found"
                        return elements, text
            elif LowLevel == path_list[2]:
                elements = EC3 + NC3
                if 'search-button' not in changed_id:
                    text = "Enter a value and press Search"
                    return elements, text
                else:
                    hg = "{}".format(hg_value)
                    if str(hg) in find_human_genes(elements):
                        edgesHG, nodesHG = edgesNnodes_HG(elements, hg)
                        el = edgesHG + nodesHG
                        text = "you searched for: " + hg
                        return el, text
                    else:
                        text = "The gene you were looking for was not found"
                        return elements, text
            elements = EC1 + NC1
            if 'search-button' not in changed_id:
                text = "Enter a value and press Search"
                return elements, text
            else:
                hg = "{}".format(hg_value)
                for i in list(range(0, len(elements))):
                    if str(hg) in find_human_genes(elements):
                        edgesHG, nodesHG = edgesNnodes_HG(elements, hg)
                        el = edgesHG + nodesHG
                        text = "you searched for: " + hg
                        return el, text
                    else:
                        text = "The gene you were looking for was not found"
                        return elements, text
        else:
            elements = EC + NC
            if 'search-button' not in changed_id:
                text = "Enter a value and press Search"
                return elements, text
            else:
                hg = "{}".format(hg_value)
                for i in list(range(0, len(elements))):
                    if str(hg) in find_human_genes(elements):
                        edgesHG, nodesHG = edgesNnodes_HG(elements, hg)
                        el = edgesHG + nodesHG
                        text = "you searched for: " + hg
                        return el, text
                    else:
                        text = "The gene you were looking for was not found"
                        return elements, text
    elif value_orth == 'zebrafish':
        if UpLow == 1:
            if LowLevel == path_list[0]:
                elements = EZ1 + NZ1
                if 'search-button' not in changed_id:
                    text = "Enter a value and press Search"
                    return elements, text
                else:
                    hg = "{}".format(hg_value)
                    if str(hg) in find_human_genes(elements):
                        edgesHG, nodesHG = edgesNnodes_HG(elements, hg)
                        el = edgesHG + nodesHG
                        text = "you searched for: " + hg
                        return el, text
                    else:
                        text = "The gene you were looking for was not found"
                        return elements, text
            elif LowLevel == path_list[1]:
                elements = EZ2 + NZ2
                if 'search-button' not in changed_id:
                    text = "Enter a value and press Search"
                    return elements, text
                else:
                    hg = "{}".format(hg_value)
                    if str(hg) in find_human_genes(elements):
                        edgesHG, nodesHG = edgesNnodes_HG(elements, hg)
                        el = edgesHG + nodesHG
                        text = "you searched for: " + hg
                        return el, text
                    else:
                        text = "The gene you were looking for was not found"
                        return elements, text
            elif LowLevel == path_list[2]:
                elements = EZ3 + NZ3
                if 'search-button' not in changed_id:
                    text = "Enter a value and press Search"
                    return elements, text
                else:
                    hg = "{}".format(hg_value)
                    if str(hg) in find_human_genes(elements):
                        edgesHG, nodesHG = edgesNnodes_HG(elements, hg)
                        el = edgesHG + nodesHG
                        text = "you searched for: " + hg
                        return el, text
                    else:
                        text = "The gene you were looking for was not found"
                        return elements, text
            elements = EZ1 + NZ1
            if 'search-button' not in changed_id:
                text = "Enter a value and press Search"
                return elements, text
            else:
                hg = "{}".format(hg_value)
                if str(hg) in find_human_genes(elements):
                    edgesHG, nodesHG = edgesNnodes_HG(elements, hg)
                    el = edgesHG + nodesHG
                    text = "you searched for: " + hg
                    return el, text
                else:
                    text = "The gene you were looking for was not found"
                    return elements, text
        else:
            elements = EZ + NZ
            if 'search-button' not in changed_id:
                text = "Enter a value and press Search"
                return elements, text
            else:
                hg = "{}".format(hg_value)
                if str(hg) in find_human_genes(elements):
                    edgesHG, nodesHG = edgesNnodes_HG(elements, hg)
                    el = edgesHG + nodesHG
                    text = "you searched for: " + hg
                    return el, text
                else:
                    text = "The gene you were looking for was not found"
                    return elements, text


# ### Running the app:
if __name__ == '__main__':
    app.run_server(debug=False)  # set false so you can load bigger networks
