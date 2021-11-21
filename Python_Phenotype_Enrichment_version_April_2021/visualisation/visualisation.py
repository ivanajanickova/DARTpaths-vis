"""
;===================================================================================================
; Title:   Visualization of enriched phenotypes
; Authors: Alexandra, Ivana, Irene, Thao
;===================================================================================================

"""
# Import modules
from time import perf_counter
import pandas as pd
import dash
import dash_cytoscape as cyto
from dash import dcc
from dash import html
import dash_bootstrap_components as dbc

import coordinates
from Python_Phenotype_Enrichment_version_April_2021 import db_retrieve

sdk_start_time = perf_counter()


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

        cy_gene = {'data': {'id': gene, 'label': gene, 'size': 4, 'fontsize': '1.5px'}, 'classes': 'blue',
                   'position': {'x': pos_x, 'y': pos_y}}

        # ortholog node
        pos_x_ortholog, pos_y_ortholog = coordinates.check_coordinates(coordinates_ort, pos_y, n_genes, 20, 45, 55, 80)
        if ortholog != "NaN":
            cy_ortholog = {'data': {'id': ortholog, 'label': ortholog, 'size': 2, 'fontsize': '1px', 'organism': organism},
            'position': {'x': pos_x_ortholog, 'y': pos_y_ortholog}}
            cy_edge = {
                'data': {'id': gene + ortholog, 'source': gene, 'target': ortholog, 'width': '0.25', 'color': '#696969'}}

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

            cy_phenotype = {'data': {'id': phenotype, 'label': phenotype, 'metadata': metadata_phenotype, 'size': 1, 'fontsize': '0.5px'}, 'classes': 'purple',
                            'position': {'x': pos_x_phenotype, 'y': pos_y_phenotype}}
            cy_edge_2 = {'data': {'id': ortholog + phenotype, 'source': ortholog, 'target': phenotype, 'width': '0.10',
                             'color': '#B8B8B8'}}

            if phenotype not in nodes_set:
                nodes_set.add(phenotype)
                nodes_list.append(cy_phenotype)

            edges_list.append(cy_edge_2)

    return nodes_list, edges_list


# call functions
gene_df = db_retrieve.select_from_enrichment_results("AHR")
gene_df, N_genes = expand_dataframe(gene_df)
# metadata1 = db_retrieve.select_from_metadata("AHR")  # load metadata
# nodes, edges = load_info_to_graph(gene_df, metadata1, N_genes)

gene_df2 = db_retrieve.select_from_enrichment_results("AmineOxidase")
gene_df_2, N_genes_2 = expand_dataframe(gene_df2)
metadata2 = db_retrieve.select_from_metadata("AmineOxidase")
nodes, edges = load_info_to_graph(gene_df_2, metadata2, N_genes_2)
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
            'background-color': '#C918F3',
            'line-color': '#C918F3'
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
            'background-color': '#18B6F3',
            'line-color': '#18B6F3'
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
                # layout={'name': 'cose'},
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
end_time = perf_counter()
print("\nElapsed time: " + str(end_time - sdk_start_time) + "s")  # get run time
if __name__ == '__main__':
    app.run_server(debug=False)  # set false so you can load bigger networks

