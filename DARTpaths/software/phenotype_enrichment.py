#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 17:57:22 2020

@author: d
"""
#command: python3 -.py R-HSA-69239 "PATHWAY_NAME"

### Reference for calculating factorial of gamma
### log factorial:
### https://www.perlmonks.org/bare/?node_id=443318
### https://www.perlmonks.org/bare/?node_id=466322

##########################
#import necessary modules#
##########################

import re
import pandas as pd
import sys #for user input
import os
import math
import glob

#this is necessary for zebrafish
pd.options.mode.chained_assignment = None

############
#the inputs#
############

### human pathway's Reactome ID that will be enriched
pathID = sys.argv[1]   ### e.g.) R-HSA-8937144

### name of the pathway.
### This is what you choose for that pathway in your folder/path
pathway_name = sys.argv[2]  ### e.g.) "Hedgehog"

### path to DARTpaths folder e.g.) $HOME/DARTpaths/
path_root=sys.argv[3]
##
#organism = sys.argv[] #the organism to map which orthologous genes and which has enriched phenotypes
##this is now disabled, with this code produce enrichment results for every species just by one run

###### ----------------------------------------------------------------------------- #############################

path_to_databases = path_root+"databases"
path_to_phenotype_database_folder = path_to_databases+"/phenotype"
path_to_ontology_database_folder = path_to_databases+"/ontology"

path_to_pathway_folder = path_root+"Proteins/"+pathway_name+"/Orthologs_"+pathway_name+"_pathway/"

Mouse_orthologs_file_path = path_to_pathway_folder +"Merged_"+pathway_name+"_Mouse_orthologs.txt"
Rat_orthologs_file_path = path_to_pathway_folder +"Merged_"+pathway_name+"_Rat_orthologs.txt"
Cele_orthologs_file_path = path_to_pathway_folder +"Merged_"+pathway_name+"_Cele_orthologs.txt"
Zebrafish_orthologs_file_path = path_to_pathway_folder +"Merged_"+pathway_name+"_Zebrafish_orthologs.txt"
Rabbit_orthologs_file_path = path_to_pathway_folder +"Merged_"+pathway_name+"_Rabbit_orthologs.txt"
Slimemould_orthologs_file_path = path_to_pathway_folder +"Merged_"+pathway_name+"_Dicty_orthologs.txt"

######  Data Resources

Reactome_file = path_to_phenotype_database_folder+"/"+"Ensembl2Reactome_All_Levels.txt"

### ENSEMBL DATABASES (orthologs) for each species ##### same as the one used for orthology mapping part

C_ele_ENSEMBL_orthology_database = path_to_databases+"/"+"Orthology_human_celegans_ensembl101_unique.txt"
Zebrafish_ENSEMBL_orthology_database = path_to_databases+"/"+"Orthology_human_zebrafish_ensembl101_unique.txt"
Mouse_ENSEMBL_orthology_database = path_to_databases+"/"+"Orthology_human_mouse_ensembl101_unique.txt"
Slime_mould_ENSEMBL_orthology_database = path_to_databases+"/"+"Orthology_human_dicty_ensembl101_unique.txt"

### Phenotype databases for each species ####

C_ele_phenotype_databse = path_to_phenotype_database_folder+"/"+"phenotype_association.WS264.wb"
Zebrafish_phenotype_databse = path_to_phenotype_database_folder+"/"+"phenoGeneCleanData_fish_2020.09.25.txt"
Mouse_phenotype_databse = path_to_phenotype_database_folder+"/"+"ALL_genotype_phenotype.csv"
Slime_mould_phenotype_databse = path_to_phenotype_database_folder+"/"+"all-mutants-ddb_g.txt"

###### Ontology databases

C_ele_ontology_data_file = path_to_ontology_database_folder+"/"+"phenotype_ontology.WS264.obo.txt"

Zebrafish_ontology_data_file = path_to_ontology_database_folder+"/"+"anatomy_item_2020.10.20.txt"
Zebrafish_GO_annotation_file = path_to_ontology_database_folder+"/"+"go_annotation.obo"

Mouse_ontology_data_file = path_to_ontology_database_folder+"/"+"MPheno_OBO.ontology.txt"

#### Write Output to File ####
#### change this to include pathway identifier and species name to keep results organized.

phenotype_enrichment_result_file = pathID + "_phenotype_enrichment_result.txt"

##### --------------------------------------------------------------------------- ###############################

#############################################################
#Opening, Extracting, and Mapping Information from Resources#
#############################################################

#open the reactome resource
def openReactome():
  reactomefile = open(Reactome_file, "r")
  return reactomefile

####load new orthologs

#load the orthologs based on the organism.
#use try and except in case there's no new orthologs
def loadNewOrtholog(organism):

  df_new_orthos_raw2 = pd.DataFrame()
  df_new_orthos = pd.DataFrame()

  if organism == "celegans":
    Cele_orthologs_file_path = path_to_pathway_folder +"Merged_"+pathway_name+"_Cele_orthologs.txt"
    try:
      df_new_orthos_raw = pd.read_csv(Cele_orthologs_file_path, delimiter=',',header = None)
      df_new_orthos_raw2 = df_new_orthos_raw.iloc[:,1].copy()
    except:
      pass

  elif organism == "zebrafish":
    Zebrafish_orthologs_file_path = path_to_pathway_folder +"Merged_"+pathway_name+"_Zebrafish_orthologs.txt"
    try:
      df_new_orthos_raw = pd.read_csv(Zebrafish_orthologs_file_path, delimiter=',',header = None)
      df_new_orthos_raw2 = df_new_orthos_raw.iloc[:,3].copy()
    except:
      pass

  elif organism == "mouse":
    Mouse_orthologs_file_path = path_to_pathway_folder +"Merged_"+pathway_name+"_Mouse_orthologs.txt"
    try:
      df_new_orthos_raw = pd.read_csv(Mouse_orthologs_file_path, delimiter=',',header = None)
      df_new_orthos_raw2 = df_new_orthos_raw.iloc[:,3].copy()
    except:
      pass

  elif organism == "slimemould":
    Slimemould_orthologs_file_path = path_to_pathway_folder +"Merged_"+pathway_name+"_Dicty_orthologs.txt"
    try:
      df_new_orthos_raw = pd.read_csv(Slimemould_orthologs_file_path, delimiter=',',header = None)
      df_new_orthos_raw2 = df_new_orthos_raw.iloc[:,1].copy()
    except:
      pass

  if df_new_orthos_raw2.empty:
    pass
  else:
    df_new_orthos = df_new_orthos_raw2.to_frame()
    df_new_orthos.rename(columns={df_new_orthos.columns[0]: 'Genes'}, inplace=True)

  return df_new_orthos


##load ID mapping database
#set the paths
#path_to_ID_DB = path_root+"/ID_mapping_database"

#set path and load in same function
#def LoadIDMapping(organism):
  #df_full_path = path_to_ID_DB+"/ENSEMBL_"+organism+"_Gene_Prot_ID_mapping.txt"
  #df_full = pd.read_csv(df_full_path,delimiter='\t',header = 'infer')
  #return df_full


#Search and select the input pathway in the reactome file
#throw exception when pathway not found
#two possible way after the exception is to ask user for another input
#or the code just stops, and the user just need to rerun the code with another pathway input
#either way in this version it's only thrown the pathway not found.
def selectPathway(reactomefile):
  pathway = []

  for line in reactomefile:
     if re.search(pathID, line):
         stripped_line = line.strip()
         line_list = stripped_line.split()
         pathway.append(line_list)
         if line == None:
           sys.stderr.write("Pathway not found") #the exception
  reactomefile.close()
  pathgenes_raw = list(zip(*pathway))[0] #select only pathway genes
  pathgenes = list(pathgenes_raw)

  return pathgenes

#Open the orthologs resource depends on the organism
#Then select the relevant column names for mapping to other database (in this case, to Phenotype database)
#Gene stable ID/name so it maps to Phenotype database
#Human gene stable ID so it maps from the reactome pathway
#I rename the files so its easier, mart means its taken from biomart
def selectForMapping(organism):

    if organism == "celegans":
      df_celegans = pd.read_csv(C_ele_ENSEMBL_orthology_database, delimiter='\t',header = 'infer')
      selectgenes= df_celegans[['Caenorhabditis elegans gene stable ID','Gene stable ID']].copy()

    elif organism == "zebrafish":
      df_zebrafish = pd.read_csv(Zebrafish_ENSEMBL_orthology_database, delimiter='\t',header = 'infer')
      selectgenes= df_zebrafish[['Zebrafish gene name','Gene stable ID']].copy()

    elif organism == "mouse":
      df_mouse = pd.read_csv(Mouse_ENSEMBL_orthology_database, delimiter='\t', header = 'infer')
      selectgenes= df_mouse[['Mouse gene name','Gene stable ID']].copy()

    elif organism == "slimemould":
      df_slime = pd.read_csv(Slime_mould_ENSEMBL_orthology_database, delimiter='\t', header='infer')
      selectgenes= df_slime[['homology_gene_stable_id','gene_stable_id ']].copy()

    selectgenes.to_csv('pd_orthofile.txt.tmp')
    orthofile = open("pd_orthofile.txt.tmp", "r")

    return orthofile

#Map the genes of human pathways to genes in the orthologs of the model organism
#basically if the orthofile have the genes from the pathway then it is ortholog
def mapOrtholog(pathgenes, orthofile):
  orthologs = open("pd_orthologs.txt.tmp", "w") #a new text file is created, this is based on the original grep when the mapping produce a text file then proceed to the next commandline step
  orthof = orthofile.readlines()
  for key in pathgenes:
    for line in orthof:
        if key in line:
          orthologs.write(line+'\n')
  orthologs.close()
  ortholog = open("pd_orthologs.txt.tmp", "r")
  return ortholog

#Small editing for the new ortholog file created
#Because in the next step it is unable to read the columns, due to empty lines in the file
#therefore this step is necessary
def removeEmptylines(ortholog):
    ortholog_clean = open("pd_ortholog_clean.txt.tmp", "w")

    for line in ortholog:
        if not line.strip(): continue  # skip the empty line
        ortholog_clean.write(line)
    ortholog_clean.close()
    df_ortholog_read = pd.read_csv("pd_ortholog_clean.txt.tmp", header=None)
    df_ortholog_raw = df_ortholog_read.iloc[:,[1,2]].copy()
    df_ortholog = df_ortholog_raw.drop_duplicates()

    return df_ortholog

#Select the genes of the organism only from the ortholog file, leaving the human gene ID
#also combine with the new orthologs
def Readgenes(df_ortholog, df_new_orthos, organism):
  orthologs_combined = pd.DataFrame()
  genes_of_interest = df_ortholog.iloc[:,0].copy()
  genes_of_interest = genes_of_interest.to_frame()
  genes_of_interest.rename(columns={genes_of_interest.columns[0]: 'Genes'}, inplace=True)

  if df_new_orthos.empty:
    orthologs_combined = genes_of_interest

  else:
    orthologs_combined = pd.merge(genes_of_interest, df_new_orthos, on=['Genes'], how='outer').copy()
    orthologs_combined.to_csv(organism+'pd_genes_combined.txt')

  return orthologs_combined

#Open and extract the phenotype data from each organism
#Only celegans dont have header name
#details/documentation on the header is in: http://geneontology.org/docs/go-annotation-file-gaf-format-2.1/
#so they follow GAF format, most-likely wont change, and consistent on their releases
#on Mouse, since there's multiple phenotypes I have to do a little editing
#to combine the phenotypes then split to different rows/entry
def Phenotypes(organism):

  if organism == "celegans":

    phcelegans = pd.read_csv(C_ele_phenotype_databse, skiprows=3, delimiter='\t', header=None)
    phcelegans2 = phcelegans[phcelegans.iloc[:,3] != 'NOT'].copy()
    phenotype = phcelegans2.iloc[:,[1,4]].copy()


  elif organism == "zebrafish":

    phzfish = pd.read_csv(Zebrafish_phenotype_databse, skiprows=1,delimiter='\t',header='infer')
    cols = ['Affected Structure or Process 1 subterm ID', 'Affected Structure or Process 1 superterm ID']
    phzfish['affected_id'] = phzfish[cols].apply(lambda row: ','.join(row.values.astype(str)), axis=1)
    phzfish2 = pd.DataFrame(phzfish.affected_id.str.split(',').values.tolist(), index=phzfish['Gene Symbol']).stack()
    phzfish2 = phzfish2.reset_index([0, 'Gene Symbol'])
    phzfish2.columns = ['Gene Symbol', 'Affected Phenotype ID']
    phzfish3 = phzfish2[(phzfish2!='nan').all(1)]
    phenotype = phzfish3.drop_duplicates()


  elif organism == "mouse":
    phmouse = pd.read_csv(Mouse_phenotype_databse, header='infer')
    cols = ['top_level_mp_term_id', 'mp_term_id']
    phmouse['all_mp'] = phmouse[cols].apply(lambda row: ','.join(row.values.astype(str)), axis=1)
    phenotype = pd.DataFrame(phmouse.all_mp.str.split(',').values.tolist(), index=phmouse.marker_symbol).stack()
    phenotype = phenotype.reset_index([0, 'marker_symbol'])
    phenotype.columns = ['marker_symbol', 'split_phen']


  elif organism == "slimemould":
    phslime = pd.read_csv(Slime_mould_phenotype_databse, delimiter='\t', header='infer')
    phenotype = pd.DataFrame(phslime['DDB_G_ID'].str.split('|').values.tolist(), index=phslime.Phenotypes).stack()
    phenotype = phenotype.reset_index([0, 'Phenotypes'])
    phenotype.columns = ['DDB_G_ID','Phenotypes']
    phenotype.columns = ['Phenotypes','DDB_G_ID']
    phenotype = pd.DataFrame(phenotype['Phenotypes'].str.split('|', expand=True).values.tolist(), index=phenotype['DDB_G_ID']).stack()
    phenotype = phenotype.reset_index([0, 'DDB_G_ID'])
    phenotype.columns = ['DDB_G_ID','Phenotypes']

  mapping = {phenotype.columns[0]:'Genes', phenotype.columns[1]: 'Phenotypes'}
  phenotype.rename(columns=mapping, inplace=True)
  phenotype = phenotype.drop_duplicates()

  return phenotype

#the original panda dataframe (before selected columns) is saved for onthology/term in later parts
#gammln, logfact, hypergeom unchanged from previous code
def gammln(xx):
    cof = [76.18009172947146, -86.50532032941677,
    24.01409824083091, -1.231739572450155,
    0.12086509738661e-2, -0.5395239384953e-5]
    y = x = xx
    tmp = x + 5.5
    lgtmp = math.log(tmp)
    c = (x + .5) * lgtmp
    tmp -= c
    ser = 1.000000000190015
    for j in range(6):
        y += 1
        ser += cof[j]/y
    return -tmp + math.log(2.5066282746310005*ser/x)

def logfact(arg):
    return gammln(arg + 1.0)

#Hypergeometric distribution
#I don't remove prof's comment here
def hypergeom(n,m,N,i):
    # There are m "bad" and n "good" balls in an urn.
    # Pick N of them. The probability of i or more successful selection

    # (m!n!N!(m+n-N)!)/(i!(n-i)!(m+i-N)!(N-i)!(m+n)!)
    loghyp1 = logfact(m)+logfact(n)+logfact(N)+logfact(m+n-N);
    loghyp2 = logfact(i)+logfact(n-i)+logfact(m+i-N)+logfact(N-i)+logfact(m+n);
    return math.exp(loghyp1 - loghyp2);

#The enrichment.
#Since pandas keep duplicate values, unlike dictionaries that only keep unique keys, thats why I have to do a lot of formatting
#to make the column values unique, hence looks a little messy
def Enrichment(organism, genes, phenotypes):
    count_phenotypes_path = 0 #this doesn't used for some reason, might delete it (consult prof)
    count_annotated_genes = 0
    #to see the annotated genes
    annotatedg = pd.merge(genes, phenotypes, on=['Genes'], how='inner').copy()
    uniquegenes = annotatedg.drop_duplicates(subset=['Genes'],keep ='first', inplace = False)
    count_annotated_genes = len(uniquegenes.index)
    #annotated genes ends here
    uniquerows = annotatedg.drop_duplicates()
    #overlap genes starts here
    uniquerow_agg = uniquerows.groupby('Phenotypes').agg(lambda x: x.tolist())
    overlapgenes = uniquerow_agg.reset_index()
    overlapgenes.rename(columns={overlapgenes.columns[0]: 'Enriched Phenotype'}, inplace=True)
    #Overlapgenes ends here
    #count the phenotypes starts here
    countphen = uniquerows.iloc[:,1].copy()
    count_phenotypes = countphen.value_counts()
    count_phenotypes = pd.DataFrame(count_phenotypes).reset_index()
    count_phenotypes.columns = ['Phenotypes', 'count']
    #count phenotypes ends here
    #to count n starts here
    uniquephenes = annotatedg.drop_duplicates(subset=['Phenotypes'],keep ='first', inplace = False)
    countphenframe = uniquephenes[['Phenotypes']].copy()
    nraw = pd.merge(countphenframe, phenotypes, on=['Phenotypes'], how='inner').copy()
    uniquenraw = nraw.drop_duplicates()
    ncountphen = uniquenraw.iloc[:,0].copy()
    count_n = ncountphen.value_counts()
    count_n = pd.DataFrame(count_n).reset_index()
    count_n.columns = ['Phenotypes', 'count']
    #count n ends here
    #to combine both n and count phenotypes for the iteration
    counts = pd.merge(count_phenotypes, count_n, on=['Phenotypes'], how='inner').copy()
    N = count_annotated_genes
    #start modify dataframe needed for m
    uniquegenesdb = phenotypes.drop_duplicates(subset=['Genes'],keep ='first', inplace = False)
    #end needed for m
    #to make the outputs
    pval = []
    nval = []
    mval = []
    ival = []
    e_phenotype = []
    #start the iteration for hypergeometric
    for tuplephen in counts.itertuples():
      i = tuplephen[2]
      n = tuplephen[3]
      m = len(uniquegenesdb.index) - n
      hyp = hypergeom(n, m, N, i)
      pval.append(hyp)
      nval.append(n)
      mval.append(m)
      ival.append(i)
      e_phenotype.append(tuplephen[1])
    df_enrichment = pd.DataFrame(list(zip(e_phenotype, pval, nval, mval, ival)), columns=['Enriched Phenotype','P-value', 'n', 'm','i'])
    df_enrichment = df_enrichment.sort_values(by ='P-value' )
    #FDR starts here
    df_enrichment['rank'] = range(1, 1+len(df_enrichment))
    M = len(counts)
    Q = 0.001
    qval = []

    for tupleBH in df_enrichment.itertuples():

      i = tupleBH[6]
      q = (i/M)*Q
      qval.append(q)
    df_enrichment['q-value'] = qval
    sigenrichment = df_enrichment[df_enrichment['P-value'] < df_enrichment['q-value']]

    if sigenrichment.empty:
      print('No enriched phenotypes!')

    return sigenrichment, overlapgenes

#Open the onthology database, except for slime mould since the phenotype database already has phenotype term
#For zebrafish there's 2 terms the ZFA and GO
def openOnthology(organism):
  onthology_zfa = pd.DataFrame()

  if organism == "celegans":
    onthology = pd.read_csv(C_ele_ontology_data_file, sep=";", names=['Ontology'])

  elif organism == "zebrafish":
    onthology_zfa_raw = pd.read_csv(Zebrafish_ontology_data_file, skiprows=1, delimiter='\t',header='infer')
    onthology_zfa = onthology_zfa_raw [['Anatomy ID','Anatomy Name']].copy()
    onthology = pd.read_csv(Zebrafish_GO_annotation_file, sep=";", names=['Ontology'])

  elif organism == "mouse":
    onthology = pd.read_csv(Mouse_ontology_data_file, sep=";", names=['Ontology'])

  return onthology, onthology_zfa

#parsing the database since most of them are obo format
def searchOnthology(organism, sigenrichment, overlapgenes, onthology, onthology_zfa):

  onthology.insert(0, 'Enriched Phenotype', onthology['Ontology'].str.extract(r'(id:\s+\S+.*)', expand=False).ffill())
  onthology.insert(1, 'Phenotype Name', onthology['Ontology'].str.extract(r'(name:\s+\S+.*)', expand=False).ffill())
  onthology['Enriched Phenotype'] = onthology['Enriched Phenotype'].str.replace(r'id:\s+', '')
  onthology['Phenotype Name'] = onthology['Phenotype Name'].str.replace(r'name:\s+', '')
  onthology_ph = onthology[['Enriched Phenotype','Phenotype Name']].copy()
  onthologyunique_raw = onthology_ph.drop_duplicates().dropna()
  onthologyunique = onthologyunique_raw.drop_duplicates(subset=['Enriched Phenotype'],keep ='last', inplace = False)

  if organism == "zebrafish":
    onthologyunique.columns = ['Enriched Phenotype', 'Affected Term']
    enrichedOnthology1 = pd.merge(onthologyunique, sigenrichment, on=['Enriched Phenotype'], how='inner').copy()
    onthology_zfa.columns = ['Enriched Phenotype', 'Affected Term']
    enrichedOnthology2 = pd.merge(onthology_zfa, sigenrichment, on=['Enriched Phenotype'], how='inner').copy()
    enrichedOnthologyComplete1 = pd.concat([enrichedOnthology1, enrichedOnthology2])
    enrichedOnthologyComplete = pd.merge(enrichedOnthologyComplete1, overlapgenes, on=['Enriched Phenotype'], how='inner').copy()

  else:
    enrichedOnthology1 = pd.merge(onthologyunique, sigenrichment, on=['Enriched Phenotype'], how='inner').copy()
    enrichedOnthologyComplete = pd.merge(enrichedOnthology1, overlapgenes, on=['Enriched Phenotype'], how='inner').copy()
  enrichedOnthologyComplete = enrichedOnthologyComplete.sort_values(by ='rank' )
  enrichedOnthologyComplete['Overlap Genes'] = [','.join(map(str, l)) for l in enrichedOnthologyComplete['Genes']]
  enrichedOnthologyFinal = enrichedOnthologyComplete.drop(columns=['Genes'])

  print(enrichedOnthologyFinal)
  return enrichedOnthologyFinal

#To run this search for the onthology terms, since slimemould isn't necessary then it passed
def onthologyStep(organism, sigenrichment, overlapgenes):
  enrichedOnthologyFinal = pd.DataFrame()

  if organism == "slimemould":
    pass

  else:
    onthology, onthology_zfa = openOnthology(organism)
    enrichedOnthologyFinal = searchOnthology(organism, sigenrichment,overlapgenes, onthology, onthology_zfa)

  return enrichedOnthologyFinal

#Since the results will be multiple it's better to place it in one folder
#thus the enrinchment results folder is created
def createResultFolder():
  directory = pathID+"_Enrichment_Results"
  path = os.path.join(parent_dir, directory)
  #if the folder is already made (since one run is in turn per organism first celegans then mouse etc) then it passed
  if os.path.exists(pathID+"_Enrichment_Results"):
    pass
  else:
    os.mkdir(path)
  #change the current directory to the folder
  os.chdir(pathID+"_Enrichment_Results")

#Add conservation information, as suggested by Prof. van Noort for shiny
#basically count the orthologs vs the total genes in that human pathway
#Add this info to the csv file with the result of enrichment
def addInfo(enrichedOnthologyFinal, genes, pathgenes, organism, sigenrichment):
  totorth = len(genes.index)
  totorth2 = str(totorth)
  pathlength = len(pathgenes)
  pathlength2 = str(pathlength)
  info = open(pathID+'_'+organism+'_Enrichment_Result.txt', 'a')
  info.write('orthologs: '+totorth2+'\n')
  info.write('totalgenes: '+pathlength2+'\n')
  info.write('organism: '+organism+'\n')
  if organism == "slimemould":
    sigenrichment.to_csv(info)
  else:
    enrichedOnthologyFinal.to_csv(info)
  info.close()

#The organism is already here!! so no need to put it in the command line anymore
organism_list = ["celegans", "zebrafish", "mouse", "slimemould"]

#run the functions that don't need the organism
reactomefile = openReactome()
pathgenes = selectPathway(reactomefile)
#this parent_dir is like path_root. make sure to run it in the pathroot so they're the same
#I think getcwd is easier than having to type it
parent_dir = os.getcwd()

#All functions regarding the enrichment that needs the organism as variable is here
def runEnrichment(organism_list):
  for organism in organism_list:
    df_new_ortholog = loadNewOrtholog(organism)
    orthofile = selectForMapping(organism)
    ortholog = mapOrtholog(pathgenes, orthofile)
    df_ortholog = removeEmptylines(ortholog)
    genes = Readgenes(df_ortholog, df_new_ortholog, organism)
    phenotypes = Phenotypes(organism)
    sigenrichment, overlapgenes = Enrichment(organism, genes, phenotypes)
    enrichedOnthologyFinal = onthologyStep(organism, sigenrichment, overlapgenes)
    createResultFolder()
    addInfo(enrichedOnthologyFinal, genes, pathgenes, organism, sigenrichment)
    #before continuing to other organism, change back the directory back to the root/parent
    os.chdir(parent_dir)

#Run the enrichment here.
runEnrichment(organism_list)

####Next section is to run the new summary conservation feature as suggested by Prof
#change first the dir to the results since we want the output to be there
os.chdir(pathID+"_Enrichment_Results")

#the summary function
#First it will import all the results of the files
#then it reads the conservation information as added before
#a new file is created with the summary conservation of the species analyzed
def runSummary():
  extension = 'txt'
  all_filenames = [i for i in glob.glob('*.{}'.format(extension))]
  df = pd.concat([pd.read_csv(f, sep=";", names=['Result']) for f in all_filenames ])
  #export to csv
  df.to_csv("combined_result.txt", index=False)
  df.insert(0, 'Orthologs', df['Result'].str.extract(r'(orthologs:\s+\S+.*)', expand=False).ffill())
  df.insert(1, 'Total Genes in Pathway', df['Result'].str.extract(r'(totalgenes:\s+\S+.*)', expand=False).ffill())
  df.insert(2, 'Species', df['Result'].str.extract(r'(organism:\s+\S+.*)', expand=False).ffill())
  df['Orthologs'] = df['Orthologs'].str.replace(r'orthologs:\s+', '')
  df['Total Genes in Pathway'] = df['Total Genes in Pathway'].str.replace(r'totalgenes:\s+', '')
  df['Species'] = df['Species'].str.replace(r'organism:\s+', '')
  dfph = df[['Orthologs','Total Genes in Pathway', 'Species']].copy()
  dfphunique = dfph.drop_duplicates().dropna()
  dfphunique2 = dfphunique.drop_duplicates(subset=['Species'],keep ='first', inplace = False)
  print(dfphunique2)
  dfphunique2.to_csv(phenotype_enrichment_result_file, index=False)

runSummary()
