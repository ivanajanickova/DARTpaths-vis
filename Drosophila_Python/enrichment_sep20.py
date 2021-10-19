#Enrichment v60
#command: python3 -.py R-HSA-69239 

##########################
#import necessary modules#
########################## 

import glob
import math
import os
import re
import sys  # for user input

import numpy as np
import pandas as pd


pd.options.mode.chained_assignment = None

############
#the inputs#
############
pathID = sys.argv[1] #the human pathway that wants to be enriched
organism = sys.argv[2] #the organism to map which orthologous genes and which has enriched phenotypes  

#############################################################
#Opening, Extracting, and Mapping Information from Resources#
#############################################################

#open the reactome resource
def openReactome():
  reactomefile = open("Ensembl2Reactome_All_Levels.txt", "r")
  return reactomefile

#load new orthologs


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
  organism = sys.argv[2]
  if organism == "celegans":
    df_celegans = pd.read_csv("Orthology_human_celegans_ensembl101_unique.txt",delimiter='\t',header = 'infer')
    selectgenes= df_celegans[['Caenorhabditis elegans gene stable ID','Gene stable ID', 'Caenorhabditis elegans homology type']].copy()
  elif organism == "zebrafish":
    df_zebrafish = pd.read_csv("Orthology_human_zebrafish_ensembl101_unique.txt", delimiter='\t',header = 'infer')
    selectgenes= df_zebrafish[['Zebrafish gene name','Gene stable ID', 'Zebrafish homology type']].copy()
  elif organism == "mouse":
    df_mouse = pd.read_csv("Orthology_human_mouse_ensembl101_unique.txt",delimiter='\t', header = 'infer')
    selectgenes= df_mouse[['Mouse gene name','Gene stable ID', 'Mouse homology type']].copy()
  elif organism == "slimemould":
    df_slime = pd.read_csv("Orthology_human_dicty_ensembl101_unique.txt",delimiter='\t', header='infer')
    selectgenes= df_slime[['homology_gene_stable_id','gene_stable_id ', 'Human homology type']].copy()
  elif organism == "dmelanogaster":
    df_dmelanogaster = pd.read_csv("Orthology_human_dmelanogaster_ensembl101_unique.txt",delimiter='\t', header='infer')
    selectgenes = df_dmelanogaster[['Drosophila melanogaster gene name','Gene stable ID', 'Drosophila melanogaster homology type']].copy()
  selectgenesfix =  selectgenes[(selectgenes.iloc[:,2] == 'ortholog_one2one') |(selectgenes.iloc[:,2] == 'ortholog_one2many') ]   
  selectgenesfix.to_csv('pd_orthofile.txt')  
  orthofile = open("pd_orthofile.txt", "r")
  return orthofile

#Map the genes of human pathways to genes in the orthologs of the model organism 
#basically if the orthofile have the genes from the pathway then it is ortholog
def mapOrtholog(pathgenes, orthofile):
  orthologs = open("pd_orthologs.txt", "w") #a new text file is created, this is based on the original grep when the mapping produce a text file then proceed to the next commandline step
  orthof = orthofile.readlines()
  for key in pathgenes:
    for line in orthof:
        if key in line:
          orthologs.write(line+'\n')
  orthologs.close()
  ortholog = open("pd_orthologs.txt", "r")  
  return ortholog

#Small editing for the new ortholog file created
#Because in the next step it is unable to read the columns, due to empty lines in the file
#therefore this step is necessary
def removeEmptylines(ortholog):
    ortholog_clean = open("pd_ortholog_clean.txt", "w")
    for line in ortholog:
        if not line.strip(): continue  # skip the empty line
        ortholog_clean.write(line)
    ortholog_clean.close()
    df_ortholog = pd.DataFrame()
    try: 
      df_ortholog_read = pd.read_csv("pd_ortholog_clean.txt", header=None)
      df_ortholog_raw = df_ortholog_read.iloc[:,[1,2]].copy()
      df_ortholog = df_ortholog_raw.drop_duplicates()
    except:
      pass  
    return df_ortholog

#Select the genes of the organism only from the ortholog file, leaving the human gene ID   
def Readgenes(df_ortholog):
  genes_of_interest = pd.DataFrame()
  if df_ortholog.empty:
    genes_of_interest['Genes'] = ""
  else:
    genes_of_interest = df_ortholog.iloc[:,0].copy()
    genes_of_interest = genes_of_interest.to_frame()
    genes_of_interest.rename(columns={genes_of_interest.columns[0]: 'Genes'}, inplace=True)
  return genes_of_interest

#Open and extract the phenotype data from each organism
#Only celegans dont have header name
#details/documentation on the header is in: http://geneontology.org/docs/go-annotation-file-gaf-format-2.1/
#so they follow GAF format, most-likely wont change, and consistent on their releases
#on Mouse, since there's multiple phenotypes I have to do a little editing
#to combine the phenotypes then split to different rows/entry
def Phenotypes(organism): 
  organism = sys.argv[2]
  phmouse = pd.DataFrame()
  phzfish = pd.DataFrame()
  if organism == "celegans":
    phcelegans = pd.read_csv("phenotype_association.WS264.wb",skiprows=3, delimiter='\t', header=None)
    phcelegans2 = phcelegans[phcelegans.iloc[:,3] != 'NOT'].copy()
    phenotype = phcelegans2.iloc[:,[1,4]].copy()
  elif organism == "zebrafish":
    phzfish = pd.read_csv("phenoGeneCleanData_fish_2020.09.25.txt",skiprows=1,delimiter='\t',header='infer')
    cols = ['Affected Structure or Process 1 subterm ID', 'Affected Structure or Process 1 superterm ID']
    phzfish['affected_id'] = phzfish[cols].apply(lambda row: ','.join(row.values.astype(str)), axis=1)
    phzfish2 = pd.DataFrame(phzfish.affected_id.str.split(',').values.tolist(), index=phzfish['Gene Symbol']).stack()
    phzfish2 = phzfish2.reset_index([0, 'Gene Symbol'])
    phzfish2.columns = ['Gene Symbol', 'Affected Phenotype ID']
    phzfish3 = phzfish2[(phzfish2!='nan').all(1)]
    phenotype = phzfish3.drop_duplicates()
  elif organism == "mouse":
    phmouse = pd.read_csv("ALL_genotype_phenotype.csv",header='infer')
    cols = ['top_level_mp_term_id', 'mp_term_id']
    phmouse['all_mp'] = phmouse[cols].apply(lambda row: ','.join(row.values.astype(str)), axis=1)
    phenotype1 = pd.DataFrame(phmouse.all_mp.str.split(',').values.tolist(), index=phmouse.marker_symbol).stack()
    phenotype1 = phenotype1.reset_index([0, 'marker_symbol'])
    phenotype1.columns = ['marker_symbol', 'split_phen']
    phenotype2 = phenotype1[(phenotype1!='nan').all(1)]
    phenotype = phenotype2.drop_duplicates()
  elif organism == "slimemould":
    phslime = pd.read_csv("all-mutants-ddb_g.txt", delimiter='\t', header='infer')
    phenotype = pd.DataFrame(phslime['DDB_G_ID'].str.split('|').values.tolist(), index=phslime.Phenotypes).stack()
    phenotype = phenotype.reset_index([0, 'Phenotypes'])
    phenotype.columns = ['DDB_G_ID','Phenotypes']
    phenotype.columns = ['Phenotypes','DDB_G_ID']
    phenotype = pd.DataFrame(phenotype['Phenotypes'].str.split('|', expand=True).values.tolist(), index=phenotype['DDB_G_ID']).stack()
    phenotype = phenotype.reset_index([0, 'DDB_G_ID'])
    phenotype.columns = ['DDB_G_ID','Phenotypes']
  elif organism == "dmelanogaster":
    col_names = ['AlleleID', 'AlleleSymbol', 'GeneID', 'GeneSymbol']
    dmelanogaster_gene2alleles = pd.read_csv("fbal_to_fbgn_fb_FB2021_01.tsv", delimiter='\t', dtype= str, names=col_names)
    dmelanogaster_alleles2phenotypes = pd.read_csv("allele_phenotypic_data_fb_2021_01.tsv", delimiter='\t', header=3)
    phdmelanogaster = pd.merge(left=dmelanogaster_gene2alleles, right=dmelanogaster_alleles2phenotypes, left_on = ["AlleleID", "AlleleSymbol"], right_on=["allele_FBal#", "##allele_symbol"], how="right")
    phenotype = phdmelanogaster.iloc[:, lambda df: [2, 6]]
  mapping = {phenotype.columns[0]:'Genes', phenotype.columns[1]: 'Phenotypes'}
  phenotype.rename(columns=mapping, inplace=True)
  phenotype = phenotype.drop_duplicates().dropna()
  #print(phenotype)
  return phenotype, phmouse, phzfish

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
def Enrichment(organism, genes, phenotypes, phmouse, phzfish):
    organism = sys.argv[2]
    count_phenotypes_path = 0 #this doesn't used for some reason, might delete it (consult prof)
    count_annotated_genes = 0
    #to see the annotated genes
    annotatedg = pd.merge(genes, phenotypes, on=['Genes'], how='inner').copy()
    uniquegenes = annotatedg.drop_duplicates(subset=['Genes'],keep ='first', inplace = False)
    count_annotated_genes = len(uniquegenes.index)
    #annotated genes ends here
    uniquerows = annotatedg.drop_duplicates()
    #print(uniquerows)
    #overlap genes starts here
    if organism == "mouse":
      musmgi = phmouse[['marker_accession_id','marker_symbol']].copy()
      musmgi.rename(columns={musmgi.columns[1]: 'Genes'}, inplace=True)
      uniquerows2 = pd.merge(uniquerows, musmgi, on=['Genes'], how='inner').copy() #lets hope this works
      uniquerows2 = uniquerows2.drop(columns=['Genes'])
      uniquerows2.rename(columns={uniquerows2.columns[1]: 'Genes'}, inplace=True)
      uniquerows3 = uniquerows2.drop_duplicates()
      uniquerow_agg = uniquerows3.groupby('Phenotypes').agg(lambda x: x.tolist())
      #print(uniquerows)
    elif organism == "zebrafish":
      zfzdb = phzfish[['Gene ID','Gene Symbol']].copy()
      zfzdb.rename(columns={zfzdb.columns[1]: 'Genes'}, inplace=True)
      uniquerows2 = pd.merge(uniquerows, zfzdb, on=['Genes'], how='inner').copy() #lets hope this works
      uniquerows2 = uniquerows2.drop(columns=['Genes'])
      uniquerows2.rename(columns={uniquerows2.columns[1]: 'Genes'}, inplace=True)
      uniquerows3 = uniquerows2.drop_duplicates()
      uniquerow_agg = uniquerows3.groupby('Phenotypes').agg(lambda x: x.tolist())
      #print(uniquerows)
    else:
      uniquerow_agg = uniquerows.groupby('Phenotypes').agg(lambda x: x.tolist())
    overlapgenes = uniquerow_agg.reset_index()    
    overlapgenes.rename(columns={overlapgenes.columns[0]: 'Enriched Phenotype'}, inplace=True)
    #Overlapgenes ends here
    #count the phenotypes starts here
    #print(annotatedg)
    #print(uniquerows)
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
    #print(counts)
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
    Q = 0.1
    qval = []
    for tupleBH in df_enrichment.itertuples():
      i = tupleBH[6]
      q = (i/M)*Q
      qval.append(q)
    df_enrichment['q-value'] = qval
    sigenrichmentraw = df_enrichment[df_enrichment['P-value'] < df_enrichment['q-value']] 
    pvals = sigenrichmentraw['P-value'] 
    maxpvals1 = pvals.max()
    maxpvals = float(maxpvals1)
    sigenrichment = df_enrichment[df_enrichment['P-value'] <= maxpvals]
    if sigenrichment.empty:
      print('No enriched phenotypes!')
    if organism == 'slimemould':
       if sigenrichment.empty:
        pass
       else: 
        sigenrichment = pd.merge(sigenrichment, overlapgenes, on=['Enriched Phenotype'], how='inner').copy() 
        sigenrichment['Overlap Genes'] = [','.join(map(str, l)) for l in sigenrichment['Genes']]
        sigenrichment = sigenrichment.drop(columns=['Genes'])
        print(sigenrichment)
    return sigenrichment, overlapgenes

def openOnthology(organism):
  organism = sys.argv[2]
  onthology_zfa = pd.DataFrame()
  if organism == "celegans":
    onthology = pd.read_csv('phenotype_ontology.WS264.obo', sep=";", names=['Ontology'])
  elif organism == "zebrafish":
    onthology_zfa_raw = pd.read_csv('anatomy_item_2020.10.14.txt', skiprows=1, delimiter='\t',header='infer')
    onthology_zfa = onthology_zfa_raw [['Anatomy ID','Anatomy Name']].copy()
    onthology = pd.read_csv('go_annotation.obo', sep=";", names=['Ontology'])
  elif organism == "mouse":
    onthology = pd.read_csv('MPheno_OBO.ontology', sep=";", names=['Ontology'])
  elif organism == "dmelanogaster":
    onthology = pd.read_csv('go-basic.obo', sep=";", names=['Ontology'])
  return onthology, onthology_zfa

def searchOnthology(organism, sigenrichment, overlapgenes, onthology, onthology_zfa):
  organism = sys.argv[2]
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
  if enrichedOnthologyFinal.empty:
    pass
  else:
    print(enrichedOnthologyFinal)
  return enrichedOnthologyFinal

def onthologyStep(organism, sigenrichment, overlapgenes):
  organism = sys.argv[2]
  enrichedOnthologyFinal = pd.DataFrame()
  if organism == "slimemould":
    pass
  else:
    onthology, onthology_zfa = openOnthology(organism)
    enrichedOnthologyFinal = searchOnthology(organism, sigenrichment,overlapgenes, onthology, onthology_zfa)
  return enrichedOnthologyFinal
  
def createResultFolder():
  directory = pathID+"_Enrichment_Results"
  path = os.path.join(parent_dir, directory) 
  if os.path.exists(pathID+"_Enrichment_Results"):
    pass
  else:
    os.mkdir(path) 
  os.chdir(pathID+"_Enrichment_Results")

def addInfo(enrichedOnthologyFinal, genes, pathgenes, organism, sigenrichment):
  organism = sys.argv[2]
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

organism_list = ["celegans", "zebrafish", "mouse", "slimemould", "dmelanogaster"]

reactomefile = openReactome()
pathgenes = selectPathway(reactomefile)
parent_dir = os.getcwd()

def runEnrichment(organism_list):
  for organism in organism_list:
    orthofile = selectForMapping(organism)
    ortholog = mapOrtholog(pathgenes, orthofile)
    df_ortholog = removeEmptylines(ortholog)
    genes = Readgenes(df_ortholog)
    phenotypes, phmouse, phzfish = Phenotypes(organism)
    #sampai sini gpp
    sigenrichment, overlapgenes = Enrichment(organism, genes, phenotypes, phmouse, phzfish)
    enrichedOnthologyFinal = onthologyStep(organism, sigenrichment, overlapgenes)
    createResultFolder()
    addInfo(enrichedOnthologyFinal, genes, pathgenes, organism, sigenrichment)
    os.chdir(parent_dir)

runEnrichment(organism_list)

os.chdir(pathID+"_Enrichment_Results")

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
  dfphunique2.to_csv(f'{pathID}_summary_conservation.txt', index=False)

runSummary()
