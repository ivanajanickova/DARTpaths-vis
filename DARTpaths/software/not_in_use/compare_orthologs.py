#!/usr/bin/python3

import sys

# Compare the orthologs from different trees and returns a file with 
# consistent and non-consistent orthologs sepperated

#Input: -current dir, -name of protein, -number of tables

path = sys.argv[1]
prot = sys.argv[2]
nr_of_tables = int(sys.argv[3])

#Read all tables in as a list of strings
tables = []
for nr in range(nr_of_tables):
    with open('{}/{}_orthologs/{}_human_orthologs.{}'.format(path,prot,prot,nr),'r') as curtable:
        tables.append(curtable.read())

#Take first string(table) as readout for a new list
lines = tables[0].split('\n\n')
human_genes = []
for line in lines:
    line.strip()
    if line == '':
        continue
    tabs = line.split('\t')
    human_genes.append(tabs[0])

# Now look through the tables per human_gene and make a temporary dict for each gene with 
# counts of its orthologs,

# at the end write those orthologs and counts to a list
final_orthos = []
for gene in human_genes:
    temp_dic = {}
    for table in tables:
        lines = table.split('\n\n')
        for line in lines:
            line.strip()
            if line == '':
                continue
            tabs = line.split('\t')
            if tabs[0] != gene:
                continue
            orthologs = tabs[1].split(', ')
            for ortho in orthologs:
                if ortho in temp_dic:
                    temp_dic[ortho] += 1
                else:
                    temp_dic[ortho] = 1
    c_ortho = []
    nc_ortho = []
    for ortho in temp_dic.keys():
        if temp_dic.get(ortho) == nr_of_tables:
            c_ortho.append(ortho)
        else:
            nc_ortho.append('{} ({})'.format(ortho,temp_dic.get(ortho)))
    final_orthos.append('{}\t{}\t{}'.format(gene,', '.join(sorted(c_ortho)),', '.join(sorted(nc_ortho))))

#Write list in file and save
with open('{}/{}_orthologs/{}_common_orthologs'.format(path,prot,prot),'w') as output_file:
    output_file.write('\n\n'.join(final_orthos))

