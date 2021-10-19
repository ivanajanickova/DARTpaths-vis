#!/usr/bin/python3

import sys
import os

#Script which uses the Notung stat-output of multiple trees to find the trees with the least duplications

#Input: -current dir, -name of protein

path=sys.argv[1]
prot=sys.argv[2]

#Look up all stat-files
tree_names = []
for f in os.listdir('{}/{}_trees'.format(path,prot)):
	if f.endswith('stats.txt'):
		tree_names.append(f)

#Check if there are multiple trees in the first plact, if not return treenum 0 by default
if tree_names == []:
	print('ERROR: No stats files found')
	sys.exit(0)
elif len(tree_names) == 1:
	sys.exit('0')

#If multiple trees exist, loop over themm, loading their stats and pulling the number of duplications from them, then putting it in a dict. 

tree_dict = {}
for tree in tree_names:
	with open('{}/{}_trees/{}'.format(path,prot,tree),'r') as stat:
		stats = stat.read()
	for line in stats.split('\n'):
		if 'Duplications' in line:
			dupes = line[16:]
			break
	[name,_] = tree.split('.stats')
	tree_nr = name[-1]
	tree_dict.update({tree_nr:dupes})

#Check dict for lowest amount of dupes, return the key which has the least
lowest_dupes = min(tree_dict, key=tree_dict.get)
sys.exit(lowest_dupes)
