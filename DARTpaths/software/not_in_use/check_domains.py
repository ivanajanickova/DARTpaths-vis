#!/usr/bin/python3
# Script which checks whether multiple domains excist on an hmm-profile and if so, runs a rigorous analysis using the characteristics of the original proteins. If the profile contains only a single domain: checks for a correct overlap of the sequence and the profile. 

import sys
import re
from myfunc import *

#Input: -protein, -minimum overlap of total profile, maximum ofset of domain-midpoint vs. original proteins, ratio for minimum count of domains vs. original proteins, ratio for maximum count of domains vs. original proteins
prot = sys.argv[1]
min_overlap = float(sys.argv[2])
max_offset = float(sys.argv[3])
min_count = float(sys.argv[4])
max_count = float(sys.argv[5])

# Open files needed
with open('{}_hmms/hmmer_search_{}fam_scores_temp2.tsv'.format(prot,prot), 'r') as fhandle:
	scores = fhandle.read()
	scores = scores.strip()
with open('{}_searches_literature/{}_reactome.fa'.format(prot,prot), 'r') as fhandle:
	original = fhandle.read()
	original = original.strip()

# First, get original proteins
original = original.strip()
lines = original.split('\n')
original_prots=[]
for line in lines:
	if '>' in line:
		[_,original_prot] = line.split('.')
		[original_prot,_,_] = original_prot.split(' ')
		original_prots.append(original_prot)

# We need to establish whether the original proteins have multiple domains, or a single one.
hard_way = False
for original_prot in original_prots:
	domains = find_match(scores,original_prot)
	# If we find a protein which has multiple domains, set hard_way and continue
	if len(domains) >= 2:
		hard_way = True
		continue

# If there is only a single domain check each score for overlap-ratio
new_score_list=[]
if not hard_way:
	for score in scores.split('\n'):
		score_tabs = score.split('\t')
		# Check if the score length is higher than profile-length times overlap ratio
		if int(score_tabs[5])-int(score_tabs[4])+1>int(score_tabs[2])*min_overlap:
			new_score_list.append(score)

# If there are multiple domains, extract each domains in a dictionary, then loop over all proteins and check whether their domains match.
# Put those proteins in a new list.
if hard_way:
	# Loop over the original proteins
	for original_prot in original_prots:
		# Find the different domains, retreive the edges of each domain in a list
		domains = find_match(scores,original_prot)
		domain_edges = []	
		for domain in domains:
			domain_tabs = domain.split('\t')
			domain_edges.append((int(domain_tabs[4]),int(domain_tabs[5])))	
		# As a single protein might have repeating patterns it can have multiple domains which match the same part of the hmm-profile. Unique domains need to be sepparated and repeating domains 			need to be grouped. 
		# Loop over these domains, put them in a dictionary: {unique_domain:counts}.
		for i,domain_edge in enumerate(domain_edges):
			# First domain starts dict
			if i==0:		
				domdict = {domain_edge:0}
			# Check if domain is already in domain dict (thus non-unique), then increase counter of that domain
			is_there, match_domain = is_domain(domdict, domain_edge, max_offset)
			if is_there:
				domdict[match_domain] += 1
			# If its not the same, create a new dict entry of and add it
			else:
				domdict.update({domain_edge:1})
		print('\nDomain dict of {} is: {}'.format(original_prot,domdict))	
		# Now we have the domain-profile dictionary, next we loop over all proteins and look for matches
		done_prots = []	
		for score in scores.split('\n'):
			score_tabs = score.split('\t')
			prot_name = score_tabs[0]
			# Check if protein is already done (multiple domains means every protein has multiple entries, we use all entries at once but the loop won't know this)
			if prot_name in done_prots:
				continue
			# If not, gather scores of protein
			prot_scores = find_match(scores,prot_name)
			# Now check if the domains match the ones of the original protein
			# First, gather these domains in a similar fashion as the original dict
			score_dict = {}
			for prot_score in prot_scores:
				prot_score_tabs = prot_score.split('\t')
				domain_edge = (int(prot_score_tabs[4]),int(prot_score_tabs[5]))
				is_match, match_domain = is_domain(domdict, domain_edge, max_offset)
				if is_match:
					# Update if it's already found
					if match_domain in score_dict:
						score_dict[match_domain] += 1
					# Add if not
					else:
						score_dict.update({match_domain:1})
	
			# Now we have a score_dict, which we need to compair with the original protein dict
			# We assume a match untill a domain is either absent, under- or over-respresented
			match = True		
			for domain, orig_count in domdict.items():
				# Check is domain is found, if not: break
				if domain not in score_dict:
					match = False
					break
				# If its present, check the count, if it does not add up: break
				prot_count = score_dict[domain]
				if (prot_count < orig_count*min_count) or (prot_count > orig_count*max_count):
					match = False
					break		
			
			# If it "survived" the comparisson, add the scores to the new file
			if match:
				new_score_list.extend(prot_scores)
			
			# Finally,put this protein in the 'done_prots' list so it is skipped from now on	
			done_prots.append(prot_name)

# If everything went right, we now have a 'new_score_list' with only proteins which match our criteria.
# However, there might be doubles, as proteins are added multiple times if they match with multiple original proteins. 
# This does not matter as after this only unique protein names are extracted.
# Save our list
with open('{}_hmms/hmmer_search_{}fam_sign_scores.tsv'.format(prot,prot), 'w') as fhandle:
	fhandle.write('\n'.join(new_score_list))
