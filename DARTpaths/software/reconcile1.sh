#!/bin/bash
set -e

# First batch script for constructing and reconciling a genetree for a protein-familie.
# Needs reconcile_scripts folder for individual tasks.
# Every task in itself has been set up so it can also be run individually.

# Input: -name of protein

prot=$1

# Add path to scripts
#PATH=$PATH:~/DARTpaths/${prot}/reconcile_scripts # original
PATH=$PATH:$HOME/Desktop/2020_0218_scripts/Mock_automatepipeline/DARTpaths/software/reconcile_scripts

#PATH=$PATH:$HOME/ncbi-blast-2.8.1+/bin
BLASTDB=$HOME/DARTpaths/software/blast_db

# Make required dirs if they do not exist yet
mkdir -p ${prot}_searches_blast
mkdir -p ${prot}_alignments
mkdir -p ${prot}_hmms

echo '----------Starting on '${prot}' protein familie----------'
nr_prot="$(grep -o \> ${prot}_searches_literature/${prot}_reactome.fa | wc -l)"
echo 'Found '${nr_prot}' proteins for blast query'

# Use blast to get similar proteins
echo '----------Executing blastp----------'
blastp_search_and_reports_multiple_queries ${prot}

# Create an alignment of the blast-proteins
echo '----------Creating alignment of blast fasta----------'
mafft_blast ${prot}

# Use this alignment to build a hmm profile
echo '----------Building hmm profile----------'
hmm_build ${prot}

# Search our own protein database (protein.sequences.usable.fa) with the hmm profile
echo '----------Searching database with hmm profile----------'
hmm_search ${prot}

# Extract the sequences with a significant score (check file for parameters)
echo '----------Extracting fasta seqs with significant score---------------'
get_sign_seqs ${prot}

# Create a new alignment using the hmm-proteins
echo '----------Creating alignment of hmm fasta----------'
mafft_hmm ${prot}

# Remove major gaps in this alignment using Gblocks
echo '----------Removing massive gaps----------'
remove_gaps ${prot}

# Next steps
echo '----------Done, next steps are:----------'
echo '1. Check alignment file in UGENE'
echo '2.1 If alignment has >4000 positions, prepare reconcile2.1-server_raxml-ng.pbs file'
echo '2.2 Copy hmm alignment and pbs-file to HPC server (check Jens/HPC_info_n_code for code)'
echo '2.3 Copy alignment file to $VSC_DATA, run "qsub reconcile2.1-server_raxml-ng.pbs'
echo '2.4 Copy results back and run reconcile2.2-servre.sh locally'
echo '3.1 If alignment has <4000 positions, run it through local batch_reconcile2.sh
