#!/bin/bash

# The script downloads the complete murine ribosomal DNA repeat,
# and stores it in file rDNA_repeat.fa
# GeneBank Accession Number BK000964
# Reference: Grozdanov et al, Genomics 82, 637 (2003)
#
# Author: Maurits Evers (maurits.evers@anu.edu.au)
# Original date: 17/11/2016
# Last changed: 29/11/2016

version=0.9

function message {
    echo [`date`] $1;
}

function usage {
    echo "downloadFastaAccFromGenBank version ${version} by Maurits Evers (maurits.evers@anu.edu.au)"
    echo "Download fasta sequence from NCBI's GenBank database."
    echo "Usage:"
    echo "  downloadFastaAccFromGenBank.sh <Acc> <fasta>"
    echo ""
    echo "  <Acc>   GenBank Accession number."
    echo "  <fasta> Fasta output file."
    echo "          Note: Will we overwritten if exists."
}

if [ $# != 2 ]; then
    usage;
    exit;
fi

acc=$1;
fasta=$2;
#genbank=${fasta%.*}".gff3"

message "Downloading fasta file"
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${acc}&rettype=fasta" > $fasta

#message "Downloading gff3 file"
#curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${acc}&rettype=gff3" > $genbank

message "Done"
