#!/bin/bash

#filename="ftp://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz";
#curl $filename > GRCh38_chr22.fa.gz

filename="ftp://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz";
curl $filename > GRCh38_MT.fa.gz
