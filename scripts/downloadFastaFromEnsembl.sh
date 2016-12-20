#!/bin/bash

# The script downloads a reference genome from Ensembl,
# and stores it in a fasta file.
#
# Author: Maurits Evers (maurits.evers@anu.edu.au)
# Original date: 07/12/2016
# Last changed: 19/12/2016

version=1.0

function message {
    echo [`date`] $1;
}

function usage {
    echo "downloadFastaFromEnsembl version ${version} by Maurits Evers (maurits.evers@anu.edu.au)"
    echo "Download fasta sequence from Ensembl."
    echo "Usage:"
    echo "  downloadFastaFromEnsembl.sh <genome> <masked> <outdir>"
    echo ""
    echo "  <genome>   Genome version; can be one the following:"
    echo "                * GRCh38, hg38"
    echo "                * GRCh37, hg19"
    echo "                * GRCm38, mm10"
    echo "                * NCBIM37, mm9"
    echo "  <masked>   unmasked, masked, soft."
    echo "  <outdir>   Output directory"
}

if [ $# != 3 ]; then
    usage;
    exit;
fi

genome=$1;
mask=$2;
refdir=$3;

file="";
if [ $genome == "mm10" ] || [ $genome == "GRCm38" ]; then
    rootpath="ftp://ftp.ensembl.org/pub/release-86/fasta/mus_musculus/dna/";
    if [ $mask == "masked" ]; then
        file="${rootpath}Mus_musculus.GRCm38.dna_rm.toplevel.fa.gz";
    elif [ $mask == "soft" ]; then
	    file="${rootpath}Mus_musculus.GRCm38.dna_sm.toplevel.fa.gz";
    else
	    file="${rootpath}Mus_musculus.GRCm38.dna.toplevel.fa.gz";
    fi
elif [ $genome == "mm9" ] || [ $genome == "NCBIM37" ]; then
    rootpath="ftp://ftp.ensembl.org/pub/release-67/fasta/mus_musculus/dna/";
    if [ $mask == "masked" ]; then
        file="${rootpath}Mus_musculus.NCBIM37.67.dna_rm.toplevel.fa.gz";
    elif [ $mask == "soft" ]; then
        echo "Soft-masked version does not exist";
        exit;
    else
        file="${rootpath}Mus_musculus.NCBIM37.67.dna.toplevel.fa.gz";
    fi
elif [ $genome == "hg38" ] || [ $genome == "GRCh38" ]; then
    rootpath="ftp://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/dna/";
    if [ $mask == "masked" ]; then
        file="${rootpath}Homo_sapiens.GRCh38.dna_rm.toplevel.fa.gz";
    elif [ $mask == "soft" ]; then
        file="${rootpath}Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz";
    else
        file="${rootpath}Homo_sapiens.GRCh38.dna.toplevel.fa.gz";
    fi
elif [ $genome == "hg19" ] || [ $genome == "GRCh37" ]; then
    rootpath="http://ftp.ensembl.org/pub/grch37/release-83/fasta/homo_sapiens/dna/";
    if [ $mask == "masked" ]; then
        file="${rootpath}Homo_sapiens.GRCh37.dna_rm.toplevel.fa.gz";
    elif [ $mask == "soft" ]; then
        file="${rootpath}Homo_sapiens.GRCh37.dna_sm.toplevel.fa.gz";
    else
        file="${rootpath}Homo_sapiens.GRCh37.dna.toplevel.fa.gz";
    fi
else
    echo "[ERROR] ${genome} is not supported.";
    exit;
fi

message "Downloading fasta file"
if [ ! -d "${refdir}/${genome}" ]; then
    mkdir ${refdir}/${genome}
fi
curl --progress-bar $file > "${refdir}/${genome}/${genome}.fa.gz";
gunzip "${refdir}/${genome}/${genome}.fa.gz";

message "Done"
