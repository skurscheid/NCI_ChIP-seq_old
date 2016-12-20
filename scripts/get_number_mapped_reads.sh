#!/bin/bash

# The script extracts the number of mapped reads from a
# samtools flagstat output file
#
# Author: Maurits Evers (maurits.evers@anu.edu.au)
# Original date: 30/11/2016
# Last changed: 30/11/2016

version=1.0

function message {
    echo [`date`] $1;
}

function usage {
    echo "get_number_mapped_reads version ${version} by Maurits Evers (maurits.evers@anu.edu.au)"
    echo "Extract number of mapped reads from samtools flagstat output file."
    echo "Usage:"
    echo "  get_number_mapped_reads.sh <file>"
    echo ""
    echo "  <file>  samtools flagstat output file"
}

if [ $# != 1 ]; then
    usage;
    exit;
fi

file=$1;

cat $1 | awk -F" " 'FNR==1{printf "%i", $1}'
