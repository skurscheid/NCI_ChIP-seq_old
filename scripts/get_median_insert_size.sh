#!/bin/bash

# The script extracts the median fragment size from a
# picard-tools CollectInsertSizeMetrics output file.
#
# Author: Maurits Evers (maurits.evers@anu.edu.au)
# Original date: 29/11/2016
# Last changed: 29/11/2016

version=1.0

function message {
    echo [`date`] $1;
}

function usage {
    echo "get_median_insert_size version ${version} by Maurits Evers (maurits.evers@anu.edu.au)"
    echo "Extract median insert size from picard-tools CollectInsertSizeMetrics text file."
    echo "Usage:"
    echo "  get_median_insert_size.sh <file>"
    echo ""
    echo "  <file>  picard-tools CollectInsertSizeMetrics output text file"
}

if [ $# != 1 ]; then
    usage;
    exit;
fi

file=$1;

# Explanation:
# grep -E /^[^#]/      Ignore lines starting with #
# -F"\t"               Separate on tabs
# FNR==2               Print only 2nd line
# printf "%i", $1      Print first column entries
cat $1 | grep -E "^[^#]" | awk -F"\t" 'FNR==2{printf "%i", $1}'
