#!/bin/bash

snakemake --dag | dot -Tpdf > dag.pdf
snakemake --rulegraph | dot -Tpdf > rulegraph.pdf

snakemake -p --configfile config.yaml \
	  --snakefile Snakefile
