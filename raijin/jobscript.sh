#!/bin/bash
# properties = {properties}

set -euo pipefail

module load sra-toolkit samtools bowtie2/2.3.1 R/3.3.0 

export TMPDIR=$PBS_JOBFS
export R_USER_LIBS="/short/rl2/R/x86_64-pc-linux-gnu-library"

{exec_job}
