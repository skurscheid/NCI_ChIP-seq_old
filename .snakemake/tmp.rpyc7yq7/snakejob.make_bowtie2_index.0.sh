#!/bin/sh
# properties = {"params": {"cmd": "bowtie2-build", "base": "ref/hg38/hg38"}, "rule": "make_bowtie2_index", "threads": 1, "input": ["ref/hg38/hg38.fa"], "local": false, "cluster": {"threads": 2, "h_vmem": "25G", "outstream": "/dev/null", "queue": "all.q,hugemem.q", "virtual_free": "24G", "errorstream": "/dev/null"}, "resources": {}, "output": ["ref/hg38/hg38.1.bt2l", "ref/hg38/hg38.2.bt2l", "ref/hg38/hg38.3.bt2l", "ref/hg38/hg38.4.bt2l"]}
cd /Users/u2528469/Work/ChIPseq_simul/snakemake && \
/Users/u2528469/miniconda3/bin/snakemake --snakefile /Users/u2528469/Work/ChIPseq_simul/snakemake/Snakefile \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files ref/hg38/hg38.fa /Users/u2528469/Work/ChIPseq_simul/.snakemake/tmp.rpyc7yq7 --latency-wait 5 \
--benchmark-repeats 1 \
 --configfile config.yaml --nocolor \
--notemp --quiet --no-hooks --nolock ref/hg38/hg38.1.bt2l ref/hg38/hg38.2.bt2l ref/hg38/hg38.3.bt2l ref/hg38/hg38.4.bt2l --force-use-threads  --allowed-rules make_bowtie2_index  && touch "/Users/u2528469/Work/ChIPseq_simul/.snakemake/tmp.rpyc7yq7/0.jobfinished" || (touch "/Users/u2528469/Work/ChIPseq_simul/.snakemake/tmp.rpyc7yq7/0.jobfailed"; exit 1)

