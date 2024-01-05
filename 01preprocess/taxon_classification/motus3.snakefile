# Snakefile for metaphaln3 classification
# works with paired end read files
from os.path import join
import sys
import snakemake
import time

# output base directory
outdir = config['outdir']
#db = config['db']
localrules: 

# function to get the sample reads from the tsv
def get_sample_reads(sample_file):
    sample_reads = {}
    paired_end = ''
    with open(sample_file) as sf:
        for l in sf.readlines():
            s = l.strip().split("\t")
            if len(s) == 1 or s[0] == 'Sample' or s[0] == '#Sample' or s[0].startswith('#'):
                continue
            sample = s[0]
            # paired end specified
            if (len(s)==3):
                reads = [s[1],s[2]]
                if paired_end != '' and not paired_end:
                    sys.exit('All samples must be paired or single ended.')
                paired_end = True
            # single end specified
            elif len(s)==2:
                reads=s[1]
                if paired_end != '' and paired_end:
                    sys.exit('All samples must be paired or single ended.')
                paired_end = False
            if sample in sample_reads:
                raise ValueError("Non-unique sample encountered!")
            sample_reads[sample] = reads
    return (sample_reads, paired_end)


# read in sample info and reads from the sample_file
sample_reads, paired_end = get_sample_reads(config['sample_file'])
sample_names = list(sample_reads.keys())
print("top 10 sample_names: {}".format(sample_names[0:10]))
print("sample numbers: {}".format(len(sample_names)))

rule all:
    input:
        expand(join(outdir, "results/{samp}.txt"), samp=sample_names),
        join(outdir, "merged_abundance_table.txt")

rule motus3:
    input:
        r1 = lambda wildcards: sample_reads[wildcards.samp][0],
        r2 = lambda wildcards: sample_reads[wildcards.samp][1],
    output:
        join(outdir, "results", "{samp}.txt")
    threads: 8
    resources: 
        mem = 64, 
        time = 24
    singularity: "docker://quay.io/biocontainers/metaphlan:3.0.13--pyhb7b1952_0"
    shell: """
    export PATH=/home1/jialh/tools/miniconda3/envs/motus3/bin:$PATH
    
    /home1/jialh/tools/miniconda3/envs/motus3/bin/motus profile \
    -f {input.r1} -r {input.r2} -t {threads} > {output}
    """

rule merge_tables:
    input:
        expand(join(outdir, "results/{samp}.txt"), samp=sample_names)
    output:
        t1 = join(outdir, "merged_abundance_table.txt")
    params:
        results_dir = join(outdir, "results")
    singularity: "docker://quay.io/biocontainers/metaphlan:3.0.13--pyhb7b1952_0"
    shell: """
        /home1/jialh/tools/miniconda3/envs/motus3/bin/motus merge -i {params.results_dir}/*.txt > {output.t1}
    """

