# Snakefile for metaphaln3 classification
# works with paired end read files
from os.path import join
import sys
import snakemake
import time
# output base directory
outdir = config['outdir']
db = config['db']
readlength = config['readlength']
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
# print(sample_names)

rule all:
    input:
        reports = expand(join(outdir,"kraken2/{samp}-kraken2-report.txt"), samp=sample_names),
        brackens = expand(join(outdir,"kraken2/{samp}-bracken-report.txt"), samp=sample_names),
        report_combined = join(outdir,"merged-kraken2-report02.txt"),
        bracken_combined =join(outdir,"merged-bracken_new-report02.txt"),
        bracken_rel =join(outdir,"merged-bracken_new-report_rel02.txt")
        #join(outdir, "merged_abundance_table_species.txt"),

##reference:https://lichenhao.netlify.app/post/2020-08-22-krakentools/
rule kraken2:
    input:
        r1 = lambda wildcards: sample_reads[wildcards.samp][0],
        r2 = lambda wildcards: sample_reads[wildcards.samp][1],
    output:
        report = join(outdir, "kraken2", "{samp}-kraken2-report.txt")
    params:
        db = db
    threads: 8
    resources: 
        mem = 64, 
        time = 24
    shell: """
    /share/home1/jialh/tools/miniconda3/bin/kraken2 \
    -db {params.db} \
    --threads {threads} \
    --report-minimizer-data --minimum-hit-groups 3 \
    --use-names --report {output.report} \
    --paired {input.r1} {input.r2} > /dev/null
    """

##reference: https://hackmd.io/@astrobiomike/kraken2-bracken-standard-build
rule bracken:
    input:
        report= join(outdir,"kraken2","{samp}-kraken2-report.txt")
    output:
        brackenout = join(outdir,"kraken2","{samp}-bracken-report.txt"),
        bracken_new = join(outdir,"kraken2","{samp}-bracken_new-report.txt")
    threads: 8
    params:
        readlength = readlength
    shell: """
        #rm -rf {input.report}
        /share/home1/jialh/brain/tools/Bracken/bracken \
        -r {params.readlength} \
        -l S -t 10 \
        -d /share/home1/jialh/brain/databases/Kraken2/db \
        -i {input.report} \
        -o {output.brackenout} \
        -w {output.bracken_new} > /dev/null
    """

rule bracken2mpa:
    input:
        report= join(outdir,"kraken2","{samp}-kraken2-report.txt"),
        bracken=join(outdir,"kraken2","{samp}-bracken_new-report.txt")
    output:
        report_mpa= join(outdir,"mpa","{samp}-kraken2-report.txt"),
        bracken_mpa=join(outdir,"mpa","{samp}-bracken_new-report.txt"),
        bracken_rel=join(outdir,"mpa","{samp}-bracken_new-report_rel.txt")
    params:
        readlength=readlength
    shell: """
    /share/home1/jialh/tools/miniconda3/bin/python \
    /share/home1/jialh/brain/tools/KrakenTools/kreport2mpa_with_name.py \
    -r {input.report} -o {output.report_mpa}
    
    /share/home1/jialh/tools/miniconda3/bin/python \
    /share/home1/jialh/brain/tools/KrakenTools/kreport2mpa_with_name.py \
    -r {input.bracken} -o {output.bracken_mpa}
    
    /share/home1/jialh/tools/miniconda3/bin/python \
    /share/home1/jialh/brain/tools/KrakenTools/kreport2mpa_with_name.py \
    --percentages \
    -r {input.bracken} -o {output.bracken_rel}
    """

rule combine_mpa:
    input:
        report = expand(join(outdir,"mpa","{samp}-kraken2-report.txt"), samp=sample_names),
        bracken = expand(join(outdir,"mpa","{samp}-bracken_new-report.txt"), samp=sample_names),
        bracken_rel = expand(join(outdir,"mpa","{samp}-bracken_new-report_rel.txt"), samp=sample_names)
    output:
        report = join(outdir, "merged-kraken2-report.txt"),
        bracken = join(outdir, "merged-bracken_new-report.txt"),
        bracken_rel = join(outdir,"merged-bracken_new-report_rel.txt")
    params:
        results_dir = join(outdir, "mpa")
    shell: """
        /share/home1/jialh/tools/miniconda3/envs/mpa/bin/python \
        /share/home1/jialh/brain/tools/KrakenTools/merge_metaphlan_tables.py \
        {params.results_dir}/*-kraken2-report.txt > {output.report}
        
        /share/home1/jialh/tools/miniconda3/envs/mpa/bin/python \
        /share/home1/jialh/brain/tools/KrakenTools/merge_metaphlan_tables.py \
        {params.results_dir}/*-bracken_new-report.txt > {output.bracken}
        
        /share/home1/jialh/tools/miniconda3/envs/mpa/bin/python \
        /share/home1/jialh/brain/tools/KrakenTools/merge_metaphlan_tables.py \
        {params.results_dir}/*-bracken_new-report_rel.txt > {output.bracken_rel}
    """

rule rename_mpa:
    input:
        report=join(outdir,"merged-kraken2-report.txt"),
        bracken=join(outdir,"merged-bracken_new-report.txt"),
        bracken_rel=join(outdir,"merged-bracken_new-report_rel.txt")
    output:
        report = join(outdir,"merged-kraken2-report02.txt"),
        bracken=join(outdir,"merged-bracken_new-report02.txt"),
        bracken_rel=join(outdir,"merged-bracken_new-report_rel02.txt")
    shell: """
    python /share/home1/jialh/brain/pipeline/2020NBTbhattlab/taxon_classification/scripts/01rename_mpa.py {input.report} {output.report}
    python /share/home1/jialh/brain/pipeline/2020NBTbhattlab/taxon_classification/scripts/01rename_mpa.py {input.bracken} {output.bracken}
    python /share/home1/jialh/brain/pipeline/2020NBTbhattlab/taxon_classification/scripts/01rename_mpa.py {input.bracken_rel} {output.bracken_rel}
    """

