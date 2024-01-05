import re,os,subprocess
from os.path import join, expanduser, abspath
import multiprocessing
import pandas as pd

################################################################################
# specify project directories
DATA_DIR    = config["raw_reads_directory"]
PROJECT_DIR = config["output_directory"]
READ_SUFFIX = config["read_specification"]
EXTENSION   = config["extension"]
# if gzipped, set this. otherwise not
gz_ext = '.gz' if EXTENSION.endswith('.gz') else ''
# print(gz_ext)

# convert PROJECT_DIR and DATA_DIR to absolute path
if PROJECT_DIR[0] == '~':
    PROJECT_DIR = expanduser(PROJECT_DIR)
PROJECT_DIR = abspath(PROJECT_DIR)
if DATA_DIR[0] == '~':
    DATA_DIR = expanduser(DATA_DIR)
DATA_DIR = abspath(DATA_DIR)

# get file names
FILES = [f for f in os.listdir(DATA_DIR) if f.endswith(EXTENSION)]
SAMPLE_PREFIX = []
for root,dirs,files in os.walk(DATA_DIR):
    for file in files:
        if file.endswith(EXTENSION):
            sample = re.split('|'.join([ a for a in READ_SUFFIX]), file)[0]
            SAMPLE_PREFIX.append(sample)
            file = os.path.join(root,file)
            FILES.append(file)
print("Top 5 files: ", FILES[0:5])
SAMPLE_PREFIX = list(set(SAMPLE_PREFIX))
#SAMPLE_PREFIX = list(set([re.split('|'.join(['_' + a for a in READ_SUFFIX]), i)[0] for i in FILES]))
print("Top 5 samples: ", SAMPLE_PREFIX[0:5])
print("Number of samples: ", len(SAMPLE_PREFIX))


# config for trim_galore: only enable start_trim if >0
start_trim = config['trim_galore']['start_trim']
end_trim = config['trim_galore']['end_trim']
if (start_trim) >0: 
    start_trim_string = '--clip_R1 {a} --clip_R2 {a}'.format(a=str(start_trim))
else: 
    start_trim_string = ''
if (end_trim) >0: 
    end_trim_string = '--three_prime_clip_R1 {a} --three_prime_clip_R2 {a}'.format(a=str(end_trim))
else: 
    end_trim_string = ''

################################################################################
localrules: assembly_meta_file, pre_multiqc, post_multiqc, cleanup
rule all:
    input:
        expand(join(PROJECT_DIR, "01_processing/05_sync/{sample}_orphans.fq.gz"), sample=SAMPLE_PREFIX),
        expand(join(PROJECT_DIR, "01_processing/00_qc_reports/pre_fastqc/{sample}{read}_fastqc.html"), sample=SAMPLE_PREFIX, read=READ_SUFFIX),
        expand(join(PROJECT_DIR, "01_processing/00_qc_reports/post_fastqc/{sample}_{read}_fastqc.html"), sample=SAMPLE_PREFIX, read=['1', '2', 'orphans']),
        join(PROJECT_DIR, "01_processing/00_qc_reports/pre_multiqc/multiqc_report.html"),
        join(PROJECT_DIR, "01_processing/00_qc_reports/post_multiqc/multiqc_report.html"),
        join(PROJECT_DIR, "01_processing/assembly_input.txt"),
        join(PROJECT_DIR, "01_processing/classification_input.txt"),
        join(PROJECT_DIR, "01_processing/readcounts.tsv"),
        join(PROJECT_DIR, "01_processing/readcounts.pdf")

################################################################################
rule pre_fastqc:
    input:  
        fwd = join(DATA_DIR, "{sample}_" + READ_SUFFIX[0] + EXTENSION),
        rev = join(DATA_DIR, "{sample}_" + READ_SUFFIX[1] + EXTENSION)
    output:
        fwd = join(PROJECT_DIR,  "01_processing/00_qc_reports/pre_fastqc/{sample}_" + READ_SUFFIX[0] + "_fastqc.html"),
        rev = join(PROJECT_DIR,  "01_processing/00_qc_reports/pre_fastqc/{sample}_" + READ_SUFFIX[1] + "_fastqc.html")
    params:
        outdir = join(PROJECT_DIR, "01_processing/00_qc_reports/pre_fastqc/")
    threads: min(4, len(READ_SUFFIX))
    resources:
            time = 6,
            mem = 32
    singularity: "docker://quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1"
    benchmark: join(PROJECT_DIR,  "01_processing/00_qc_reports/pre_fastqc/{sample}_time.txt")
    shell: """
        mkdir -p {params.outdir}
        fastqc {input} --outdir {params.outdir} --threads {threads}
    """

rule pre_multiqc:
    input:
        expand(join(PROJECT_DIR,  "01_processing/00_qc_reports/pre_fastqc/{sample}{read}_fastqc.html"), sample=SAMPLE_PREFIX, read=READ_SUFFIX)
    output: join(PROJECT_DIR,  "01_processing/00_qc_reports/pre_multiqc/multiqc_report.html")
    params:
        indir = join(PROJECT_DIR,  "01_processing/00_qc_reports/pre_fastqc"),
        outdir = join(PROJECT_DIR,  "01_processing/00_qc_reports/pre_multiqc/")
    singularity: "docker://quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0"
    shell: """
        /share/home1/jialh/tools/miniconda3/envs/GNNBin/bin/multiqc --force {params.indir} -o {params.outdir}
    """

###############################################################################
rule deduplicate:
    input:
        fwd = join(DATA_DIR, "{sample}_" + READ_SUFFIX[0] + EXTENSION),
        rev = join(DATA_DIR, "{sample}_" + READ_SUFFIX[1] + EXTENSION)
    output:
        fwd = join(PROJECT_DIR, "01_processing/01_dedup/{sample}_1.fq.gz"),
        rev = join(PROJECT_DIR, "01_processing/01_dedup/{sample}_2.fq.gz"),
    params:
        tmp_fwd = '{sample}_R1.fastq.gz',
        tmp_rev = '{sample}_R2.fastq.gz',
        outdir = join(PROJECT_DIR, "01_processing/01_dedup/")
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        mem = lambda wildcards, attempt: attempt * 16, 
        time = 24
    singularity: "docker://dzs74/htstream"
    benchmark: join(PROJECT_DIR, "01_processing/01_dedup/{sample}_time.txt")
    shell: """
        mkdir -p {params.outdir} && cd {params.outdir}
        hts_SuperDeduper -1 {input.fwd} -2 {input.rev} -f {wildcards.sample} -F
        mv {params.tmp_fwd} {output.fwd}
        mv {params.tmp_rev} {output.rev}
    """
###############################################################################

################################################################################
rule trim_galore:
    input:
        fwd = rules.deduplicate.output.fwd,
        rev = rules.deduplicate.output.rev,
    output:
        fwd = join(PROJECT_DIR, "01_processing/02_trimmed/{sample}_1_val_1.fq.gz"),
        rev = join(PROJECT_DIR, "01_processing/02_trimmed/{sample}_2_val_2.fq.gz"),
        orp = join(PROJECT_DIR, "01_processing/02_trimmed/{sample}_unpaired.fq.gz")
    threads: 4
    resources:
        mem=32,
        time=lambda wildcards, attempt: attempt * 24
    params:
        orp_fwd = join(PROJECT_DIR, "01_processing/02_trimmed/{sample}_1_unpaired_1.fq.gz"),
        orp_rev = join(PROJECT_DIR, "01_processing/02_trimmed/{sample}_2_unpaired_2.fq.gz"),
        q_min   = config['trim_galore']['quality'],
        left    = config['trim_galore']['start_trim'],
        min_len = config['trim_galore']['min_read_length'],
        outdir  = join(PROJECT_DIR, "01_processing/02_trimmed/"),
        gz_output = str(gz_ext == '.gz').lower()
    singularity: "docker://quay.io/biocontainers/trim-galore:0.6.7--hdfd78af_0"
    benchmark: join(PROJECT_DIR, "01_processing/02_trimmed/{sample}_time.txt")
    shell: """
        mkdir -p {params.outdir}
        trim_galore \
        --path_to_cutadapt /share/home1/jialh/tools/miniconda3/bin/cutadapt \
        --quality {params.q_min} \
            --length {params.min_len} \
            --output_dir {params.outdir} \
            --paired {input.fwd} {input.rev} \
            --retain_unpaired \
            --cores {threads} \
            {start_trim_string} \
            {end_trim_string}

        #merge unpaired and gzip
        zcat -f {params.orp_fwd} {params.orp_rev} | /share/home1/jialh/tools/miniconda3/bin/pigz -b 32 -p {threads} > {output.orp}
        # delete intermediate files
        rm {params.orp_fwd} {params.orp_rev}
    """

################################################################################
rule rm_host_reads:
    input:
        index_amb = config['rm_host_reads']['host_genome'] + '.amb',
        index_ann = config['rm_host_reads']['host_genome'] + '.ann',
        index_bwt = config['rm_host_reads']['host_genome'] + '.bwt',
        index_pac = config['rm_host_reads']['host_genome'] + '.pac',
        index_sa = config['rm_host_reads']['host_genome'] + '.sa',
        fwd       = rules.trim_galore.output.fwd,
        rev       = rules.trim_galore.output.rev,
        orp       = rules.trim_galore.output.orp
    output:
        unmapped_1 = join(PROJECT_DIR, "01_processing/05_sync/{sample}_1.fq.gz"),
        unmapped_2 = join(PROJECT_DIR, "01_processing/05_sync/{sample}_2.fq.gz"),
        unmapped_singletons = join(PROJECT_DIR, "01_processing/05_sync/{sample}_orphans.fq.gz"),
    params:
        bwa_index_base = join(config['rm_host_reads']['host_genome']),
        singelton_temp_1 = join(PROJECT_DIR, "01_processing/04_host_align/{sample}_rmHost_singletons1.fq.gz"),
        singelton_temp_2 = join(PROJECT_DIR, "01_processing/04_host_align/{sample}_rmHost_singletons2.fq.gz"),
    threads: 8
    resources:
        mem_mb=32000,
        mem=32,
        time=24
    singularity: "shub://bsiranosian/bens_1337_workflows:align"
    # conda: "envs/align.yaml"
    benchmark: join(PROJECT_DIR, "01_processing/04_host_align/{sample}_time.txt")
    shell: """
        mkdir -p {PROJECT_DIR}/01_processing/04_host_align/
        # if an index needs to be built, use bwa index ref.fa
        # run on paired reads
        /share/home1/jialh/tools/miniconda3/bin/bwa mem -t {threads} {params.bwa_index_base} {input.fwd} {input.rev} | \
            /share/home1/jialh/tools/samtools-1.9/samtools fastq -@ {threads} -t -T BX -f 4 -1 {output.unmapped_1} -2 {output.unmapped_2} -s {params.singelton_temp_1} -
        # run on unpaired reads
        /share/home1/jialh/tools/miniconda3/bin/bwa mem -t {threads} {params.bwa_index_base} {input.orp} | \
            /share/home1/jialh/tools/samtools-1.9/samtools fastq -@ {threads} -t -T BX -f 4 - > {params.singelton_temp_2}
        # combine singletons 
        zcat -f {params.singelton_temp_1} {params.singelton_temp_2} | /share/home1/jialh/tools/miniconda3/bin/pigz -p {threads} > {output.unmapped_singletons}
        #rm {params.singelton_temp_1} {params.singelton_temp_2}
    """

################################################################################
rule post_fastqc:
    input:  join(PROJECT_DIR, "01_processing/05_sync/{sample}_{read}.fq.gz")
    output: join(PROJECT_DIR,  "01_processing/00_qc_reports/post_fastqc/{sample}_{read}_fastqc.html"),
    params:
        outdir = join(PROJECT_DIR, "01_processing/00_qc_reports/post_fastqc/")
    threads: 4
    resources:
            time = 6,
            mem = 32
    singularity: "docker://quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1"
    benchmark: join(PROJECT_DIR, "01_processing/00_qc_reports/post_fastqc/{sample}_{read}_time.txt")
    shell: """
        mkdir -p {params.outdir}
        if [ -z $(gzip -cd {input} | head -c1) ]; then
            echo EMPTY!
            touch {output}
        else
            fastqc {input} -f fastq --outdir {params.outdir} -t {threads}
        fi
    """

rule post_multiqc:
    input: expand(join(PROJECT_DIR,  "01_processing/00_qc_reports/post_fastqc/{sample}_{read}_fastqc.html"), sample=SAMPLE_PREFIX, read=['1', '2', 'orphans']),
    output: join(PROJECT_DIR,  "01_processing/00_qc_reports/post_multiqc/multiqc_report.html")
    params:
        indir = join(PROJECT_DIR,  "01_processing/00_qc_reports/post_fastqc"),
        outdir = join(PROJECT_DIR,  "01_processing/00_qc_reports/post_multiqc/")
    singularity: "docker://quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0"
    shell: """
        /share/home1/jialh/tools/miniconda3/envs/GNNBin/bin/multiqc --force {params.indir} -o {params.outdir}
    """

################################################################################
rule assembly_meta_file:
    input: expand(join(PROJECT_DIR, "01_processing/05_sync/{sample}_orphans.fq.gz"), sample=SAMPLE_PREFIX)
    output: join(PROJECT_DIR, "01_processing/assembly_input.txt")
    run:
        outfile = str(output)
        if (os.path.exists(outfile)):
            os.remove(outfile)
        with open(outfile, 'w') as outf:
            outf.writelines(['# Sample\tReads1.fq[.gz][,Reads2.fq[.gz][,orphans.fq[.gz]]]\n'])
            for sample in SAMPLE_PREFIX:
                outline = [sample, ','.join([
                join(PROJECT_DIR, "01_processing/05_sync/" + sample + "_1.fq.gz"),
                join(PROJECT_DIR, "01_processing/05_sync/" + sample + "_2.fq.gz"),
                join(PROJECT_DIR, "01_processing/05_sync/" + sample + "_orphans.fq.gz")])]
                outf.writelines('\t'.join(outline) + '\n')

################################################################################
rule classification_meta_file:
    input: expand(join(PROJECT_DIR, "01_processing/05_sync/{sample}_orphans.fq.gz"), sample=SAMPLE_PREFIX)
    output: join(PROJECT_DIR, "01_processing/classification_input.txt")
    run:
        outfile = str(output)
        if (os.path.exists(outfile)):
            os.remove(outfile)
        with open(outfile, 'w') as outf:
            outf.writelines(['# Sample\tr1\tr2\n'])
            for sample in SAMPLE_PREFIX:
                outline = [sample, '\t'.join([
                join(PROJECT_DIR, "01_processing/05_sync/" + sample + "_1.fq.gz"),
                join(PROJECT_DIR, "01_processing/05_sync/" + sample + "_2.fq.gz")])]
                outf.writelines('\t'.join(outline) + '\n')

################################################################################
def file_len(fname, outfile):
    read_count = 0
    read_len = 0
    #print("outfile: {}".format(outfile))
    if os.path.exists(outfile) and os.path.getsize(outfile)>0:
        read_stat_df = pd.read_table(outfile, header=0, index_col=0)
        read_count = read_stat_df.loc[fname, 'num_seqs']
        read_len = read_stat_df.loc[fname, 'sum_len']
    else:
        os.system("/home1/jialh/anaconda3/bin/seqkit stat -T -j 4 {} > {}".format(fname, outfile))
        #print("outfile:", outfile)
        if os.path.exists(outfile):
            read_stat_df = pd.read_table(outfile, header=0, index_col=0)
            print("read_stat_df: ", read_stat_df)
            read_count = read_stat_df.loc[fname, 'num_seqs']
            read_len = read_stat_df.loc[fname, 'sum_len']
    return read_count, read_len

def fastq_stat(sample, DATA_DIR, PROJECT_DIR, READ_SUFFIX, EXTENSION, gz_ext):
    raw_file = join(DATA_DIR, sample + READ_SUFFIX[0] + EXTENSION)
    dedup_file = join(PROJECT_DIR, "01_processing/01_dedup/" + sample + "_1.fq.gz")
    trimmed_file = join(PROJECT_DIR, "01_processing/02_trimmed/" + sample + "_1_val_1.fq" + gz_ext)
    rmhost_file = join(PROJECT_DIR, "01_processing/05_sync/" + sample + "_1.fq.gz")
    orphans_file = join(PROJECT_DIR, "01_processing/05_sync/" + sample + "_orphans.fq.gz")

    raw_stat = join(PROJECT_DIR,  "01_processing/00_rawstat/" + sample + READ_SUFFIX[0] + EXTENSION) + '.stat.txt'
    dedup_stat = dedup_file + '.stat.txt'
    trimmed_stat = trimmed_file + '.stat.txt'
    rmhost_stat = rmhost_file + '.stat.txt'
    orphans_stat = orphans_file + '.stat.txt'

    raw_count, raw_len = file_len(raw_file, raw_stat)
    dedup_count, dedup_len = file_len(dedup_file, dedup_stat)
    trimmed_count, trimmed_len = file_len(trimmed_file, trimmed_stat)
    rmhost_count, rmhost_len = file_len(rmhost_file, rmhost_stat)
    orphans_count, orphans_len = file_len(orphans_file, orphans_stat)

    dedup_count_frac = round(dedup_count / float(raw_count), 3)
    trimmed_count_frac = round(trimmed_count / float(raw_count), 3)
    rmhost_count_frac = round(rmhost_count / float(raw_count), 3)
    orphans_count_frac = round(orphans_count / float(raw_count), 3)

    dedup_len_frac = round(dedup_len / float(raw_len), 3)
    trimmed_len_frac = round(trimmed_len / float(raw_len), 3)
    rmhost_len_frac = round(rmhost_len / float(raw_len), 3)
    orphans_len_frac = round(orphans_len / float(raw_len), 3)

    line = '\t'.join([sample, str(raw_count),  str(raw_len),
                      str(dedup_count), str(dedup_count_frac), str(dedup_len), str(dedup_len_frac),
                      str(trimmed_count), str(trimmed_count_frac), str(trimmed_len), str(trimmed_len_frac),
                      str(rmhost_count), str(rmhost_count_frac), str(rmhost_len), str(rmhost_len_frac),
                      str(orphans_count), str(orphans_count_frac), str(orphans_len), str(orphans_len_frac)])
    print(line)
    return(line)

rule readcounts:
    input:
        raw = expand(join(DATA_DIR, "{sample}" + READ_SUFFIX[0] + EXTENSION), sample=SAMPLE_PREFIX),
        dedup = expand(join(PROJECT_DIR, "01_processing/01_dedup/{sample}_1.fq.gz"), sample=SAMPLE_PREFIX),
        trimmed = expand(join(PROJECT_DIR, "01_processing/02_trimmed/{sample}_1_val_1.fq.gz"), sample=SAMPLE_PREFIX),
        rmhost = expand(join(PROJECT_DIR, "01_processing/05_sync/{sample}_1.fq.gz"), sample=SAMPLE_PREFIX),
        orphans = expand(join(PROJECT_DIR, "01_processing/05_sync/{sample}_orphans.fq.gz"), sample=SAMPLE_PREFIX)
    output:
        join(PROJECT_DIR, "01_processing/readcounts.tsv")
    resources:
        time = 24
    run:
        '''
        outfile = str(output)
        if (os.path.exists(outfile)):
            os.remove(outfile)
    
        pool = multiprocessing.Pool(processes=16)
        result_list = []
        for sample in SAMPLE_PREFIX:
            #print("readcounts:", sample)
            sample = sample.strip()
            print("readcounts:", sample)
            #run(sample, DATA_DIR, PROJECT_DIR, READ_SUFFIX, EXTENSION, gz_ext)
            result_list.append(pool.apply_async(func=run, args=(sample,DATA_DIR, PROJECT_DIR, READ_SUFFIX, EXTENSION, gz_ext)))
        pool.close()
        pool.join()
    
        with open(outfile, 'w') as outf, open(logfile, "w") as log:
            outf.writelines('Sample\traw_count\traw_len\t'+
                            'dedup_count\tdedup_count_frac\tdedup_len\tdedup_len_frac\t' +
                            'trimmed_count\ttrimmed_count_frac\ttrimmed_len\ttrimmed_len_frac\t'+
                            'host_removed_count\thost_removed_count_frac\thost_removed_len\thost_removed_len_frac\t'+
                            'orphan_count\torphan_count_frac\torphan_len\torphan_len_frac\n')
            for result in result_list:
                outf.write(result.get())
        '''

rule readcounts_graph:
    input:
        rules.readcounts.output
    output:
        join(PROJECT_DIR, "01_processing/readcounts.pdf")
    singularity: "shub://bhattlab/bhattlab_workflows:plotting"
    script:
        "/home1/jialh/anaconda3/envs/R420/bin/Rscript scripts/plot_readcounts.R"

################################################################################
rule cleanup:
    input: expand(join(PROJECT_DIR, "01_processing/05_sync/{sample}_orphans.fq.gz"), sample=SAMPLE_PREFIX)
    output: join(PROJECT_DIR, "cleaned")
    params:
        rmdir_1 = join(PROJECT_DIR, '01_processing/01_dedup'),
        rmdir_2 = join(PROJECT_DIR, '01_processing/02_trimmed'),
    shell: """
        rm -f {params.rmdir_1}/*.fq.gz
        rm -f {params.rmdir_2}/*.fq.gz
        touch {output}
    """