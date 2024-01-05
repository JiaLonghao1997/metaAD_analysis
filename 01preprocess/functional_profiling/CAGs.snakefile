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


# function to get the sample contigs from the tsv
def get_sample_contigs(sample_file):
    sample_contigs = {}
    with open(sample_file) as lines:
        for line in lines:
            sample = line.strip().split(",")[0]
            contig = line.strip().split(",")[1]
            sample_contigs[sample] = contig
    return sample_contigs


# read in sample info and reads from the sample_file
sample_contigs = get_sample_contigs(config['sample_contigs'])
sample_names = list(sample_contigs.keys())
print("top 10 sample_names: {}".format(sample_names[0:10]))
print("sample numbers: {}".format(len(sample_names)))

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
sample_reads, paired_end = get_sample_reads(config['sample_reads'])
#sample_names = list(sample_reads.keys())

rule all:
    input:
        genes = expand(join(outdir,"genes","{samp}_gene.fna"), samp=sample_names),
        proteins = expand(join(outdir,"genes","{samp}_protein.faa"), samp=sample_names),
        cdhit_rep_seq = join(outdir, "cdhit_rep_seq.fna"),
        annotations= join(outdir,"cdhit_rep_seq.emapper.annotations"),
        gene_abundance= expand(join(outdir,"coverm","{samp}_gene_abundance.txt"), samp=sample_names),
        merge_GOs= join(outdir,"merge_GOs.txt"),
        merge_KOs=join(outdir,"merge_KOs.txt"),
        merge_pathways=join(outdir,"merge_pathways.txt"),
        merge_modules = join(outdir, "merge_modules.txt")

rule filter_500bp:
    input:
        contigs = lambda wildcards: sample_contigs[wildcards.samp]
    output:
        contigs_500bp = join(outdir,"contigs_500bp","{samp}_contigs.fa")
    params:
        sample = lambda wildcards: wildcards.samp
    shell:
        """
        python /home1/jialh/brain/pipeline/2020NBTbhattlab/CAGs/scripts/filter_contigs_500bp.py \
        {params.sample} {input.contigs} 500 {output.contigs_500bp}
        """

rule prodigal:
    input:
        contigs_500bp = join(outdir,"contigs_500bp","{samp}_contigs.fa")
    output:
        gene = join(outdir,"genes","{samp}_gene.fna"),
        protein = join(outdir, "genes", "{samp}_protein.faa"),
        gff = join(outdir, "genes", "{samp}_gene.gff"),
        stat = join(outdir, "genes", "{samp}_stat.txt")
    shell: """
    /home1/jialh/anaconda3/envs/MAG/bin/prodigal -i {input.contigs_500bp} \
    -d {output.gene} -a {output.protein} -o {output.gff} \
    -f gff -p meta -s {output.stat} -q
    """

rule merge_gene_fna:
    input:
        genes = expand(join(outdir,"genes","{samp}_gene.fna"), samp=sample_names)
    output:
        merge_genes = join(outdir,"merge_gene.fna")
    params:
        gene_dir=join(outdir,"genes")
    shell:
        """
        /usr/bin/cat {params.gene_dir}/*_gene.fna > {output.merge_genes}
        """

##reference: https://github.com/UriNeri/RVMT/blob/main/Clustering/DoubleClustering.sh#L20
rule linclust:
    input:
        merge_genes = join(outdir,"merge_gene.fna")
    output:
        linclust_genes = join(outdir, "linclust_rep_seq.fasta")
    params:
        out_prefix = join(outdir,"linclust"),
        out_temp = join(outdir,"linclust_temp")
    threads: 32
    shell:
        """
        /home1/jialh/anaconda3/envs/MAG/bin/mmseqs easy-linclust \
        --min-seq-id 0.95 -c 0.9 --cov-mode 1 --threads {threads} --split-memory-limit 16G \
        {input.merge_genes} {params.out_prefix} {params.out_temp}
        """

##https://github.com/jiaonall/CRC-multi-kingdom/blob/main/1_raw%20sequence%20process.sh
##multi_process: https://openmp.llvm.org/
##Please install cdhit from source(https://github.com/weizhongli/cdhit).
##`conda install -c agbiome cdhit` ===> Fail multi_preprocessing.
rule cdhit:
    input:
        linclust_genes = join(outdir, "linclust_rep_seq.fasta")
    output:
        representative_genes = join(outdir, "cdhit_rep_seq.fna")
    threads: 32
    shell:
        """
        /home1/jialh/brain/tools/cd-hit-v4.8.1-2019-0228/cd-hit-est \
        -i {input.linclust_genes} -o {output.representative_genes} \
        -aS 0.9 -c 0.95 -G 0 -g 0 -T {threads} -M 0
        """

#
rule EggNOG_mapper_search:
    input:
        representative_genes = join(outdir,"cdhit_rep_seq.fna"),
        EggNOGdb = config['EggNOGdb']
    output:
        hits = join(outdir,"cdhit_rep_seq.emapper.hits"),
        seed_orthologs = join(outdir,"cdhit_rep_seq.emapper.seed_orthologs")
    params:
        out_prefix = join(outdir,"cdhit_rep_seq")
    threads: 32
    shell:
        """
        alias diamond=/home1/jialh/anaconda3/envs/MAG/bin/diamond
        /home1/jialh/anaconda3/envs/MAG/bin/python \
        /home1/jialh/anaconda3/envs/MAG/bin/emapper.py \
        -m diamond --no_annot --no_file_comments --itype CDS --cpu {threads} \
        --data_dir {input.EggNOGdb} -i {input.representative_genes} -o {params.out_prefix}
        """

rule EggNOG_mapper_annotate:
    input:
        seed_orthologs = join(outdir,"cdhit_rep_seq.emapper.seed_orthologs"),
        EggNOGdb= config['EggNOGdb']
    output:
        annotations = join(outdir,"cdhit_rep_seq.emapper.annotations")
    params:
        out_prefix = join(outdir,"cdhit_rep_seq")
    threads: 32
    shell:
        """
        /home1/jialh/anaconda3/envs/MAG/bin/python \
        /home1/jialh/anaconda3/envs/MAG/bin/emapper.py \
        --annotate_hits_table {input.seed_orthologs}  \
        --no_file_comments --dbmem --cpu {threads} \
        --data_dir {input.EggNOGdb} \
        -o {params.out_prefix}
        """

##reference: https://github.com/deng-lab/viroprofiler/blob/main/modules/local/abundance.nf
##[WARNING] For a multi-part index, no @SQ lines will be outputted. Please use --split-prefix. https://github.com/lh3/minimap2/issues/301
##node: One CPUs contains two cores. If threads=4, %CPU will be around 200%.
rule mapping_reads_to_contigs:
    input:
        r1 = lambda wildcards: sample_reads[wildcards.samp][0],
        r2 = lambda wildcards: sample_reads[wildcards.samp][1],
        representative_genes = join(outdir,"cdhit_rep_seq.fna")
    output:
        bam = join(outdir, "bam", "{samp}.bam")
    params:
        tempdir = join(outdir, "bam", "{samp}_temp_bam")
    threads: 4
    shell:"""
        /home1/jialh/tools/anaconda3/envs/coverm/bin/minimap2 -t {threads} \
        -ax sr --split-prefix {params.tempdir} {input.representative_genes} {input.r1} {input.r2} | \
        /home1/jialh/tools/anaconda3/bin/samtools sort --threads {threads} -o {output.bam}
    """

rule CoverM:
    input:
        bam = join(outdir, "bam", "{samp}.bam")
    output:
        gene_abundance = join(outdir, "coverm", "{samp}_gene_abundance.txt")
    params:
        temp_dir = join(outdir, "coverm", "{samp}_coverm_temp")
    threads: 4
    conda: "/home1/jialh/brain/pipeline/2020NBTbhattlab/CAGs/scripts/coverm.yaml"
    shell: """
        /home1/jialh/tools/anaconda3/envs/coverm/bin/coverm contig \
        --bam-files {input.bam} -t {threads} --min-read-percent-identity 0.95 \
        --output-file {output.gene_abundance}
    """

##reference: https://github.com/NBISweden/nbis-meta/blob/main/workflow/scripts/eggnog-parser.py
rule eggNOG_parse:
    input:
        KEGGdb = config['KEGGdb'],
        annotations = join(outdir,"cdhit_rep_seq.emapper.annotations")
    output:
        EggNOG_cazy = join(outdir, 'EggNOGdb', "cazy.parsed.tsv"),
        EggNOG_enzymes= join(outdir,'EggNOGdb',"enzymes.parsed.tsv"),
        EggNOG_gos= join(outdir,'EggNOGdb',"gos.parsed.tsv"),
        EggNOG_kos= join(outdir,'EggNOGdb',"kos.parsed.tsv"),
        EggNOG_modules= join(outdir,'EggNOGdb',"modules.parsed.tsv"),
        EggNOG_pathways= join(outdir,'EggNOGdb',"pathways.parsed.tsv"),
        EggNOG_tc= join(outdir,'EggNOGdb',"tc.parsed.tsv")
    params:
        EggNOGdb_dir = join(outdir,"EggNOGdb")
    shell:
        """
        python /home1/jialh/brain/pipeline/2020NBTbhattlab/CAGs/scripts/eggnog-parser.py parse \
        --map_go {input.KEGGdb} {input.annotations} {params.EggNOGdb_dir}
        """

rule eggNOG_quantify:
    input:
        EggNOG_gos = join(outdir,'EggNOGdb',"gos.parsed.tsv"),
        EggNOG_kos = join(outdir,'EggNOGdb',"kos.parsed.tsv"),
        EggNOG_pathways = join(outdir,'EggNOGdb',"pathways.parsed.tsv"),
        EggNOG_modules= join(outdir,'EggNOGdb',"pathways.parsed.tsv"),
        gene_abundance = join(outdir,"coverm", "{samp}_gene_abundance.txt")
    output:
        GOs = join(outdir,"functional_annotation", "{samp}_GOs.txt"),
        KOs = join(outdir,"functional_annotation", "{samp}_KOs.txt"),
        pathways = join(outdir,"functional_annotation", "{samp}_pathways.txt"),
        modules = join(outdir,"functional_annotation", "{samp}_modules.txt")
    shell:
        """
        python /home1/jialh/brain/pipeline/2020NBTbhattlab/CAGs/scripts/eggnog-parser.py quantify \
        {input.gene_abundance} {input.EggNOG_gos} {output.GOs}
        
        python /home1/jialh/brain/pipeline/2020NBTbhattlab/CAGs/scripts/eggnog-parser.py quantify \
        {input.gene_abundance} {input.EggNOG_kos} {output.KOs}
        
        python /home1/jialh/brain/pipeline/2020NBTbhattlab/CAGs/scripts/eggnog-parser.py quantify \
        {input.gene_abundance} {input.EggNOG_pathways} {output.pathways}
        
        python /home1/jialh/brain/pipeline/2020NBTbhattlab/CAGs/scripts/eggnog-parser.py quantify \
        {input.gene_abundance} {input.EggNOG_modules} {output.modules}
        """

rule eggNOG_merge:
    input:
        GOs = expand(join(outdir,"functional_annotation","{samp}_GOs.txt"), samp=sample_names),
        KOs = expand(join(outdir,"functional_annotation","{samp}_KOs.txt"), samp=sample_names),
        pathways = expand(join(outdir,"functional_annotation","{samp}_pathways.txt"), samp=sample_names),
        modules= expand(join(outdir,"functional_annotation","{samp}_modules.txt"),samp=sample_names)
    output:
        merge_GOs = join(outdir,"merge_GOs.txt"),
        merge_KOs= join(outdir,"merge_KOs.txt"),
        merge_pathways= join(outdir,"merge_pathways.txt"),
        merge_modules= join(outdir,"merge_modules.txt")
    params:
        functional_annotation = join(outdir,"functional_annotation")
    shell:
        """
        python /home1/jialh/brain/pipeline/2020NBTbhattlab/CAGs/scripts/eggnog-parser.py merge \
        {params.functional_annotation}/*_GOs.txt {output.merge_GOs}
        
        python /home1/jialh/brain/pipeline/2020NBTbhattlab/CAGs/scripts/eggnog-parser.py merge \
        {params.functional_annotation}/*_KOs.txt {output.merge_KOs}
        
        python /home1/jialh/brain/pipeline/2020NBTbhattlab/CAGs/scripts/eggnog-parser.py merge \
        {params.functional_annotation}/*_pathways.txt {output.merge_pathways}
        
        python /home1/jialh/brain/pipeline/2020NBTbhattlab/CAGs/scripts/eggnog-parser.py merge \
        {params.functional_annotation}/*_modules.txt {output.merge_modules}
        """
