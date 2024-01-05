###assembly
configfile=$1

export PATH=/home1/jialh/tools/miniconda3/bin:$PATH

/home1/jialh/tools/miniconda3/bin/python \
/home1/jialh/tools/miniconda3/bin//snakemake \
-s /home1/jialh/metaHiC/pipelines/2020NBTbhattlab/assembly/assembly.snakefile \
--configfile ${configfile} \
--restart-times 0 --keep-going --cores 32 --latency-wait 100 --rerun-incomplete #--unlock

