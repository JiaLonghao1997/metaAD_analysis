dataset=$1
configfile=$2

export PATH="/home1/jialh/tools/miniconda3/envs/GNNBin/bin:$PATH"
export PATH="/home1/jialh/metaHiC/tools/biotools/HTStream_v1.3.0:$PATH"
export PATH="/home1/jialh/metaHiC/tools/biotools/TrimGalore-0.6.7:$PATH"
alias pigz=/home1/jialh/tools/miniconda3/bin/pigz

if [[ $dataset == "Gut" || $dataset == "Water" ||  $dataset == "Yeast" ]]
then
/home1/jialh/tools/miniconda3/bin/python \
/home1/jialh/tools/miniconda3/bin/snakemake \
-s /home1/jialh/metaHiC/pipelines/2020NBTbhattlab/preprocessing/preprocessing.snakefile \
--configfile ${configfile} \
--restart-times 0 --keep-going --cores 48 --latency-wait 100 --rerun-incomplete #--unlock
elif [[ $dataset == "GIS20" || $dataset == "CAMISIM" || $dataset == "RealData" ]]
then
/home1/jialh/tools/miniconda3/bin/python \
/home1/jialh/tools/miniconda3/bin/snakemake \
-s /home1/jialh/metaHiC/pipelines/2020NBTbhattlab/preprocessing/preprocessing_R1R2.snakefile \
--configfile ${configfile} \
--restart-times 0 --keep-going --cores 48 --latency-wait 100 --rerun-incomplete #--unlock
fi