import pandas as pd
import os

sample_names = snakemake@params[['sample_names']]
assembly_dir = snakemake@params[['assembly_dir']]

merge_df = pd.DataFrame()
for sample in sample_names:
    print("deal with sample: {}".format(sample))
    infile = os.path.join(assembly_dir, sample, 'quast/report.tsv')
    if not os.path.exists(infile):
        print(('**********Cant find {}!***********'.format(infile)))
    else:
        quast_df = pd.read_table(infile, header=0, index_column=0)
        if merge_df.shape[1] == 0:
            merge_df = quast_df
        else:
            merge_df = pd.merge(merge_df, quast_df, left_index=True, right_index=True)

merge_df.columns = merge_df.columns.str.rstrip('.contig')
merge_df = merge_df.T
merge_df.to_csv(snakemake@output[[1]], sep = '\t', index=True, index_label="Sample")