import re,os,subprocess
import sys
from os.path import join, expanduser, abspath
import multiprocessing
import pandas as pd


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

def run(sample, DATA_DIR, PROJECT_DIR, READ_SUFFIX, EXTENSION, gz_ext):
    raw_file = join(DATA_DIR, sample + '_' + READ_SUFFIX[0] + EXTENSION)
    dedup_file = join(PROJECT_DIR, "01_processing/01_dedup/" + sample + "_1.fq.gz")
    trimmed_file = join(PROJECT_DIR, "01_processing/02_trimmed/" + sample + "_1_val_1.fq" + gz_ext)
    rmhost_file = join(PROJECT_DIR, "01_processing/05_sync/" + sample + "_1.fq.gz")
    orphans_file = join(PROJECT_DIR, "01_processing/05_sync/" + sample + "_orphans.fq.gz")

    raw_stat = join(PROJECT_DIR,  "01_processing/00_rawstat/" + sample + '_' + READ_SUFFIX[0] + EXTENSION) + '.stat.txt'
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
    #print(line)
    #outf.writelines(line + '\n')



DATA_DIR= sys.argv[1]
PROJECT_DIR= sys.argv[2]
SAMPLE_PREFIX = sys.argv[3]
output= sys.argv[4]
READ_SUFFIX = sys.argv[5]
EXTENSION = sys.argv[6]
gz_ext = sys.argv[7]

outfile = str(output)
if (os.path.exists(outfile)):
    os.remove(outfile)

pool = multiprocessing.Pool(processes=16)
result_list = []
#with open("/share/home1/jialh/brain/metaAD/00metadata/runid.head5.list") as SAMPLE_PREFIX:
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
                    'host_removed_count\thost_removed_count_frac\thost_removed_len\thost_removed_len_frac'+
                    'orphan_count\torphan_count_frac\torphan_len\torphan_len_frac\n')
    for result in result_list:
        outf.write(result.get())