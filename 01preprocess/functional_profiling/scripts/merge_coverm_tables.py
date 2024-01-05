#!/usr/bin/env python3
##copy from: /home1/jialh/tools/miniconda3/envs/mpa/lib/python3.7/site-packages/metaphlan/utils
import argparse
import os
import sys
import re
import pandas as pd
from itertools import takewhile

def merge( aaastrIn, ostm ):
    """
    Outputs the table join of the given pre-split string collection.

    :param  aaastrIn:   One or more split lines from which data are read.
    :type   aaastrIn:   collection of collections of string collections
    :param  iCol:       Data column in which IDs are matched (zero-indexed).
    :type   iCol:       int
    :param  ostm:       Output stream to which matched rows are written.
    :type   ostm:       output stream
    """

    listmpaVersion = set()
    merged_tables = pd.DataFrame()

    for f in aaastrIn:
        iIn = pd.read_csv(f, sep='\t', header=0, index_col=0).fillna('')
        #iIn = iIn.drop(['Mean'])
        #         # print("iIn.iloc[0:5, 0:5]", iIn.iloc[0:5, 0:5])
        #         # print("iIn.shape: ", iIn.shape)
        #         # index_col = [0]
        #         # iIn = iIn.set_index(iIn.columns[index_col].to_list())
        if merged_tables.empty:
            merged_tables = iIn
        else:
            merged_tables = pd.merge(iIn, merged_tables, how='outer', left_index=True, right_index=True )
    if listmpaVersion:
        ostm.write(list(listmpaVersion)[0]+'\n')
    merged_tables.fillna('0').reset_index().to_csv(ostm, index=False, sep = '\t', index_label="clade_name")

argp = argparse.ArgumentParser( prog = "merge_metaphlan_tables.py",
    description = """Performs a table join on one or more metaphlan output files.""")
argp.add_argument( "aistms",    metavar = "input.txt", nargs = "+",
    help = "One or more tab-delimited text tables to join" )
argp.add_argument( '-o',    metavar = "output.txt", nargs = 1,
    help = "Name of output file in which joined tables are saved" )

__doc__ = "::\n\n\t" + argp.format_help( ).replace( "\n", "\n\t" )

argp.usage = argp.format_usage()[7:]+"\n\n\tPlease make sure to supply file paths to the files to combine. If combining 3 files (Table1.txt, Table2.txt, and Table3.txt) the call should be:\n\n\t\tpython merge_metaphlan_tables.py Table1.txt Table2.txt Table3.txt > output.txt\n\n\tA wildcard to indicate all .txt files that start with Table can be used as follows:\n\n\t\tpython merge_metaphlan_tables.py Table*.txt > output.txt"


def main( ):
    args = argp.parse_args( )
    if args.o is None:
        merge(args.aistms, sys.stdout)
    else:
        with open(args.o[0], 'w') as fout:
            merge(args.aistms, fout)

if __name__ == '__main__':
    main()