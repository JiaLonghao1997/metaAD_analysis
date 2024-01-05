import os
import multiprocessing


def run(species):
    os.system("/share/apps/R/4.2.0/lib64/R/bin/Rscript " +
              "/home1/jialh/brain/01meta/multikingdom/04mediation/01metaAD_mediationAnalysis_parallel.R " +
              "{} ".format(species))



# with open("Exposures.txt") as input:
#     for line in input:
#         species = line.strip()
#         print("deal with {}".format(species))
#         run(species)

pool = multiprocessing.Pool(processes=16)
with open("Exposures_unfinished.txtaa") as input:
    for line in input:
        species = line.strip()
        print("deal with {}".format(species))
        pool.apply_async(func=run, args=(species, ))
pool.close()
pool.join()