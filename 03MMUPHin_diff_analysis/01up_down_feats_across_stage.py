import pandas as pd
import os

taxon_types = ["KOs", 'GMMs', 'GBMs', 'pathways']
for taxon_type in taxon_types:
    # taxon_type = "KOs"
    taxons = ["Archaea", "Bacteria", "Fungi", "Viruses"]
    stages = ["SCS", "SCD", "MCI", "AD"]

    workdir = os.path.join("/home1/jialh/brain/01meta/multikingdom/03MMUPHin_diff", taxon_type)
    infile = os.path.join(workdir, "MMUPHin_diff_merge.csv")
    outfile = os.path.join(workdir, "01up_down_stage.csv")
    data = pd.read_csv(infile, header=0, index_col=0)
    ####
    data_filter = data.loc[data.isna().sum(axis=1)<=2, ]
    data_filter.fillna(0, inplace=True)
    print("早期或者晚期都显著上升: ")
    species_list = data_filter.loc[(((data_filter[['SCS', 'SCD']]>=0).sum(axis=1)>=2) & ((data_filter[['MCI', 'AD']]>=0).sum(axis=1)>=2)) & ((data_filter[['MCI', 'AD']]>0).sum(axis=1)>=1), ].index.values.tolist()
    # for species in species_list:
    #     print(species)
    print(data_filter.loc[(((data_filter[['SCS', 'SCD']]>=0).sum(axis=1)>=2) & ((data_filter[['MCI', 'AD']]>=0).sum(axis=1)>=2)), ].index.values.tolist())
    print("只在早期显著上升:")
    print(data_filter.loc[(((data_filter[['SCS', 'SCD']]>=0).sum(axis=1)>=2) & ((data_filter[['MCI', 'AD']]>0).sum(axis=1)==0)), ].index.values.tolist())

    print("只在后期显著上升:")
    print(data_filter.loc[(((data_filter[['SCS', 'SCD']]>0).sum(axis=1)==0) & ((data_filter[['MCI', 'AD']]>=0).sum(axis=1)>=2)), ].index.values.tolist())

    with open(outfile, "w") as output:
        output.write("stage,taxon,change,count\n")
        for stage in stages:
            if taxon_type == "species":
                for taxon in taxons:
                    up_count = data.loc[(data.index.str.contains("k__{}".format(taxon))) & (data[stage]>0), ].shape[0]
                    down_count = data.loc[(data.index.str.contains("k__{}".format(taxon))) & (data[stage] < 0),].shape[0]
                    output.write("{},{},Elevation,{}\n".format(stage, taxon, up_count))
                    output.write("{},{},Depletion,{}\n".format(stage, taxon, down_count))
            else:
                up_count = data.loc[(data[stage]>0), ].shape[0]
                down_count = data.loc[(data[stage] < 0), ].shape[0]
                output.write("{},{},Elevation,{}\n".format(stage, taxon_type, up_count))
                output.write("{},{},Depletion,{}\n".format(stage, taxon_type, down_count))


###
#  8 Actinobacteria
#  7 Ascomycota
# 19 Bacteroidota
#  1 Crenarchaeota
# 18 Euryarchaeota
# 47 Firmicutes
#  6 Proteobacteria
# 22 Uroviricota

outfile_phlyum = os.path.join(workdir, "01up_down_stage_phylum.csv")
taxons = ["Actinobacteria", "Ascomycota", "Bacteroidota", "Crenarchaeota",
          "Euryarchaeota", "Firmicutes", "Proteobacteria", "Uroviricota"]
data = pd.read_csv(infile, header=0, index_col=0)
with open(outfile_phlyum, "w") as output:
    output.write("stage,taxon,change,count\n")
    for taxon in taxons:
        for stage in stages:
            up_count = data.loc[(data.index.str.contains("p__{}".format(taxon))) & (data[stage]>0), ].shape[0]
            down_count = data.loc[(data.index.str.contains("p__{}".format(taxon))) & (data[stage] < 0),].shape[0]
            output.write("{},{},Elevation,{}\n".format(stage, taxon, up_count))
            output.write("{},{},Depletion,{}\n".format(stage, taxon, down_count))

