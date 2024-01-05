import pandas as pd
import os

diff_dir = "/home1/jialh/brain/01meta/multikingdom/03MMUPHin_diff/KOs"
diff_KOs = read.csv(os.path.join(diff_dir, "MMUPHin_features_merge.csv"), header=0)

profile_dir = "/home1/jialh/brain/01meta/multikingdom/00profile"
KOs_desc = read.csv(os.path.join(profile_dir, "KEGG_KOs_to_pathways_metabolism.csv"), header=0)