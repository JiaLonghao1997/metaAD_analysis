import pandas as pd
import os
import numpy as np
import sys
sys.path.append("/home1/jialh/brain/01meta/multikingdom/scripts/utils")
from prepare_data import data_multikingdom, data_multikingdom_transfer


def merge_features(features_df, feat_type):
    feature_list = []
    if feat_type in ['metabolites', 'species', 'KOs', 'pathways', 'GMMs', 'GBMs']:
        X, feature_list = data_multikingdom(features_df, feat_type, feature=feat_type)
    elif feat_type == 'ABFV+KOs+pathways':
        ABFV_X, ABFV_feature_list = merge_features(features_df, feat_type="ABFV")
        KOs_X, KOs_feature_list = data_multikingdom(features_df, "KOs", feature="KOs")
        pathways_X, pathways_feature_list = data_multikingdom(features_df, "pathways", feature="pathways")
        X = np.concatenate((ABFV_X, KOs_X, pathways_X), axis=1)
        feature_list = ABFV_feature_list, KOs_feature_list, pathways_feature_list
    else:
        A , A_feats = data_multikingdom(features_df, 'A', feature='species')
        B, B_feats = data_multikingdom(features_df, 'B', feature='species')
        F, F_feats = data_multikingdom(features_df, 'F', feature='species')
        V, V_feats = data_multikingdom(features_df, 'V', feature='species')
        X = None
        if feat_type == 'A':
            X = A
            feature_list = A_feats
        elif feat_type == "B":
            X = B
            feature_list = B_feats
        elif feat_type == "F":
            X = F
            feature_list = F_feats
        elif feat_type == "V":
            X = V
            feature_list = V_feats

        ##two kingdom: 6
        elif feat_type == 'AB':
            X = np.concatenate((A, B), axis=1)
            feature_list = A_feats + B_feats
        elif feat_type == "AF":
            X = np.concatenate((A, F), axis=1)
            feature_list = A_feats + F_feats
        elif feat_type == "AV":
            X = np.concatenate((A, V), axis=1)
            feature_list = A_feats + V_feats
        elif feat_type == "BF":
            X = np.concatenate((B, F), axis=1)
            feature_list = B_feats + F_feats
        elif feat_type == "BV":
            X = np.concatenate((B, V), axis=1)
            feature_list = B_feats + V_feats
        elif feat_type == "FV":
            X = np.concatenate((F, V), axis=1)
            feature_list = F_feats + V_feats
        ##three kingdom: 4
        elif feat_type == "ABF":
            X = np.concatenate((A, B, F), axis=1)
            feature_list = A_feats + B_feats + F_feats
        elif feat_type == "ABV":
            X = np.concatenate((A, B, V), axis=1)
            feature_list = A_feats + B_feats + V_feats
        elif feat_type == "AFV":
            X = np.concatenate((A, F, V), axis=1)
            feature_list = A_feats + F_feats + V_feats
        elif feat_type == "BFV":
            X = np.concatenate((B, F, V), axis=1)
            feature_list = B_feats + F_feats + V_feats
        ##four kingdom: 1
        elif feat_type == "ABFV":
            X = np.concatenate((A, B, F, V), axis=1)
            feature_list = A_feats + B_feats + F_feats + V_feats

    return X, feature_list

def merge_multiple_kingdoms(train_features_df, test_features_df, feat_type):
    train_A, test_A, A_feats = data_multikingdom_transfer(train_features_df, test_features_df, 'A', feature='species')
    train_B, test_B, B_feats = data_multikingdom_transfer(train_features_df, test_features_df, 'B', feature='species')
    train_F, test_F, F_feats = data_multikingdom_transfer(train_features_df, test_features_df, 'F', feature='species')
    train_V, test_V, V_feats = data_multikingdom_transfer(train_features_df, test_features_df, 'V', feature='species')
    train_X = None
    test_X = None
    feature_list = []
    if feat_type == 'A':
        train_X = train_A
        test_X = test_A
        feature_list = A_feats
    elif feat_type == "B":
        train_X = train_B
        test_X = test_B
        feature_list = B_feats
    elif feat_type == "F":
        train_X = train_F
        test_X = test_F
        feature_list = F_feats
    elif feat_type == "V":
        train_X = train_V
        test_X = test_V
        feature_list = V_feats

    ##two kingdom: 6
    elif feat_type == 'AB':
        train_X = np.concatenate((train_A, train_B), axis=1)
        test_X = np.concatenate((test_A, test_B), axis=1)
        feature_list = A_feats + B_feats
    elif feat_type == "AF":
        train_X = np.concatenate((train_A, train_F), axis=1)
        test_X = np.concatenate((test_A, test_F), axis=1)
        feature_list = A_feats + F_feats
    elif feat_type == "AV":
        train_X = np.concatenate((train_A, train_V), axis=1)
        test_X = np.concatenate((test_A, test_V), axis=1)
        feature_list = A_feats + V_feats
    elif feat_type == "BF":
        train_X = np.concatenate((train_B, train_F), axis=1)
        test_X = np.concatenate((test_B, test_F), axis=1)
        feature_list = B_feats + F_feats
    elif feat_type == "BV":
        train_X = np.concatenate((train_B, train_V), axis=1)
        test_X = np.concatenate((test_B, test_V), axis=1)
        feature_list = B_feats + V_feats
    elif feat_type == "FV":
        train_X = np.concatenate((train_F, train_V), axis=1)
        test_X = np.concatenate((test_F, test_V), axis=1)
        feature_list = F_feats + V_feats
    ##three kingdom: 4
    elif feat_type == "ABF":
        train_X = np.concatenate((train_A, train_B, train_F), axis=1)
        test_X = np.concatenate((test_A, test_B, test_F), axis=1)
        feature_list = A_feats + B_feats + F_feats
    elif feat_type == "ABV":
        train_X = np.concatenate((train_A, train_B, train_V), axis=1)
        test_X = np.concatenate((test_A, test_B, test_V), axis=1)
        feature_list = A_feats + B_feats + V_feats
    elif feat_type == "AFV":
        train_X = np.concatenate((train_A, train_F, train_V), axis=1)
        test_X = np.concatenate((test_A, test_F, test_V), axis=1)
        feature_list = A_feats + F_feats + V_feats
    elif feat_type == "BFV":
        train_X = np.concatenate((train_B, train_F, train_V), axis=1)
        test_X = np.concatenate((test_B, test_F, test_V), axis=1)
        feature_list = B_feats + F_feats + V_feats
    ##four kingdom: 1
    elif feat_type == "ABFV":
        train_X = np.concatenate((train_A, train_B, train_F, train_V), axis=1)
        test_X = np.concatenate((test_A, test_B, test_F, test_V), axis=1)
        feature_list = A_feats + B_feats + F_feats + V_feats
    return train_X, test_X, feature_list

#merge_features_transfer
#data_multikingdom_transfer(train_data, test_data, feat_type, feature="species"):
def merge_features_transfer(train_features_df, test_features_df, feat_type):
    #feature_list = []
    if feat_type == "all":
        train_X, test_X, feature_list = data_multikingdom_transfer(train_features_df, test_features_df, 'all', feature='species')
    elif feat_type in ['metabolites', 'metabolites+ConQuR']:
        train_X, test_X, feature_list = data_multikingdom_transfer(train_features_df, test_features_df, 'metabolites', feature='metabolites')
    elif feat_type in ["pathways", "pathways+ConQuR"]:
        train_X, test_X, feature_list = data_multikingdom_transfer(train_features_df, test_features_df, 'pathways', feature='pathways')
    elif feat_type in ["KOs", "KOs+ConQuR"]:
        train_X, test_X, feature_list = data_multikingdom_transfer(train_features_df, test_features_df, 'KOs', feature='KOs')
    elif feat_type in ["GMMs", "GMMs+ConQuR"]:
        train_X, test_X, feature_list = data_multikingdom_transfer(train_features_df, test_features_df, 'GMMs', feature='GMMs')
    elif feat_type in ["GBMs", "GBMs+ConQuR"]:
        train_X, test_X, feature_list = data_multikingdom_transfer(train_features_df, test_features_df, 'GBMs', feature='GBMs')
    elif feat_type in ['ABFV+KOs+pathways', 'ABFV+KOs+pathways+ConQuR']:  ##ABFV+KOs+pathways
        print("deal with {}".format(feat_type))
        train_ABFV_X, test_ABFV_X, ABFV_feature_list = merge_multiple_kingdoms(train_features_df, test_features_df, feat_type="ABFV")
        train_KOs_X, test_KOs_X, KOs_feature_list = data_multikingdom_transfer(train_features_df, test_features_df, 'KOs',  feature="KOs")
        train_pathways_X, test_pathways_X, pathways_feature_list = data_multikingdom_transfer(train_features_df, test_features_df, 'pathways', feature="pathways")
        print("ABFV+KOs+pathways: merge {} species, {} KOs, {} pathways. ".format(len(ABFV_feature_list),
                                                                                  len(KOs_feature_list),
                                                                                  len(pathways_feature_list)))
        train_X = np.concatenate((train_ABFV_X, train_KOs_X, train_pathways_X), axis=1)
        test_X = np.concatenate((test_ABFV_X, test_KOs_X, test_pathways_X), axis=1)
        feature_list = ABFV_feature_list+KOs_feature_list+pathways_feature_list
    else:
        train_X, test_X, feature_list = merge_multiple_kingdoms(train_features_df, test_features_df, feat_type)
    return train_X, test_X, feature_list