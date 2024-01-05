import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold

# def data_preprocess(data, min_abundance, min_prevalance,  feature="function"):
#     """
#     :param data: row with samples and columns with features; relative abundance
#     :param feature: function profile or species profile
#     :return:
#     """
#     # convert to relative abundance
#     data_rel = data.div(data.sum(axis=1), axis=0)  ##row is samples, column is features.
#     data = data_rel
#     #data['sum'] = data.sum(axis=1)
#     if feature == "gf":
#         func_cutoff = min_abundance
#         func_pesudo = 1e-4 * min_abundance
#         filter_data = data.loc[:, (data > func_cutoff).mean(axis=0) >= min_prevalance]
#         feature_list = filter_data.columns.values.tolist()
#         print("filtered featrues: " + str(data.shape[1] - filter_data.shape[1]))
#         # log10-transformed
#         log_data = np.log10(filter_data.astype(np.float64) + func_pesudo)
#         # z-score normalization
#         scale = StandardScaler()
#         z_data = scale.fit_transform(log_data)
#         return z_data,feature_list
#
#     elif feature == "species":
#         tax_cutoff = min_abundance
#         feat_pesudo = 1e-2 * min_abundance
#         ##axis=1是跨行，axis=0是跨列。
#         #data_ = data.loc[:, data.mean(axis=0)>1e-4]
#         filter_data = data.loc[:, (data > tax_cutoff).mean(axis=0) >= min_prevalance] ##至少在3个样本中，丰度超过10.
#         feature_list = filter_data.columns.values.tolist()
#         print("{} - {} = {} features." .format(data.shape[1], data.shape[1]-filter_data.shape[1], filter_data.shape[1]))
#         # log10-transformed
#         log_data = np.log10(filter_data.astype(np.float64) + feat_pesudo)
#         # z-score normalization
#         scale = StandardScaler()
#         z_data = scale.fit_transform(log_data)
#         return z_data, feature_list
#     elif feature == "02diff_species":
#         feat_pesudo = 1e-6
#         filter_data = data
#         feature_list = filter_data.columns.values.tolist()
#         print(
#             "{} - {} = {} features.".format(data.shape[1], data.shape[1] - filter_data.shape[1], filter_data.shape[1]))
#         # log10-transformed
#         log_data = np.log10(filter_data.astype(np.float64) + feat_pesudo)
#         # z-score normalization
#         scale = StandardScaler()
#         z_data = scale.fit_transform(log_data)
#         return z_data, feature_list
#
# ##
def data_multikingdom(data, feat_type, feature="species"):
    """
    :param data: row with samples and columns with features; relative abundance
    :param feature: function profile or species profile
    :return:
    """
    feat_type_dic = {"A":"k__Archaea", "B":"k__Bacteria", "F":"k__Fungi", "V":"k__Viruses"}
    if feat_type in feat_type_dic.keys():
        feat_type_str = feat_type_dic[feat_type]
        data = data.loc[:, data.columns.str.contains(feat_type_str)]
    elif feat_type in ["all", 'metabolites', 'pathways', 'KOs', 'GMMs', 'GBMs']:
        if feat_type == "KOs":
            data = data.loc[:, data.columns.str.startswith('K')]
        elif feat_type == 'pathways':
            data = data.loc[:, data.columns.str.startswith('map')]
        elif feat_type == 'metabolites':
            data = data.loc[:, data.columns.str.startswith(('LI', 'ME', 'MW'))]
        elif feat_type == "GMMs":
            data = data.loc[:, data.columns.str.startswith(('MF'))]
        elif feat_type == "GBMs":
            data = data.loc[:, data.columns.str.startswith(('MGB'))]

    # convert to relative abundance
    data_rel = data.div(data.sum(axis=1), axis=0)  ##row is samples, column is features.
    data = data_rel
    # data['sum'] = data.sum(axis=1)
    if (feature == "KOs") | (feature == "GMMs") | (feature == "GBMs"):
        min_abundance = 1e-6
        min_prevalence = 0.1
        feat_cutoff = min_abundance
        feat_pesudo = 0.01 * min_abundance
        #var_threshold = 0.001
    elif feature == "pathways":
        min_abundance = 1e-6
        min_prevalence = 0.1
        feat_cutoff = min_abundance
        feat_pesudo = 0.01 * min_abundance
        #var_threshold = 0.001
    elif feature == "species":
        min_abundance = 1e-4
        min_prevalence = 0.1
        feat_cutoff = min_abundance
        feat_pesudo = 0.01 * min_abundance
        #var_threshold = 0.001
    elif feature == "metabolites":
        min_abundance = 1e-4
        min_prevalence = 0.1
        feat_cutoff = min_abundance
        feat_pesudo = 0.01 * min_abundance
        #var_threshold = 0.001
    ##axis=1是跨行，axis=0是跨列。
    # data_ = data.loc[:, data.mean(axis=0)>1e-4]
    filter_data = data.loc[:, (data > feat_cutoff).mean(axis=0) >= min_prevalence]  ##至少在3个样本中，丰度超过10.
    feature_list = filter_data.columns.values.tolist()
    # print("{} - {} = {} features.".format(data.shape[1], data.shape[1] - filter_data.shape[1],
    #                                       filter_data.shape[1]))
    # log10-transformed
    log_data = np.log10(filter_data.astype(np.float64) + feat_pesudo)
    # z-score normalization
    scale = StandardScaler()
    z_data = scale.fit_transform(log_data)
    return z_data, feature_list


def data_multikingdom_transfer(train_data, test_data, feat_type, feature="species"):
    """
    :param data: row with samples and columns with features; relative abundance
    :param feature: function profile or species profile
    :return:
    """
    feat_type_dic = {"A":"k__Archaea", "B":"k__Bacteria", "F":"k__Fungi", "V":"k__Viruses"}
    if feat_type in feat_type_dic.keys():
        feat_type_str = feat_type_dic[feat_type]
        train_data = train_data.loc[:, train_data.columns.str.contains(feat_type_str)]
        test_data = test_data.loc[:, test_data.columns.str.contains(feat_type_str)]
    elif feat_type in ["all", 'metabolites', 'pathways', 'KOs', 'GMMs', 'GBMs']:
        if feat_type == "KOs":
            train_data = train_data.loc[:, train_data.columns.str.startswith('K')]
            test_data = test_data.loc[:, test_data.columns.str.startswith('K')]
        elif feat_type == 'pathways':
            train_data = train_data.loc[:, train_data.columns.str.startswith('map')]
            test_data = test_data.loc[:, test_data.columns.str.startswith('map')]
        elif feat_type == 'GMMs':
            train_data = train_data.loc[:, train_data.columns.str.startswith('MF')]
            test_data = test_data.loc[:, test_data.columns.str.startswith('MF')]
        elif feat_type == 'GBMs':
            train_data = train_data.loc[:, train_data.columns.str.startswith('MGB')]
            test_data = test_data.loc[:, test_data.columns.str.startswith('MGB')]
        elif feat_type == 'metabolites':
            train_data = train_data.loc[:, train_data.columns.str.startswith(('LI', 'ME', 'MW'))]
            test_data = test_data.loc[:, test_data.columns.str.startswith(('LI', 'ME', 'MW'))]
    # convert to relative abundance
    train_data_rel = train_data.div(train_data.sum(axis=1), axis=0)  ##row is samples, column is features.
    train_data = train_data_rel

    ##*****************test_data******************##
    # convert to relative abundance
    test_data_rel = test_data.div(test_data.sum(axis=1), axis=0)  ##row is samples, column is features.
    test_data = test_data_rel

    # data['sum'] = data.sum(axis=1)
    if feature in ["KOs", 'pathways', 'GMMs', 'GBMs']:
        min_abundance = 1e-6
        min_prevalence = 0.1
        feat_cutoff = min_abundance
        feat_pesudo = 0.01 * min_abundance
        var_threshold = 0.001
    elif feature == "species":
        min_abundance = 1e-4
        min_prevalence = 0.1
        feat_cutoff = min_abundance
        feat_pesudo = 0.01 * min_abundance
        var_threshold = 0.001
    elif feature == "metabolites":
        min_abundance = 1e-4
        min_prevalence = 0.1
        feat_cutoff = min_abundance
        feat_pesudo = 0.01 * min_abundance
        var_threshold = 0.001
    ##axis=1是跨行，axis=0是跨列。
    # data_ = data.loc[:, data.mean(axis=0)>1e-4]
    train_data_filter = train_data.loc[:, (train_data > feat_cutoff).mean(axis=0) >= min_prevalence]  ##至少在3个样本中，丰度超过10.
    train_feature_list = train_data_filter.columns.values.tolist()

    test_data_filter = test_data.loc[:, train_feature_list]
    # print("{} - {} = {} features.".format(data.shape[1], data.shape[1] - filter_data.shape[1],
    #                                       filter_data.shape[1]))
    # log10-transformed
    train_log_data = np.log10(train_data_filter.astype(np.float64) + feat_pesudo)
    test_log_data = np.log10(test_data_filter.astype(np.float64) + feat_pesudo)

    #z-score normalization
    #The test set must use identical scaling to the training set.
    #https://datascience.stackexchange.com/questions/39932/feature-scaling-both-training-and-test-data
    scale = StandardScaler()
    train_z_data = scale.fit_transform(train_log_data)
    test_z_data = scale.transform(test_log_data)
    ##
    selector = VarianceThreshold(threshold=var_threshold)
    train_z_data = selector.fit_transform(train_z_data)
    test_z_data = selector.transform(test_z_data)
    mask = selector.get_support()
    train_feature_list = [b for a, b in zip(mask, train_feature_list) if a]
    return train_z_data, test_z_data, train_feature_list