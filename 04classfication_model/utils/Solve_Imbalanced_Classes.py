from imblearn.over_sampling import RandomOverSampler
from imblearn.over_sampling import SMOTE, ADASYN
from imblearn.under_sampling import RandomUnderSampler
from imblearn.under_sampling import NearMiss
from collections import Counter

def imblearn(X, y, stage, study, sampler):
    X_resampled = None
    y_resampled = None
    if sampler == "None":
        X_resampled, y_resampled = X, y
    elif sampler == "RandomOverSampler":
        ros = RandomOverSampler(random_state=0)
        X_resampled, y_resampled = ros.fit_resample(X, y)
    elif sampler == "SMOTEOverSampler":
        X_resampled, y_resampled = SMOTE().fit_resample(X, y)
    elif sampler == "ADASYNOverSampler":
        X_resampled, y_resampled = ADASYN().fit_resample(X, y)
    elif sampler == "RandomUnderSampler":
        rus = RandomUnderSampler(random_state=0)
        X_resampled, y_resampled = rus.fit_resample(X, y)
    elif sampler == "NearMissUnderSampler":
        nm1 = NearMiss(version=1)
        X_resampled, y_resampled = nm1.fit_resample(X, y)
    ##
    count = sorted(Counter(y_resampled).items())
    print("stage: {}, study: {}, sampler: {}, count: {}".format(stage, study, sampler, count))
    return X_resampled, y_resampled
