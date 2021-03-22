# import packages
import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from operator import itemgetter
from sklearn.model_selection import train_test_split
from sklearn.model_selection import LeaveOneOut, GridSearchCV, KFold, StratifiedKFold
from sklearn.metrics import confusion_matrix
from irf import irf_utils, irf_jupyter_utils
from irf.ensemble import RandomForestClassifierWithWeights


# load data
with open('C:/Users/Mischa/Documents/Uni Masters/Module 5- complex/sample_sheet.csv', 'r') as f:
    idx = np.array([r for r in csv.reader(f, delimiter = ',')])

with open('C:/Users/Mischa/Documents/Uni Masters/Module 5- complex/polar_pos_pqn_imputed_glog.csv', 'r') as f:
    dmx = np.array([r for r in csv.reader(f, delimiter = ',')])
    

# prep
X = dmx[1:,1:].astype(float).T
y = np.array(itemgetter(*idx[1:,1])({'Control': 0, '1x': 0, '10x': 1}))

#X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.3)

# trying a better version of testtrain split
kf = KFold(n_splits=5, random_state=15, shuffle=True)


for count_k,(train_index, test_index) in enumerate(kf.split(X)):
    X_train = X[train_index]
    X_test  = X[test_index]
    y_train = y[train_index]
    y_test  = y[test_index]
    irfres = irf_utils.run_iRF(
        X_train = X_train, X_test = X_test,
        y_train = y_train, y_test = y_test,
        rf = RandomForestClassifierWithWeights(n_estimators = 30),
        K = 10, # number of iteration
        B = 30, # The number of bootstrap samples
        M = 20, # number of trees (RIT) to build
        max_depth = 5,
)



rf_weights, K_iter_rf_data, rf_bootstrap_output, rit_bootstrap_output, stability_score = irfres


# feature importance
fids = dmx[1:,0]

iteration = 'rf_iter5'
impt = K_iter_rf_data[iteration]['feature_importances']
impt_std = K_iter_rf_data[iteration]['feature_importances_std']
impt_rank_idx = K_iter_rf_data[iteration]['feature_importances_rank_idx']

impt, impt_std = impt[impt_rank_idx], impt_std[impt_rank_idx]
impt, impt_std, impt_rank_idx = impt[impt > 0], impt_std[impt > 0], impt_rank_idx[impt > 0]

xids = np.arange(impt.shape[0])
plt.figure()
plt.bar(
    xids, impt, yerr = impt_std, align = 'center'
)
plt.xticks(xids, fids[impt_rank_idx], rotation = 'vertical')
plt.xlim([-1, xids.shape[0]])
plt.show()


for f,i in zip(fids[impt_rank_idx],impt): 
    print(f'feature [{f}] with importance {i:.3f}')

    
# interaction stability scores
for i,s in sorted(stability_score.items(), key = lambda x: -x[1]):
    print(f'interaction [{i}] with stability score {s:.3f}')

irf_jupyter_utils._get_histogram(
    {k: v for k,v in stability_score.items() if v > 0}, 
    sort = True,
)

feature_importance = []
imp = []
for f,i in zip(fids[impt_rank_idx],impt): 
    feature_importance.append(f)
    imp.append(f'{i:.3f}')

feature_stab = []
stability = []  
# interaction stability scores
for i,s in sorted(stability_score.items(), key = lambda x: -x[1]):
    feature_stab.append(i)
    stability.append(f'{s:.3f}')

imp_stab_df = pd.DataFrame(list(zip(feature_importance, imp, feature_stab, stability)))
imp_stab_df.columns = ['feature_importance', 'imp', 'interaction_stab', 'stability']
imp_stab_df.to_csv(r'C:\Users\Mischa\Documents\Uni Masters\Module 5- complex\METABOLOMEimp_stab_rna.csv', index=False)