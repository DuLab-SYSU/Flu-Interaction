# -*- encoding: utf-8 -*-
'''
@File    :   2_3-Model-influenza_H1N1B.py
@Time    :   2024/4/5 19:51:54
@Author  :   DuLab
@Desc    :   Classification model for influenza H1N1 and B
'''

# here put the import lib

import matplotlib.pyplot as plt
plt.switch_backend('agg')
from sklearn.feature_selection import RFECV
import warnings
warnings.filterwarnings('ignore')
from matplotlib import pyplot
import sklearn
import xgboost as xgb
from sklearn.metrics import classification_report,roc_curve, auc
from sklearn.model_selection import KFold,cross_val_predict,cross_val_score,GridSearchCV
from xgboost import plot_importance, XGBClassifier
from xgboost import plot_tree
from imblearn.combine import SMOTETomek
import pandas as pd
from collections import Counter
from featexp import *
import numpy as np
import shap
import os
shap.initjs()
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.feature_selection import SelectFromModel
from sklearn.linear_model import LogisticRegression

from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap



df=pd.read_csv('./data/2_3-Model-H1B.csv',index_col=u'Country Name')
y=df['H1B'].values
df=df.drop(['H1B', 'H1N1 cluster', 'B betweenness', 'H3N2 betweenness', 'Aging'],axis=1)

df=df[['GDP level', 'CPI', 'Children', 'Adults', 'Sex ratio', 'Density', 'Flights', 'Chronic diseases', 'Respiratory diseases', 'Forest area', 'H3N2 cluster', 'H1N1 trunk', 'B gc', 'B pi', 'Duration', 'H3N2 infection', 'B infection', 'Flu vaccination']]

print(df.columns)
X=df.values
smote_tomek = SMOTETomek(random_state=0)
X_resampled, y_resampled = smote_tomek.fit_resample(X, y)

cv = KFold(n_splits=5)
# model=XGBClassifier(gamma=0, learning_rate=0.01, max_depth=3, min_child_weight=1, n_estimators=100)
model=XGBClassifier()
model.fit(X_resampled, y_resampled)
print(model)
accuracy = cross_val_score(model, X_resampled, y_resampled, scoring='accuracy', cv=5)
precision = cross_val_score(model, X_resampled, y_resampled, scoring='precision', cv=5)
recall = cross_val_score(model,X_resampled, y_resampled, scoring='recall', cv=5)
f1_score = cross_val_score(model,X_resampled, y_resampled, scoring='f1',cv=5)
auc = cross_val_score(model, X_resampled, y_resampled, scoring='roc_auc', cv=5)
print("Accuracy:",accuracy.mean())
print("Precision:",precision.mean())
print("Recall:",recall.mean())
print("F1_score:",f1_score.mean())
print("AUC:",auc.mean())


# prediction
""" y_pred = model.predict(df)
predictions = [round(value) for value in y_pred]
dfpred=pd.DataFrame()
dfpred=pd.concat([dfpred, pd.DataFrame(y,columns=['y'])],axis=1)
dfpred=pd.concat([dfpred, pd.DataFrame(y_pred,columns=['y_pred'])],axis=1)
#print(dfpred)
dfpred.to_csv('./H1N1B/H1N1Bprediction.csv') """


# #######################################################################################################
# shap value
# summaryplot
""" explainer = shap.TreeExplainer(model)
shap_values = explainer.shap_values(df)  
shap.summary_plot(shap_values, X,feature_names=list(df.columns),max_display=18,plot_size=(5.5,8),show=False,plot_type="violin")
plt.rcParams["pdf.fonttype"]=42
plt.rcParams["ps.fonttype"]=42
plt.savefig("./H1N1B/fig3H1N1Bshap.pdf",bbox_inches='tight',dpi=600,format='pdf')
plt.close() """ 

#shap value abs
""" a=pd.DataFrame(shap_values).apply(lambda x:abs(x)).mean()
# print(a)
data2=pd.Series(a,name='shapvalue')
data1=pd.DataFrame()
data1['varible']=df.columns
data1=pd.concat([data1,data2],axis=1)
# print(data1)
data1.to_csv(r'./H1N1B/H1N1Bshapvalue.csv', sep=',', header=True,index=False) """

# shap interaction values
""" shap_interaction_values = explainer.shap_interaction_values(df)
interaction = pd.DataFrame()
feature_names=list(df.columns)
new_feature_names = []
for c1 in feature_names:
    for c2 in feature_names:
        if c1 == c2:
            new_feature_names.append(c1)
        else:
            new_feature_names.append(c1 + "-" + c2)
for i in range(18):
    for j in range(18):
        data=shap_interaction_values[:, i, j]
        interaction_i=pd.DataFrame([data])
        interaction=interaction.append(interaction_i)
interaction=interaction.T
interaction.columns=new_feature_names
interaction.index=df.index
interaction.to_csv(r'./H1N1B/H1N1Binteractionvalue.csv', sep=',', header=True,index=True)

interaction1=interaction.abs()
#print(interaction1)
mean_interaction=interaction1.mean()
# print(mean_interaction)
mean_interaction.to_csv(r'./H1N1B/H1N1Binteractionvaluemean.csv', sep=',', header=True,index=True)

interactiononly=interaction[['Respiratory diseases-B gc','Density-Forest area','CPI-Density','B infection-Density','Density-Children','B gc-CPI','B gc-Sex ratio','Flights-Density','CPI-B infection','Duration-B gc']]
dfonly=df[['Respiratory diseases','B gc','Density','Forest area','CPI','B infection','Children','Sex ratio','Flights','Children','Duration']]
interaction3d=interactiononly.merge(dfonly, how='inner', left_index=True, right_index=True)
interaction3d.to_csv(r'./H1N1B/H1N1Binteraction3d.csv', sep=',', header=True,index=True) """


# shap individual country
""" shap_values2 = explainer(df)
plt.figure(figsize=(12,6),frameon=True)
shap.plots.bar(shap_values2[9], max_display=10, show_data=True)
plt.rcParams["pdf.fonttype"]=42
plt.rcParams["ps.fonttype"]=42
plt.savefig('./H1N1B/singlecountry/H1N1Bcanada0.pdf',bbox_inches = 'tight',dpi=300,format='pdf')
plt.close() """

# shap dependence plot
""" shapsingle = pd.DataFrame()
feature_names1=list(df.columns)
new_feature_names = []
for c1 in feature_names1:
    new_feature_names.append(c1+'_shapvalue')
for i in range(18):
    data=shap_values[:, i]
    single_i=pd.DataFrame([data])
    shapsingle=shapsingle.append(single_i)
shapsingle=shapsingle.T
shapsingle.columns=new_feature_names
shapsingle.index=df.index
shapsingle.to_csv(r'./H1N1B/H1N1Bshapsinglevalue.csv', sep=',', header=True,index=True)
shapsingleonly=shapsingle
dfonly=df
singlesandian=shapsingleonly.merge(dfonly, how='inner', left_index=True, right_index=True)
singlesandian.to_csv(r'./H1N1B/H1N1Bsinglesandian.csv', sep=',', header=True,index=True) """


# #######################################################################################################
# #Feature selection
# model=XGBClassifier()
# rfecv = RFECV(estimator=model,step=1,verbose=1, cv=KFold(n_splits=5, random_state=1,shuffle=True),scoring='roc_auc',min_features_to_select=10)
# rfecv.fit(X_resampled, y_resampled)
# print("Number of features : %s" % rfecv.n_features_)
# print("Ranking of features : %s" % rfecv.ranking_)
# print("support of features : %s" % rfecv.support_)
# flist=df.columns.tolist()
# ranklist=rfecv.ranking_.tolist()
# c = {'a': flist,'b': ranklist}
# data1 = pd.DataFrame(c)
# data2=data1[data1['b'].isin([1])]
# featurelist=data2['a'].tolist()
# print(featurelist)

#######################################################################################################
#GridSearch parameters
#parameters= [{'learning_rate':[0.01,0.05,0.1,0.2,0.3],'n_estimators':[50,70,100,200,300],'max_depth':[i for i in range(3,10,1)],'min_child_weight':[i for i in range(1,6,1)],'gamma':[0,0.05,0.1,0.3,0.5,0.7,0.9,1]}]
# parameters= [{'learning_rate':[0.01,0.05,0.1,0.2,0.3],'n_estimators':[50,100,200,300,400],'max_depth':[i for i in range(3,10,1)],'min_child_weight':[i for i in range(1,6,1)],'gamma':[0,0.05,0.1,0.3,0.5,0.7,0.9,1]}]
# gsearch = GridSearchCV(model, param_grid=parameters, scoring='roc_auc', cv=5)
# gsearch.fit(X_resampled, y_resampled)
# print(gsearch.cv_results_['mean_test_score'],gsearch.cv_results_['params'])
# print(gsearch.best_params_)
# print("Best score: %0.3f" % gsearch.best_score_)
