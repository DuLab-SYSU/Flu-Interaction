# -*- encoding: utf-8 -*-
'''
@File    :   2_2-Model-influenza_H3N2B.py
@Time    :   2024/4/5 19:51:20
@Author  :   DuLab
@Desc    :   Classification model for influenza H3N2 and B
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



df=pd.read_csv('./data/2_2-Model-H3B.csv',index_col=u'Country Name')
y=df['H3B'].values
df=df.drop(['H3B','H1N1 cluster', 'B betweenness', 'H3N2 betweenness', 'Aging'],axis=1)
#f1
df=df[['GDP growth', 'GDP level', 'CPI', 'Population', 'Children', 'Adults', 'Urban Pop', 'Density', 'Flights', 'Chronic diseases', 'Respiratory diseases', 'Smoking', 'Pet', 'Humidity', 'Forest area', 'H1N1 betweenness', 'H3N2 trunk', 'B trunk', 'H1N1 Fst', 'B gc', 'H3N2 gc', 'B pi', 'H1N1 pi', 'Beds', 'H3N2 infection', 'B infection', 'Flu vaccination']]

print(df.columns)
X=df.values
smote_tomek = SMOTETomek(random_state=0)
X_resampled, y_resampled = smote_tomek.fit_resample(X, y)
cv = KFold(n_splits=5)
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
print(dfpred)
dfpred.to_csv('./H3N2B/H3N2Bprediction.csv') """

# #######################################################################################################
# shap value
# summaryplot
explainer = shap.TreeExplainer(model)
shap_values = explainer.shap_values(df)  
shap.summary_plot(shap_values, X,feature_names=list(df.columns),max_display=18,plot_size=(5.5,8),show=False,plot_type="violin")
plt.rcParams["pdf.fonttype"]=42
plt.rcParams["ps.fonttype"]=42
plt.savefig("./H3N2B/figH3N2Bshap1.pdf",bbox_inches='tight',dpi=600,format='pdf')
plt.close() 

#shap value abs
""" a=pd.DataFrame(shap_values).apply(lambda x:abs(x)).mean()
# print(a)
data2=pd.Series(a,name='shapvalue')
data1=pd.DataFrame()
data1['varible']=df.columns
data1=pd.concat([data1,data2],axis=1)
# print(data1)
data1.to_csv(r'./H3N2B/H3N2Bshapvalue.csv', sep=',', header=True,index=False)
 """
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
for i in range(27):
    for j in range(27):
        data=shap_interaction_values[:, i, j]
        interaction_i=pd.DataFrame([data])
        interaction=interaction.append(interaction_i)
interaction=interaction.T
interaction.columns=new_feature_names
interaction.index=df.index
interaction.to_csv(r'./H3N2B/H3N2Binteractionvalue.csv', sep=',', header=True,index=True)

interaction1=interaction.abs()
#print(interaction1)
mean_interaction=interaction1.mean()
# print(mean_interaction)
mean_interaction.to_csv(r'./H3N2B/H3N2Binteractionvaluemean.csv', sep=',', header=True,index=True)

interactiononly=interaction[['Population-B trunk','Population-H3N2 gc','Urban Pop-B trunk','Flights-B trunk','Adults-Urban Pop','Population-Density','Population-Urban Pop','Population-Flights','Population-Children','H3N2 gc-Beds']]
dfonly=df[['Population','B trunk','H3N2 gc','Urban Pop','Flights','Adults','Density','Children','Beds']]
interaction3d=interactiononly.merge(dfonly, how='inner', left_index=True, right_index=True)
interaction3d.to_csv(r'./H3N2B/H3N2Binteraction3d.csv', sep=',', header=True,index=True) """


# shap individual country
""" shap_values2 = explainer(df)
plt.figure(figsize=(4,3),frameon=True)
shap.plots.bar(shap_values2[37], max_display=10, show_data=True)
plt.rcParams["pdf.fonttype"]=42
plt.rcParams["ps.fonttype"]=42
plt.savefig('./H3N2B/singlecountry/H3N2BNorway0.pdf',bbox_inches = 'tight',dpi=300,format='pdf')
plt.close() """

# shap dependence plot
""" shapsingle = pd.DataFrame()
feature_names1=list(df.columns)
new_feature_names = []
for c1 in feature_names1:
    new_feature_names.append(c1+'_shapvalue')
for i in range(27):
    data=shap_values[:, i]
    single_i=pd.DataFrame([data])
    shapsingle=shapsingle.append(single_i)
shapsingle=shapsingle.T
shapsingle.columns=new_feature_names
shapsingle.index=df.index
shapsingle.to_csv(r'./H3N2B/H3N2Bshapsinglevalue.csv', sep=',', header=True,index=True)
shapsingleonly=shapsingle
dfonly=df
singlesandian=shapsingleonly.merge(dfonly, how='inner', left_index=True, right_index=True)
singlesandian.to_csv(r'./H3N2B/H3N2Bsinglesandian.csv', sep=',', header=True,index=True) """



# #######################################################################################################
# #Feature selection
# model=XGBClassifier()
# rfecv = RFECV(estimator=model,step=1,verbose=1, cv=KFold(n_splits=5, random_state=1,shuffle=True),scoring='f1',min_features_to_select=10)
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
