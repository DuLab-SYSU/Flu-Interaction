# -*- encoding: utf-8 -*-
'''
@File    :   2_1-Model-influenza_A-B.py
@Time    :   2024/4/5 19:40:00
@Author  :   DuLab
@Desc    :   Classification model for influenza A and B
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



df=pd.read_csv('./data/2_1-Model-AB.csv',index_col=u'Country Name')
y=df['AB'].values
df=df.drop(['AB','H1N1 cluster', 'B betweenness', 'H3N2 betweenness', 'Aging'],axis=1)
#f1
df=df[['GDP level', 'CPI', 'Population', 'Children', 'Urban Pop', 'Sex ratio', 'Density', 'Flights', 'Smoking', 'Pet', 'Humidity', 'Forest area', 'Zone', 'B degree', 'H1N1 degree', 'H1N1 trunk', 'H1N1 Fst', 'B Fst', 'H3N2 Fst', 'H1N1 gc', 'H3N2 gc', 'Health expenditure', 'H1N1 infection', 'B infection', 'Flu vaccination']]


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
dfpred.to_csv('./AB/ABprediction.csv') """


# #######################################################################################################
# shap value
# summaryplot
explainer = shap.TreeExplainer(model)
shap_values = explainer.shap_values(df)

shap.summary_plot(shap_values, X,feature_names=list(df.columns),max_display=18,plot_size=(5.5,8),show=False,plot_type="violin")
plt.rcParams["pdf.fonttype"]=42
plt.rcParams["ps.fonttype"]=42
plt.savefig("./AB/fig3ABshap.pdf",bbox_inches='tight',dpi=600,format='pdf')
plt.close() 

#shap value abs
""" a=pd.DataFrame(shap_values).apply(lambda x:abs(x)).mean()
# print(a)
data2=pd.Series(a,name='shapvalue')
data1=pd.DataFrame()
data1['varible']=df.columns
data1=pd.concat([data1,data2],axis=1)
# print(data1)
data1.to_csv(r'./AB/ABshapvalue.csv', sep=',', header=True,index=False) """

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
for i in range(25):
    for j in range(25):
        data=shap_interaction_values[:, i, j]
        interaction_i=pd.DataFrame([data])
        interaction=interaction.append(interaction_i)
interaction=interaction.T
interaction.columns=new_feature_names
interaction.index=df.index
interaction.to_csv(r'./AB/ABinteractionvalue.csv', sep=',', header=True,index=True)

interaction1=interaction.abs()
#print(interaction1)
mean_interaction=interaction1.mean()
# print(mean_interaction)
mean_interaction.to_csv(r'./AB/ABinteractionvaluemean.csv', sep=',', header=True,index=True)

interactiononly=interaction[['Zone-B Fst','Smoking-Density','Population-GDP level','B infection-Density','Forest area-Sex ratio','Population-CPI','Population-Children','Population-H1N1 Fst','Population-H3N2 Fst','Zone-GDP level']]
dfonly=df[['Zone','B Fst','Smoking','Density','Population','GDP level','B infection','Forest area','Sex ratio','CPI','Children','H1N1 Fst','H3N2 Fst']]
interaction3d=interactiononly.merge(dfonly, how='inner', left_index=True, right_index=True)
interaction3d.to_csv(r'./AB/ABinteraction3d.csv', sep=',', header=True,index=True) """


# shap individual country
""" shap_values2 = explainer(df)
plt.figure(figsize=(4,3),frameon=True)
shap.plots.bar(shap_values2[54], max_display=10, show_data=True)
plt.rcParams["pdf.fonttype"]=42
plt.rcParams["ps.fonttype"]=42
plt.savefig('./AB/singlecountry/ABUSA1.pdf',bbox_inches = 'tight',dpi=300,format='pdf')
plt.close()  """

# shap dependence plot
""" shapsingle = pd.DataFrame()
feature_names1=list(df.columns)
new_feature_names = []
for c1 in feature_names1:
    new_feature_names.append(c1+'_shapvalue')
for i in range(25):
    data=shap_values[:, i]
    single_i=pd.DataFrame([data])
    shapsingle=shapsingle.append(single_i)
shapsingle=shapsingle.T
shapsingle.columns=new_feature_names
shapsingle.index=df.index
shapsingle.to_csv(r'./AB/ABshapsinglevalue.csv', sep=',', header=True,index=True)
shapsingleonly=shapsingle
dfonly=df
singlesandian=shapsingleonly.merge(dfonly, how='inner', left_index=True, right_index=True)
singlesandian.to_csv(r'./AB/ABsinglesandian.csv', sep=',', header=True,index=True)
 """



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
