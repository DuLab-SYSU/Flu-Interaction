# -*- encoding: utf-8 -*-
'''
@File    :   1_1-Correlation_analysis.py
@Time    :   2024/4/5 16:12:42
@Author  :   DuLab
@Desc    :   Calculate the correlation coefficient between influenza types/subtypes A and B across 55 countries, using the positivity rate per influenza season as the unit.
'''

# here put the import lib

from pandas.core.frame import DataFrame
import pandas as pd
import numpy as np
import seaborn as sns
import scipy.stats as stats
import matplotlib.pyplot as plt

def corr(lis0):
    data = pd.DataFrame(lis0, columns=['AH3','INF_B','ANOTSUBTYPED','SPEC_PROCESSED_NB','INF_A','AH1'])
    #print(data)
    f1 = lambda x: (x / (data['AH3'] + data['AH1'])) * data['ANOTSUBTYPED'] + x
    data1 = data[['AH3', 'AH1']].apply(f1)
    data1['B'] = data['INF_B']
    data1['SPEC_PROCESSED_NB'] = data['SPEC_PROCESSED_NB']
    data1['A']=data['INF_A']
    #data1['A'] = data['AH3']+data['AH1']+data['ANOTSUBTYPED']
    #data1['A'] = data1['AH3']+data1['AH1']
    print(data1)

    f2 = lambda x: x / data1['SPEC_PROCESSED_NB']
    data4 = data1.apply(f2)
    data4 = data4.drop('SPEC_PROCESSED_NB', axis=1)
    print(data4)
    if len(data4) > 2:
        k1_corr = data4.corr('spearman')
        #print(k1_corr)
        i=stats.spearmanr(data4['AH3'],data4['AH1'])
        j=stats.spearmanr(data4['AH3'], data4['B'])
        q=stats.spearmanr(data4['B'], data4['AH1'])
        t=stats.spearmanr(data4['A'], data4['B'])
        namelist.append(name)
        cA_B.append(t[0])
        pA_B.append(t[1])
        cH3_B.append(j[0])
        pH3_B.append(j[1])
        cH1_B.append(q[0])
        pH1_B.append(q[1])
        cH3_H1.append(i[0])
        pH3_H1.append(i[1])


namelist=[]
cA_B=[]
pA_B=[]
cH3_B=[]
pH3_B=[]
cH1_B=[]
pH1_B=[]
cH3_H1=[]
pH3_H1=[]

df=pd.read_csv("./data/1_1-Flu_data.csv")
df=df.replace({'Iran (Islamic Republic of)':'Iran','United Kingdom of Great Britain and Northern Ireland':'UK','United States of America':'America','Republic of Korea':'South Korea',
'Russian Federation':'Russia'})
df['edate']=pd.to_datetime(df.edate)
df=df[['Country','edate','AH1N12009','AH1','AH3','INF_B','ANOTSUBTYPED','SPEC_PROCESSED_NB','INF_A']]
df=df.fillna(0)
# 10-5
northtemp=['America','Austria', 'Belgium','Bulgaria','Canada','Croatia','Denmark','Estonia', 'Finland',  'France',
    'Georgia', 'Germany','Greece','Hungary','Iceland','Iran','Ireland','Israel','Italy','Kazakhstan',
    'Latvia','Luxembourg','Mongolia','Netherlands','Norway', 'Pakistan','Poland','Portugal', 'Ukraine', 
    'Russia', 'Romania','South Korea','Sweden',  'Slovenia','Serbia','Spain','Turkey','UK']
# 4-11
southtemp=['Argentina','Chile','New Zealand']
# according to different countries
tropical=['Algeria','Australia', 'Bangladesh','Brazil','Cambodia','Colombia','China', 'Egypt', 'India', 'Indonesia', 
    'Mexico', 'South Africa',  'Paraguay', 'Thailand']
group1=df.groupby('Country')

for name,group in group1:
    print(name)
    k1 = group.groupby(by=['edate'])[['AH1N12009','AH1','AH3','INF_B','ANOTSUBTYPED','SPEC_PROCESSED_NB','INF_A']].sum()
    k1['A1H1'] = k1['AH1N12009'] + k1['AH1']
    k2 = k1.drop(['AH1N12009', 'AH1'], axis=1)
    k2['AH1'] = k2['A1H1']
    k2 = k2.drop(['A1H1'], axis=1)
    #print('k2',k2)
    if name in northtemp:
    #if name == 'Estonia':
        p = pd.Period('2010/10', 'm')
        #print(p)
        lis0 = []
        while p.year <= k2.index[k2.shape[0] - 1].year:
            k3 = k2[p.strftime('%Y/%m'):(p + 7).strftime('%Y/%m')]
            p += 12
            print(k3)
            F = (k3.sum()).AH3 + (k3.sum()).INF_B + (k3.sum()).AH1 + (k3.sum()).ANOTSUBTYPED
            if F > 30 and float(F) < float((k3.sum()).SPEC_PROCESSED_NB) and (k3.sum()).AH3 + (k3.sum()).AH1 != 0:
                if (k3.sum()).ANOTSUBTYPED / ((k3.sum()).AH3 + (k3.sum()).AH1) < 2 or (k3.sum()).ANOTSUBTYPED / \
                        ((k3.sum()).AH3 + (k3.sum()).AH1) == 2:
                    lis0.append(list(k3.sum()))
                if 2 < (k3.sum()).ANOTSUBTYPED / ((k3.sum()).AH3 + (k3.sum()).AH1) < 5 and (k3.sum()).AH3 != 0 and (
                        k3.sum()).AH1 != 0:
                    lis0.append(list(k3.sum()))
        corr(lis0)
    if name in southtemp:
        p = pd.Period('2010/04', 'm')
        #print(p)
        lis0 = []
        while p.year <= k2.index[k2.shape[0] - 1].year-1:
            
            k3 = k2[p.strftime('%Y/%m'):(p + 7).strftime('%Y/%m')]
            p += 12
            #print(k3)
            F = (k3.sum()).AH3 + (k3.sum()).INF_B + (k3.sum()).AH1 + (k3.sum()).ANOTSUBTYPED
            if F > 30 and float(F) < float((k3.sum()).SPEC_PROCESSED_NB) and (k3.sum()).AH3 + (k3.sum()).AH1 != 0:
                if (k3.sum()).ANOTSUBTYPED / ((k3.sum()).AH3 + (k3.sum()).AH1) < 2 or (k3.sum()).ANOTSUBTYPED / \
                        ((k3.sum()).AH3 + (k3.sum()).AH1) == 2:
                    lis0.append(list(k3.sum()))
                if 2 < (k3.sum()).ANOTSUBTYPED / ((k3.sum()).AH3 + (k3.sum()).AH1) < 5 and (k3.sum()).AH3 != 0 and (
                        k3.sum()).AH1 != 0:
                    lis0.append(list(k3.sum()))
        corr(lis0)
    
    if name == 'Algeria':
        p = pd.Period('2010/11', 'm')
        #print(p)
        lis0 = []
        while p.year <= k2.index[k2.shape[0] - 1].year:
            k3 = k2[p.strftime('%Y/%m'):(p + 4).strftime('%Y/%m')]
            p += 12
            #print(k3)
            F = (k3.sum()).AH3 + (k3.sum()).INF_B + (k3.sum()).AH1 + (k3.sum()).ANOTSUBTYPED
            if F > 30 and float(F) < float((k3.sum()).SPEC_PROCESSED_NB) and (k3.sum()).AH3 + (k3.sum()).AH1 != 0:
                if (k3.sum()).ANOTSUBTYPED / ((k3.sum()).AH3 + (k3.sum()).AH1) < 2 or (k3.sum()).ANOTSUBTYPED / \
                        ((k3.sum()).AH3 + (k3.sum()).AH1) == 2:
                    lis0.append(list(k3.sum()))
                if 2 < (k3.sum()).ANOTSUBTYPED / ((k3.sum()).AH3 + (k3.sum()).AH1) < 5 and (k3.sum()).AH3 != 0 and (
                        k3.sum()).AH1 != 0:
                    lis0.append(list(k3.sum()))
        corr(lis0)
    if name == 'Australia':
        p = pd.Period('2010/06', 'm')
        #print(p)
        lis0 = []
        #print(k2.index[k2.shape[0] - 1].year)
        while p.year <= k2.index[k2.shape[0] - 1].year:
            k3 = k2[p.strftime('%Y/%m'):(p + 4).strftime('%Y/%m')]
            p += 12
            #print(k3)
            F = (k3.sum()).AH3 + (k3.sum()).INF_B + (k3.sum()).AH1 + (k3.sum()).ANOTSUBTYPED
            if F > 30 and float(F) < float((k3.sum()).SPEC_PROCESSED_NB) and (k3.sum()).AH3 + (k3.sum()).AH1 != 0:
                if (k3.sum()).ANOTSUBTYPED / ((k3.sum()).AH3 + (k3.sum()).AH1) < 2 or (k3.sum()).ANOTSUBTYPED / \
                        ((k3.sum()).AH3 + (k3.sum()).AH1) == 2:
                    lis0.append(list(k3.sum()))
                if 2 < (k3.sum()).ANOTSUBTYPED / ((k3.sum()).AH3 + (k3.sum()).AH1) < 5 and (k3.sum()).AH3 != 0 and (
                        k3.sum()).AH1 != 0:
                    lis0.append(list(k3.sum()))
        corr(lis0)
    if name == 'Bangladesh':
        p = pd.Period('2010/03', 'm')
        #print(p)
        lis0 = []
        while p.year <= k2.index[k2.shape[0] - 1].year-1:
            k3 = k2[p.strftime('%Y/%m'):(p + 7).strftime('%Y/%m')]
            p += 12
            #print(k3)
            F = (k3.sum()).AH3 + (k3.sum()).INF_B + (k3.sum()).AH1 + (k3.sum()).ANOTSUBTYPED
            if F > 30 and float(F) < float((k3.sum()).SPEC_PROCESSED_NB) and (k3.sum()).AH3 + (k3.sum()).AH1 != 0:
                if (k3.sum()).ANOTSUBTYPED / ((k3.sum()).AH3 + (k3.sum()).AH1) < 2 or (k3.sum()).ANOTSUBTYPED / \
                        ((k3.sum()).AH3 + (k3.sum()).AH1) == 2:
                    lis0.append(list(k3.sum()))
                if 2 < (k3.sum()).ANOTSUBTYPED / ((k3.sum()).AH3 + (k3.sum()).AH1) < 5 and (k3.sum()).AH3 != 0 and (
                        k3.sum()).AH1 != 0:
                    lis0.append(list(k3.sum()))
        #print(lis0)
        corr(lis0)
    if name == 'Brazil':
        p = pd.Period('2010/03', 'm')
        #print(p)
        lis0 = []
        while p.year <= k2.index[k2.shape[0] - 1].year-1:
            k3 = k2[p.strftime('%Y/%m'):(p + 5).strftime('%Y/%m')]
            p += 12
            #print(k3)
            F = (k3.sum()).AH3 + (k3.sum()).INF_B + (k3.sum()).AH1 + (k3.sum()).ANOTSUBTYPED
            if F > 30 and float(F) < float((k3.sum()).SPEC_PROCESSED_NB) and (k3.sum()).AH3 + (k3.sum()).AH1 != 0:
                if (k3.sum()).ANOTSUBTYPED / ((k3.sum()).AH3 + (k3.sum()).AH1) < 2 or (k3.sum()).ANOTSUBTYPED / \
                        ((k3.sum()).AH3 + (k3.sum()).AH1) == 2:
                    lis0.append(list(k3.sum()))
                if 2 < (k3.sum()).ANOTSUBTYPED / ((k3.sum()).AH3 + (k3.sum()).AH1) < 5 and (k3.sum()).AH3 != 0 and (
                        k3.sum()).AH1 != 0:
                    lis0.append(list(k3.sum()))
        #print(lis0)
        corr(lis0)
    if name == 'Cambodia':
        p = pd.Period('2010/07', 'm')
        #print(p)
        lis0 = []
        while p.year <= k2.index[k2.shape[0] - 1].year-1:
            k3 = k2[p.strftime('%Y/%m'):(p + 5).strftime('%Y/%m')]
            p += 12
            #print(k3)
            F = (k3.sum()).AH3 + (k3.sum()).INF_B + (k3.sum()).AH1 + (k3.sum()).ANOTSUBTYPED
            if F > 30 and float(F) < float((k3.sum()).SPEC_PROCESSED_NB) and (k3.sum()).AH3 + (k3.sum()).AH1 != 0:
                if (k3.sum()).ANOTSUBTYPED / ((k3.sum()).AH3 + (k3.sum()).AH1) < 2 or (k3.sum()).ANOTSUBTYPED / \
                        ((k3.sum()).AH3 + (k3.sum()).AH1) == 2:
                    lis0.append(list(k3.sum()))
                if 2 < (k3.sum()).ANOTSUBTYPED / ((k3.sum()).AH3 + (k3.sum()).AH1) < 5 and (k3.sum()).AH3 != 0 and (
                        k3.sum()).AH1 != 0:
                    lis0.append(list(k3.sum()))
        #print(lis0)
        corr(lis0)
    if name == 'Colombia':
        p = pd.Period('2010/04', 'm')
        #print(p)
        lis0 = []
        while p.year <= k2.index[k2.shape[0] - 1].year-1:
            k3 = k2[p.strftime('%Y/%m'):(p + 5).strftime('%Y/%m')]
            p += 12
            #print(k3)
            F = (k3.sum()).AH3 + (k3.sum()).INF_B + (k3.sum()).AH1 + (k3.sum()).ANOTSUBTYPED
            if F > 30 and float(F) < float((k3.sum()).SPEC_PROCESSED_NB) and (k3.sum()).AH3 + (k3.sum()).AH1 != 0:
                if (k3.sum()).ANOTSUBTYPED / ((k3.sum()).AH3 + (k3.sum()).AH1) < 2 or (k3.sum()).ANOTSUBTYPED / \
                        ((k3.sum()).AH3 + (k3.sum()).AH1) == 2:
                    lis0.append(list(k3.sum()))
                if 2 < (k3.sum()).ANOTSUBTYPED / ((k3.sum()).AH3 + (k3.sum()).AH1) < 5 and (k3.sum()).AH3 != 0 and (
                        k3.sum()).AH1 != 0:
                    lis0.append(list(k3.sum()))
        #print(lis0)
        corr(lis0)

    if name == 'Egypt':
        p = pd.Period('2010/12', 'm')
        #print(p)
        lis0 = []
        while p.year <= k2.index[k2.shape[0] - 1].year:
            k3 = k2[p.strftime('%Y/%m'):(p + 4).strftime('%Y/%m')]
            p += 12
            #print(k3)
            F = (k3.sum()).AH3 + (k3.sum()).INF_B + (k3.sum()).AH1 + (k3.sum()).ANOTSUBTYPED
            if F > 30 and float(F) < float((k3.sum()).SPEC_PROCESSED_NB) and (k3.sum()).AH3 + (k3.sum()).AH1 != 0:
                if (k3.sum()).ANOTSUBTYPED / ((k3.sum()).AH3 + (k3.sum()).AH1) < 2 or (k3.sum()).ANOTSUBTYPED / \
                        ((k3.sum()).AH3 + (k3.sum()).AH1) == 2:
                    lis0.append(list(k3.sum()))
                if 2 < (k3.sum()).ANOTSUBTYPED / ((k3.sum()).AH3 + (k3.sum()).AH1) < 5 and (k3.sum()).AH3 != 0 and (
                        k3.sum()).AH1 != 0:
                    lis0.append(list(k3.sum()))
        #print(lis0)
        corr(lis0)
    if name == 'Mexico':
        p = pd.Period('2010/11', 'm')
        #print(p)
        lis0 = []
        while p.year <= k2.index[k2.shape[0] - 1].year:
            k3 = k2[p.strftime('%Y/%m'):(p + 5).strftime('%Y/%m')]
            p += 12
            #print(k3)
            F = (k3.sum()).AH3 + (k3.sum()).INF_B + (k3.sum()).AH1 + (k3.sum()).ANOTSUBTYPED
            if F > 30 and float(F) < float((k3.sum()).SPEC_PROCESSED_NB) and (k3.sum()).AH3 + (k3.sum()).AH1 != 0:
                if (k3.sum()).ANOTSUBTYPED / ((k3.sum()).AH3 + (k3.sum()).AH1) < 2 or (k3.sum()).ANOTSUBTYPED / \
                        ((k3.sum()).AH3 + (k3.sum()).AH1) == 2:
                    lis0.append(list(k3.sum()))
                if 2 < (k3.sum()).ANOTSUBTYPED / ((k3.sum()).AH3 + (k3.sum()).AH1) < 5 and (k3.sum()).AH3 != 0 and (
                        k3.sum()).AH1 != 0:
                    lis0.append(list(k3.sum()))
        #print(lis0)
        corr(lis0)
    if name == 'Indonesia':
        p = pd.Period('2010/09', 'm')
        #print(p)
        lis0 = []
        while p.year <= k2.index[k2.shape[0] - 1].year:
            k3 = k2[p.strftime('%Y/%m'):(p + 9).strftime('%Y/%m')]
            p += 12
            #print(k3)
            F = (k3.sum()).AH3 + (k3.sum()).INF_B + (k3.sum()).AH1 + (k3.sum()).ANOTSUBTYPED
            if F > 30 and float(F) < float((k3.sum()).SPEC_PROCESSED_NB) and (k3.sum()).AH3 + (k3.sum()).AH1 != 0:
                if (k3.sum()).ANOTSUBTYPED / ((k3.sum()).AH3 + (k3.sum()).AH1) < 2 or (k3.sum()).ANOTSUBTYPED / \
                        ((k3.sum()).AH3 + (k3.sum()).AH1) == 2:
                    lis0.append(list(k3.sum()))
                if 2 < (k3.sum()).ANOTSUBTYPED / ((k3.sum()).AH3 + (k3.sum()).AH1) < 5 and (k3.sum()).AH3 != 0 and (
                        k3.sum()).AH1 != 0:
                    lis0.append(list(k3.sum()))
        #print(lis0)
        corr(lis0)
    if name == 'Paraguay':
        p = pd.Period('2010/04', 'm')
        #print(p)
        lis0 = []
        while p.year <= k2.index[k2.shape[0] - 1].year-1:
            k3 = k2[p.strftime('%Y/%m'):(p + 6).strftime('%Y/%m')]
            p += 12
            #print(k3)
            F = (k3.sum()).AH3 + (k3.sum()).INF_B + (k3.sum()).AH1 + (k3.sum()).ANOTSUBTYPED
            if F > 30 and float(F) < float((k3.sum()).SPEC_PROCESSED_NB) and (k3.sum()).AH3 + (k3.sum()).AH1 != 0:
                if (k3.sum()).ANOTSUBTYPED / ((k3.sum()).AH3 + (k3.sum()).AH1) < 2 or (k3.sum()).ANOTSUBTYPED / \
                        ((k3.sum()).AH3 + (k3.sum()).AH1) == 2:
                    lis0.append(list(k3.sum()))
                if 2 < (k3.sum()).ANOTSUBTYPED / ((k3.sum()).AH3 + (k3.sum()).AH1) < 5 and (k3.sum()).AH3 != 0 and (
                        k3.sum()).AH1 != 0:
                    lis0.append(list(k3.sum()))
        #print(lis0)
        corr(lis0)
    if name == 'South Africa':
        p = pd.Period('2010/05', 'm')
        #print(p)
        lis0 = []
        while p.year <= k2.index[k2.shape[0] - 1].year-1:
            k3 = k2[p.strftime('%Y/%m'):(p + 5).strftime('%Y/%m')]
            p += 12
            #print(k3)
            F = (k3.sum()).AH3 + (k3.sum()).INF_B + (k3.sum()).AH1 + (k3.sum()).ANOTSUBTYPED
            if F > 30 and float(F) < float((k3.sum()).SPEC_PROCESSED_NB) and (k3.sum()).AH3 + (k3.sum()).AH1 != 0:
                if (k3.sum()).ANOTSUBTYPED / ((k3.sum()).AH3 + (k3.sum()).AH1) < 2 or (k3.sum()).ANOTSUBTYPED / \
                        ((k3.sum()).AH3 + (k3.sum()).AH1) == 2:
                    lis0.append(list(k3.sum()))
                if 2 < (k3.sum()).ANOTSUBTYPED / ((k3.sum()).AH3 + (k3.sum()).AH1) < 5 and (k3.sum()).AH3 != 0 and (
                        k3.sum()).AH1 != 0:
                    lis0.append(list(k3.sum()))
        #print(lis0)
        corr(lis0)

    if name == 'Thailand':
        p = pd.Period('2010/06', 'm')
        #print(p)
        lis0 = []
        while p.year <= k2.index[k2.shape[0] - 1].year:
            #print(k2)
            k3 = k2[p.strftime('%Y/%m'):(p + 5).strftime('%Y/%m')]
            k31=k2[(p+6).strftime('%Y/%m'):(p + 11).strftime('%Y/%m')]
            #print(k3)
            p += 12
            F = (k3.sum()).AH3 + (k3.sum()).INF_B + (k3.sum()).AH1 + (k3.sum()).ANOTSUBTYPED
            if F > 30 and float(F) < float((k3.sum()).SPEC_PROCESSED_NB) and (k3.sum()).AH3 + (k3.sum()).AH1 != 0:
                if (k3.sum()).ANOTSUBTYPED / ((k3.sum()).AH3 + (k3.sum()).AH1) < 2 or (k3.sum()).ANOTSUBTYPED / \
                        ((k3.sum()).AH3 + (k3.sum()).AH1) == 2:
                    lis0.append(list(k3.sum()))
                if 2 < (k3.sum()).ANOTSUBTYPED / ((k3.sum()).AH3 + (k3.sum()).AH1) < 5 and (k3.sum()).AH3 != 0 and (
                        k3.sum()).AH1 != 0:
                    lis0.append(list(k3.sum()))
            
            #print(k31)
            F1 = (k31.sum()).AH3 + (k31.sum()).INF_B + (k31.sum()).AH1 + (k31.sum()).ANOTSUBTYPED
            if F1 > 30 and float(F1) < float((k3.sum()).SPEC_PROCESSED_NB) and (k31.sum()).AH3 + (k31.sum()).AH1 != 0:
                if (k31.sum()).ANOTSUBTYPED / ((k31.sum()).AH3 + (k31.sum()).AH1) < 2 or (k31.sum()).ANOTSUBTYPED / \
                        ((k31.sum()).AH3 + (k31.sum()).AH1) == 2:
                    lis0.append(list(k31.sum()))
                if 2 < (k31.sum()).ANOTSUBTYPED / ((k31.sum()).AH3 + (k31.sum()).AH1) < 5 and (k31.sum()).AH3 != 0 and (
                        k31.sum()).AH1 != 0:
                    lis0.append(list(k31.sum()))
        corr(lis0)
    if name == 'India':
        p = pd.Period('2010/02', 'm')
        #print(p)
        lis0 = []
        while p.year <= k2.index[k2.shape[0] - 1].year:
            #print(k2)
            k3 = k2[p.strftime('%Y/%m'):(p + 2).strftime('%Y/%m')]
            k31=k2[(p+5).strftime('%Y/%m'):(p + 10).strftime('%Y/%m')]
            #k3=pd.concat([k30,k31],axis=0)
            #print(k3)
            p += 12
            F = (k3.sum()).AH3 + (k3.sum()).INF_B + (k3.sum()).AH1 + (k3.sum()).ANOTSUBTYPED
            if F > 30 and float(F) < float((k3.sum()).SPEC_PROCESSED_NB) and (k3.sum()).AH3 + (k3.sum()).AH1 != 0:
                if (k3.sum()).ANOTSUBTYPED / ((k3.sum()).AH3 + (k3.sum()).AH1) < 2 or (k3.sum()).ANOTSUBTYPED / \
                        ((k3.sum()).AH3 + (k3.sum()).AH1) == 2:
                    lis0.append(list(k3.sum()))
                if 2 < (k3.sum()).ANOTSUBTYPED / ((k3.sum()).AH3 + (k3.sum()).AH1) < 5 and (k3.sum()).AH3 != 0 and (
                        k3.sum()).AH1 != 0:
                    lis0.append(list(k3.sum()))
            F1 = (k31.sum()).AH3 + (k31.sum()).INF_B + (k31.sum()).AH1 + (k31.sum()).ANOTSUBTYPED
            if F1 > 30 and float(F1) < float((k3.sum()).SPEC_PROCESSED_NB) and (k31.sum()).AH3 + (k31.sum()).AH1 != 0:
                if (k31.sum()).ANOTSUBTYPED / ((k31.sum()).AH3 + (k31.sum()).AH1) < 2 or (k31.sum()).ANOTSUBTYPED / \
                        ((k31.sum()).AH3 + (k31.sum()).AH1) == 2:
                    lis0.append(list(k31.sum()))
                if 2 < (k31.sum()).ANOTSUBTYPED / ((k31.sum()).AH3 + (k31.sum()).AH1) < 5 and (k31.sum()).AH3 != 0 and (
                        k31.sum()).AH1 != 0:
                    lis0.append(list(k31.sum()))
        corr(lis0)
    if name == 'China':
        p = pd.Period('2010/05', 'm')
        #print(p)
        lis0 = []
        while p.year <= k2.index[k2.shape[0] - 1].year-1:
            #print(k2)
            k3 = k2[p.strftime('%Y/%m'):(p + 5).strftime('%Y/%m')]
            k31=k2[(p+6).strftime('%Y/%m'):(p + 11).strftime('%Y/%m')]
            #k3=pd.concat([k30,k31],axis=0)
            print(k3)
            print(k31)
            p += 12
            F = (k3.sum()).AH3 + (k3.sum()).INF_B + (k3.sum()).AH1 + (k3.sum()).ANOTSUBTYPED
            if F > 30 and float(F) < float((k3.sum()).SPEC_PROCESSED_NB) and (k3.sum()).AH3 + (k3.sum()).AH1 != 0:
                if (k3.sum()).ANOTSUBTYPED / ((k3.sum()).AH3 + (k3.sum()).AH1) < 2 or (k3.sum()).ANOTSUBTYPED / \
                        ((k3.sum()).AH3 + (k3.sum()).AH1) == 2:
                    lis0.append(list(k3.sum()))
                if 2 < (k3.sum()).ANOTSUBTYPED / ((k3.sum()).AH3 + (k3.sum()).AH1) < 5 and (k3.sum()).AH3 != 0 and (
                        k3.sum()).AH1 != 0:
                    lis0.append(list(k3.sum()))
            F1 = (k31.sum()).AH3 + (k31.sum()).INF_B + (k31.sum()).AH1 + (k31.sum()).ANOTSUBTYPED
            if F1 > 30 and float(F1) < float((k3.sum()).SPEC_PROCESSED_NB) and (k31.sum()).AH3 + (k31.sum()).AH1 != 0:
                if (k31.sum()).ANOTSUBTYPED / ((k31.sum()).AH3 + (k31.sum()).AH1) < 2 or (k31.sum()).ANOTSUBTYPED / \
                        ((k31.sum()).AH3 + (k31.sum()).AH1) == 2:
                    lis0.append(list(k31.sum()))
                if 2 < (k31.sum()).ANOTSUBTYPED / ((k31.sum()).AH3 + (k31.sum()).AH1) < 5 and (k31.sum()).AH3 != 0 and (
                        k31.sum()).AH1 != 0:
                    lis0.append(list(k31.sum()))
        corr(lis0)


c = {'Country': namelist,'A-B': cA_B, 'A-B(p)': pA_B, 'AH3-AH1': cH3_H1, 'AH1-AH3(p)': pH3_H1,
'AH3-B': cH3_B,'AH3-B(p)': pH3_B,'AH1-B': cH1_B, 'AH1-B(p)': pH1_B}
data5 = DataFrame(c)
print(data5)
data5.to_csv('./data/country_corr.csv')
