# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 18:19:03 2024

@author: Nate Whitsett
"""

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import simpson
from scipy.stats import poisson
font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 12,
        }
font_small = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 12,
        }
data_T1 = pd.read_csv('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Injection_Test_T1.csv', index_col=None)
data_T2 = pd.read_csv('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Injection_Test_T2.csv', index_col=None)
data_T1.sort_values(by=['Sigma'], ascending = True, inplace=True, ignore_index=True)
data_T2.sort_values(by=['Sigma'], ascending = True, inplace=True, ignore_index=True)
Amp = data_T2['Amplitude']
FWHM = data_T2['FWHM']
Error = data_T2['Error']
Bool = data_T2['Accepted?']
False_Pos = data_T2['False_Pos?']
sigma = data_T2['Sigma']

Amp1 = data_T1['Amplitude']
FWHM1 = data_T1['FWHM']
Error1 = data_T1['Error']
Bool1 = data_T1['Accepted?']
False_Pos1 = data_T1['False_Pos?']
sigma1 = data_T1['Sigma']
FN = 0
TP = 0
FP = 0
Amp_Sum = []
val = 0
x = []
y = []
x1 = []
y1 = []
count = 0
x_err_l1 = []
x_err_h1 = []
y_err_l1 = []
y_err_h1 = []

sigma_set = [2, 2.5, 3, 3.5, 4, 5]
counts = 0
for bins in range(len(sigma_set)):
    for index, num in enumerate(sigma1):
        if sigma1[index] == sigma_set[counts]:
            if Bool1[index] == 1 and False_Pos1[index] == 0:
                TP += 1
                count += 1
            if Bool1[index] == 0 and False_Pos1[index] == 0:
                FN += 1
                count += 1
            if Bool1[index] == 0 and False_Pos1[index] == 1:
                FP += 1
                count +=1
    
    precision_list = []
    recall_list = []
    precision = TP/(TP+FP)
    recall = TP/(TP+FN)
    for trials in range(10000):
        TP_t = poisson.rvs(TP)
        FP_t = poisson.rvs(FP)
        FN_t = poisson.rvs(FN)
        recall_list.append(((TP_t)/(TP_t+FN_t)))
        precision_list.append((TP_t)/(TP_t+FP_t))
    precision_list = np.sort(np.array(precision_list))
    recall_list = np.sort(np.array(recall_list))
    y_err_l1.append(np.abs(precision - precision_list[500]))
    y_err_h1.append(np.abs(precision - precision_list[9500]))
    x_err_l1.append(np.abs(recall - recall_list[500]))
    x_err_h1.append(np.abs(recall - recall_list[9500]))
    x1.append((TP)/(TP+FN))
    y1.append((TP)/(TP+FP))
    counts += 1
    TP = 0
    FP = 0
    FN = 0


# upper = 100
# lower = 94
x_err_l = []
x_err_h = []
y_err_l = []
y_err_h = []
counts = 0
for bins in range(len(sigma_set)):
    for index, num in enumerate(sigma):
        if sigma[index] == sigma_set[counts]:
            if Bool[index] == 1 and False_Pos[index] == 0:
                TP += 1
                count += 1
            if Bool[index] == 0 and False_Pos[index] == 0:
                FN += 1
                count += 1
            if Bool[index] == 0 and False_Pos[index] == 1:
                FP += 1
                count +=1
    precision_list = []
    recall_list = []
    precision = TP/(TP+FP)
    recall = TP/(TP+FN)
    for trials in range(10000):
        TP_t = poisson.rvs(TP)
        FP_t = poisson.rvs(FP)
        FN_t = poisson.rvs(FN)
        recall_list.append((TP_t)/(TP_t+FN_t))
        precision_list.append((TP_t)/(TP_t+FP_t))
        
    precision_list = np.sort(np.array(precision_list))
    recall_list = np.sort(np.array(recall_list))
    y_err_l.append(np.abs(precision-precision_list[500]))
    y_err_h.append(np.abs(precision - precision_list[9500]))
    x_err_l.append(np.abs(recall-recall_list[500]))
    x_err_h.append(np.abs(recall-recall_list[9500]))
    x.append((TP)/(TP+FN))
    y.append((TP)/(TP+FP))
    counts +=1 
    TP = 0
    FP = 0
    FN = 0
# print(count)
x_err_h1[0] = 0
x_err_h[0] = 0
x_err1 = [x_err_l1, x_err_h1]
y_err1 = [y_err_l1, y_err_h1]
x_err = [x_err_l, x_err_h]
y_err = [y_err_l, y_err_h]
# print(simpson(y, x=x))
# plt.xlim(0, 1.1)
x = x/np.max(x)
x1 = x1/np.max(x1)
fig, ax = plt.subplots()
ax.set_ylim(0,1.1)
ax.set_xlabel('Recall', fontdict=font)
ax.set_ylabel('Precision', fontdict=font)
ax.text(0.15, 0.2, 'Tier 1 AUC: ' + str(round(-simpson(y1, x=x1)+0.1*0.3, 2)), fontdict=font_small)
ax.text(0.15, 0.125, 'Tier 2 AUC: ' + str(round(-simpson(y, x=x)+0.066*.56, 2)), fontdict=font_small)
# plt.plot(x,y, linestyle='--', marker='s',alpha=0.5, label='Tier2')
# plt.plot(x1,y1, linestyle='--', marker='s',alpha=0.5, label='Tier1')
ax.errorbar(x, y, xerr=x_err, yerr=y_err, linestyle='--', marker='s', alpha=0.5, label='Tier 2')
ax.errorbar(x1, y1, xerr=x_err1, yerr=y_err1, linestyle='--', marker='s', alpha=0.5, label = 'Tier 1')
ax.set_xticklabels(ax.get_xticklabels(), fontdict = font)
ax.set_yticklabels(ax.get_yticklabels(), fontdict = font)
ax.legend(prop={"family":"serif"})
plt.savefig('C:/Users/Nate Whitsett/OneDrive - Washington University in St. Louis/Grad School/Fall 2023/Research/Injection Tests/Precision_Recall/Precision_Recall_Curve.png', dpi=600, bbox_inches="tight")
plt.show()