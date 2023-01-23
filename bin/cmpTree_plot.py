#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import argparse
import seaborn as sns


parser = argparse.ArgumentParser(description='''You did not specify any parameters.''')
parser.add_argument('-p', "--PB_dist",type=str, help='PB_dist')
parser.add_argument('-g', "--Gubb_dist",type=str, help='Gubb_dist')
parser.add_argument('-m', "--CFML_dist",type=str, help='CFML_dist')
parser.add_argument('-prf', "--PB_wrf",type=str, help='wrf_dist')
parser.add_argument('-grf', "--Gubb_wrf",type=str, help='wrf_dist')
parser.add_argument('-mrf', "--CFML_wrf",type=str, help='wrf_dist')
args = parser.parse_args()



dist_PBtree = args.PB_dist
dist_Gubb = args.Gubb_dist
dist_CFML = args.CFML_dist
wrf_PBtree = args.PB_wrf
wrf_Gubb = args.Gubb_wrf
wrf_CFML = args.CFML_wrf

f = open(dist_PBtree, "r")
df = pd.read_csv(f, sep=';', names=['PB'], header=None)
w = open(wrf_PBtree, "r")
wdf = pd.read_csv(w, sep=';', names=['PB'], header=None)

f1 = open(dist_Gubb, "r")
df1 = pd.read_csv(f1, sep=';', names=['Gubbins'], header=None)
w1 = open(wrf_Gubb, "r")
wdf1 = pd.read_csv(w1, sep=';', names=['Gubbins'], header=None)

f2 = open(dist_CFML, "r")
df2 = pd.read_csv(f2, sep=';', names=['CFML'], header=None)
w2 = open(wrf_CFML, "r")
wdf2 = pd.read_csv(w2, sep=';', names=['CFML'], header=None)


all_dist = pd.concat([df2,df1,df], axis=1)
all_wrf = pd.concat([wdf2,wdf1,wdf], axis=1)

fig2 = plt.figure(figsize=(10,5))
ax1 = fig2.add_subplot(1, 2, 1)
ax1 = sns.boxplot(data=all_dist )
ax1 = sns.stripplot(data=all_dist,  jitter=True, dodge=True, marker='o', color=".1")
ax1.set_title('euclidean_distance' , fontsize=9)

ax2 = fig2.add_subplot(1, 2, 2)
ax2 = sns.boxplot(data=all_wrf)
ax2 = sns.stripplot(data=all_wrf,  jitter=True, dodge=True, marker='o', color=".1")
ax2.set_title('RFWeighted' , fontsize=9)


plt.savefig("Dist_summary.jpeg")