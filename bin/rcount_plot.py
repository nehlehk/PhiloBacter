#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import argparse
import seaborn as sns
import numpy as np



parser = argparse.ArgumentParser(description='''You did not specify any parameters.''')
parser.add_argument('-b', "--Baci", type=str, help='rcount_BaciSim')
parser.add_argument('-t', "--PB_two", type=str, help='rcount_philobacter_two')
parser.add_argument('-g', "--Gubbins", type=str, help='rcount_gubbins')
parser.add_argument('-c', "--CFML", type=str, help='rcount_CFML')
args = parser.parse_args()


rcount_baci = args.Baci
rcount_PB_two = args.PB_two
rcount_Gubbins = args.Gubbins
rcount_CFML = args.CFML


# rcount_baci = '/home/nehleh/PhiloBacteria/Summary_Results/rcount_baci.csv'
# rcount_PB_two = '/home/nehleh/PhiloBacteria/Summary_Results/rcount_PB_two.csv'
# rcount_Gubbins = '/home/nehleh/PhiloBacteria/Summary_Results/rcount_Gubbins.csv'
# rcount_CFML = '/home/nehleh/PhiloBacteria/Summary_Results/rcount_CFML.csv'


df = pd.read_csv(open(rcount_PB_two, "r"), sep=';', names=['PB_two'], header=None)
df1 = pd.read_csv(open(rcount_Gubbins, "r"), sep=';', names=['Gubbins'], header=None)
df2 = pd.read_csv(open(rcount_CFML, "r"), sep=';', names=['CFML'], header=None)
df3 = pd.read_csv(open(rcount_baci, "r"), sep=';', names=['True'], header=None)


one = abs(df3['True'] - df['PB_two'])
two = abs(df3['True'] - df1['Gubbins'])
three = abs(df3['True'] - df2['CFML'])
new = pd.concat([one,two,three],axis = 1 ,ignore_index=True , names=['PB_two','Gubbins','CFML'])
new = new.set_axis(['PhiloBacter','Gubbins','CFML'], axis=1, inplace=False)
print(new)


fig = plt.figure(figsize=(15,5))
ax = fig.add_subplot(1, 3, 1)
ax = sns.scatterplot(x=df3['True'] , y= df['PB_two'])
z = np.polyfit(df3['True'], df['PB_two'], 1)
p = np.poly1d(z)
plt.plot(df3['True'],p(df3['True']),"r--")
ax.set_title('PhiloBacter' , fontsize=9)
plt.xlabel('Number of simulated recombinations')
plt.ylabel('Number of infered recombinations')


ax = fig.add_subplot(1, 3, 2)
ax = sns.scatterplot(x=df3['True'] , y= df1['Gubbins'])
z = np.polyfit(df3['True'], df1['Gubbins'], 1)
p = np.poly1d(z)
plt.plot(df3['True'],p(df3['True']),"r--")
ax.set_title('Gubbins' , fontsize=9)
plt.xlabel('Number of simulated recombinations')
plt.ylabel('Number of infered recombinations')


ax = fig.add_subplot(1, 3, 3)
ax = sns.scatterplot(x=df3['True'] , y= df2['CFML'])
z = np.polyfit(df3['True'], df2['CFML'], 1)
p = np.poly1d(z)
plt.plot(df3['True'],p(df3['True']),"r--")
ax.set_title('CFML' , fontsize=9)
plt.xlabel('Number of simulated recombinations')
plt.ylabel('Number of infered recombinations')
plt.savefig("rcount_scatter.jpeg")


fig2 = plt.figure(figsize=(5,5))
ax1 = fig2.add_subplot(1, 1, 1)
ax1 = sns.boxplot( data=new )
ax1 = sns.stripplot(data=new,  jitter=True, dodge=True, marker='o', color=".1")
ax1.set_title('Comparison' , fontsize=9)
plt.ylabel('number of recombination events')
plt.savefig("rcount_boxplot.jpeg")
# plt.show()
