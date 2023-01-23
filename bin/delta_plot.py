#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import argparse
import seaborn as sns
import numpy as np



parser = argparse.ArgumentParser(description='''You did not specify any parameters.''')
parser.add_argument('-b', "--Baci", type=str, help='delta_BaciSim')
parser.add_argument('-t', "--PB_two", type=str, help='delta_philobacter_two')
parser.add_argument('-g', "--Gubbins", type=str, help='delta_gubbins')
parser.add_argument('-c', "--CFML", type=str, help='delta_CFML')
args = parser.parse_args()


delta_baci = args.Baci
delta_PB_two = args.PB_two
delta_Gubbins = args.Gubbins
delta_CFML = args.CFML


# delta_baci = '/home/nehleh/PhiloBacteria/Summary_Results/delta_baci.csv'
# delta_PB_two = '/home/nehleh/PhiloBacteria/Summary_Results/delta_PB_two.csv'
# delta_Gubbins = '/home/nehleh/PhiloBacteria/Summary_Results/delta_Gubbins.csv'
# delta_CFML = '/home/nehleh/PhiloBacteria/Summary_Results/delta_CFML.csv'


df = pd.read_csv(open(delta_PB_two, "r"))
df1 = pd.read_csv(open(delta_Gubbins, "r"))
df2 = pd.read_csv(open(delta_CFML, "r"))
df3 = pd.read_csv(open(delta_baci, "r"))


fig = plt.figure(figsize=(15,5))
ax = fig.add_subplot(1, 3, 1)
ax = sns.scatterplot(x=df3['len'] , y= df['len'])
ax.set_title('PhiloBacter' , fontsize=9)
plt.xlabel('Truth - length of recombination')
plt.ylabel('Inferred - length of recombination')


ax2 = fig.add_subplot(1, 3, 2)
ax2 = sns.scatterplot(x=df3['len'] , y= df1['len'])
ax2.set_title('Gubbins' , fontsize=9)
plt.xlabel('Truth - length of recombination')
plt.ylabel('Inferred - length of recombination')


ax3 = fig.add_subplot(1, 3, 3)
ax3 = sns.scatterplot(x=df3['len'] , y= df2['length'])
ax3.set_title('CFML' , fontsize=9)
plt.xlabel('Truth - length of recombination')
plt.ylabel('Inferred - length of recombination')
plt.savefig("delta_scatter.jpeg")



# plt.show()
