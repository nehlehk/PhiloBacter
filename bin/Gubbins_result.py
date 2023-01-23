#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from dendropy import Tree
import dendropy
from sklearn.metrics import mean_squared_error
from statsmodels.genmod.families.links import cloglog
from utility import *
from BCBio import GFF
from dendropy.calculate import treecompare


def set_index(tree,dna):
    for node in tree.postorder_node_iter():
      node.index = -1
      node.annotations.add_bound_attribute("index")

    s = len(tree.leaf_nodes())
    for node in tree.postorder_node_iter():
      if not node.is_leaf():
          node.index = s
          node.label = str(node.index)
          s += 1
      else:
          # node.index = int(node.taxon.label)
          # node.label = str(node.index)
          for idx, name in enumerate(dna):
              if str(name) == str(node.taxon):
                  node.index = idx
                  break
# **********************************************************************************************************************
def Gubbins_recombination(gubbins_log,gubbins_tree,nodes_number,alignment_len,clonal_tree):
    starts = []
    ends = []
    desc = []
    mrca = []
    mrca_clonal = []
    length = []
    gubb_handle = open(gubbins_log)
    for rec in GFF.parse(gubb_handle):
        for feature in rec.features:
            starts.append(feature.location.start)
            ends.append(feature.location.end)
            length.append(feature.location.end - feature.location.start)
            d = str(feature.qualifiers['taxa'][0])
            kids = [int(s) for s in d.split() if s.isdigit()]
            # print(kids)
            desc.append(kids)
            mrca.append(my_mrca(gubbins_tree, kids))
            mrca_clonal.append(my_mrca(clonal_tree, kids))
    gubb_handle.close()

    all_data = {'nodes': desc, 'start': starts, 'end': ends , 'mrca' :mrca , 'mrca_clonal' :mrca_clonal , 'len':length}
    # all_data = {'nodes': desc, 'start': starts, 'end': ends, 'mrca': mrca, 'len': length}
    df = pd.DataFrame(all_data)
    # print((df))
    Gubb_recom_count = len(df)

    GubbData = np.zeros((alignment_len, nodes_number))
    rmseData = np.zeros((alignment_len, tips_num))
    for i in range(len(df)):
        s = int(df['start'][i])
        e = int(df['end'][i])
        # node = int(df['mrca'][i])
        node = int(df['mrca_clonal'][i])
        # print((node))
        if node >= tips_num:
            internal_nodes = df['nodes'][i]
            for j in range(len(internal_nodes)):
                inode = int(internal_nodes[j])
                rmseData[s:e, inode] = 1
        else:
            rmseData[s:e, node] = 1
        GubbData[s:e, node] = 1

    return GubbData,rmseData,df,Gubb_recom_count
# **********************************************************************************************************************
def Gubbins_resultFig(gubbins_tree,GubbData,tips_num,nodes_number,df):
    fig = plt.figure(figsize=(tips_num + 9, tips_num / 2))
    color = ['red', 'green', 'purple', 'blue', 'black']
    for i in range(nodes_number):
        ax = fig.add_subplot(nodes_number, 1, i + 1)
        if i >= tips_num:
            # d = (df.loc[df['mrca'] == i]['nodes'])
            d = (df.loc[df['mrca_clonal'] == i]['nodes'])
            d = d.drop_duplicates()
            if len(d.values) > 0 :
                txt = str(d.values[0])
            else:
                desc = set()
                # txt = str(give_descendents(gubbins_tree, i, desc, tips_num))
                txt = str(give_descendents(clonal_tree, i, desc, tips_num))
                # print(txt)
            # print(i, txt)
            ax.plot(GubbData[:, i], label=str(i) + ' is mrca:' + txt, color=color[i % 5])
        else:
            ax.plot(GubbData[:, i], label=i, color=color[i % 5])
        ax.legend(bbox_to_anchor=(0.045, 1.5), prop={'size': 10})
        ax.set_frame_on(False)
        ax.axis('off')
    ax.axis('on')
    ax.set_yticklabels([])
    # plt.show()
    plt.savefig("Gubbins_Recombination.jpeg")
# **********************************************************************************************************************
def rescale_gubbtree(gubbins_tree,gubb_csv,alignment_len):
    gubbins_tree.scale_edges(1/alignment_len)
    g_tree = gubbins_tree.as_string(schema="newick")
    g = g_tree.replace('\n',"")
    myfile = open('./Gubbinstree_rescale.tree', 'w')
    myfile.write(g)
    myfile.close()
# **********************************************************************************************************************
if __name__ == "__main__":

    # clonal_path = '/home/nehleh/PhiloBacteria/Results/num_3/num_3_Clonaltree.tree'
    # genomefile = '/home/nehleh/PhiloBacteria/Results/num_3/num_3_recom_1_Wholegenome_3_1.fasta'
    # baciSimLog = '/home/nehleh/PhiloBacteria/Results/num_3/num_3_recom_1_BaciSim_Log.txt'
    # gubbins_log = '/home/nehleh/PhiloBacteria/Results/num_3/num_3_recom_1_gubbins.recombination_predictions.gff'
    # gubb_tree = '/home/nehleh/PhiloBacteria/Results/num_3/num_3_recom_1_gubbins.node_labelled.final_tree.tre'
    # gubb_csv = '/home/nehleh/PhiloBacteria/Results/num_3/num_3_recom_1_gubbins.per_branch_statistics.csv'

    parser = argparse.ArgumentParser(description='''You did not specify any parameters.''')
    parser.add_argument('-cl', "--clonaltreeFile", type=str,  help='clona tree from BaciSim')
    parser.add_argument('-a', "--alignmentFile", type=str,  help='fasta file')
    parser.add_argument('-rl', "--recomlogFile", type=str, help='BaciSim recombination log file')
    parser.add_argument('-gl', "--gubblogFile", type=str, help='Gubbins Log File')
    parser.add_argument('-gt', "--gubbtreefile", type=str, help='Gubbins tree File')
    parser.add_argument('-gs', "--gubbcsvfile", type=str, help='Gubbins per_branch_statistics csv File')
    parser.add_argument('-sim', "--simulation", type=int, default=1, help='1 for the simulation data and 0 for emprical sequence')

    args = parser.parse_args()
    genomefile = args.alignmentFile
    gubbins_log = args.gubblogFile
    gubb_tree = args.gubbtreefile
    gubb_csv = args.gubbcsvfile
    simulation = args.simulation

    tns = dendropy.TaxonNamespace()
    alignment = dendropy.DnaCharacterMatrix.get(file=open(genomefile), schema="fasta")
    gubbins_tree = Tree.get_from_path(gubb_tree, 'newick', taxon_namespace=tns)
    set_index(gubbins_tree,alignment)
    nodes_num_g = len(gubbins_tree.nodes())
    tips_num = len(alignment)
    alignment_len = alignment.sequence_size
    rescale_gubbtree(gubbins_tree, gubb_csv, alignment_len)


    if simulation > 0:
        gubbins_tree.reroot_at_midpoint()
        clonal_path = args.clonaltreeFile
        clonal_tree = Tree.get_from_path(clonal_path, 'newick', taxon_namespace=tns)
        set_index(clonal_tree,alignment)
        Gubb_euclidean_distance = treecompare.euclidean_distance(clonal_tree, gubbins_tree, edge_weight_attr="length")
        write_value(Gubb_euclidean_distance, 'Gubb_dist.csv')
        Gubb_WRF_distance = treecompare.weighted_robinson_foulds_distance(clonal_tree, gubbins_tree , edge_weight_attr="length")
        write_value(Gubb_WRF_distance, 'Gubb_WRF_distance.csv')
    if simulation == 1:
        baciSimLog = args.recomlogFile
        nodes_num_c = len(clonal_tree.nodes())
        realData, rmse_real = real_recombination(baciSimLog, clonal_tree, nodes_num_c, alignment_len, tips_num)
        GubbData, rmse_Gubb, df, Gubb_recom_count = Gubbins_recombination(gubbins_log, gubbins_tree, nodes_num_g,alignment_len,clonal_tree)
        Gubbins_resultFig(gubbins_tree, GubbData, tips_num, nodes_num_g, df)
        rmse_real_Gubb = mean_squared_error(rmse_real, rmse_Gubb, squared=False)
        write_value(rmse_real_Gubb, 'RMSE_Gubbins.csv')
        write_value(Gubb_recom_count, 'Gubb_rcount.csv')
        df[['len']].to_csv('Gubbins_delta.csv', index=False)









