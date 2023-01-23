#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from dendropy import Tree
import dendropy
from sklearn.metrics import mean_squared_error
from utility import *
from dendropy.calculate import treecompare



def set_index(tree,dna):
    for node in tree.postorder_node_iter():
      node.index = -1
      node.annotations.add_bound_attribute("index")

    s = len(tree.leaf_nodes())
    for node in tree.postorder_node_iter():
      if not node.is_leaf():
          node.index = s
          node.label = str(node.label)
          s += 1
      else:
          # node.index = int(node.taxon.label)
          # node.label = str(node.index)
          for idx, name in enumerate(dna):
              if str(name) == str(node.taxon):
                  node.index = idx
                  break
# **********************************************************************************************************************
def give_descendents_CFML(tree,node_label,result):
    if "NODE" in str(node_label):
        internal_recom_node = tree.find_node_with_label(node_label)
        # internal_recom_node = tree.find_node_with_label(node_label[5:])
        children = internal_recom_node.child_nodes()
        for n in range(len(children)):
            if not (children[n].taxon is None):
                r_node= (children[n].taxon.label)
            else:
                r_node = children[n].label
            if "NODE" in str(r_node):
                give_descendents_CFML(tree,r_node,result)
            else:
                result.add(r_node)
    return result
# **********************************************************************************************************************
def CFML_resultFig(tree,CFMLData):
    fig = plt.figure(figsize=(tips_num + 9, tips_num / 2))
    color = ['red', 'green', 'purple', 'blue', 'black']
    taxa = CFMLData.shape[1]
    for i in range(taxa):
        ax = fig.add_subplot(taxa, 1, i + 1)
        if i >= tips_num:
            node = int(my_label[my_index.index(i)])
            label = str('NODE '+ str(node))
            desc = set()
            d = give_descendents_CFML(tree, label, desc)
            ax.plot(CFMLData[:, i], label=label + ' is mrca:' + str(d), color=color[i % 5])
        else:
            ax.plot(CFMLData[:, i], label=i, color=color[i % 5])
        ax.legend(bbox_to_anchor=(0.045, 1.5), prop={'size': 10})
        # ax.plot(CFMLData[:, i],label=i, color=color[i % 5])
        # ax.legend(bbox_to_anchor=(0.04, 1.33) ,prop={'size':10} )
        ax.set_frame_on(False)
        ax.axis('off')
    ax.axis('on')
    ax.set_yticklabels([])
    # plt.show()
    plt.savefig("CFML_Recombination.jpeg")
# **********************************************************************************************************************
def set_to_list(set):
    taxon_label = []
    for elm in set:
        taxon_label.append(""+str(elm)+"")
    return taxon_label
    # **********************************************************************************************************************
def CFML_recombination(CFML_recomLog,cfml_tree,tips_num):
    CFMLData = np.zeros((alignment_len, nodes_number))
    rmseData = np.zeros((alignment_len, tips_num))
    df = pd.read_csv(CFML_recomLog, sep='\t', engine='python')
    CFML_recom_count=len(df)
    for i in range(len(df)):
        s = df['Beg'][i]
        e = df['End'][i]
        node = df['Node'][i]
        if "NODE_" in str(node):
            temp_node = (node[5:])
            node = my_index[(my_label.index(temp_node))]
            desc = set()
            d = give_descendents_CFML(cfml_tree,"NODE "+str(temp_node),desc)
            for elm in d:
                rmseData[s:e, int(elm)] = 1
        else:
            node = my_index[my_label.index(str(node))]
            rmseData[s:e, int(node)] = 1
        CFMLData[s:e, int(node)] = 1


    return CFMLData,rmseData,CFML_recom_count,df
# **********************************************************************************************************************


if __name__ == "__main__":
    # clonal_path = '/home/nehleh/PhiloBacteria/Results/num_1/num_1_Clonaltree.tree'
    # genomefile = '/home/nehleh/PhiloBacteria/Results/num_1/num_1_nu_0.05_Rlen_500_Rrate_0.01_Wholegenome.fasta'
    # baciSimLog = '/home/nehleh/PhiloBacteria/Results/num_1/num_1_nu_0.05_Rlen_500_Rrate_0.01_BaciSim_Log.txt'
    # cfml_log = '/home/nehleh/PhiloBacteria/Results/num_1/num_1_nu_0.05_Rlen_500_Rrate_0.01_CFML.importation_status.txt'
    # CFML_tree = '/home/nehleh/PhiloBacteria/Results/num_1/num_1_nu_0.05_Rlen_500_Rrate_0.01_CFML.labelled_tree.newick'

    parser = argparse.ArgumentParser(description='''You did not specify any parameters.''')
    parser.add_argument('-cl', "--clonaltreeFile", type=str, help='tree')
    parser.add_argument('-a', "--alignmentFile", type=str, help='fasta file')
    parser.add_argument('-cfl', "--cfmllogFile", type=str, help='cfmlFile')
    parser.add_argument('-cft', "--cfmltreefile", type=str, help='cfmltreefile')
    parser.add_argument('-rl', "--recomlogFile", type=str, help='BaciSim recombination log file')
    parser.add_argument('-sim', "--simulation", type=int, default= 1 , help='1 for the simulation data and 0 for emprical sequence')
    args = parser.parse_args()

    CFML_tree = args.cfmltreefile
    cfml_log = args.cfmllogFile
    genomefile = args.alignmentFile
    simulation = args.simulation

    # simulation = 0

    tns = dendropy.TaxonNamespace()
    cfml_tree = Tree.get_from_path(CFML_tree,'newick',taxon_namespace=tns)
    nodes_number = len(cfml_tree.nodes())
    alignment = dendropy.DnaCharacterMatrix.get(file=open(genomefile), schema="fasta")
    tips_num = len(alignment)
    alignment_len = alignment.sequence_size
    set_index(cfml_tree,alignment)


    my_label = []
    my_index = []
    for node in cfml_tree.postorder_node_iter():
        if "NODE" in str(node.label):
            my_label.append((node.label[5:]))
        else:
            my_label.append((node.taxon.label))
        my_index.append(node.index)


    CFMLData,rmse_CFML,CFML_recom_count,df = CFML_recombination(cfml_log,cfml_tree,tips_num)
    CFML_resultFig(cfml_tree, CFMLData)
    write_value(CFML_recom_count,'CFML_rcount.csv')
    (df['End']-df['Beg']).to_csv('CFML_delta.csv', index=False , header = ['length'])



    if simulation > 0:
        cfml_tree.reroot_at_midpoint(update_bipartitions=True)
        clonal_path = args.clonaltreeFile
        clonal_tree = Tree.get_from_path(clonal_path, 'newick', taxon_namespace=tns)
        set_index(clonal_tree,alignment)
        CFML_euclidean_distance = treecompare.euclidean_distance(clonal_tree, cfml_tree , edge_weight_attr="length")
        write_value(CFML_euclidean_distance, 'CFML_dist.csv')
        CFML_WRF_distance = treecompare.weighted_robinson_foulds_distance(clonal_tree, cfml_tree , edge_weight_attr="length")
        write_value(CFML_WRF_distance, 'CFML_WRF_distance.csv')
    if simulation == 1:
        baciSimLog = args.recomlogFile
        nodes_number_c = len(clonal_tree.nodes())
        realData, rmse_real = real_recombination(baciSimLog, clonal_tree, nodes_number_c, alignment_len, tips_num)
        rmse_real_CFML = mean_squared_error(rmse_real, rmse_CFML, squared=False)
        write_value(rmse_real_CFML, 'RMSE_CFML.csv')


