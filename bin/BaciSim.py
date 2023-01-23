#!/usr/bin/env python


from utility import *
import argparse
from dendropy import Tree
import numpy as np
import dendropy
import random
import matplotlib.pyplot as plt
import pandas as pd
import csv
import os.path




def set_index(tree):
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
          node.index = int(node.taxon.label)
          node.label = str(node.index)
# **********************************************************************************************************************
def give_recom_num(tree,recom_rate,alignment_len):
    # Poisson( tree.sum() * rel_recomb_rate_per_site * alignment_length)
    recom_num =  np.random.poisson(tree.length() * recom_rate * alignment_len)
    return recom_num
# **********************************************************************************************************************
def make_nodes_weight(tree,status):
    set_index(tree)
    node_labels = []
    node_weight = []
    for node in tree.postorder_node_iter():
      if not node==tree.seed_node:
        if status == 1:
          node_labels.append(node.index)
          node_weight.append(node.edge_length /tree.length())
        elif status == 0:
          if node.is_leaf():
            node_labels.append(node.index)
            node_weight.append(node.edge_length /tree.length())
        elif status == 2:
          if node.is_internal():
            node_labels.append(node.index)
            node_weight.append(node.edge_length /tree.length())

    return node_labels,node_weight
# **********************************************************************************************************************
def write_sim_nu(nu,node):
    headers = ['nu','node']
    with open('./Sim_nu.csv', mode='a') as sim_nu:
        nu_writer = csv.writer(sim_nu, delimiter=';', quotechar='"', quoting=csv.QUOTE_MINIMAL )
        file_is_empty = os.stat('./Sim_nu.csv').st_size == 0
        if file_is_empty:
            nu_writer.writerow(headers) # file doesn't exist yet, write a header
        nu_writer.writerow([nu,node])
# **********************************************************************************************************************
def recom_on_alignment(recom_num,recom_len,alignment_len,clonal_tree,node_labels,node_weight,nu_ex,taxa,threshold_len):
    if recom_num != 0:
        ideal_gap = int(alignment_len / recom_num)
    else:
        ideal_gap = 0
    my_trees = []
    nodes = []
    starts = []
    ends = []
    recomlens = []
    rand_nu = []
    edge_len = []
    for recom_id in range(int(recom_num)):
      starting_falg = False
      tree = Tree.get_from_string(clonal_tree,schema='newick')
      set_index(tree)
      random_tip = np.random.choice(node_labels,1,node_weight)
      start_pos = random.randint(ideal_gap * recom_id, ideal_gap * (recom_id+1) )
      r_len = random.randint(recom_len, recom_len + threshold_len)
      if (start_pos + r_len > alignment_len):
          r_len = alignment_len - start_pos
      if r_len > 0:
          nodes.append(int(random_tip))
          starts.append(start_pos)
          recomlens.append(r_len)
          ends.append(start_pos + r_len)
          recom_node = tree.find_node_with_label(str(int(random_tip)))
          edge_len.append(recom_node.edge_length)
          # recom_tree= ex_recom_maker(tree,recom_node,nu_ex,taxa) # make external recombination
          recom_tree, my_nu = ex_recom_maker(tree, recom_node, nu_ex,taxa)
          my_trees.append(recom_tree)
          rand_nu.append(my_nu)



    all_data = {'nodes':nodes , 'start':starts , 'end':ends, 'len':recomlens , 'tree':my_trees , 'nu' : rand_nu , 'edge_len':edge_len}
    df = pd.DataFrame(all_data)


    gf = df.value_counts(['nodes','edge_len','nu']).reset_index(name='count')
    gf.to_csv('./Recombination_count.csv', sep='\t', header=True)
    df[['nodes','start','end' ,'len','nu']].to_csv('./Recom_stat.csv', sep=',', header=True)

    return df,all_data
# **********************************************************************************************************************
def ex_recom_maker(tree ,node ,nu ,taxa):
    # rand_nu = np.random.normal(nu,0.007)
    rand_nu = nu
    write_sim_nu(rand_nu, node)
    co_recom = rand_nu/2
    new_tree = dendropy.Tree(taxon_namespace=taxa)
    parent = node.parent_node

    # changeing in topology when recombiantion is larger than tree.max_distance_from_root()
    if (co_recom + node.edge_length + node.distance_from_tip()) >= tree.max_distance_from_root():
      if (node.edge_length is None):
        node.edge.length = 0
      if node.is_leaf():
        external_len = (co_recom + node.edge_length) - tree.max_distance_from_root()
        tree.prune_subtree(node)
        other_nodes = tree.seed_node
        new_tree.seed_node.add_child(other_nodes)
        new_tree.seed_node.add_child(node)
        other_nodes.edge_length = co_recom + node.edge_length - other_nodes.distance_from_tip()
        node.edge_length = co_recom + node.edge_length
      else:
        tree.prune_subtree(node)
        other_nodes = tree.seed_node
        new_tree.seed_node.add_child(other_nodes)
        new_tree.seed_node.add_child(node)
        node.edge_length = co_recom + node.edge_length
        other_nodes.edge_length = node.distance_from_tip() + node.edge_length - other_nodes.distance_from_tip()
    #*************************************************************************
    # topology does not change in this case:
    elif  node.is_leaf() and ((co_recom + node.edge_length) < parent.distance_from_tip()) and (parent != tree.seed_node):
        parent.edge.length = parent.distance_from_tip() - (co_recom + node.edge_length)
        node.edge.length = node.edge.length + co_recom
        sister = node.sister_nodes()
        sister[0].edge.length = sister[0].edge.length + co_recom
        new_tree = tree
    elif ((co_recom + node.edge_length) < parent.distance_from_tip()) and (parent == tree.seed_node):
        node.edge.length = node.edge.length + co_recom
        sister = node.sister_nodes()
        sister[0].edge.length = sister[0].edge.length + co_recom
        new_tree = tree
    #*************************************************************************
    # changing in topology to make recombination tree:
    elif  ((co_recom + node.edge_length) > parent.distance_from_tip())  and ((co_recom + node.edge_length + node.distance_from_tip()) < tree.max_distance_from_root()):
        ancestor = []
        recom_length = co_recom + node.edge_length
        totip = node.distance_from_tip()
        for id,tmp_node in enumerate(node.ancestor_iter()):
            ancestor.append(tmp_node)
            if recom_length + totip < tmp_node.distance_from_tip() :
                attached_node = tmp_node
                attached_id = id
                break

        sister = node.sister_nodes()
        grandparent = parent.parent_node
        relocated_nodes = ancestor[attached_id-1]  # relocated node is the adjacent node of recombinant node
        if parent == relocated_nodes :
            node.edge_length = co_recom + node.edge_length
            sister[0].edge_length = sister[0].edge_length + co_recom
            parent.edge_length = parent.edge_length - co_recom
            new_tree = tree
        else:
            parent.remove_child(node)         # the original recombinant node was removed to reinsert in the other side
            attached_node.remove_child(relocated_nodes) # relocated node was removed to reinsert in to right side
            newborn = dendropy.datamodel.treemodel.Node()  # newborn in the new mrca of recombinant node and its sister
            newborn.edge_length = attached_node.distance_from_tip() - (recom_length + totip)
            node.edge_length = recom_length
            newborn.add_child(node)
            relocated_nodes.edge_length = relocated_nodes.edge_length - newborn.edge_length
            if (not grandparent is None):
                grandparent.remove_child(parent)
                grandparent.add_child(sister[0])
                sister[0].edge_length = sister[0].edge_length + parent.edge_length
            newborn.add_child(relocated_nodes)
            # sister[0].edge_length = sister[0].edge_length + parent.edge_length
            attached_node.add_child(newborn)
            new_tree.seed_node = tree.seed_node
    #*************************************************************************
    return new_tree.as_string(schema="newick"),rand_nu
#***********************************************************************************************************************
def recom_repeat(recom_num, recom_len, alignment_len, clonal_tree, nu_ex, taxa, r_nodes):
    if recom_num != 0:
        ideal_gap = int(alignment_len / recom_num)
    else:
        ideal_gap = 0
    my_trees = []
    nodes = []
    starts = []
    ends = []
    recomlens = []
    rand_nu = []
    edge_len = []
    for random_tip in range(r_nodes):
        for recom_id in range(int(recom_num)):
            starting_falg = False
            tree = Tree.get_from_string(clonal_tree, schema='newick')
            set_index(tree)
            while not starting_falg:
                # random_tip = np.random.choice(node_labels,1,node_weight)
                start_pos = random.randint(ideal_gap * recom_id, ideal_gap * (recom_id + 1))
                r_len = recom_len
                if (start_pos + r_len <= alignment_len):
                    nodes.append(int(random_tip))
                    starts.append(start_pos)
                    recomlens.append(r_len)
                    ends.append(start_pos + r_len)
                    recom_node = tree.find_node_with_label(str(int(random_tip)))
                    edge_len.append(recom_node.edge_length)
                    # recom_tree= ex_recom_maker(tree,recom_node,nu_ex,taxa) # make external recombination
                    recom_tree, my_nu = ex_recom_maker(tree, recom_node, nu_ex, taxa)
                    my_trees.append(recom_tree)
                    rand_nu.append(my_nu)
                    starting_falg = True

    all_data = {'nodes': nodes, 'start': starts, 'end': ends, 'len': recomlens, 'tree': my_trees, 'nu': rand_nu, 'edge_len': edge_len}
    df = pd.DataFrame(all_data)



    gf = df.value_counts(['nodes', 'edge_len', 'nu']).reset_index(name='count')
    gf.to_csv('./Recombination_count.csv', sep='\t', header=True)
    df.to_csv('./Recombination_Log.txt', sep='\t', header=True)

    return df, all_data
# **********************************************************************************************************************
def make_recom_fig(all_data,alignment_len,nodes_number,tips_num,clonal_tree,recom_num):
    output = np.zeros((alignment_len, nodes_number))
    for i in range(nodes_number):
      for j in range(recom_num):
        if int(all_data['nodes'][j]) == i:
          s = int(all_data['start'][j])
          e = int(all_data['end'][j])
          output[s:e,i] = 1

    fig = plt.figure(figsize=(tips_num+9,tips_num/2))
    color=['red', 'green' ,'purple', 'blue','black']
    clonaltree = Tree.get_from_string(clonal_tree,schema='newick')
    set_index(clonaltree)
    for i in range(nodes_number):
      ax = fig.add_subplot(nodes_number,1,i+1)
      if i >= tips_num:
        desc = set()
        d = give_descendents(clonaltree,i,desc,tips_num)
        ax.plot(output[:,i] ,label= str(i)+' is mrca:'+ str(d) ,color = color[i%5])
      else:
        ax.plot(output[:,i] ,label= i ,color = color[i%5])
      ax.legend( bbox_to_anchor=(0.045, 1.5) ,prop={'size':10} )
      ax.set_frame_on(False)
      ax.axis('off')

    ax.axis('on')
    ax.set_yticklabels([])
    plt.savefig("./BaciSim_Recombination.jpeg")
    # plt.show()
    return  output
# **********************************************************************************************************************
def simple_merge(recomtrees , recomnodes):
    clonaltree = Tree.get_from_string(clonal_tree, schema='newick')
    set_index(clonaltree)

    equ = np.zeros((len(recomtrees), 4))
    for treeid in range(len(recomtrees)):
        rtree = Tree.get_from_string(recomtrees[treeid], schema='newick')
        set_index(rtree)
        equ[treeid, 0] = recomnodes[treeid]
        equ[treeid, 1:3] = give_equivalent_node(rtree)
        equ[treeid, 3] = treeid

    for i in range(len(equ)):
      main_node = clonaltree.find_node_with_label(str(int(equ[i][0])))
      main_node.edge_length = equ[i][2]

    return clonaltree.as_string(schema="newick")
# **********************************************************************************************************************
def generate_final_report(df,alignment_len,clonal_tree,tips_num):
    endpoints = df[['start', 'end']].stack().sort_values().reset_index(drop=True)
    intervals = pd.DataFrame({'start': endpoints.shift().fillna(0), 'end': endpoints}).astype(int)
    # construct the list of intervals from the endpoints
    intervals['intv'] = [pd.Interval(a, b) for a, b in zip(intervals.start, intervals.end)]
    # these are the original intervals
    orig_invt = pd.arrays.IntervalArray([pd.Interval(a, b) for a, b in zip(df.start, df.end)])
    # walk through the intervals and compute the intersections
    intervals['total'] = intervals.intv.apply(lambda x: orig_invt.overlaps(x).sum())
    bounds = np.unique(df[['start', 'end']])
    if 0 not in bounds: bounds = np.insert(bounds, 0, 0)
    end = alignment_len
    bounds = np.append(bounds, end)
    total = []
    interval = []
    stat = []
    final_len = []
    nodes = []
    r_trees = []
    final_tree = []
    clonaltree = Tree.get_from_string(clonal_tree, schema='newick')
    set_index(clonaltree)
    children = []
    for i in range(len(bounds) - 1):
        # Find which intervals fit
        ix = (df['start'] <= bounds[i]) & (df['end'] >= bounds[i + 1])
        final_len.append(bounds[i + 1] - bounds[i])
        total.append(np.sum(ix))

        interval.append(df[ix].values.tolist())

        temp_node = []
        temp_tree = []
        kids = []
        temp = df[ix].values.tolist()
        for j in range(len(temp)):
            temp_node.append(int(temp[j][0]))
            if int(temp[j][0]) >= tips_num:
                desc = set()
                d = give_descendents(clonaltree, int(temp[j][0]), desc,tips_num)
                kids.append(d)
            else:
                kids.append("")
            temp_tree.append(temp[j][4])
        nodes.append(temp_node)
        r_trees.append(temp_tree)
        children.append(kids)

        if (np.sum(ix) == 0):
            stat.append('clonal')
            final_tree.append(clonal_tree)
        elif (np.sum(ix) == 1):
            stat.append("recom")
            final_tree.append(r_trees[i])
        else:
            stat.append("overlap")
            final_tree.append(simple_merge(r_trees[i], nodes[i]))

    final = pd.DataFrame({'start': bounds[:-1], 'end': bounds[1:] ,'nodes':nodes, 'descendants':children,  'len':final_len , 'status':stat ,'final_tree': final_tree ,'total': total ,'tree':r_trees })
    final[['nodes', 'start', 'end', 'len', 'descendants', 'status']].to_csv('./BaciSim_Log.txt', sep='\t', header=True)

    # delete row with len zero, it made problem in seq-gen
    delete_row = final[final['len'] == 0].index
    final = final.drop(delete_row)

    myfile = open('./BaciSimTrees.tree', 'w')
    total = 0
    for id in range(len(final)):
        if len(final['final_tree'][id]) == 1:
            end_tree = remove_internal_labels(final['final_tree'][id][0],tips_num)
        else:
            end_tree = remove_internal_labels(final['final_tree'][id],tips_num)
        tmp = "[" + str(final['len'][id]) + "]" + " " + end_tree
        myfile.write("%s" % tmp)

    myfile.close()

    return final
# **********************************************************************************************************************
# simulate an alignment of length seq_len from a tree
def simulate(tree, seq_len, alphabet='ACTG'):
    rates = [0.975070 ,4.088451 ,0.991465 ,0.640018 ,3.840919 ]
    pi = [0.2184,0.2606,0.3265,0.1945]
    GTR_sample = GTR_model(rates, pi)
    edge_lengths = np.empty(len(tree.taxon_namespace) * 2 - 2)
    for node in tree.postorder_node_iter():
        if node != tree.seed_node:
            edge_lengths[node.index] = node.edge_length
    edge_lengths /= np.sum(edge_lengths)

    mats = GTR_sample.p_t(np.expand_dims(edge_lengths, axis=1))

    alignment = [None] * (len(tree.taxon_namespace) * 2 - 1)
    alignment[tree.seed_node.index] = np.random.choice(np.arange(4), seq_len, p=pi)
    for node in tree.preorder_node_iter():
        if node != tree.seed_node:
            node.edge_length = edge_lengths[node.index]
            sequence = []
            for i in range(seq_len):
                probs = mats[node.index][alignment[node.parent_node.index][i],]
                sequence.append(np.random.choice(np.arange(len(alphabet)), 1, p=probs))
            alignment[node.index] = np.concatenate(sequence)
    alignment = {taxon.label: ''.join(list(map(lambda x: alphabet[x], seq))) for taxon, seq in  zip(tree.taxon_namespace, alignment[:len(tree.taxon_namespace)])}
    return alignment
# **********************************************************************************************************************

if __name__ == "__main__":

    path = os.getcwd()

    parser=argparse.ArgumentParser(
        description='''You did not specify any parameters. You must at least specify the number of chromosomes sampled and the sequence length. ''',
        epilog="""All's well that ends well.""")
    parser.add_argument('-cl', "--clonaltree", type=str , default=path+'/Clonaltree.tree' , help='clonaltree')
    parser.add_argument('-n', "--tips_number", type=int, default=10 , help='Sets the number of isolates (default is 10)')
    parser.add_argument('-g', "--alignment_len", type=int, default=100000 , help='Sets the number and lengths of fragments of genetic material (default is 5000)')
    parser.add_argument('-l', "--recom_len", type=int, default=500, help='Sets the average length of an external recombinant interval, (default is 500)')
    parser.add_argument('-r', "--recom_rate",type=float, default=0.005, help='Sets the site-specific rate of external (between species) recombination, (default is 0.05)')
    parser.add_argument('-nu',"--nu" ,  type=float, default=0.01, help='nu')
    parser.add_argument('-s',"--status" ,  type=int, default=1, help='0 is just leaves, 1 is for both internal nodes and leaves and 2 is just internal nodes')
    parser.add_argument('-f', "--fixed", type=int, default=1, help='0 for fixed number and fixed len of recombination and 1 for normal/random way making recombination events.')
    parser.add_argument('-e', "--each_recom", type=int, default=1, help='each_recom')



    # Read arguments from command line
    args = parser.parse_args()

    clonaltree = args.clonaltree
    tips_num = args.tips_number
    alignment_len = args.alignment_len
    recom_len = args.recom_len
    recom_rate = args.recom_rate
    nu_ex = args.nu
    status = args.status
    fixed_status = args.fixed
    threshold_len = 50
    each_recom = args.each_recom


    taxa = make_taxa(tips_num)
    tree = Tree.get_from_path(clonaltree, 'newick')
    set_index(tree)
    clonal_tree = tree.as_string(schema="newick")
    # print(tree.as_ascii_plot(show_internal_node_labels=True))
    nodes_num = len(tree.nodes())

    if status == 0:
        r_nodes = tips_num
    else:
        r_nodes = nodes_num - 1

    if fixed_status == 0:
        df,all_data = recom_repeat(each_recom, recom_len, alignment_len, clonal_tree, nu_ex, taxa,r_nodes)
        recom_num = len(df)
    else:
        recom_num = give_recom_num(tree, recom_rate, alignment_len)
        node_labels, node_weight = make_nodes_weight(tree, status)
        df, all_data = recom_on_alignment(recom_num,recom_len,alignment_len,clonal_tree,node_labels,node_weight,nu_ex,taxa,threshold_len)

    output = make_recom_fig(all_data,alignment_len, nodes_num, tips_num, clonal_tree,recom_num)
    final_report = generate_final_report(df, alignment_len, clonal_tree, tips_num)

    tree.deroot()
    unroot_clonaltree = tree.as_string(schema="newick")
    unroot_clonaltree = unroot_clonaltree.replace('\n',"")
    myfile = open('./unroot_Clonaltree.tree', 'w')
    myfile.write(unroot_clonaltree)
    myfile.close()