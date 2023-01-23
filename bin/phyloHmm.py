#!/usr/bin/env python
import concurrent.futures
import numpy as np
from utility import *
import pandas as pd
import dendropy
import argparse
import csv
import json
import scipy.optimize as spo
import time
from scipy.optimize import Bounds
import operator
import itertools
import Calculation
from dendropy.calculate import treecompare



class GTR_model:
    def __init__(self, rates, pi):
        self.rates = rates
        self.pi = pi
    #     ========================================================================
    def get_pi(self):
        return self.pi
    #     ========================================================================
    def p_matrix(self , br_length):
        p = Calculation.p_matrix(self.rates,self.pi,br_length)
        return p
    def p_t(self , br_length):
        p = Calculation.p_t(self.rates,self.pi,br_length)
        return p
#***********************************************************************************************************************
def compute_logprob_phylo_bw(X,recom_trees,model,tip_partial,alignment_len):
    n, dim = X.shape
    result = np.zeros((n, len(recom_trees)))
    for tree_id, item in enumerate(recom_trees):
        state_tree = dendropy.Tree.get(data=item, schema="newick")
        set_index(state_tree,alignment)
        persite_ll, partial = computelikelihood_mixture_bw(state_tree, alignment_len, tip_partial, model, tips_num)
        result[:, tree_id] = persite_ll

    return result
# **********************************************************************************************************************
def set_tips_partial(column,tips_num):
    partial = np.zeros(((alignment_len, tips_num, 4)))
    for tip in range(tips_num):
      for site in range(alignment_len):
        dna = column[site]
        i = give_index(str(dna[tip]))
        if i >= 0 :
            partial[site][tip][i] = 1
        else:
            partial[site][tip][0:4] = 1
    return partial
# **********************************************************************************************************************
def computelikelihood_mixture(tree,alignment_len,tip_partial,model,tips_num):
    partial = np.zeros(((alignment_len,(2 * tips_num) -1, 4)))
    partial[:,0:tips_num,:] = tip_partial
    persite_ll = np.zeros(alignment_len)
    partial_new =  np.rollaxis(partial, 2, 0)
    for node in tree.postorder_node_iter():
        if not node.is_leaf():
            children = node.child_nodes()
            partial_new[..., node.index] = np.dot(model.p_matrix(children[0].edge_length), partial_new[..., children[0].index])
            for i in range(1, len(children)):
                partial_new[..., node.index] *= np.dot(model.p_matrix(children[i].edge_length), partial_new[..., children[i].index])


    persite_ll = np.log(model.get_pi() @ partial_new[..., tree.seed_node.index])

    return persite_ll, partial
# **********************************************************************************************************************
def computelikelihood_mixture_bw(tree,alignment_len,tip_partial,model,tips_num):
    partial = np.zeros(((alignment_len,(2 * tips_num) - 1,4)))
    partial[:,0:tips_num,:] = tip_partial
    persite_ll = np.zeros(alignment_len)
    partial_new =  np.rollaxis(partial,2,0)
    for node in tree.postorder_node_iter():
        if not node.is_leaf():
            children = node.child_nodes()
            partial_new[..., node.index] = np.dot(model.p_matrix(children[0].edge_length), partial_new[..., children[0].index])
            for i in range(1, len(children)):
                partial_new[..., node.index] *= np.dot(model.p_matrix(children[i].edge_length), partial_new[..., children[i].index])

    persite_ll = (model.get_pi() @ partial_new[..., tree.seed_node.index])

    return persite_ll, partial
# **********************************************************************************************************************
def tree_evolver_rerooted(tree ,node ,nu):
    co_recom = nu
    if (node.edge_length is None):
       node.edge.length = 0
    node.edge.length = node.edge.length + co_recom
    recombination_tree = tree.as_string(schema="newick")
    return recombination_tree
# **********************************************************************************************************************
def recom_maker(r_tree,index,nu):
    filter_fn = lambda n: hasattr(n, 'index') and n.index == index
    recombined_node = r_tree.find_node(filter_fn=filter_fn)
    return tree_evolver_rerooted(r_tree, recombined_node, nu)
# **********************************************************************************************************************
def update_mixture_partial_PSE(column,node,tipdata,posterior_recom,nu):
  for site in range(alignment_len):
    dna = column[site]
    my_number = give_index(dna[node.index])
    epsilon = nu * posterior_recom[site]
    # print(epsilon)
    for i in range(4):
        if i == my_number:
            tipdata[site,node.index,i] = 1 - (epsilon)
        else:
            tipdata[site,node.index,i] = (epsilon)/3
  return tipdata
# **********************************************************************************************************************
def my_best_nu(X,tree_path,clonal,target_node,tipdata,p_trans,p_start):
    def fn(nu):
        temptree = Tree.get_from_path(tree_path, 'newick')
        set_index(temptree, alignment)
        r_trees = recom_maker(temptree, target_node.index, nu)
        emission = compute_logprob_phylo_bw(X, [clonal, r_trees], GTR_sample, tipdata, alignment_len)
        alpha , c = forward(X, p_trans, emission, p_start)
        s = np.sum(np.log(c))
        return s

    result = spo.minimize_scalar(fn, method="bounded", bounds=(0.0, 0.5) , options={'disp': 1})
    return result.x
# *********************************************************************************************************************
def write_best_nu(best_nu,outputname):
    with open(outputname, mode='w') as bestnu_file:
        nu_writer = csv.writer(bestnu_file, delimiter=';', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        nu_writer.writerow(best_nu)

# **********************************************************************************************************************
def phylohmm(tree,alignment_len,column,p_start,p_trans,tips_num):
    mytree = []
    tipdata = Calculation.set_tips_partial(column,tips_num,alignment_len)
    posterior1 = [None] * (nodes_number-1)
    best_nu = [None] * (nodes_number-1)
    # make per_site likelihood and per_site partials for all nodes including tips and internal nodes
    persite_ll, partial = Calculation.computelikelihood_mixture(tree,GTR_sample, alignment_len, tips_num,tipdata)

    # each node play the role of target nodes
    for id_tree, target_node in enumerate(tree.postorder_node_iter()):
        if target_node != tree.seed_node:
            # print(target_node)
            recombination_trees = []
            mytree.append(Tree.get_from_path(tree_path, 'newick'))
            set_index(mytree[id_tree],alignment)
            #     step 1 --- make hmm input
            #  take the partials of the target node as input of hmm
            X = partial[:,target_node.index,:]


            #     step 2 --- make recombination tree
            my_nu = np.arange(0.0001, 0.1, 0.01)
            recombination_trees.append(mytree[id_tree].as_string(schema="newick"))
            # find the best nu for target branch based on the maximizing score value of hmm
            nu = my_best_nu(X,tree_path,recombination_trees[0],target_node,tipdata,p_trans,p_start)
            best_nu[target_node.index] = nu
            # print("best_nu:",nu)

            # make recombination tree for target node using best nu
            recombination_trees.append(recom_maker(mytree[id_tree],target_node.index,nu))
            emission = Calculation.compute_logprob_phylo_bw(X, recombination_trees, GTR_sample, tipdata, alignment_len,alignment,tips_num)
            # emission = compute_logprob_phylo_bw(X, recombination_trees, GTR_sample, tipdata, alignment_len)
            best_trans , p = Calculation.baum_welch(X, p_trans, emission, p_start, n_iter=1)
            p = p.T
            p_trans_nu0 = np.array([[1, 0],
                                    [1, 0]])
            if nu <= my_nu[0]: # if the best nu is smaller than a threshold, we consider there is no recombination on that branch
                trans = p_trans_nu0
            else:
                trans = best_trans

            R_over_theta = -(np.log(trans[0][0]) - target_node.edge_length)
            if trans[1][1] == 0:
                delta = 0
            else:
                delta = -1/np.log(trans[1][1])


            posterior1[target_node.index] = p[:, 1]
            # Update tip partials based on the posterior probability
            if target_node.is_leaf():
                Calculation.update_mixture_partial_PSE(column, target_node,tipdata,p[:, 1],nu,alignment_len)

    np.set_printoptions(threshold=np.inf)
    myposterior1 = np.array(posterior1, dtype='double')
    recom_prob = pd.DataFrame({'posterior_rec': pd.Series(list(myposterior1))})
    write_best_nu(best_nu,'PB_nu_two.txt')

    return tipdata,recom_prob,myposterior1,best_nu
# **********************************************************************************************************************
def recom_output(pb_tree,recom_prob,tips_num,threshold,nodes_number):
    output = np.zeros((alignment_len,nodes_number))
    rmsedata = np.zeros((alignment_len,tips_num))
    for i in range(len(recom_prob)):
        if (i < tips_num):
            for j in range(alignment_len):
                if (float(recom_prob['posterior_rec'][i][j]) >= threshold):
                    output[j, i] = 1
                    rmsedata[j, i] = 1
        else:
            for k in range(alignment_len):
                if (float(recom_prob['posterior_rec'][i][k]) >= threshold):
                    # output[k, recom_prob['recom_nodes'][i]] = 1
                    output[k, i] = 1
                    desc = set()
                    d = give_descendents(pb_tree, i, desc, tips_num)
                    # print(d)
                    for elm in d:
                        if simulation == 0:
                            node = pb_tree.find_node_with_taxon_label(elm)
                            elm = node.index
                        rmsedata[k, int(elm)] = 1
    return output,rmsedata
# **********************************************************************************************************************
def recom_resultFig_dm(pb_tree,recom_prob,tips_num,mixtureProb,outputname,nodes_number):
    clonaltree = Tree.get_from_path(pb_tree, 'newick')
    set_index(clonaltree, alignment)
    output,rmsedata = recom_output(clonaltree,recom_prob,tips_num,mixtureProb,nodes_number)
    fig = plt.figure(figsize=(tips_num + 9, tips_num / 2))
    color = ['red', 'green', 'purple', 'blue', 'black']
    for i in range(nodes_number):
        ax = fig.add_subplot(nodes_number, 1, i + 1)
        if i >= tips_num:
            desc = set()
            d = give_descendents(clonaltree, i, desc,tips_num)
            ax.plot(output[:, i], label=str(i) + ' is mrca:' + str(d), color=color[i % 5])
        else:
            ax.plot(output[:, i], label=give_taxon(clonaltree, i), color=color[i % 5])
        ax.legend(bbox_to_anchor=(0.045, 1.5), prop={'size': 10})
        ax.set_frame_on(False)
        ax.axis('off')

    ax.axis('on')
    ax.set_yticklabels([])
    plt.text(0.045,2.7, r'', fontsize=14,
             horizontalalignment='left', verticalalignment='top',
             bbox=dict(boxstyle="sawtooth", ec=(1., 0.5, 0.5), fc=(1., 0.8, 0.8),)
             )
    plt.savefig(outputname)
    return output,rmsedata
# **********************************************************************************************************************
def phyloHMM_Log(c_tree,output,outputname):
    nodes = []
    starts = []
    ends = []
    recomlens = []

    for i in range(output.shape[1]):
        non_zeros = [[i for i, value in it] for key, it in itertools.groupby(enumerate(output[:, i]), key=operator.itemgetter(1)) if key != 0]
        for j in range(len(non_zeros)):
            if i < tips_num:
                n = give_taxon(c_tree, i)
            else:
                n = str(i)
            nodes.append(n)
            starts.append(non_zeros[j][0])
            ends.append(non_zeros[j][len(non_zeros[j]) - 1])
            recomlens.append(non_zeros[j][len(non_zeros[j]) - 1] - non_zeros[j][0])

    all_data = {'nodes': nodes, 'start': starts, 'end': ends, 'len': recomlens}
    df = pd.DataFrame(all_data)
    df = df.sort_values(by=['nodes'], ascending=[True])
    df.to_csv(outputname, sep='\t', header=True , index = False)

    return df
# **********************************************************************************************************************
def forward(X, trans, emission, initial_distribution):
    alpha = np.zeros((X.shape[0], trans.shape[0]))
    c = np.zeros(X.shape[0])
    alpha[0, :] = initial_distribution * emission[0]
    c[0] = 1/alpha[0].sum(axis=0)
    alpha[0] *= c[0]
    for t in range(1, X.shape[0]):
        for j in range(trans.shape[0]):
            alpha[t, j] = alpha[t - 1].dot(trans[:, j]) * emission[t, j]
        c[t] = 1/alpha[t].sum(axis=0)
        alpha[t] *= c[t]

    return alpha , c
# **********************************************************************************************************************
def backward(X, trans, emission,c):
    beta = np.zeros((X.shape[0], trans.shape[0]))
    # setting beta(T) = 1
    beta[X.shape[0] - 1] = np.ones((trans.shape[0]))

    # Loop in backward way from T-1 to
    for t in range(X.shape[0] - 2, -1, -1):
        for j in range(trans.shape[0]):
            beta[t, j] = (beta[t + 1] * emission[t + 1, :]).dot(trans[j, :])
        beta[t] *= c[t]

    return beta
# **********************************************************************************************************************
def baum_welch(X, trans, emission, initial_distribution, n_iter=1):
    M = trans.shape[0]
    T = len(X)

    for n in range(n_iter):
        alpha,c  = forward(X, trans, emission, initial_distribution)
        beta = backward(X, trans, emission,c)

        gamma = np.zeros((M, T - 1))
        xi = np.zeros((M, M, T - 1))
        for t in range(T - 1):
            gamma[:, t] = (alpha[t, :] * beta[t, :]) / np.dot(alpha[t, :], beta[t, :])
            denominator = np.dot(np.dot(alpha[t, :].T, trans) * emission[t + 1].T, beta[t + 1, :])
            for i in range(M):
                numerator = alpha[t, i] * trans[i, :] * emission[t + 1].T * beta[t + 1, :].T
                xi[i, :, t] = numerator / denominator

        trans = np.sum(xi, 2) / np.sum(gamma, axis=1).reshape((-1, 1))

        # Add additional T'th element in gamma
        gamma = np.hstack((gamma , (alpha[T-1, :] * beta[T-1, :] / np.dot(alpha[T-1, :] , beta[T-1, :])).reshape((-1, 1)) ))



    return trans,gamma
# **********************************************************************************************************************
def computelikelihood(brs, nu, posterior):
    # brs: [B] B: number of branches
    # nu: [B] B: number of branches
    # posterior: [B,N] N: number of sites
    # P [B,4,4]

    partial = np.empty((4,nodes_number, alignment_len))
    partial[:,0:tips_num,:] = tipdata.T
    P = GTR_sample.p_t(brs)
    P_nu = GTR_sample.p_t(brs + nu)
    for node in tree.postorder_node_iter():
        if not node.is_leaf():
            children = node.child_nodes()
            if children[0].is_leaf():
                partial[:, node.index, :] = np.dot(P[children[0].index], partial[:,children[0].index, :])
            else:
                post = np.expand_dims(np.expand_dims(posterior[children[0].index], -1), -1)
                P_mixture = ((1 - post) * P[children[0].index] + post * P_nu[children[0].index])
                partial[:, node.index, :] = np.squeeze(P_mixture @ np.expand_dims(partial[:, children[0].index, :].T, -1), -1).T
            for i in range(1, len(children)):
                if children[i].is_leaf():
                    partial[:, node.index, :] *= np.dot(P[children[i].index], partial[:,children[i].index, :])
                else:
                    post = np.expand_dims(np.expand_dims(posterior[children[i].index], -1), -1)
                    P_mixture = ((1 - post) * P[children[i].index] + post * P_nu[children[i].index])
                    partial[:, node.index, :] *= np.squeeze(P_mixture @ np.expand_dims(partial[:, children[i].index, :].T, -1), -1).T

    persite_ll = np.log(GTR_sample.get_pi() @ partial[:, tree.seed_node.index, :])
    return persite_ll.sum()
# **********************************************************************************************************************
def computelikelihood_normaltips(brs, nu, posterior):
    # brs: [B] B: number of branches
    # nu: [B] B: number of branches
    # posterior: [B,N] N: number of sites
    # P [B,4,4]
    tipdata = set_tips_partial(column, tips_num)
    partial = np.empty((4,nodes_number, alignment_len))
    partial[:,0:tips_num,:] = tipdata.T
    P = GTR_sample.p_t(brs)
    P_nu = GTR_sample.p_t(brs + nu)
    for node in tree.postorder_node_iter():
        if not node.is_leaf():
            children = node.child_nodes()
            post = np.expand_dims(np.expand_dims(posterior[children[0].index], -1), -1)
            P_mixture = ((1 - post) * P[children[0].index] + post * P_nu[children[0].index])
            partial[:, node.index, :] = np.squeeze(P_mixture @ np.expand_dims(partial[:, children[0].index, :].T, -1), -1).T
            for i in range(1, len(children)):
                post = np.expand_dims(np.expand_dims(posterior[children[i].index], -1), -1)
                P_mixture = ((1 - post) * P[children[i].index] + post * P_nu[children[i].index])
                partial[:, node.index, :] *= np.squeeze(P_mixture @ np.expand_dims(partial[:, children[i].index, :].T, -1), -1).T

    persite_ll = np.log(GTR_sample.get_pi() @ partial[:, tree.seed_node.index, :])

    return persite_ll.sum()
# **********************************************************************************************************************


if __name__ == "__main__":

    # path = '/home/nehleh/Desktop/sisters/mutiple_sisters'
    # tree_path = path+'/num_1_RAxML_bestTree.tree'
    # genomefile = path+'/num_1_wholegenome_1.fasta'
    # baciSimLog = path+'/BaciSim_Log.txt'
    # clonal_path = path+'/clonaltree.tree'
    # json_path = '/home/nehleh/PhiloBacteria/bin/template/GTR_temp_partial.json'

    # path = '/home/nehleh/PhiloBacteria/Results/num_1'
    # tree_path = '/home/nehleh/PhiloBacteria/Emprical_Result/num_1/num_1_nu_1_Rlen_1_Rrate_1_RAxML_bestTree.tree'
    # clonal_path = path+'/num_1_Clonaltree.tree'
    # genomefile = '/home/nehleh/PhiloBacteria/hpc/0_RealData/SE_2018-20_outbreak_ten.fst'
    # baciSimLog = path+'/num_1_nu_0.05_Rlen_500_Rrate_0.01_BaciSim_Log.txt'
    # baciSimStat = path +'/num_1_nu_0.05_Rlen_500_Rrate_0.01_Recom_stat.csv'
    # json_path = '/home/nehleh/PhiloBacteria/bin/template/GTR_temp_partial.json'


    parser = argparse.ArgumentParser(description='''You did not specify any parameters.''')
    parser.add_argument('-t', "--raxmltree", type=str, help='tree')
    parser.add_argument('-a', "--alignmentFile", type=str, help='fasta file')
    parser.add_argument('-cl', "--clonaltreeFile", type=str, help='clonaltreeFile tree from BaciSim')
    parser.add_argument('-rl', "--recomlogFile", type=str, help='BaciSim recombination log file')
    parser.add_argument('-p', "--threshold", type=float, default=0.5, help='threshold')
    parser.add_argument('-f', "--frequencies", type=list, default= [0.2184,0.2606,0.3265,0.1946],help='frequencies')
    parser.add_argument('-r', "--rates", type=list, default= [0.975070 ,4.088451 ,0.991465 ,0.640018 ,3.840919 ], help='rates')
    parser.add_argument('-js', "--jsonFile", type=str, default='/home/nehleh/PhiloBacteria/bin/template/GTR_temp_partial.json', help='jsonFile')
    parser.add_argument('-sim', "--simulation", type=int, default=1, help='1 for the simulation data and 0 for emprical sequence')
    parser.add_argument('-rs', "--recomstat", type=str, help='recomstat')
    args = parser.parse_args()

    tree_path = args.raxmltree
    genomefile = args.alignmentFile
    pi = np.asarray(args.frequencies)
    rates = np.asarray(args.rates)
    threshold = args.threshold
    simulation = args.simulation
    # json_path = args.jsonFile
    # ============================================ setting parameters ==================================================
    tree = Tree.get_from_path(tree_path,schema='newick')
    alignment = dendropy.DnaCharacterMatrix.get(file=open(genomefile),schema="fasta")
    nodes_number = len(tree.nodes())
    tips_num = len(alignment)
    alignment_len = alignment.sequence_size
    GTR_sample = GTR_model(rates, pi)
    column = get_DNA_fromAlignment(alignment)
    set_index(tree, alignment)
    # print(tree.as_ascii_plot(show_internal_node_labels=True))


    p_start = np.array([0.9, 0.1])
    p_trans = np.array([[1-(1/alignment_len), 1/alignment_len],
                        [1/alignment_len, 1-(1/alignment_len)]])


    brs = [None]*nodes_number
    for _, node in enumerate(tree.postorder_node_iter()):
        brs[node.index] = node.edge_length
    brs.pop()


    tipdata, recom_prob, posterior, best_nu = phylohmm(tree, alignment_len,column,p_start,p_trans, tips_num)

    # start = time.time()
    tipdata = Calculation.set_tips_partial(column,tips_num,alignment_len)
    initial_guess = np.asarray(brs)
    bounds = [[1.e-10, tree.max_distance_from_root()]]*len(brs)
    result = spo.minimize(lambda x: -Calculation.computelikelihood_normaltips(tree,GTR_sample,x,best_nu,posterior,alignment_len,tips_num,nodes_number,tipdata), initial_guess, method='TNC', bounds=bounds)
    new_edges = np.asarray(result.x)
    for node in tree.postorder_node_iter():
        if node != tree.seed_node:
            node.edge_length = new_edges[node.index]
    tree.write(path="PhiloBacter.tree", schema="newick")
    # end = time.time()
    # print("cython_LL:",end - start)


    phyloHMMData2, rmse_PB2 = recom_resultFig_dm("PhiloBacter.tree", recom_prob, tips_num, threshold,'PB_Recom_two.jpeg', nodes_number)
    phyloHMM_log = phyloHMM_Log(tree, phyloHMMData2, 'PB_Log_two.txt')
    recom_count_pb2 = len(phyloHMM_log)
    write_value(len(phyloHMM_log), 'PB_rcount_two.csv')
    phyloHMM_log[['len']].to_csv('PB_delta_two.csv', index=False)
    mean_dalta = phyloHMM_log[['len']].mean()


    if simulation > 0:
        tns = dendropy.TaxonNamespace()
        clonal_path = args.clonaltreeFile
        clonal_tree = Tree.get_from_path(clonal_path, 'newick', taxon_namespace=tns)
        set_index(clonal_tree, alignment)
        PBtree = Tree.get(path='./PhiloBacter.tree',schema='newick',taxon_namespace=tns)
        PBtree.reroot_at_midpoint(update_bipartitions=True)
        PB_euclidean_distance = treecompare.euclidean_distance(clonal_tree,PBtree,edge_weight_attr="length")
        write_value(PB_euclidean_distance, 'PB_dist.csv')
        PB_WRF_distance = treecompare.weighted_robinson_foulds_distance(clonal_tree,PBtree,edge_weight_attr="length")
        write_value(PB_WRF_distance, 'PB_WRF_distance.csv')
    if simulation == 1 :
        baciSimLog = args.recomlogFile
        baciSimStat = args.recomstat
        nodes_number_c = len(clonal_tree.nodes())
        realData, rmse_real = real_recombination(baciSimLog, clonal_tree, nodes_number_c, alignment_len, tips_num)
        recom_stat = pd.read_csv(open(baciSimStat, "r"), sep=',')
        num_recom_real = len(recom_stat)
        write_value(len(recom_stat), 'baci_rcount.csv')
        recom_stat[['len']].to_csv('baci_delta.csv', index=False)
        rmse_real_philo2 = mean_squared_error(rmse_real, rmse_PB2, squared=False)
        write_value(rmse_real_philo2, 'RMSE_PB_two.csv')
