from __future__ import print_function
import cython
import math
import numpy as np
cimport numpy as np
import dendropy
import scipy.optimize as spo
from utility import *


cpdef p_matrix(np.ndarray[np.double_t, ndim=1] rates , np.ndarray[np.double_t, ndim=1] pi , double br_length):
    cdef:
      np.ndarray[np.double_t, ndim=2] p,freq,q,sqrtPi,sqrtPiInv,exchang,s
      np.ndarray[np.double_t, ndim=1] fun
      np.double_t mu
      int f = 1

    p = np.zeros((4, 4))
    mu = 0
    freq = np.zeros((4, 4))
    q = np.zeros((4, 4))
    sqrtPi = np.zeros((4, 4))
    sqrtPiInv = np.zeros((4, 4))
    exchang = np.zeros((4, 4))
    s = np.zeros((4, 4))
    fun = np.zeros(1)
    a, b, c, d, e = rates
    f = 1

    freq = np.diag(pi)
    sqrtPi = np.diag(np.sqrt(pi))
    sqrtPiInv = np.diag(1.0 / np.sqrt(pi))
    mu = 1 / (2 * ((a * pi[0] * pi[1]) + (b * pi[0] * pi[2]) + (c * pi[0] * pi[3]) + (d * pi[1] * pi[2]) + (
            e * pi[1] * pi[3]) + (pi[2] * pi[3])))
    exchang[0][1] = exchang[1][0] = a
    exchang[0][2] = exchang[2][0] = b
    exchang[0][3] = exchang[3][0] = c
    exchang[1][2] = exchang[2][1] = d
    exchang[1][3] = exchang[3][1] = e
    exchang[2][3] = exchang[3][2] = f

    q = np.multiply(np.dot(exchang, freq), mu)

    cdef int i
    cdef np.ndarray[np.double_t, ndim=1] eigval
    cdef np.ndarray[np.double_t, ndim=2] eigvec
    cdef np.ndarray[np.double_t, ndim=2] eigvec_inv
    cdef np.ndarray[np.double_t, ndim=2] left
    cdef np.ndarray[np.double_t, ndim=2] right


    for i in range(4):
        q[i][i] = -sum(q[i][0:4])

    s = np.dot(sqrtPi, np.dot(q, sqrtPiInv))
    eigval, eigvec = la.eig(s)
    eigvec_inv = la.inv(eigvec)
    left = np.dot(sqrtPiInv, eigvec)
    right = np.dot(eigvec_inv, sqrtPi)
    p = np.dot(left, np.dot(np.diag(np.exp(eigval * br_length)), right))
    return p


cpdef p_t(np.ndarray[np.double_t, ndim=1] rates,np.ndarray[np.double_t, ndim=1] pi,np.ndarray[np.double_t, ndim=1] br_length):
    cdef np.ndarray[np.double_t, ndim=2] blens = np.expand_dims(br_length, -1)
    rates = np.concatenate((rates, np.array([1.0])))
    cdef np.ndarray[np.double_t, ndim=2] exchang = np.zeros((4, 4))
    cdef np.ndarray[np.double_t, ndim=2] sqrtPi = np.diag(np.sqrt(pi))
    cdef np.ndarray[np.double_t, ndim=2] sqrtPiInv = np.diag(1.0 / np.sqrt(pi))
    iu = np.triu_indices(4, 1)
    exchang[np.triu_indices(4, 1)] = exchang[iu[1], iu[0]] = rates
    exchang = np.dot(exchang, np.diag(pi))
    exchang[range(4), range(4)] = -exchang.sum(-1)
    exchang = -exchang/np.sum(exchang.diagonal()*pi)
    s = np.dot(sqrtPi, np.dot(exchang, sqrtPiInv))


    cdef np.ndarray[np.double_t, ndim=1] eigval
    cdef np.ndarray[np.double_t, ndim=2] eigvec
    cdef np.ndarray[np.double_t, ndim=2] eigvec_inv
    cdef np.ndarray[np.double_t, ndim=2] left
    cdef np.ndarray[np.double_t, ndim=2] right
    eigval, eigvec = la.eig(s)
    eigvec_inv = la.inv(eigvec)

    left = np.dot(sqrtPiInv, eigvec)
    right = np.dot(eigvec_inv, sqrtPi)
    return left @ (np.expand_dims(np.exp(eigval * blens), 1)*np.eye(4)) @ right


cdef give_index(c):
    if c == "A":
        return 0
    elif c == "C":
        return 1
    elif c == "G":
        return 2
    elif c == "T":
        return 3
    elif c == "-":
        return -1
    elif c == "N":
        return -1


cdef set_index(tree, dna):
    cdef unsigned int idx
    cdef unsigned int sequence_count = len(dna)

    for node in tree.postorder_node_iter():
        node.index = -1
        node.annotations.add_bound_attribute("index")

    cdef unsigned int s = sequence_count
    for node in tree.postorder_node_iter():
        if not node.is_leaf():
            node.index = s
            node.label = str(node.index)
            s += 1
        else:
            for idx, name in enumerate(dna):
                 if str(name).replace("_"," ") == str(node.taxon):
                    node.index = idx
                    node.label = str(node.index)
                    break


cpdef  set_tips_partial(column,unsigned int tips_num,unsigned int alignment_len):
    cdef size_t site, tip
    cdef int i
    cdef np.ndarray[np.double_t, ndim=3] partial
    partial = np.zeros((alignment_len,tips_num, 4),dtype=np.double)
    for tip in range(tips_num):
      for site in range(alignment_len):
        i = give_index(column[site][tip])
        if i >= 0 :
            partial[site][tip][i] = 1
        else:
            partial[site][tip][0:4] = 1
    return partial


cpdef computelikelihood_mixture(tree,model, unsigned int alignment_len,unsigned int tips_num,np.ndarray[np.double_t, ndim=3] tip_partial):
    partial = np.zeros((alignment_len,(2 * tips_num) -1, 4))
    cdef np.ndarray[np.double_t, ndim=1] persite_ll = np.zeros(alignment_len)
    partial[:,0:tips_num,:] = tip_partial
    partial_new =  np.rollaxis(partial, 2, 0)
    for node in tree.postorder_node_iter():
        if not node.is_leaf():
            children = node.child_nodes()
            partial_new[:,:, node.index] = np.dot(model.p_matrix(children[0].edge_length), partial_new[:,:, children[0].index])
            for i in range(1, len(children)):
                partial_new[:,:, node.index] *= np.dot(model.p_matrix(children[i].edge_length), partial_new[:,:, children[i].index])

    persite_ll = np.log(model.get_pi() @ partial_new[:,:, tree.seed_node.index])
    return persite_ll,partial


cdef computelikelihood_mixture_bw(tree,unsigned int alignment_len,np.ndarray[np.double_t, ndim=3] tip_partial,model,unsigned int tips_num):
    cdef np.ndarray[np.double_t, ndim=3] partial = np.zeros((alignment_len,(2 * tips_num) - 1,4))
    partial[:,0:tips_num,:] = tip_partial
    cdef np.ndarray[np.double_t, ndim=1] persite_ll = np.zeros(alignment_len)
    partial_new =  np.rollaxis(partial,2,0)
    for node in tree.postorder_node_iter():
        if not node.is_leaf():
            children = node.child_nodes()
            partial_new[:,:, node.index] = np.dot(model.p_matrix(children[0].edge_length), partial_new[:,:, children[0].index])
            for i in range(1, len(children)):
                partial_new[:,:, node.index] *= np.dot(model.p_matrix(children[i].edge_length), partial_new[:,:, children[i].index])

    persite_ll = (model.get_pi() @ partial_new[:,:, tree.seed_node.index])
    return persite_ll, partial


cpdef compute_logprob_phylo_bw(np.ndarray[np.double_t, ndim=2] X,recom_trees,model,np.ndarray[np.double_t, ndim=3] tip_partial,unsigned int alignment_len,alignment,unsigned int tips_num):
    cdef size_t d0 = X.shape[0]
    cdef size_t d1 = X.shape[1]
    cdef size_t r = len(recom_trees)
    cdef size_t n = d0
    cdef size_t dim = d1
    result = np.zeros((n,r))
    cdef size_t tree_id
    for tree_id, item in enumerate(recom_trees):
        state_tree = dendropy.Tree.get(data=item, schema="newick")
        set_index(state_tree,alignment)
        persite_ll, partial = computelikelihood_mixture_bw(state_tree, alignment_len, tip_partial, model, tips_num)
        result[:, tree_id] = persite_ll
    return result


cdef forward(np.ndarray[np.double_t, ndim=2] X, np.ndarray[np.double_t, ndim=2] trans,np.ndarray[np.double_t, ndim=2] emission,np.ndarray[np.double_t, ndim=1] initial_distribution):
    cdef size_t d0 = X.shape[0]
    cdef size_t t0 = trans.shape[0]
    cdef size_t t,j
    alpha = np.zeros((d0, t0))
    cdef np.ndarray[np.double_t, ndim=1] c = np.zeros(d0)
    alpha[0, :] = initial_distribution * emission[0]
    c[0] = 1/alpha[0].sum(axis=0)
    alpha[0, :] *= c[0]
    for t in range(1, d0):
        for j in range(t0):
            alpha[t, j] = alpha[t - 1].dot(trans[:, j]) * emission[t, j]
        c[t] = 1/alpha[t].sum(axis=0)
        alpha[t] *= c[t]

    return alpha , c


cdef backward(np.ndarray[np.double_t, ndim=2] X, np.ndarray[np.double_t, ndim=2] trans, np.ndarray[np.double_t, ndim=2] emission,np.ndarray[np.double_t, ndim=1] c):
    cdef size_t t,j
    cdef size_t d0 = X.shape[0]
    cdef size_t t0 = trans.shape[0]
    beta = np.zeros((d0, t0))
    beta[d0 - 1] = np.ones(t0)

    # Loop in backward way from T-1 to
    for t in range(d0-2,-1,-1):
        for j in range(t0):
            beta[t, j] = (beta[t + 1] * emission[t + 1, :]).dot(trans[j, :])
        beta[t] *= c[t]

    return beta


cpdef baum_welch(np.ndarray[np.double_t, ndim=2] X, np.ndarray[np.double_t, ndim=2] trans,np.ndarray[np.double_t, ndim=2] emission,np.ndarray[np.double_t, ndim=1] initial_distribution, unsigned int n_iter=1):
    cdef size_t M = trans.shape[0]
    cdef size_t T = len(X)
    cdef size_t n,t,i
    cdef np.ndarray[np.double_t, ndim=2] gamma
    cdef np.ndarray[np.double_t, ndim=3] xi

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


cpdef update_mixture_partial_PSE(column,node,np.ndarray[np.double_t, ndim=3] tipdata,np.ndarray[np.double_t, ndim=1] posterior_recom,double nu,unsigned int alignment_len):
  cdef size_t site,my_number,i
  cdef double epsilon
  for site in range(alignment_len):
    dna = column[site]
    my_number = give_index(dna[node.index])
    epsilon = nu * posterior_recom[site]
    for i in range(4):
        if i == my_number:
            tipdata[site,node.index,i] = 1 - (epsilon)
        else:
            tipdata[site,node.index,i] = (epsilon)/3
  return tipdata


cpdef computelikelihood_normaltips(tree,GTR_sample,brs,nu,np.ndarray[np.double_t, ndim=2] posterior, int alignment_len,int tips_num ,int nodes_number,np.ndarray[np.double_t, ndim=3] tipdata):
    # brs: [B] B: number of branches
    # nu: [B] B: number of branches
    # posterior: [B,N] N: number of sites
    # P [B,4,4]
    cdef:
        np.ndarray[np.double_t, ndim=3] partial = np.empty((4,nodes_number, alignment_len))
        np.ndarray[np.double_t, ndim=3] P
        np.ndarray[np.double_t, ndim=3] P_nu
        np.ndarray[np.double_t, ndim=3] post
        np.ndarray[np.double_t, ndim=3] P_mixture
        int i

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
