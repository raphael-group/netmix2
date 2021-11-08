import csv
import pandas as pd
import numpy as np
from scipy.stats import norm
import networkx as nx

def read_genes_from_file(genes_file):
    L = []
    with open(genes_file) as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        for row in reader:
            if row[0] != "gene":
                L.append(row[0])
    return L


def write_list_to_file(filename, l):
    with open(filename, 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        for i in range(len(l)):
            tsv_writer.writerow([ str(l[i]) ])

def load_network(edge_list_file, verbosity=0):
    if verbosity > 0:
        print('loading network')
    df_el = pd.read_csv(edge_list_file, header = None, sep = "\t")
    node_list = sorted(set(df_el[0]) | set(df_el[1]))
    node_set = set(node_list)
    num_nodes = len(node_list)

    # load edges
    A_network=np.zeros((num_nodes, num_nodes))

    with open(edge_list_file, 'r') as f:
        arrs = [l.rstrip().split("\t") for l in f if not l.startswith("#")]
        for row in arrs:
            i = node_list.index(row[0])
            j = node_list.index(row[1])
            A_network[i,j] = 1
            A_network[j,i] = 1
    num_edges = int(np.sum(A_network)/2)
    
    if verbosity > 0:
        print("Number of nodes: {}, Number of edges: {}".format(num_nodes, num_edges))
    
    return (node_list, A_network)

            
def load_pvalues(pvalues_file, node_list, verbosity=0):
    if verbosity > 0:
        print('loading genescores')
    node_set = set(node_list)
    pval_list = np.zeros(len(node_set))
    
    with open(pvalues_file) as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        for row in reader:
            gene=row[0]
            if gene in node_set:
                ind = node_list.index(gene)
                pval=float(row[1])                
                pval_list[ind]=pval
    return pval_list
    
def restrict_to_genes_in_network(pvals_list, node_list, A_network, verbosity=0):
    # find nodes in the network that has scores
    nodes_with_pval_in_network_indices = np.array([ind for ind,val in enumerate(pvals_list) if val > 0])
    
    node_list = np.array(node_list)[nodes_with_pval_in_network_indices]
    A_network = A_network[np.ix_(nodes_with_pval_in_network_indices, nodes_with_pval_in_network_indices)]
    pvals_list = pvals_list[nodes_with_pval_in_network_indices]
    
    # find largest connected component
    G=nx.Graph(A_network)
    lcc=list(max(list(nx.connected_components(G)),key=lambda x: len(x)))
    
    node_list = np.array(node_list)[lcc]
    A_network=A_network[np.ix_(lcc,lcc)]
    pvals_list = pvals_list[lcc]
    
    n = len(pvals_list)
    
    if verbosity > 0:
        print('number of nodes in G: {}'.format(n))
        print('number of edges in G: {}'.format(int(np.sum(A_network)/2)))
    
    return (pvals_list, node_list, A_network)


def compute_zscores(pvalues):
    zscores=-1*norm.ppf(pvalues)

    return zscores

def post_process_zscores(zscores):
    # correct for -np.inf from genes with P=1
    new_zscores = []
    for i in range(len(zscores)):
        zscore = zscores[i]
        if zscore != -np.inf:
            new_zscores.append(zscore)
    
    min_zscore = min(new_zscores)
    print(min_zscore)
    
    # swap -inf with next smallest z-scores
    for i in range(len(zscores)):
        if zscores[i] == -np.inf:
            zscores[i] = min_zscore
    return zscores


def correct_nans_from_locfdr(r_locfdr, scores, nulltype_name="mlest", verbosity=0):
    # correct the nans in the locfdr
    resps=1-r_locfdr['fdr']
    mu = r_locfdr['fp0']["delta"][nulltype_name]
    nan_count = 0
    nan_count2 = 0
    original_nonnull_count = 0
    nonnull_count = 0
    
    for ind,t in enumerate(resps):
        if np.isnan(t):
            if scores[ind] > 4:
                resps[ind]=1
            else:
                resps[ind]=0
                nan_count2 += 1
        elif t > 0 and scores[ind] < mu: # mu vs 0
            resps[ind]=0
            if t > 0.5:
                nan_count += 1
        
        if resps[ind] > 0.5:
            nonnull_count += 1
        if not np.isnan(t) and t > 0.5:
            original_nonnull_count += 1

    nonnull_count_locfdr = int((1-r_locfdr['fp0']["p0"]["mlest"])*len(scores))
    if verbosity > 0:
        print("Estimate size of altered subnetwork: {}".format(nonnull_count))
    return nonnull_count




def compute_ppr_kernel(A_network, verbosity = 0):
    if verbosity > 0:
        print('computing the PPR kernel')
    r = 0.4
    degs_network=np.sum(A_network,0)
    num_nodes = len(A_network)
    
    D_network = np.diag(1.0/degs_network)
    P_network = A_network.dot(D_network)
    
    # PPR kernel
    PPR_mat = r * np.linalg.inv(np.eye(num_nodes) - (1-r)*P_network)

    # Compute rowsums of PPR mat
    PPR_mat_rowsums = np.sum(PPR_mat, axis=1)
    
    # Similarity matrix
    PPR_sim_mat=np.minimum(PPR_mat,PPR_mat.T)
    #print(PPR_sim_mat.shape)
    
    # return similarity matrix, rowsums
    return PPR_sim_mat, PPR_mat_rowsums
