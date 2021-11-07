import numpy as np

from collections import defaultdict

from gurobipy import *

import scipy as sp
import csv
import argparse

from time import time
import os

import sys
from locfdr import locfdr

from src.common import *
from src.netmix2 import *

####################################
# Parse arguments
def get_parser():
    description = ''
    parser = argparse.ArgumentParser(description=description, fromfile_prefix_chars='@')
    parser.add_argument('-el', '--edge_list', required=False, type=str, help='edge list file')
    parser.add_argument('-gs', '--gene_score', required=False, type=str, help='gene-to-score file')
    parser.add_argument('-st', '--score_type', required=False, type=str, default='p', help='type of gene scores ("p" for pvalue or "z" for zscore')
    parser.add_argument('-d', '--delta', required=False, type=float, help='delta')
    parser.add_argument('-num_edges', '--num_edges', required=False, type=int, default=175000, help='Number of edges in G_delta for computing delta')
    parser.add_argument('-density', '--density', required=False, type=float, default=0.05, help='Target edge density (quadratic)')
    parser.add_argument('-time_limit', '--time_limit', required=False, type=int, default=12, help='Time limit (hours) for netmix function')
    parser.add_argument('-v', '--verbose', required=False, type=int, default=0, help='Print the progress of running NetMix2 (0, 1, or 2)')
    parser.add_argument('-o', '--output', required=False, type=str, help='directory for netmix results')
    return parser

###########################################################
def run(args):
    ###########################################################
    # load network
    (node_list, A_network, degs_network, n_network, num_edges_network) = load_network(args.edge_list, args.verbose)
    ###########################################################
    # load gene scores (as z-scores)
    gene_scores = load_genescores(args.gene_score, args.score_type, args.verbose)
    ###########################################################
    # restric to genes in the network
    (gene_scores_network, node_list_network_gwas, A_network_gwas, degs_network_gwas, n_network_gwas, nodes_in_gwas) = restrict_to_genes_in_network(p_val_list, node_list, A_network, args.verbose)
    ###########################################################
    # Compute PPR matrix
    (PPR_sim_mat_gwas, PPR_mat_rowsums_gwas, PPR_mat_network_gwas) = compute_ppr_kernel(A_network_gwas, degs_network_gwas, n_network_gwas, args.verbose)
    ###########################################################
    # fit local fdr
    nulltype_ind = 1
    nulltype_name = "mlest"
    r_locfdr=locfdr(gene_scores_network, nulltype=1, plot=0)
    # print('{} estimated size of altered subnetwork: {}'.format(nulltype_name, (1-r_locfdr['fp0']["p0"][nulltype_name])*len(gene_scores_network)))
    
    ###########################################################
    # correct the nans in the locfdr
    _, nonnull_count = correct_nans_from_locfdr(r_locfdr, gene_scores_network, nulltype_name)
    print("k: {}".format(nonnull_count))
    ###########################################################
    # compute G_delta
    # threshold similarity matrix
    sim_mat = PPR_sim_mat_gwas
    
    sim_mat_nodiag=sim_mat-np.diag(np.diag(sim_mat))
    sim_mat_nodiag_sorted=np.sort(sim_mat_nodiag[np.triu_indices(sim_mat_nodiag.shape[0],1)])[::-1]
    
    if args.delta:
        delta=args.delta # delta for string network
    else:
        delta=sim_mat_nodiag_sorted[args.num_edges]

    sim_mat_delta=1*(sim_mat_nodiag > delta)
    degs_PPR_gwas = np.sum(sim_mat_delta,0)
    num_edges_delta = int(np.sum(sim_mat_delta)/2)

    print("delta: {}".format(delta))
    print("kernel: {}".format(args.kernel))
    print('number of edges in G: {}'.format(np.sum(A_network_gwas)/2))
    print('number of edges in G_delta: {}'.format(np.sum(sim_mat_delta)/2))
    print('number of edges in G_delta/number of edges in G: {}'.format( (np.sum(sim_mat_delta)/2)/(np.sum(A_network_gwas)/2) ))
    
    ###########################################################
    write_list_to_file(os.path.join(args.output, 'node_list_network_gwas.tsv'), node_list_network_gwas)
    
    ###########################################################
    # other parameters for netmix
    s=nonnull_count
    alpha=s/len(gene_scores_network)
    target_edge_density = args.density
    rho=target_edge_density*(s-1)

    print('edge density: {}'.format(target_edge_density))
    print('rho: {}'.format(rho))
    print('s: {}'.format(s))
    print('alpha: {}'.format(alpha))
    
    ###########################################################
    # RUN NETMIX
    print('running netmix')
    print("time_limit", args.time_limit)

    est_subnetwork_zscore = netmix_edgedense(sim_mat_delta, rho, gene_scores_network, alpha, edge_dense_linear=True, time_limit=3600*args.time_limit)
    
    ###########################################################
    ######### print solution
    
    solution_size = len(est_subnetwork_zscore)
    est_subnetwork_genes = node_list_network_gwas[est_subnetwork_zscore] if solution_size>0 else []

    solution_network = sim_mat_delta[np.ix_(est_subnetwork_zscore,est_subnetwork_zscore)] if solution_size>0 else None
    num_edges_in_solution = sum(sum(solution_network))/2 if solution_size>0 else 0
    solution_network_density = num_edges_in_solution/(solution_size*(solution_size-1)/2) if solution_size>0 else 0

    print("Number of vertices in A: {}".format(solution_size))
    print("Density of A: {}".format(solution_network_density))
    ###########################################################
    # write solution
    if args.output:
        write_list_to_file(os.path.join(args.output, 'netmix_subnetwork.tsv'), est_subnetwork_genes)
    
if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))
