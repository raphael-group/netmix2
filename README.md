# NetMix2

NetMix2 is an algorithm for identifying altered subnetworks from a wide range of subnetwork families, 
including the **propagation family** which approximates the subnetworks ranked highly by network propagation. 

This repository includes instructions for installation and tutorials using example data for NetMix2.  
*This README is work in progress.*

## Algorithm
--------------

<img src="doc/netmix2_overview.png" width="600">

The main goal of NetMix2 is to identify an altered subnetwork $A$ in a subnetwork family $\mathcal{S}$ from a graph $G=(V,E)$ with gene scores $X_v$ for all $v \in V$.  
NetMix2 consists of two main steps:
1. Estimate the size of altered subnetwork.
2. Identify the altered subnetwork from a subnetwork family.

For step 1, NetMix2 uses local false discovery (local FDR) method to estimate the altered distribution and the background distribution from a set of gene scores. Local FDR uses a semi-parametric model that has weaker assumptions than the parametric Gaussian mixture model and thus can flexibly model a wide range of gene score distributions.  
For step 2, NetMix2 identifies the altered subnetwork with size estimated from step 1 and largest total vertex scores from an input subnetwork family. By default, NetMix2 uses the propagation family which approximates the altered subnetworks found by network propagation, thereby unifying the principles of network propagation with altered subnetwork identification.

## Setup
--------------

Setting up NetMix2 requires the following steps:

### Download

Download NetMix2 using the following command. This command clones the NetMix2 repository from Github.

`git clone https://github.com/raphael-group/netmix2.git`

### Requirements

NetMix2 is written in Python 3 and requires several dependencies listed below.  
We recommend `virtualenv` or `conda` for managing the required dependencies.

- Python (3.6)
- NumPy
- SciPy
- Matpotlib
- pandas
- statsmodels
- locfdr-python (v0.1a) - available at https://github.com/leekgroup/locfdr-python
- Gurobi

### Testing NetMix2
NetMix2 using the *propagation family* can be executed on example data using the following command.
```
python run_netmix2.py -el data/edge_list.tsv -gs data/gene_scores.tsv -o results/example_output.tsv
```
Detailed instructions for running NetMix2 including the input file format and command-line options are described below. 


## Usage
----------

NetMix2 uses the *propagation family* by default.  
For this subnetwork family, NetMix2 constructs a similarity threshold graph $G_\delta$ using the similarity matrix $M \in \mathbb{R}^{|V|\times|V|}$ where $M_{v, w}$ is the Personalized PageRank from vertex $v$ to vertex $w$.  
The propagation family $\mathcal{M}_{\delta,p}$ is then equal to the edge-dense family $\mathcal{E}_{G_\delta,p}$ for the similarity threshold graph $G_\delta = (V, E_\delta)$ which has edge $(v, w)$ if $M_{v,w} \geq \delta$ and $M_{w,v} \geq \delta$.

Instructions for using other subnetwork families are described in Additional Infromation.

### Input

NetMix2 requires two tab-separated text file - an edge list for interaction network $G$ and a gene scores file $X$. 

The following example demonstrates a network with three vertices `A`, `B`, and `C` that have gene scores ($p$-values) of `0.1`, `0.5`, and `0.9`, respectively.


#### Edge list file
Each row in this file corresponds to an edge in the network.

    A    C
    B    C

#### Gene-to-score file
Each line in this file associates a node with a score:

    A    0.1
    B    0.5
    C    0.9

#### Parameters for the propagation family
In addition to the files above, running NetMix2 using the propagation family requires two family-specific parameters:
- Similarity threshold $\delta$. Alternatively, users can choose the number of edges $|E_\delta|$ for the similarity threshold graph $G_\delta$.
- The minimum edge density $p$ of the altered subnetwork in $G_\delta$.

Below are the command line options for NetMix2.

### NetMix2 (*propagation family*) command line options
| Flag | Name | Description | Default Value |
| --- | --- | --- | --- |
| -el | edge_list | Edge list file | None |
| -gs | gene_scores | Gene-to-score file | None |
| -d | delta | The similarity threshold $\delta$ for $G_\delta$ | 175,000 |
| -ne | num_edges | The number of edges in $G_\delta$ | 175,000 |
| -p | density | The minimum edge density $p$ of the altered subnetwork in $G_\delta$ | 0.05 |
| -t | time_limit | Time limit for running the Gurobi solver | 12 (hours) | 
| -o | output | Directory for the NetMix2 output | None


### Output
NetMix2 outputs a list of vertices corresponding to the altered subnetwork $\hat{A}_{NetMix2}$. Each line in the output file is a vertex:

    B
    C


## Additional information
--------
### Tutorial
A tutorial with step-by-step instructions for NetMix2 is available in the Jupyter notebook.

### Running NetMix2 using other subnetwork families
In addition to the propagation family, users can also choose to run NetMix2 using one of the following subnetwork families:

- Connected family $\mathcal{C}_G$

Execution command:  
` python run_netmix2_connected -el [EDGE_LIST] -gs [GENE_SCORES] (-o [OUTDIR])`

- Edge-dense family $\mathcal{E}_{G,p}$

Execution command:  
` python run_netmix2_edge_dense -el [EDGE_LIST] -gs [GENE_SCORES] -p [MINIMUM_EDGE_DENSITY] (-o [OUTDIR])`

- Cut family $\mathcal{T}_{G,\rho}$

Execution command:  
` python run_netmix2_connected -el [EDGE_LIST] -gs [GENE_SCORES] -rho [MAXIMUM_CUTSIZE] (-o [OUTDIR])`


Please refer to the NetMix2 manuscript for defition of each subnetwork family.

## Contacts
----------
NetMix2 has been developed by members of the research group of prof. Ben Raphael at Princeton University.
For any related question, please email Uthsav Chitra (uchitra@princeton.edu) or Tyler Park (typark@princeton.edu).


### License
See `LICENSE` for license information.
### Citation
The NetMix2 manuscript is currently under review and will be available soon.