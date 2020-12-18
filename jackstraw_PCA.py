import numpy as np
import pandas as pd
import seaborn as sns
import spams
import os,glob
from pdb import set_trace as bp
import scanpy as sc
import anndata
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from tqdm import tqdm
from statsmodels.distributions.empirical_distribution import ECDF
import statsmodels.stats.multitest as stats_mt

# Calculate null distributions for each PC
def calc_PCA_null(X_mat, n_cells=20000, n_genes=150, n_reps=100, n_comps=50):

  # Can pre-allocate null_distribs for speed

	# Generate null distributions for each PC
	for rep_i in tqdm(range(n_reps)):

		# Randomly select n_cells and n_genes from each cell from X_mat
		cell_inds = np.random.choice(np.arange(X_mat.shape[0]), size=n_cells, replace=False)
		gene_inds = np.random.choice(np.arange(X_mat.shape[1]), size=n_genes, replace=False)
		X_mat_downsample = X_mat[cell_inds, :][:, gene_inds]

		pca = PCA(n_components=n_comps)
		pca.fit(X_mat_downsample.T)
		X_pca = pca.fit_transform(X_mat_downsample.T)

		if rep_i == 0:
			null_distribs = X_pca.T
		else:
			null_distribs = np.hstack((null_distribs, X_pca.T))

	return null_distribs
 
 # Calculate p values by gene given full PCA and null distributions
 def calc_gene_pvals(PCA_comps, null_distribs):

	assert(PCA_comps.shape[0] == null_distribs.shape[0])
	num_comps = PCA_comps.shape[0]

	p_vals = np.ones(PCA_comps.shape)
	for comp_i in range(num_comps):
		curr_ECDF = ECDF(null_distribs[comp_i, :])

		neg_inds = np.where(PCA_comps[comp_i, :] < 0)[0]
		pos_inds = np.where(PCA_comps[comp_i, :] >= 0)[0]

		p_vals[comp_i, neg_inds] = curr_ECDF(PCA_comps[comp_i, neg_inds])
		p_vals[comp_i, pos_inds] = 1-curr_ECDF(PCA_comps[comp_i, pos_inds])

	return p_vals
