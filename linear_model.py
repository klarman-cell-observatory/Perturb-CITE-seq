import glob,os
import numpy as np
import scipy.stats as sp_stats
import matplotlib.pyplot as plt
from sklearn.linear_model import ElasticNet
from statsmodels.distributions.empirical_distribution import ECDF
import argparse
import fastcluster
from scipy.cluster.hierarchy import dendrogram
import seaborn as sns
from statsmodels.stats.multitest import multipletests
from sklearn.cluster import KMeans
from adjustText import adjust_text

## Fitting and clustering

# Fit linear model with elastic net regularization
# Documentation: https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.ElasticNet.html
def fit_lm(X, y, l1_ratio=0.5, alpha=0.0005, max_iter=10000, z_score=False):

	lmfit = ElasticNet(precompute=True, l1_ratio=l1_ratio, alpha=alpha, max_iter=max_iter)

	if z_score:
		y = sp_stats.zscore(y, axis=0)

	lmfit.fit(X, y)

	return lmfit.coef_

# Calculate and cluster correlation matrix
def calc_corr_matrix(coefs):

	corr_mat = pd.DataFrame(coefs).corr().to_numpy()
	corr_mat[np.isnan(corr_mat)] = 0

	row_linkage = fastcluster.linkage(corr_mat, method='complete')
	col_linkage = fastcluster.linkage(corr_mat.T, method='complete')
	row_dend = dendrogram(row_linkage, no_plot=True)
	col_dend = dendrogram(col_linkage, no_plot=True)
	row_inds = np.array(row_dend['leaves'])
	col_inds = np.array(col_dend['leaves'])

	assert(np.all(row_inds==col_inds)) # Symmetric matrix

	return corr_mat, row_inds, col_inds

# Sweep number of clusters for k-means algorithm to determine "elbow"
# Documentation: https://scikit-learn.org/stable/modules/generated/sklearn.cluster.KMeans.html
def calculate_inertia(mat, max_clusters=20):

	clusters_to_test = np.arange(max_clusters)+1

	all_distances = np.zeros(clusters_to_test.size)

	for cluster_i, num_clusters in enumerate(clusters_to_test):
		km = KMeans(n_clusters=num_clusters)
		km.fit(mat)
		all_distances[cluster_i] = km.inertia_

	return clusters_to_test, all_distances

# Cluster matrix using k-means clustering
# Documentation: https://scikit-learn.org/stable/modules/generated/sklearn.cluster.KMeans.html
def cluster_by_kmeans(mat, num_clusters):

	km = KMeans(n_clusters=num_clusters, random_state=3)
	km.fit(mat)

	ret_order = np.array([], dtype=np.int)
	for clust_i in np.arange(num_clusters):
		curr_clust_inds = np.where(clust_i == km.labels_)[0]
		ret_order = np.append(ret_order, curr_clust_inds.astype(np.int))

	assert(ret_order.size == km.labels_.size)

	return ret_order, km.labels_[ret_order]

## Permutation test

# Calculate null distributions for given covariate (specified by "cov_ind" in X matrix)
def shuffle_and_fit(X, y, cov_ind, num_iters=1000):

	# Can pre-allocate "all_nulls" for speed
	
	for iter in np.arange(num_iters):

		X_shuffled = X.copy()
		X_shuffled[:, cov_ind] = np.random.permutation(X_shuffled[:, cov_ind])

		lm_coefs = fit_lm(X_shuffled, y)

		if iter == 0:
			all_nulls = lm_coefs[:, cov_ind].flatten()
		else:
			all_nulls = np.append(all_nulls, lm_coefs[:, cov_ind].flatten())

		del X_shuffled

	return all_nulls

# Calculate empirical p values for a given covariate (specified by "cov_ind") given regulatory matrix and null distributions
def calc_p_vals(beta_mat, null_distrib, cov_ind):

	curr_coeffs = beta_mat[:, cov_ind].flatten()
	curr_ECDF = ECDF(null_distrib)

	p_vals = np.ones(curr_coeffs.size)

	neg_inds = np.where(curr_coeffs < 0)[0]
	pos_inds = np.where(curr_coeffs >= 0)[0]

	p_vals[neg_inds] = curr_ECDF(curr_coeffs[neg_inds])
	p_vals[pos_inds] = 1-curr_ECDF(curr_coeffs[pos_inds])

	return p_vals

## Helper

# Return a design matrix from an array containing sgRNA assignments
def design_from_arrs(sgRNA_names, sgRNA_assignments):

	design_mat = np.zeros((sgRNA_assignments.size, sgRNA_names.size))

	for cell_i in range(sgRNA_assignments.size):
		sgRNA_ind = np.where(sgRNA_assignments[cell_i] == sgRNA_names)[0]
		assert(sgRNA_ind.size == 1)
		design_mat[cell_i, sgRNA_ind] = 1

	return design_mat

## Plotting

# Remove sparse rows and columns from a matrix for improved visualization
def remove_sparse_rows_and_columns(mat, row_names, col_names, approx_zero = 0.02, row_sparse_thresh = 0.7, col_sparse_thresh = 0.7):

	row_sparsity = np.array([ (np.where( np.abs(mat[curr_row, :]) <= approx_zero)[0].size / col_names.size) for curr_row in range(row_names.size)])
	good_rows = np.where(row_sparsity < row_sparse_thresh)[0]
	filtered_rows = row_names[good_rows]

	col_sparsity = np.array([ (np.where( np.abs(mat[:, curr_col]) <= approx_zero)[0].size / row_names.size) for curr_col in range(col_names.size)])
	good_cols = np.where(col_sparsity < col_sparse_thresh)[0]
	filtered_cols = col_names[good_cols]

	filtered_mat = mat[good_rows, :][:, good_cols]

	return filtered_mat, filtered_rows, filtered_cols

# Plot cluster map
def plot_clustermap(plot_arr, row_inds, col_inds, row_names, col_names, save_pn,
					remove_sparse_rows=False, row_sparse_thresh=0.8, clim=1,
					cmap_name='seismic', DPI=800, xlabel='', ylabel='',
					label_ticks=False, fig_size_mult=1):

	fig, ax = plt.subplots(figsize=[fig_size_mult*6.4, fig_size_mult*4.8])
	filtered_plot_mat = plot_arr[row_inds, :][:, col_inds]

	if remove_sparse_rows:
		row_sparsity = np.array([ (np.where(filtered_plot_mat[curr_row, :] == 0)[0].size / col_inds.size) for curr_row in range(row_inds.size)])
		good_rows = np.where(row_sparsity < row_sparse_thresh)[0]
		filtered_plot_mat = filtered_plot_mat[good_rows, :]
		row_inds = row_inds[good_rows]

	ax_ret = ax.imshow(filtered_plot_mat, aspect=(filtered_plot_mat.shape[1]/filtered_plot_mat.shape[0]), cmap=cmap_name, clim=[-clim, clim])
	cbar = fig.colorbar(ax_ret, ax=ax, extend='both')
	cbar.minorticks_on()

	if label_ticks:
		ax.set_yticks(np.arange(row_inds.size))
		ax.set_yticklabels(row_names[row_inds], fontsize=1)

		ax.set_xticks(np.arange(col_inds.size))
		ax.set_xticklabels(col_names[col_inds], fontsize=1, rotation=90)
	else:
		ax.set_yticks([])
		ax.set_xticks([])

	ax.set_ylabel(ylabel)
	ax.set_xlabel(xlabel)

	ax.tick_params(width=0.1)
	ax.set_ylim([filtered_plot_mat.shape[0]-0.5, -0.5])

	fig.savefig(save_pn, bbox_inches='tight', dpi=DPI)
	plt.close()
