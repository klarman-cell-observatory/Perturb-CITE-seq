import numpy as np
import pandas as pd
import seaborn as sns
import glob, os
import h5py
from pdb import set_trace as bp
import matplotlib.pyplot as plt
import matplotlib
import scipy.stats as sp_stats
import scanpy as sc
import anndata
from tqdm import tqdm
from statsmodels.stats.multitest import multipletests
import pickle
from adjustText import adjust_text

# Volcano plot given effect size and p-value of each perturbation on program score
def plot_volcano(eff_size, p_vals, dot_labels, fig_save_pn, title='',
				 lfc_thresh=1, p_thresh=0.01, DPI=800, label_fontsize=2,
				 dot_size=1, FDR_correction=False, alpha_val=0.05,
				 FDR_method='fdr_bh', fig_size_mult=1, adjust_labels=True,
				 draw_arrows=False, arrow_style='-', arrow_color='red',
				 shade_program=False, program_inds=None, program_color='r',
				 program_name='', legend=True, arrow_width=0.5):

	fig, ax = plt.subplots(figsize=[fig_size_mult*6.4, fig_size_mult*4.8])

	if shade_program:
		non_prog_inds = np.delete(np.arange(dot_labels.size), program_inds)
		ax.scatter(eff_size[non_prog_inds], -np.log10(p_vals[non_prog_inds]), s=dot_size, c='b')
		ax.scatter(eff_size[program_inds], -np.log10(p_vals[program_inds]), s=dot_size, c=program_color, label=program_name)
	else:
		ax.scatter(eff_size, -np.log10(p_vals), s=dot_size)

	ax.set_xlabel('LFC'), ax.set_ylabel(u'-log\u2081\u2080(p)')

	if FDR_correction:
		label_inds = np.where(multipletests(p_vals, alpha=alpha_val, method=FDR_method)[0])[0]
	else:
		label_inds = np.where( (np.abs(eff_size) > lfc_thresh) & (p_vals < p_thresh) )[0]

	if label_inds.size > 0:
		xlabels = eff_size[label_inds]
		ylabels = -np.log10(p_vals)[label_inds]
		label_text = dot_labels[label_inds]
		text_objs = [ax.text(xlabels[ind], ylabels[ind], label_text[ind], fontsize=label_fontsize, fontweight='bold') for ind in np.arange(label_inds.size)]

		if adjust_labels:
			if draw_arrows:
				adjust_text(text_objs, arrowprops=dict(arrowstyle=arrow_style, color=arrow_color, lw=arrow_width))
			else:
				adjust_text(text_objs)

	ax.set_title(title)
	if shade_program and legend:
		ax.legend()
	fig.savefig(fig_save_pn, bbox_inches='tight', dpi=DPI)
	plt.close()
