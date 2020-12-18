import glob,os
import numpy as np
import scipy.stats as sp_stats
import matplotlib.pyplot as plt
from sklearn.linear_model import ElasticNet
from statsmodels.distributions.empirical_distribution import ECDF
import argparse

def fit_lm(X, y, l1_ratio=0.5, alpha=0.0005, max_iter=10000, z_score=False):

	lmfit = ElasticNet(precompute=True, l1_ratio=l1_ratio, alpha=alpha, max_iter=max_iter)

	if z_score:
		y = sp_stats.zscore(y, axis=0)

	lmfit.fit(X, y)

	return lmfit.coef_
 
 
