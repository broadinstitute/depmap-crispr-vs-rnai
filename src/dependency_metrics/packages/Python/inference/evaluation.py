from __future__ import print_function
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
from sklearn.metrics import roc_auc_score

def ssmd(positive_controls, negative_controls):
	return (positive_controls.mean()- negative_controls.mean()) / (positive_controls.var() + negative_controls.var())

def roc_auc(positive_controls, negative_controls):
	return roc_auc_score(
		[0]*len(positive_controls) + [1] * len(negative_controls), 
		list(positive_controls) + list(negative_controls)
		)

def overlap(positive_controls, negative_controls, bw_positive='scott', bw_negative='scott',
		gridsize=10000
	):
	est_pos = gaussian_kde(positive_controls, bw_method=bw_positive)
	est_neg = gaussian_kde(negative_controls, bw_method=bw_negative)
	points = np.linspace(
		min([min(positive_controls), min(negative_controls)]), 
		max([max(positive_controls), max(negative_controls)]),
		gridsize
		)
	overlap = (points[1] - points[0]) * np.sum(np.maximum(est_pos(points), est_neg(points)))
