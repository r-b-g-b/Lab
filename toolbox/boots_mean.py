import numpy as np
from numpy.linalg import lstsq
from scipy.stats import norm


def run(x, r):
	'''
	run bias-corrected, accelerated bootstrap test for mean value
	Applied Regression Analysis and Generalized Linear Models, Second Edition, John Fox
	ch. 21 - Bootstrapping Regression Models
	'''
	
	theta_hat = x.mean()
	theta_star = boots(x, r)
	theta_star.sort()
	
	Z = calc_Z(theta_hat, theta_star)
	
	A = calc_A(x)
	
	ci = calc_CI(A, Z, r, 0.05)
	
	CI = np.array([theta_star[ci[0]], theta_star[ci[1]]])

	return theta_hat, CI

def boots(x, r):
	
	n = x.shape[0]
	theta_star = np.zeros(r, dtype = np.float32)
	for i in range(r):
		ix = np.random.random_integers(0, n-1, n)
		theta_star[i] = x[ix].mean()
	
	return theta_star
	
def calc_CI(A, Z, r, alpha):
	'''
	A	:	from calc_A
	Z	:	from calc_Z
	alpha	:	confidence bound
	'''
	z = norm.isf(alpha/2.)
	A1_ = Z + (Z - z) / (1. - A*(Z - z))
	A1 = norm.cdf(A1_)
	
	A2_ = Z + (Z + z) / (1. - A*(Z + z))
	A2 = norm.cdf(A2_)
	
	lo = np.int32(A1*r)
	up = np.int32(A2*r)
	
	ci = np.array((lo, up))
	
	return ci

def calc_Z(theta_hat, theta_star):
	'''
	Calculates the bias of the bootstrap distribution. This is simply
	the proportion of bootstrapped samples that are less than the 
	parameter's point estimate
	theta_hat	:	the parameter estimate
	theta_star	:	vector of bootstrapped parameter estimates
	'''
	
	r = np.float32(theta_star.shape[0])
	Z_ = (theta_star<theta_hat).sum() / r
	Z = norm.ppf(Z_)	
	
	return Z
	
def calc_A(x):
	'''
	theta_bar is the average of jackknifed (leave-one-out) parameter estimates
	'''

	n = x.shape[0]
	theta_jack = np.zeros(n)
	for i in range(n):
		ix = np.concatenate((np.arange(i), np.arange(i+1, n)))
		theta_jack[i] = x[ix].mean()
	
	theta_bar = theta_jack.mean()
	A_num = ((theta_bar - theta_jack)**3.).sum()
	A_denom = 6.*((((theta_bar - theta_jack)**2.).sum())**(3/2.))
	A = A_num / A_denom
	
	return A
