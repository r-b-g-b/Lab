import numpy as np
from numpy.linalg import lstsq
from scipy.stats import norm
import matplotlib.pyplot as plt


class boots_lr():
	
	def __init__(self, X, y, r, alpha = 0.05):

		'''
		run bias-corrected, accelerated bootstrap test for mean value
		Applied Regression Analysis and Generalized Linear Models, Second Edition, John Fox
		ch. 21 - Bootstrapping Regression Models
		'''
		
		self.X = X
		self.y = y
		self.alpha = alpha
		self.r = r
		self.nsamp = y.size
		self.params = {}
	
	def fit(self):
		
		beta1_hat, beta0_hat = lstsq(self.X, self.y)[0]
		
		beta1_star, beta0_star = self.boots()
		beta1_star.sort()
		beta0_star.sort()
	
		Z_beta1 = self.calc_Z(beta1_hat, beta1_star)
		Z_beta0 = self.calc_Z(beta0_hat, beta0_star)
	
		A_beta1, A_beta0 = self.calc_A()
	
		ci_beta1 = self.calc_CI(A_beta1, Z_beta1)
		ci_beta0 = self.calc_CI(A_beta0, Z_beta0)
	
		CI_beta1 = np.array((beta1_star[ci_beta1[0]], beta1_star[ci_beta1[1]]))
		CI_beta0 = np.array((beta0_star[ci_beta0[0]], beta0_star[ci_beta0[1]]))
	
		self.B = np.array((beta0_hat, beta1_hat))
		self.B_int = np.vstack((CI_beta0, CI_beta1))
		self.B_boot = np.vstack((beta0_star, beta1_star)).T
	
	def plot(self, ax = None, color = 'b'):
	
		if ax is None:
			fig = plt.figure();
			ax = fig.add_subplot(111);
	
		ax.scatter(self.X[:, 0], self.y, c = color)
		xs_ = ax.get_xlim()
		xs = np.linspace(0, xs_[1], 20)
		ys = self.B[0] + self.B[1] * xs
		ax.plot(xs, ys, color = color, lw = 2, ls = '-')
	
		# create a matrix with every possible regression line from the bootstrapping analysis
		linesmat = np.zeros((self.r, xs.size))
		for i in range(self.r):
			linesmat[i, :] = self.B_boot[i][0] + self.B_boot[i][1]*xs
	
		# using the bootstrapped data, find (1-(alpha/2))*100% and (alpha/2)*100% ranked values to make the CI boundary
		ci_band = np.zeros((2, xs.size));
		n_max = np.round((1-(self.alpha/2))*self.r+0.05); # nearest rank for (1-(alpha/2))*100%
		n_min = np.round((self.alpha/2)*self.r+0.05); # nearest rank for (alpha/2)*100%
		for i in range(linesmat.shape[1]):
			tmp = np.sort(linesmat[:,i])
			ci_band[0, i] = tmp[n_max]
			ci_band[1, i] = tmp[n_min]

		ax.plot(xs, ci_band[0, :], color = color, lw = 2, ls = '--'); # lower bound
		ax.plot(xs, ci_band[1, :], color = color, lw = 2, ls = '--'); # upper bound
		plt.show();


	def boots(self):
	
		beta1 = np.zeros(self.r)
		beta0 = np.zeros(self.r)
		for i in range(self.r):
			ix = np.random.random_integers(0, self.nsamp-1, self.nsamp)
			q = lstsq(self.X[ix, :], self.y[ix])
			beta1[i], beta0[i] = q[0]
	
		return beta1, beta0
	
	def calc_CI(self, A, Z):
		'''
		A	:	from calc_A
		Z	:	from calc_Z
		alpha	:	confidence bound
		'''
		z = norm.isf(self.alpha/2.)
		A1_ = Z + (Z - z) / (1. - A*(Z - z))
		A1 = norm.cdf(A1_)
	
		A2_ = Z + (Z + z) / (1. - A*(Z + z))
		A2 = norm.cdf(A2_)
	
		lo = np.int32(A1*self.r)
		up = np.int32(A2*self.r)
	
		ci = np.array((lo, up))
	
		return ci

	@staticmethod
	def calc_Z(theta_hat, theta_star):
		'''
		theta_hat	:	the parameter estimate
		theta_star	:	vector of bootstrapped parameter estimates
		'''
	
		Z_ = (theta_star<theta_hat).mean() + (theta_star==theta_hat).mean()/2.
		Z = norm.ppf(Z_)
	
		return Z
	
	def calc_A(self):
		'''
		theta_bar is the average of jackknifed (leave-one-out) parameter estimates
		'''

		beta1_jack = np.zeros(self.nsamp)
		beta0_jack = np.zeros(self.nsamp)
		for i in range(self.nsamp):
			ix = np.concatenate((np.arange(i), np.arange(i+1, self.nsamp)))
			beta1_jack[i], beta0_jack[i] = lstsq(self.X[ix, :], self.y[ix])[0]
	
		beta1_bar = beta1_jack.mean()
		beta0_bar = beta0_jack.mean()
	
		score = (beta1_bar - beta1_jack)
		A_num = (score**3).sum()
		A_denom = 6*((score**2).sum())**(3/2.)
		A_beta1 = A_num / A_denom
	
		A_num = ((beta0_bar - beta0_jack)**3).sum()
		A_denom = 6*(((beta0_bar - beta0_jack)**2).sum())**(3/2.)
		A_beta0 = A_num / A_denom
	
		return A_beta1, A_beta0

