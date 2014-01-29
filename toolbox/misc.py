import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import pandas as pd
import os, h5py
import scipy.stats as st
import itertools
import Tkinter
import tkFileDialog
import datetime

def get_subplot_grid(n):

	if n<4:
		nrows = n; ncols = 1
	else:
		nrows = ncols = np.ceil(np.sqrt(n))

	return nrows, ncols

def str2date(s, delimiter = '_', format = 'MMDDYYYY'):
	'''
	Input:
		s: a string representation of a date MO_DA_YR
	Output:
		date: datetime formatted date
	'''
	if delimiter=='':
		if format == 'YYYYMMDD':
			yr = s[:4]; mo = s[4:6]; da = s[6:]
	
	else:
		if format == 'MMDDYYYY':
			mo, da, yr = s.split(delimiter)
		elif format == 'YYYYMMDD':
			yr, mo, da = s.split(delimiter)
	

	return datetime.date(int(yr), int(mo), int(da))

def pd_sem(ds):

	n = len(ds)
	ds_mean = np.mean(ds)
	runningsum = 0
	for i in ds:
		runningsum = runningsum + (i-ds_mean)**2
	std = np.sqrt(runningsum / (n-1))
	sem = std / np.sqrt(n)
	return sem

def pd_sem_2d(ds):

	n = len(ds)
	nsamp = len(ds[ds.index[0]])
	ds_mean = np.mean(ds)
	runningsum = np.zeros(nsamp)
	for i in ds:
		runningsum = runningsum + (i-ds_mean)**2
	std = np.sqrt(runningsum / (n-1))
	sem = std / np.sqrt(n)
	return sem

def get_first(df):
    x = df.iloc[0]
    return x

def samexaxis(ax, xlim = None):
	
	minx = +np.inf
	maxx = -np.inf
	
	if xlim is None:
		for ax_ in ax:
			xlim = ax_.get_xlim()
			minx = np.min([minx, xlim[0]])
			maxx = np.max([maxx, xlim[1]])
	else:
		(minx, maxx) = xlim
		
	[a.set_xlim([minx, maxx]) for a in ax]
	
def sameyaxis(ax, ylim = None):

	minx = +np.inf
	maxx = -np.inf
	
	if ylim is None:
		for ax_ in ax:
			xlim = ax_.get_ylim()
			minx = np.min([minx, xlim[0]])
			maxx = np.max([maxx, xlim[1]])
	else:
		(minx, maxx) = ylim
		
	[a.set_ylim([minx, maxx]) for a in ax]
	
def closest(vec, x, log = False):
	
	if log:
		orig_vec = vec.copy()
		vec = np.log2(orig_vec)
		x = np.log2(x)
	differ = np.abs(vec - x)
	min_differ = differ.min()
	ix = differ.argmin()
	if log:
		val = orig_vec[ix]
		error = vec[ix] - x
	else:
		val = vec[ix]
		error = val - x


	return val, ix, error
	
def vline(ax, xloc, **kwargs):
	
	if len(xloc) == 1:
		xloc = [xloc]
	
	for ax_ in ax:
		for xloc_ in xloc:
			ax_.axvline(xloc_, **kwargs)

def hline(ax, yloc, **kwargs):
		
	if len(yloc) == 1:
		yloc = [yloc]
	
	for ax_ in ax:
		for yloc_ in yloc:
			ax_.axhline(yloc_, **kwargs)

def nansem(x, axis = None):
	
	return st.sem(x[~np.isnan(x)], axis = axis)
		
def isbetween(x, xmin, xmax):
	x = np.array(x)
	ix = np.logical_and(xmin <= x, x < xmax)
	return ix

def octave_diff(f1, f2):
	
	return np.log2(f1)-np.log2(f2)
	
def bar2(y, x = None, yerr = None, width = None, ax = None, groups = None, grouplabels = None):
	
	if x is None:
		x = np.arange(y.shape[0])
		
	if width is None:
		width = y.shape[0]
		
	x = np.asarray(x)
	
	ax, fig = axis_check(ax)

	ndim = len(y.shape)
	
	ax.bar(x, y[:, 0], yerr = yerr[:, 0], ecolor = 'k')
	ax.bar(x+width, y[:, 1], yerr = yerr[:, 1], color = 'g', ecolor = 'k')

	if groups is not None:
		labels = []
		for i in itertools.product(groups[0], groups[1]):
			labels.append(str(i[0]) + '.' + str(i[1]))
		
		ax.set_xticks(np.arange(np.prod(y.shape))+0.5)
		ax.set_xticklabels(labels)
	
	
	return ax
	
def jitter(x, amp):
	'''
	returns a Gaussian jittered version of the input. Only good for 1-D vectors.
	'''
	return x + amp*np.random.randn(np.asarray(x).size)
	
def order(x):
	x_ = np.empty_like(x)
	ux = np.unique(x)

	for i, ux_ in enumerate(ux):
		x_[x==ux_] = i
		
	return x_
	
def plot_spread_y(x):
	
	npts, nobs = x.shape
	y_ = 2*np.abs(x).max(1).mean()
	y = np.arange(0, y_*nobs, y_)
	return x + np.tile(y, (npts, 1))

def errorfill(x, y, yerr = None, ax = None, color = None, err_type = 'sem', **kwargs):

	ax, fig = axis_check(ax)

	if color is None:
		color = ax._get_lines.color_cycle.next()

	if len(y.shape)==1:
		y_mean = y
	else:
		y_mean = y.mean(0)

	if yerr is None:
		if err_type == 'sem':
			yerr = 0.5*st.sem(y, 0)
		elif err_type == 'std':
			yerr = 0.5*st.std(y, 0)

	x_hi = y_mean + yerr
	x_lo = y_mean - yerr
	l = ax.plot(x, y_mean, color = color, **kwargs)
	l_up = ax.plot(x, y_mean+yerr, color = color, alpha = 0.2)
	l_lo = ax.plot(x, y_mean-yerr, color = color, alpha = 0.2)
	ax.fill_between(x, x_hi, x_lo, alpha = 0.2, color = color)
	
	return ax

def bin(x, bins, ordered_x = False):
	
	# returns a 'rounded down' version of data (i.e. the cloeset bin value below the raw value)
	bins = np.asarray(bins)
	nbins = bins.size - 1
	x_bin = np.empty(x.size) * np.nan
	for i in range(nbins):
		x_bin[isbetween(x, bins[i], bins[i+1])] = i
	if not ordered_x:
		x_bin = np.asarray([bins[x_] for x_ in x_bin])
		
	return x_bin
		
def slopeintercept(m, b, ax, **kwargs):
	(xi, xf) = ax.get_xlim()
	yi = m*xi + b
	yf = m*xf + b
	ax.plot([xi, xf], [yi, yf], **kwargs)

def scatterfit(x, y, ax = None, **kwargs):
	
	if ax is None:
		fig = plt.figure();
		ax = fig.add_subplot(111);
		
	ax.scatter(x, y)
	(slope, inter) = np.polyfit(x, y, 1)
	slopeintercept(slope, inter, ax, color = 'r')
	
def axis_check(ax, **kwargs):
	'''
	Returns ax, fig
	'''
	if ax is None:
		fig = plt.figure();
		ax = fig.add_subplot(111, **kwargs);
	else:
		fig = ax.get_figure();
	
	return ax, fig
	
def target(x, y, ax, color = 'r'):
	
	ax.plot(x, y, marker = 'o', color = color, ms = 12)
	ax.plot(x, y, marker = 'o', color = 'w', ms = 10)
	ax.plot(x, y, marker = 'o', color = color, ms = 6)
	ax.plot(x, y, marker = 'o', color = 'w', ms = 3)
	
def find_on_off(x):

	x = np.concatenate((np.array([0]), x, np.array([x[-1]+10])))
	d = np.diff(x)
	ix = (d>10).nonzero()[0]
	offset = x[ix][1:]+1000
	onset = x[ix+1][:-1]-1000

	return onset, offset

def find_consec_middles(Q):
	'''
	Input:
		Q	:	list of integers

	Output
		Q	: 	same list except culled down to not contain consecutive ranges
				(new list only contains the middle number of that range of consecutive numbers)
	'''
	X = [True]
	while sum(X) > 0:
		L = [l+1 == Q[i+1] for i, l in enumerate(Q[:-1])] # check left
		L.extend([False])
		R = [l-1 == Q[i] for i, l in enumerate(Q[1:])] # check right
		R.insert(0, False)
		M = [l and r for (l, r) in zip(L, R)] # middle values
		L2 = [l if not m else False for (l, m) in zip(L, M)]
		R2 = [r if not m else False for (r, m) in zip(R, M)]
		Z = [False if not l else l and r for (l, r) in zip(L2[:-1], R2[1:])] # no LR pairs
		Z.extend([False])
		Y = [l ^ r for (l, r) in zip(L, R)] # left or right
		X = [y if not z else False for (y, z) in zip(Y[:-1], Z)] # what to remove (L/R but not LR pairs)
		X.extend([Y[-1]])
		Q = [q for (q, x) in zip(Q, X) if not x]

	return Q

def find_edges(Q):
	dQ = np.diff(Q)
	return (dQ == 1).nonzero()[0], (dQ == -1).nonzero()[0]
	
def find_peaks(Q, pplot = False):
	dQ = np.diff(Q)
	w = dQ[:-1] * dQ[1:]
	peaks_x = np.vstack((w<0, dQ[:-1]>0)).all(0).nonzero()[0]+1
	peaks_y = Q[peaks_x]
	
	if pplot:
		fig = plt.figure()
		ax = fig.add_subplot(111)
	
		ax.plot(Q, 'b')
		ax.plot(dQ, 'r')
		ax.plot(w*1E12, 'g')
		ax.plot(peaks_x, peaks_y, '*m')
	
	return peaks_x, peaks_y
	
def cos2ramp(x, ramplen):
	
	ramplen = np.int32(ramplen)
	rampoff = np.cos(np.pi/(2*ramplen) * np.arange(ramplen))**2
	rampon = rampoff[::-1]
	x[:ramplen] = x[:ramplen] * rampon
	x[-ramplen:] = x[-ramplen:] * rampoff
	
	return x
	
	return x2
	
def medfilt1(x=None,L=None):

	'''
	a simple median filter for 1d numpy arrays.

	performs a discrete one-dimensional median filter with window
	length L to input vector x. produces a vector the same size 
	as x. boundaries handled by shrinking L at edges; no data
	outside of x used in producing the median filtered output.
	(upon error or exception, returns None.)

	inputs:
		x, Python 1d list or tuple or Numpy array
		L, median filter window length
	output:
		xout, Numpy 1d array of median filtered result; same size as x

	bdj, 5-jun-2009
	'''

	# input checks and adjustments 
	try:
		N = len(x)
		if N < 2:
			print 'Error: input sequence too short: length =',N
			return None
		elif L < 2:
			print 'Error: input filter window length too short: L =',L
			return None
		elif L > N:
			print 'Error: input filter window length too long: L = %d, len(x) = %d'%(L,N)
			return None
	except:
		print 'Exception: input data must be a sequence'
		return None

	xin = np.array(x)
	if xin.ndim != 1:
		print 'Error: input sequence has to be 1d: ndim =',xin.ndim
		return None

	xout = np.zeros(xin.size)

	# ensure L is odd integer so median requires no interpolation
	L = int(L)
	if L%2 == 0: # if even, make odd
		L += 1
	else: # already odd
		pass
	Lwing = (L-1)/2

	for i,xi in enumerate(xin):

		# left boundary (Lwing terms)
		if i < Lwing:
			xout[i] = np.median(xin[0:i+Lwing+1]) # (0 to i+Lwing)

		# right boundary (Lwing terms)
		elif i >= N - Lwing:
			xout[i] = np.median(xin[i-Lwing:N]) # (i-Lwing to N-1)
   
		# middle (N - 2*Lwing terms; input vector and filter window overlap completely)
		else:
			xout[i] = np.median(xin[i-Lwing:i+Lwing+1]) # (i-Lwing to i+Lwing)

	return xout	
	

def get_path():

	root = Tkinter.Tk()
	root.withdraw()
	
	path = tkFileDialog.askopenfilename()
	root.quit()
	
	return path


def export_h5(f, fname = 'OUT'):

	fout = open(os.path.join('/Users/robert/Desktop/', fname+'.txt'), 'w')
	header = ['gen', 'exp', 'sess']
	for i, sess in enumerate(f.keys()):
		fsess = f[sess]
		gen, exp, _ = sess.split('_')
		for j, unit in enumerate(fsess.keys()):
			funit = fsess[unit]
			newline = [gen, exp, sess]
			for m, field in enumerate(np.sort(funit.keys())):
				ffield = funit[field]
				fieldsizes = np.array(ffield.shape)
				if fieldsizes.size == 0:
					if i==0 and j==0:
						header.append(field)
					newline.append('%4.4f' % ffield.value)
				elif fieldsizes.size == 1 and fieldsizes[0]<10:
					for n in range(fieldsizes[0]):
						if i==0 and j==0:
							header.append('%s-%2.2u' % (field, n))
						newline.append('%4.4f' % ffield.value[n])
							
			if i==0 and j==0:
				fout.write('\t'.join(header) + '\n')	
			fout.write('\t'.join(newline) + '\n')
	
	fout.close()
			

def export_DB(DB, fname = 'DB'):

	nunits = DB.shape[0]
	txt = np.empty((nunits+1, 0), dtype = 'str')

	for name in DB.dtype.names:
		if len(DB[name].shape) == 1: # save out 1D DB entries
			txt = np.hstack((txt, np.vstack((name, DB[name][:, np.newaxis]))))
		elif len(DB[name].shape) == 2: # save out 2D DB entries
			if DB[name].shape[1] < 10:
				for i in range(DB[name].shape[1]):
					headername = '%s%i' % (name, i)
					txt = np.hstack((txt, np.vstack((headername, DB[name][:, i, np.newaxis]))))

	np.savetxt(os.path.join('/Users/robert/Desktop/', fname+'.txt'), txt, fmt = '%s')


def export_DB2(data, fname = 'DB'):
	
	savepath = os.path.join('/Users/robert/Desktop', fname+'.txt')

	assert not os.path.exists(savepath)

	f = open(savepath, 'w')
	line = ''
	dat = data[0]
	for key_ in dat.dtype.descr:
		key = key_[0]
		newfield = dat[key]
		if newfield.size==1:
			newfield = [newfield]
		for i, newfield_ in enumerate(newfield):
			header = key
			if len(newfield)>1:
				header = '%s%2.2i' % (header, i+1)
			line += '%s\t' % header
	line += '\n'
	f.write(line)
	
	for dat in data:
		line = ''
		for key_ in dat.dtype.descr:
			key = key_[0]
			newfield = dat[key]
			if newfield.size==1:
				newfield = [newfield]
			for newfield_ in newfield:
				if not type(newfield_) in [str, np.string_]:
					newfield_ = '%.5e' % newfield_
				line += '%s\t' % newfield_
		line += '\n'
		f.write(line)

	f.close()
		
def make_combinations(x):
	ndim = x.shape[1]
	nlevels = np.empty(ndim, dtype = np.int32)
	ulevels = []
	for i in range(ndim):
		ulevels.append(list(np.unique(x[:, i])))
		nlevels[i] = len(ulevels[-1])
	
	for i in itertools.product(*ulevels):
		print i
		
def np_to_pd(x):
	
	df = pd.DataFrame(index = np.arange(x.size))
	for i in x.dtype.descr:
		if len(i)==2:
			name, dtype = i
			df[name] = x[name]
		elif len(i)==3:
			name, dtype, shape = i
			if len(shape)==1 and shape[0]<100:
				for j in xrange(shape[0]):
					df['%s%2.2i' % (name, j+1)] = x[name][:, j]
		
	return df
	
	
def get_rast_stimparams_from_df(df_, studyprefix = 'rr', unitprefix = 'RR', basedir = '/Volumes/BOB_SAGET/Fmr1_RR'):
	
	path = os.path.join(basedir, 'Sessions', '%s_%s_%s_%s' % (studyprefix, df_['gen'], df_['exp'], df_['sess']), 'fileconversion', '%s%3.3i.h5' % (unitprefix, df_['unit']))
	if os.path.exists(path):
		print 'Loading %s' % path
		f = h5py.File(path, 'r')
	
		rast = f['rast'].value
		stimparams = f['stimID'].value

		f.close()
	else:
		print 'Could not find %s' % path
		return None, None
		
	return rast, stimparams

def add_unity_line(ax = None, **kwargs):
	if ax is None:
		ax = plt.gca()
	
	xlim = ax.get_xlim()
	ylim = ax.get_ylim()
	minlim = min((xlim[0], ylim[0]))
	maxlim = max((xlim[1], ylim[1]))

	ax.plot([minlim, maxlim], [minlim, maxlim], **kwargs)
	
	ax.set_xlim(xlim)
	ax.set_ylim(ylim)
	
	plt.draw()	
	
	
	
def plot_matrix(x, N = 10, ax = None, **kwargs):
	'''
	plots the mean trace of a matrix and transparent sample traces
	'''
	ax, fig = axis_check(ax)
	
	l, = ax.plot(x.mean(1), lw = 2, **kwargs)
	color = l.get_color() 
	ax.plot(x[:, np.linspace(0, x.shape[1]-1, N).astype(int)], alpha = 0.1, color = color, **kwargs)
	
	return ax
	
	
def pd_errorbar(y, yerr, ax = None, **kwds):
	
	ax, fig = axis_check(ax)

	colors = kwds.pop('color', 'brgyk')

	ax_pos = np.arange(len(y)) + 0.25
	tickoffset = 0.375

	K = y.columns.size
	rects = []
	labels = []
	for i, label in enumerate(y.columns):
		y_ = y[label]
		yerr_ = yerr[label]
		kwds['color'] = colors[i % len(colors)]
		rect = ax.bar(ax_pos + i*0.75/K, y_, 0.75/K, label = label, **kwds)
		ax.errorbar(ax_pos + (i*0.75/K) + (0.75/K)/2, y_, yerr_, fmt = None, ecolor = 'k')
		rects.append(rect)
		labels.append(label)
	
	patches = [r[0] for r in rects]
	ax.legend(patches, labels, loc = 'best')
	ax.set_xticks(ax_pos + tickoffset)
	ax.set_xticklabels(y.index.values, rotation = 90)
	plt.show();		

	return ax

def objectarray2floatarray(x):
	x_ = np.empty((len(x), len(x[x.index[0]])))
	for i, (k, v) in enumerate(x.iterkv()):
		x_[i, :] = v
	return x_

def facet_wrap(df, df_err = None, fig = None):
	
	if fig is None:
		fig = plt.figure();

	index = zip(*df.index)
	nindex = len(index)
	ulevels = np.unique(index[-1])
	nlevels = ulevels.size


	ax = []
	for i, lev in enumerate(ulevels):
		ax.append(fig.add_subplot(nlevels, 1, i+1))
		df.xs(lev, level=nindex-1).transpose().plot(kind = 'line', ax = ax[-1])
		# df.xs(lev, level=nindex-1)
		# df_err.xs(lev, level=nindex-1)
		# pd_errorbar(df, df_err, ax = ax[-1])

	sameyaxis(ax)








