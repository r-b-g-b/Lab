import misc; # reload(misc);
import h5py
import itertools
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
import RF; reload(RF);

	
def vR_dist(psth_1, psth_2, tau = 0.03):
	
	# psth_1 = calc_psth(spktimes_1, nbins)
	# psth_2 = calc_psth(spktimes_2, nbins)
	
	psth_smoo_1 = exp_smoo(psth_1, tau = tau)
	psth_smoo_2 = exp_smoo(psth_2, tau = tau)
	
	return (1./(1000*tau)) * ((psth_smoo_1 - psth_smoo_2)**2).sum()
	
def hamming_smoo(psth, windlen = 10):
	
	wind = np.hamming(windlen)
	wind = wind / wind.sum()
	
	psth_smoo = np.convolve(psth, wind, 'same')
	return psth_smoo
	
def exp_smoo(psth, tau = 0.01):
	
	t_wind = np.arange(0, 0.030, 0.001)
	wind = np.exp(-t_wind/tau)
	wind = wind / wind.sum()

	psth_smoo = np.convolve(psth, wind)
	psth_smoo = psth_smoo[0:psth.size]

	return psth_smoo

def calc_psth(x, nbins = None, ntrials = None, normed = False):

	if len(x.shape) == 1:
		psth = calc_psth_from_spktimes(x, nbins, ntrials = ntrials, normed = False)
	elif len(x.shape) == 2:
		psth = calc_psth_from_rast(x, normed = False)
	
	return psth

def calc_psth_from_rast(rast, normed = False):
	
	psth = rast.sum(0)
	if normed:
		psth = 1000 * psth / np.float(rast.shape[1])
	
	return psth
	
def calc_psth_from_spktimes(spktimes, nbins, ntrials = None, normed = False):
 	'''
	Output:
		psth : spike counts if normed == False, firing rate (spk/s) if normed == True
	'''
	
	spktimes_ms = np.array(spktimes * 1000, dtype = np.int32)
	psth = np.zeros(nbins)
	if spktimes_ms.size == 1:
		spktimes_ms = [spktimes_ms]
	for spk in spktimes_ms:
		psth[spk] = psth[spk] + 1

	if normed:
		psth = 1000 * psth / np.float32(ntrials)
		
	return psth

def calc_rast(spktimes, spktrials, ntrials, nbins):

	if np.isnan(spktimes[0]):
		print 'No spikes!'
		
	else:
		spktimes = 1000 * spktimes.round(3)
		rast = np.zeros((ntrials, nbins))
		for spktime, spktrial in zip(spktimes, spktrials):
			rast[spktrial, spktime] += 1
		
		return rast

def calc_psth_by_stim(rast, stimparams, bins = 0.001):

	nstimparams = stimparams.shape[1]
	
	usp = []
	for i in range(nstimparams):	
		usp.append(list(np.unique(stimparams[:, i])))

	nparamlevels = np.empty(nstimparams, dtype = np.int32)
	for i in range(nstimparams):
		nparamlevels[i] = len(usp[i])

	ntrials_per_stim = np.zeros(nparamlevels)

	'''
	compute
	nbins	:	the number of bins
	bins	:	times (seconds) for each bin
	'''
	dur_ms = rast.shape[1] # number of milliseconds
	t_ms = np.arange(dur_ms) # time indices in ms
	if type(bins) in (int, np.int32):
		nbins = bins
		bins = np.linspace(0, dur_ms, nbins) / 1000
		bindur = bins[1]-bins[0]
	elif type(bins) in (float, np.float):
		bindur = bins
		bins = np.arange(0, (dur_ms/1000.)+bindur, bindur)
		nbins = bins.size-1
	elif type(bins) is np.ndarray:
		nbins = bins.size-1
		bindur = bins[1]-bins[0]
	

	assert np.unique(ntrials_per_stim).size == 1

	ntrials = np.int32(ntrials_per_stim[0])
	
	psth_shape = np.hstack((nparamlevels, nbins))
	psth = np.zeros(psth_shape)
	
	combinations = []
	combinations_ix = []
	for i in itertools.product(*usp):
		combinations.append(i)
		combinations_ix_ = []
		for j, i_ in enumerate(i):
			q = (np.array(usp[j])==i_).nonzero()[0][0]
			combinations_ix_.append(q)
		combinations_ix.append(combinations_ix_)
		
	for m, n in zip(combinations, combinations_ix):
		ix = RF.get_trials(stimparams, m)
		ntrials = ix.size
		spktimes_ = rast2spktimes(rast[ix, :])
		psth_, edges_ = np.histogram(spktimes_, bins = bins)
		psth_ = (1./bindur) * (psth_.astype(float)/ntrials)
		psth[tuple(n)] = psth_
	
	psth = psth.squeeze()		
	return psth, usp


def calc_rast_by_stim(rast, stimparams, bins = None):

	
	usp1 = np.unique(stimparams[:, 0])
	usp2 = np.unique(stimparams[:, 1])
		
	nstimparams = stimparams.shape[1]
	ustimparams = [np.unique(stimparams[:, 0]), np.unique(stimparams[:, 1])]
	nparamlevels = np.array([usp1.size, usp2.size])

	nbins = rast.shape[1]
	
	ntrials_per_stim = np.zeros((nparamlevels[0], nparamlevels[1]))
	for i in range(nparamlevels[0]):
		for j in range(nparamlevels[1]):
			ntrials_per_stim[i, j] = RF.get_trials(stimparams, np.array([usp1[i], usp2[j]])).size
	
	if np.unique(ntrials_per_stim).size > 1:
		print 'Different numbers of trials per stimulus!'
	elif np.unique(ntrials_per_stim).size == 1:
		ntrials = np.int32(ntrials_per_stim[0])
	
	rast2 = np.zeros((ntrials, nparamlevels[0], nparamlevels[1], nbins))
	for i in range(nparamlevels[0]):
		for j_ in range(nparamlevels[1]):
			
			ix = RF.get_trials(stimparams, (usp1[i], usp2[j_]))
			rast2[:, i, j_, :] = rast[ix, :]
				
	
	return rast2, ustimparams


def zscore(psth, base_range):
	base_psth = psth[base_range[0]:base_range[1]]
	base_mean = base_psth.mean()
	zpsth = psth - base_mean

def calc_on_off(psth, stim_on = 50, max_duration = 50, pplot = False):
	'''
	Calculates psth onset and offset defined as:
		onset : the first millisecond after stimulus onset in which the firing rate is above the half-max value
		offset : the first millisecond after the response peak time in which the firing rate is below the half-max value (or provided maxmimum offset)
	'''
	
	base_mean = psth[:stim_on].mean()
	on_psth = psth[stim_on:stim_on+25] # psth for first 25 ms of response
	resp_peak = on_psth.max() # highest value in on_psth
	resp_peaktime = on_psth.argmax() # peak time in on_psth
	thresh = ((resp_peak-base_mean) / 2.) + base_mean
	
	# calculate response ON
	onset = np.float((psth[stim_on:] > thresh).nonzero()[0][0])
	
	# calculate response OFF
	offset = np.float((psth[stim_on+resp_peaktime:] < thresh).nonzero()[0][0] + stim_on + resp_peaktime)
	offset = np.min([offset, onset + max_duration]) # check if calculated offset is greater than maximum offset
	
	if pplot:
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.plot(psth)
		misc.hline([ax], [base_mean, thresh], color = 'g')
		misc.vline([ax], np.array([onset, offset]) + stim_on, color = 'r', ls = '--')

	return onset, offset
	
	

def pca(spkwav, pcs = None, npcs = 2, visible = False, ax = None):
	
	# z-score
	x_mean = np.tile(spkwav.mean(1), (spkwav.shape[1], 1)).T
	x_std = np.tile(spkwav.std(1), (spkwav.shape[1], 1)).T
	
	x_z = (spkwav - x_mean) / x_std
	
	# if principal components are not provided, calculate PCs on the data set
	if pcs is None:
		C = np.cov(x_z)
		w, v = np.linalg.eig(C)

		pcs = v[:, :npcs]

	x_pc = np.dot(x_z.T, pcs)
	
	if visible:
		if ax is None:
			fig = plt.figure();
			ax1 = fig.add_subplot(211);
			ax1.plot(x_pc[:, 0], x_pc[:, 1], '.', ms = 3)
			ax2 = fig.add_subplot(212);
			ax2.plot(misc.plot_spread_y(pcs))
	
	return x_pc, pcs

def plot_spktimes(spktimes, yloc = 0, ax = None, **kwargs):
	if ax is None:
		fig = plt.figure();
		ax = fig.add_subplot(111);
	
	ax.plot(spktimes, np.ones(spktimes.size)*yloc, linestyle = '', marker = '|', markersize = 40, **kwargs)
	

def plot_sorted_raster(rast, stimparams, lineheight = None, ax = None):
	
	if ax is None:
		fig = plt.figure();
		ax = fig.add_subplot(111);
	
	sort_map = stimparams[:, 0].argsort()
	rast_sort = rast[sort_map]
	
	if not lineheight is None:
		ulevels = np.unique(stimparams[:, 0])
		nlevels = ulevels.size
		levelsize = [(stimparams[:, 0]==i).sum() for i in ulevels]
		lineheight = np.array(levelsize).cumsum()[:-1]-0.5

	plot_raster(rast_sort, lineheight = lineheight, ax = ax)
	
	return ax

def plot_raster(rast, lineheight = None, ax = None):
	
	ax, _ = misc.axis_check(ax)

	ntrials, nbins = rast.shape
	spk = rast.nonzero()
	spktrials = spk[0]
	spktimes = spk[1]/1000.
	
	ax.plot(spktimes, spktrials, marker = '|', linestyle = '', ms = 1)

	if not lineheight is None:
		[ax.axhline(i, color = 'r', lw = 0.5) for i in lineheight]

	ax.set_xlim([0, nbins/1000.])
	ax.set_ylim([0, ntrials])
	
	return ax


# def rast2spktimes(rast, ix_times):
# 	
# 	psth = rast.sum(0)
# 	spktimes = np.empty(psth.sum())
# 	i = 0
# 	for ix_time, bin in zip(ix_times, psth):
# 		spktimes[i:i+bin] = ix_time
# 		i += bin
# 	
# 	return spktimes

def make_spktimes_matrix(rast_shape):
	
	ntrials, npts = rast_shape
	spktimes_mat = np.tile(np.arange(0, npts), (ntrials, 1))
	return spktimes_mat

def rast2spktimes(rast):
	
	spktimes_mat = make_spktimes_matrix(rast.shape)
	# ok a little weird, but you have to add 0.1 ms to avoid edge effects when you bin on millisecond bin edges, just trust me
	spktimes = 0.0001+spktimes_mat[rast.astype(np.bool)] / 1000.
	return spktimes
