import numpy as np
from scipy.stats import circmean, circstd
import matplotlib.pyplot as plt
import glob
import os
import RF; reload(RF)
import Spikes; reload(Spikes)
import misc; reload(misc)
from itertools import combinations
import matplotlib.mpl as mpl

def calc_rrtf_all(rast, stimparams, freq, urrs, npips = 6, onset = 0.05, norm = True):
	'''
	Takes the full block raster and stimparams, finds the response to the first tone by filtering the rast to include only responses to the given frequency (but for all repetition rates), then it filters the rast to only the given freq/rr pair and passes that	to calc_rrtf, along with the first tone response for all repetition rates

	Input:
		rast : full block raster
		stimparams : full block stimulus parameters
		freq : in Hz, the frequency played to this unit
		urrs : sorted list (lo-hi) of repetition rates played
		npips : list of the number of pips at each rate, or scalar if each rate had the same number of pips
		onset : in seconds, the onset time for the first pip
		norm : for the response to each pip, subtract the pre-pip response

	Output:
		rrtf : n-length vector, where n is the number of repetition rates
	'''
	if type(npips) is int:
		npips = [npips]*len(urrs)

	nrrs = urrs.size
	rrtf = np.empty(nrrs) * np.nan

	# get raster subset for this frequency
	ix = stimparams[:, 0] == freq

	# response onset and offset (5 - 25 ms after stimulus onset)
	resp_start = np.int32((0.005 + onset) * 1000)
	resp_end = resp_start + 20
	nspks_1st = (rast[ix, resp_start:resp_end]).mean() # spikes per millisecond

	# normalize by pre-pip baseline (-20 - 0 ms before stimulus onset)
	if norm:
		pip_end_pre = resp_start - 5
		pip_start_pre = pip_end_pre - 20
		nspks_1st = nspks_1st - (rast[ix, pip_start_pre: pip_end_pre]).mean()

	# loop through repetition rates, get raster subset, and calculate RRTF
	for i, (rr, npip) in enumerate(zip(urrs, npips)):
		rast_ = rast[RF.get_trials(stimparams, np.array([freq, rr])), :]
		rrtf[i] = calc_rrtf(rast_, rr, nspks_1st, npips = npip, norm = norm)
	
	return rrtf
	
def calc_rrtf(rast_, rr, nspks_1st, npips = 6., onset = 0.05, norm = True):
	'''
	rast_
	rr
	nspks_1st
	npips
	onset (s)
	norm
	If norm is True, than we make an attempt to calculate the following of a neuron while controlling for baseline increases / non-phase locked increases. Simply considers the "response" to a tone as the post-tone firing rate minus the pre-tone firing rate.		
	'''

	
	pip_start = (0.005 + onset + np.arange(npips) / rr) * 1000
	pip_end = pip_start + 20

	if norm:
		pip_end_pre = pip_start - 5
		pip_start_pre = pip_end_pre - 20

	nspks = np.empty(npips-1) * np.nan
	for i in np.arange(1, npips):
		nspks[i-1] = (rast_[:, pip_start[i] : pip_end[i]]).mean()
		
		if norm:
			nspks_pre = (rast_[:, pip_start_pre[i] : pip_end_pre[i]]).mean()
			nspks[i-1] = nspks[i-1] - nspks_pre
			
	return nspks.mean() / nspks_1st

def calc_rrtf_lfp_all(lfp, lfp_t, stimparams, freq, rrs, onset = 0.05):

	nrrs = rrs.size
	rrtf_lfp = np.empty(nrrs) * np.nan
	trial_ix = stimparams[:, 0] == freq
	pip_start = 0.005 + onset
	pip_end = pip_start + 0.02
	time_ix = np.vstack((pip_start<lfp_t, lfp_t<pip_end)).all(0)

	lfp_mag_1st = (lfp[trial_ix, :][:, time_ix]).mean(0).min()
	print lfp_mag_1st
	for i, rr in enumerate(rrs):
		lfp_ = lfp[RF.get_trials(stimparams, np.array([freq, rr])), :]
		rrtf_lfp[i] = calc_rrtf_lfp(lfp_, lfp_t, rr, lfp_mag_1st)
	
	return rrtf_lfp
	
def calc_rrtf_lfp(lfp_, lfp_t, rr, lfp_mag_1st, npips = 6., onset = 0.05):
	
	pip_start = (0.005 + onset + np.arange(npips) / rr)
	pip_end = pip_start + 0.020

	lfp_mag = np.empty(npips-1) * np.nan
	for i in np.arange(1, npips):
		ix = np.vstack((pip_start[i]<lfp_t, lfp_t<pip_end[i])).all(0)
		lfp_mag[i-1] = (lfp_[:, ix]).mean(0).min()
					
	return lfp_mag.mean() / lfp_mag_1st
	
def calc_fft_following(psth_smoo, rr, npips, remfirst = True):

	if remfirst:
		pip_start = 1.
	else:
		pip_start = 0.
		
	trial_start = np.int32(1000 * pip_start/rr)
	trial_end = np.int32(1000 * npips/rr)
	
	psth_smoo = psth_smoo[trial_start:trial_end]

	P = np.abs(np.fft.fftshift(np.fft.fft(psth_smoo)))
	P_freq = np.fft.fftshift(np.fft.fftfreq(psth_smoo.size, d = 0.001))
	[_, ix] = misc.closest(P_freq, rr)
	power = P[ix[0]] / psth_smoo.sum()

	return power

def calc_vR_dist_intertrial(spktimes, spktrials, trials, rr, nbins, npips = 6.):

	ntrials = trials.size
	dist = []
	# compare all pairwise combinations of trials
	for m in combinations(range(ntrials), 2):
		spktimes_1 = spktimes[spktrials == trials[m[0]]]
		spktimes_2 = spktimes[spktrials == trials[m[1]]]
	
		# remove response to first pip
		train_start = 0.05 + 1./rr
		train_end = 0.05 + 6./rr
		spktimes_1 = spktimes_1[misc.isbetween(spktimes_1, train_start, train_end)]
		spktimes_2 = spktimes_2[misc.isbetween(spktimes_2, train_start, train_end)]
		duration = npips / rr
		dist.append(Spikes.vR_dist(spktimes_1, spktimes_2, nbins) / duration)
	
	return np.mean(dist)

def calc_vR_dist_intertrain(spktimes_, spktrials_, trials, rr, nbins, npips = 6):
	'''
	Calculates the reliability of firing to each PIP in the train
	Input:
		spktimes_	:	spike times for the trials involved
		spktrials_	:	spike trials for the above spike times
		stimparams	:	stimulus parameters for the trials
	'''

	dist = []
	# loop through each trial
	for t, trial in enumerate(trials):
		spktimes__ = spktimes_[spktrials_ == trial]
		# compare all pairwise combinations of pips
		for m in combinations(range(1, npips), 2):
			train_start_1 = 0.05 + (m[0] / rr)
			train_end_1 = 0.05 + ((m[0]+1) / rr)
			train_start_2 = 0.05 + (m[1] / rr)
			train_end_2 = 0.05 + ((m[1]+1) / rr)
			duration = npips / rr
			spktimes_1 = spktimes__[misc.isbetween(spktimes__, train_start_1, train_end_1)] - train_start_1
			spktimes_2 = spktimes__[misc.isbetween(spktimes__, train_start_2, train_end_2)] - train_start_2
			dist.append(Spikes.vR_dist(spktimes_1, spktimes_2, nbins) / duration)
	
	return np.mean(dist)

def calc_r_from_psth(psth, rr, npips, onset, base_mean = None, remove_first = False):

	onset_bin = np.int32(onset*1000)
	train_start = onset_bin
	if remove_first:
		train_start = train_start + 1./rr
	train_end = onset_bin + np.int32(1000 * np.float32(npips) / rr)	
	psth_ = psth[train_start : train_end]
	
	if base_mean is not None:
		psth_ = psth_ - base_mean

	alpha = np.arange(0, psth_.size, 1)
	alpha = alpha/1000. * 2*np.pi * rr
	
	r = np.sum(psth_ * np.exp(1j * alpha))
	
	if base_mean is not None:
		r = np.abs(r) / np.sum(base_mean)
	else:
		r = np.abs(r) / np.sum(psth_)
	
	return r

def calc_vs_psth_indiv(spktimes_, rr, npips, onset):
	
	vs = np.zeros(npips)
	# spktimes_ = spktimes_.round(3)*1000
	pipstarts = np.arange(npips, dtype = 'float32')/rr
	cyclelength = 1./rr
	for i in range(npips):
		spktime_start = onset+pipstarts[i]
		train_end = onset+pipstarts[i]+cyclelength
		spktimes_chunk = spktimes_[np.logical_and(spktime_start<spktimes_, spktimes_<train_end)]

		alpha = np.exp(1j * 2*np.pi * rr * spktimes_chunk)
		# np.cos(2 * np.pi * rr * spktimes_chunk) + 1j * np.sin(2 * np.pi * rr * spktimes_chunk)
		vs[i] = np.abs(alpha.mean())

	return vs

def calc_vs_all(rast, stimparams, ufreqs, urrs, npips = 6, onset = 0.05):
	
	ufreqs = np.asarray(ufreqs)
	urrs = np.asarray(urrs)
	nfreqs = ufreqs.size
	nrrs = urrs.size
	vs = np.empty((nfreqs, nrrs))*np.nan
	vs_p = np.empty((nfreqs, nrrs))*np.nan
	for f in range(nfreqs):
		for r in range(nrrs):
			ix = RF.get_trials(stimparams, [ufreqs[f], urrs[r]])
			rast_ = rast[ix, :]
			vs[f, r], vs_p[f, r] = calc_vs(rast_, urrs[r], npips, onset)
		
	return vs, vs_p

def calc_vs(rast_, rr, npips, onset = 0.05, analysistype = 'mean', remove_first = False):

	ix_times = np.arange(rast_.shape[1]) / 1000.
	spktimes_ = Spikes.rast2spktimes(rast_, ix_times)
	spktimes_ = spktimes_ - onset
	train_end = npips / rr # end at last pip
	if remove_first:
		train_start = 1./rr
	else:
		train_start = 0.
	spktimes_chunk = spktimes_[np.logical_and(train_start<spktimes_, spktimes_<train_end)]

	v = np.exp(1j * 2*np.pi * rr * spktimes_chunk)
	# v = np.cos(2 * np.pi * rr * spktimes_chunk) + 1j * np.sin(2 * np.pi * rr * spktimes_chunk)

	if analysistype == 'mean':
		vs = np.abs(v.mean()) #/ duration
	if analysistype == 'sum':
		vs = np.abs(v.sum())

	alpha = np.angle(v)
	p = rayleigh_test(np.angle(v))
	
	return vs, p

def calc_vs_indiv(spktimes_, rr, npips, onset):
	
	vs = np.zeros(npips)
	# spktimes_ = spktimes_.round(3)*1000
	pipstarts = np.arange(npips, dtype = 'float32')/rr
	cyclelength = 1./rr
	for i in range(npips):
		spktime_start = onset+pipstarts[i]
		train_end = onset+pipstarts[i]+cyclelength
		spktimes_chunk = spktimes_[np.logical_and(spktime_start<spktimes_, spktimes_<train_end)]

		alpha = np.cos(2 * np.pi * rr * spktimes_chunk) + 1j * np.sin(2 * np.pi * rr * spktimes_chunk)
		vs[i] = np.abs(alpha.mean())

	return vs

def circ_psth_all(rast, stimparams, freq, npips, onset = 0.05, bins = 20, color = 'b', remove_first = False, axs = None):
	'''
	Input:
	Output:
		r :  the mean vector length for each repetition rate
		V : the summed vector length for each repetition rate
		theta : the mean vector angle for each repetition rate
	'''
	urrs = np.unique(stimparams[:, 1])
	nrrs = urrs.size
	ix = RF.get_trials(stimparams, (freq, np.nan))
	rast_ = rast[ix, :]
	stimparams_ = stimparams[ix, :]
	r = []; V = []; theta = []
	
	for i in xrange(nrrs):
		ix = RF.get_trials(stimparams_, (np.nan, urrs[i]))
		r_, V_, theta_ = circ_psth(rast_[ix, :], urrs[i], npips[i], onset = onset, bins = bins, color = color, remove_first = remove_first, ax = axs[i])
		r.append(r_); V.append(V_); theta.append(theta_)
	
	misc.sameyaxis(axs)

	return np.array(r), np.array(V), np.array(theta)

def circ_psth(rast_, rr, npips, onset = 0.05, bins = 20, color = 'b', remove_first = False, ax = None):
	'''
	Input:
		rast_ : subset of block raster corresponding to only the stimuli you want plotted
		rr : the repetition rate
		npips : number of pips presented at this rate
		onset : in seconds, onset time of first pip
		bins : number of bins to include
		remove_first : exclude the respose to the first pip
	'''
	
	spktimes_ = Spikes.rast2spktimes(rast_)
	period = 1. / rr
	train_end = float(npips) / rr
	spktimes_chunk = spktimes_ - onset
	spktimes_chunk = spktimes_chunk[spktimes_chunk > 0]
	spktimes_chunk = spktimes_chunk[spktimes_chunk < train_end]
	if remove_first:
		spktimes_chunk = spktimes_chunk[spktimes_chunk > period]
	spktimes_pol = spktimes_chunk * 2*np.pi * rr # convert to radians
	spktimes_pol = np.mod(spktimes_pol, 2*np.pi) # force range 0-2pi
	nspks = spktimes_pol.size

	r, V, theta = circ_stat(spktimes_pol)
	
	ax, fig = misc.axis_check(ax, polar = True)
	if spktimes_pol.size > 0:
		ax.hist(spktimes_pol, bins = bins, color = color, histtype = 'stepfilled')
		ax.plot([0, theta], [0, V], 'r.-')
		ax.plot([0, 0], [0, V*np.cos(theta)], 'r.-')

	return r, V, theta

def circ_stat(alpha):
	'''
	Takes alpha, angles in radians
	Returns:
		r : the mean vector length (0-1)
		theta : the mean vector angle
		V : the summed vector length
	'''
	n = alpha.size
	z = np.exp(1j*alpha)
	z_ = z.mean() # mean vector
	r = np.abs(z_) # mean vector strength
	V = np.abs(z.sum()) # summed vector strength
	theta = np.angle(z_) # mean vector angle

	return r, V, theta

def rayleigh_test(alpha):

	n = alpha.size
	r, _, _ = circ_stat(alpha)
	# compute Rayleigh's R (equ. 27.1)
	R = n*r

	# compute Rayleigh's z (equ. 27.2)
	z = R**2 / n;

	# compute p value using approxation in Zar, p. 617
	p = np.exp(np.sqrt(1+4*n+4*(n**2-R**2))-(1+2*n));

	return p

def aligned_psth_separate_all(rast, stimparams, freq, npips, onset = 0.05, axs = None):
	'''
	Input:
		rast : full block raster
		stimparams : full block stimulus parameters
		freq : in Hz, the frequency played to this neuron
		npips : number of pips for each repetition rate
		onset : in ms, onset time of first pip
		axs : a list of axis to which the output will be displayed
	'''
	urrs = np.unique(stimparams[:, 1])
	nrrs = urrs.size
	if axs is None:
		axs = [None]*nrrs

	aligned_psths = []
	for rr, npip, ax in zip(urrs, npips, axs):
		ix = RF.get_trials(stimparams, (freq, rr))
		psth = 1000*rast[ix, :].mean(0)
		aligned_psths.append(aligned_psth_separate(psth, rr, npip, onset = onset, ax = ax))

	return aligned_psths

def aligned_psth_separate(psth, rr, npips, onset = 0.050, wind_offset = -0.020, ax = None, nshow = 6):

	ccycle = plt.rcParams['axes.color_cycle']

	period = 1./rr # period in ms
	nbins_period = np.round(1000*period).astype(int)

	wind_start = onset+wind_offset

	show_ix = np.round(np.linspace(0, npips-1, nshow)).astype(int)
	aligned_psth = np.empty((nshow, nbins_period))
	for i, ix in enumerate(show_ix):
		start_ix = np.round(1000*(wind_start+ix*period)).astype(int)
		stop_ix = start_ix+nbins_period
		aligned_psth[i, :] = psth[start_ix : stop_ix]

	if not ax is None:
		for i in xrange(nshow):
			color = ccycle[(nshow-(i+1))]
			lw = 7.25-(1.25*i)
			alpha = (2/7.)+(i/7.)
			ax.plot(Spikes.exp_smoo(aligned_psth[-(i+1), :]), color = color, lw = lw, alpha = alpha)

		ax.set_xlim((1000*(wind_start-0.010), 150))

	return aligned_psth

def aligned_psth(spktimes_, rr, npips, onset, visible = False):
	
	if spktimes_.size > 0:
		cyclelength = 1./rr

		spktimes_shift = spktimes_ - onset
		spktimes_shift = spktimes_shift[spktimes_shift > 0]
		spktimes_shift = np.round(spktimes_shift * 1000)
		spktimes_aligned = np.mod(spktimes_shift, cyclelength*1000)
		spktimes_aligned = spktimes_aligned / 1000.

		x = ax.hist(spktimes_aligned, cyclelength * 1000, visible = visible)
		apsth = x[0]

	else:
		apsth = np.nan
	
	return apsth

def plot_rrtf(t, psth, rr, npips, onset = 0.05, duration = 0.025, ax = None):

	ax, fig = misc.axis_check(ax)

	ax.plot(t, Spikes.exp_smoo(psth))
	plot_tone_pips(rr, npips, onset = onset, duration = duration, ax = ax, color = 'r')


def tone_times(rr, npips, onset):
	'''
	Calcualtes tone onset times given a repetition rate and number of pips
	Input:
		rr : the repetition rate
		npips : the number of pips
		onset : first pip onset time in seconds
	'''
	return np.arange(npips, dtype = 'float32')/rr + onset
	
def plot_tone_pips(rr, npips, onset, duration, ax = None, **kwargs):
	
	if ax is None:
		ax = plt.gca()
	
	ylim = ax.get_ylim()
	yrange = ylim[1] - ylim[0]
	rect_height = 0.1 * yrange
	
	piptimes = tone_times(rr, npips, onset)
	for i in range(npips):
		rect = mpl.patches.Rectangle((piptimes[i], 0), duration, rect_height, zorder = 10, **kwargs)
		ax.add_patch(rect)
		# ax.plot(np.array([piptimes[i], piptimes[i]+duration]), [0, 0], lw = 6, **kwargs)

def plot_stim_psth(stim_psth, stim_params):
	
	fig = plt.figure();
	nfreqs, nrrs, nbins = stim_psth.shape
	i = 0
	rrs = stim_params[1]
	for r in xrange(nrrs):
		for f in xrange(nfreqs):
			i += 1
			ax = fig.add_subplot(nrrs, nfreqs, i)
			ax.plot(Spikes.exp_smoo(stim_psth[f, r, :]))
			plot_tone_pips(rrs[r], rrs[r], 0.05, 0.025)

	plt.show();