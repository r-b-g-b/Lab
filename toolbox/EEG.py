import numpy as np
import os, glob, h5py, itertools
from matplotlib import pyplot as plt
from scipy.signal import butter, filtfilt
import Spikes, RF
import misc
from scipy.io import loadmat
from matplotlib.mlab import specgram

studydir = '/Volumes/BOB_SAGET/TNFalpha/tinnitus/'

freqbands = dict(delta=[1,4], theta=[4,8], beta=[13,30], gamma=[30,70])

def fileconvert_all(experiment = 'unilateralstim', epoch = 'Pre'):
	experiments = ['awakeeeg', 'unilateralstim']
	epochs = ['Pre', 'Post']
	for epoch in epochs:
		for experiment in experiments:
			sesspaths = glob.glob(os.path.join(studydir, experiment, 'Sessions', epoch, '*'))
			for sesspath in sesspaths:

				convertpath = os.path.join(sesspath, 'fileconversion')
				if not os.path.exists(convertpath):
					print 'Making fileconversion folder: %s' % convertpath
					os.mkdir(convertpath)

				matpaths = glob.glob(os.path.join(sesspath, 'data', '[A-Za-z]*.mat'))	
				for matpath in matpaths:
					fileconvert(matpath)

def fileconvert(matpath):

	outpath = matpath.replace('data', 'fileconversion').replace('.mat', '.h5')
	if not os.path.exists(outpath):
		print matpath
		tmpfile = loadmat(matpath)
		Data0 = tmpfile['Data']
		u_ = h5py.File(outpath, 'w')
		
		lfp, stimID = export_unit(Data0)
		# save out to file
		u_.create_dataset('lfp', data = lfp, compression = 'gzip')
		u_.create_dataset('stimID', data = stimID)
		u_.close()


def export_unit(Data0):
	
	# number of time bins to include in the LFP array
	nlfpsamp = 0
	for tt, trial in enumerate(Data0['trial'][0][0][0]):

		thislfpsamp = trial['LFP'].shape[1]
		if thislfpsamp>nlfpsamp:
			nlfpsamp = thislfpsamp
	
	ntrials = Data0['trial'][0][0][0].size # find number of trials
	nstimID = Data0['trial'][0][0][0][0]['Epoch_Value'][0].size
	
	# initialze LFP
	lfp = np.ndarray((2, ntrials, nlfpsamp), dtype = 'float32')
	# initialize frequency and attenuation IDs
	stimID = np.ndarray((ntrials, nstimID), dtype = 'float32')

	for tt in range(ntrials):
		trial = Data0['trial'][0][0][0][tt]

		thisstimID = np.float32(trial['Epoch_Value'][0])
		# get the LFP for this trial and pad it with nans so it can fit in a matrix (since some of the trials have +/-1 data point for LFP)
		for cc in range(2):
			lfpchannel = trial['LFP'][cc]
			lfp[cc, tt, :len(lfpchannel)] = lfpchannel
			
		# add to Epoch_Value
		stimID[tt, :] = thisstimID

	remID = np.array([0., 0.])
	trial_mask = RF.make_trial_mask(stimID, remID)
	lfp = lfp[:, ~trial_mask, :]
	stimID = stimID[~trial_mask, :]

	return lfp, stimID

def add_psd(experiment = 'awakeeeg', epoch = 'Pre'):
	duration = 6.
	fpaths = glob.glob(os.path.join(studydir, experiment, 'Sessions', epoch,  'fileconversion', '*.h5'))
	for fpath in fpaths:
		print fpath
		f = h5py.File(fpath)
		S = []
		nchan, ntrials, nsamp = f['lfp'].shape
		Fs = nsamp / duration
		for j in range(nchan):
			S.append([])
			for k in range(ntrials):
				F, S_, Serr = multi_taper_psd(f['lfp'][j, k, :], Fs = Fs)
				S[-1].append(S_)
		S = np.asarray(S)
		f.create_dataset('psd', data = S, compression = 'gzip')
		f.create_dataset('F', data = F)
		f.close()

def add_spectrogram(experiment = 'awakeeeg', epoch = 'Pre'):
	
	duration = 6.
	fpaths = glob.glob(os.path.join(studydir, experiment, 'Sessions', epoch,  'fileconversion', '*.h5'))
	for fpath in fpaths:
		print fpath
		f = h5py.File(fpath)
		S = []
		nchan, ntrials, nsamp = f['lfp'].shape
		Fs = nsamp / duration
		for j in range(nchan):
			S.append([])
			for k in range(ntrials):
				S_, F, T = specgram(f['lfp'][j, k, :], Fs = Fs, NFFT = 128, noverlap = 64)
				S[-1].append(S_)
		S = np.asarray(S)
		g = f.create_group('spectrogram')
		g.create_dataset('S', data = S, compression = 'gzip')
		g.create_dataset('F', data = F)
		g.create_dataset('T', data = T)
		f.close()

def remove_60hz(lfp, fs = 384.384384384):

	fs_nyq = fs / 2.

	b, a = butter(5, [59./fs_nyq, 61./fs_nyq], btype = 'bandstop', output = 'ba')
	lfp_filt = filtfilt(b, a, lfp)

	return lfp_filt


def calc_power_spectrum(lfp):
	
	npts, ntrials = lfp.shape
	Fs = 384.
	X = np.zeros()
	x = np.zeros((129, 7, ntrials))
	for i in range(ntrials):
		x[:, :, i] = specgram(lfp[:127, i], Fs = Fs)[0]

def calc_lfp_by_stim(rast, stimparams):

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
	'''
	dur_ms = rast.shape[1] # number of milliseconds
	t_ms = np.arange(dur_ms) # time indices in ms
	nbins = rast.shape[-1]
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
		lfp_ = rast[ix, :]
		psth[tuple(n)] = lfp_.sum(0)
				
	return psth, usp

def zero_60Hz(df, fieldname = 'FFTabs'):

	ix = np.logical_and(fft_freq>58, fft_freq<62)
	for i, ds in df.iterrows():
		ds['FFTabs'][ix] = 0.

	return df

def calc_freqband_power(df, freqband, bandname, fft_freq):

	pwr = []
	ix = np.logical_and(fft_freq>freqband[0], fft_freq<freqband[1])
	for i, ds in df.iterrows():
		pwr_ = ds['FFTabs'][ix].sum()
		pwr.append(pwr_)

	df[bandname] = pwr

	return df

def calc_freqband_power_all(df, freqbands, fft_freq):

	for (bandname, freqband) in freqbands.iteritems():
		df = calc_freqband_power(df, freqband, bandname, fft_freq)

	return df