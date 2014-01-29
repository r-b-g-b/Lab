import os, glob, h5py
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import EEG, Spikes, voc
import misc

studydir = '/Volumes/BOB_SAGET/TNFalpha/tinnitus/awakeeeg'
stim_duration = 2. # 3. before 20130606
trial_duration = 3. # 4.5 before 20130606
subjIDs = ['Red102', 'Green102', 'Blue102', 'Red101']
epochs = ['Pre', 'Post']
hemis = 'rl'
stimtypes = ['noise', 'silent']
lfp_fs = 384.384384384

f = h5py.File('/Volumes/BOB_SAGET/TNFalpha/tinnitus/awakeeeg/Sessions/Pre/tnfaKO_awakeeeg_2013/fileconversion/Blue102_b04_silent.h5')
fft_freq = f['fft_freq'].value
f.close()

tmp = np.loadtxt('/Volumes/BOB_SAGET/TNFalpha/tinnitus/lesion_notes.txt', 'S').tolist()
lesionear = dict()
for (k, v) in tmp:
	lesionear[k] = v

def make_single_trial_sheets():

	outdir = os.path.join(studydir, 'Sheets', 'Singletrial')
	for subjID in subjIDs:
		d = pd.HDFStore(os.path.join(studydir, 'Combined', '%s_silent.h5' % subjID))
		df = d['df'][::10]
		d.close()

		fig = plt.figure();
		ax_lfp = fig.add_subplot(211);
		ax_fft = fig.add_subplot(212);
		l1, = ax_lfp.plot(np.zeros(len(df.lfp[0])))
		l3, = ax_lfp.plot(np.zeros(len(df.lfpfilt[0])), color = 'r')
		l2, = ax_fft.plot(fft_freq, np.zeros(len(df.FFT[0])))
		ax_lfp.axhline(0.0005, color = 'r', ls = '--')
		ax_lfp.axhline(-0.0005, color = 'r', ls = '--')
		ax_lfp.set_ylim([-0.001, 0.001])
		ax_lfp.set_xlim([0, len(df.lfp[0])])
		ax_fft.set_ylim([0, 0.04])

		df['FFTabs'] = np.abs(df['FFT'])
		
		for i, (ix, ds) in enumerate(df.iterrows()):
			l1.set_ydata(ds['lfpfilt'])
			l3.set_ydata(ds['lfp'])
			l2.set_ydata(ds['FFTabs'])

			fig.savefig(os.path.join(outdir, '%s_trial%u.png' % (ds['sess'], i)))

		plt.close(fig)

def pre_post_freqratio_comparison():

	for subjID in subjIDs:
		d = pd.HDFStore(os.path.join(studydir, 'Combined', '%s_silent.h5' % subjID))
		df = d['df']
		d.close()

		df = mark_good(df)
		df = df[df.good]
		df['FFTabs'] = np.abs(df.FFT)

		lowhighbands = dict(low=[1,8], high=[30,70])
		df = EEG.calc_freqband_power_all(df, lowhighbands, fft_freq)

		df['highlow'] = df.high / df.low
		epochgp = df.groupby(['hemi', 'epoch'])
		epochmean = epochgp.highlow.apply(np.mean)
		epocherr = epochgp.highlow.apply(np.std)

		colors = dict(Pre = 'b', Post = 'r')
		axs = dict(l = 1, r = 2)
		fig = plt.figure()
		for hemi in hemis:
			ax = fig.add_subplot(1, 2, axs[hemi])
			ax.set_title('%s hemisphere' % hemi)
			ax.bar([0, 1], epochmean[hemi][::-1], yerr = epocherr[hemi])
			ax.set_xticks([0.4, 1.4])
			ax.set_xticklabels(epochmean[hemi].index[::-1])

		fig.suptitle('%s (%s ear lesion)' % (subjID, lesionear[subjID]))
		fig.savefig(os.path.join(studydir, 'Analysis', '%s_freqratio.png' % subjID))
		fig.suptitle('')
		fig.clf();



def pd_std(ds):
	'''
	Calculates 
	'''
	n = len(ds)
	nsamp = len(ds[ds.index[0]])
	ds_mean = np.mean(ds)
	runningsum = np.zeros(nsamp)
	for i in ds:
		runningsum = runningsum + (i-ds_mean)**2
	std = np.sqrt(runningsum / (n-1))
	return std

def pd_sem(ds):

	n = len(ds)
	nsamp = len(ds[ds.index[0]])
	ds_mean = np.mean(ds)
	runningsum = np.zeros(nsamp)
	for i in ds:
		runningsum = runningsum + (i-ds_mean)**2
	std = np.sqrt(runningsum / (n-1))
	sem = std / np.sqrt(n)
	return sem


def post_lesion_silent_comparison2():

	fig = plt.figure()
	for subjID in subjIDs:
		d = pd.HDFStore(os.path.join(studydir, 'Combined', '%s_silent.h5' % subjID))
		df = d['df']
		d.close()

		df = df[df.good]

		epochgp = df.groupby(['epoch', 'hemi'])
		df['FFTabs'] = np.abs(df.FFT)
		epochmean = epochgp.FFTabs.apply(np.mean)
		epocherr = epochgp.FFTabs.apply(pd_std)

		colors = dict(Pre = 'b', Post = 'r')
		axs = dict(l = 1, r = 2)

		for (k, v) in epochmean.iterkv():
			ax = fig.add_subplot(1, 2, axs[k[1]])
			ax.set_title('%s hemisphere' % k[1])
			misc.errorfill(fft_freq, np.log10(v), epocherr[k], color = colors[k[0]], ax = ax)

		fig.suptitle('%s (%s ear lesion)' % (subjID, lesionear[subjID]))
		fig.savefig(os.path.join(studydir, 'Analysis', '%s_silent_pre_post.png' % subjID))
		fig.suptitle('')
		fig.clf()

def post_lesion_silent_comparison():


	f = h5py.File('/Volumes/BOB_SAGET/TNFalpha/tinnitus/awakeeeg/Sessions/Pre/tnfaKO_awakeeeg_2013/fileconversion/Blue102_b04_silent.h5')
	fft_freq = f['fft_freq'].value
	f.close()

	tmp = np.loadtxt('/Volumes/BOB_SAGET/TNFalpha/tinnitus/lesion_notes.txt', 'S').tolist()
	lesionear = dict()
	for (k, v) in tmp:
		lesionear[k] = v

	fig = plt.figure()
	for subjID in subjIDs:
		d = pd.HDFStore(os.path.join(studydir, 'Combined', '%s_silent.h5' % subjID))
		df = d['df']
		d.close()

		df = df[df.good]

		df['FFTabs'] = np.abs(df.FFT)
		sessgp = df.groupby(['sess', 'hemi'])
		epochmean = sessgp.agg(dict(epoch = misc.get_first, hemi = misc.get_first))
		fftmean = sessgp.FFTabs.apply(np.mean)

		fftmean = pd.DataFrame(dict(FFTabs = fftmean), index = fftmean.index)
	
		sessdf = pd.concat([fftmean, epochmean], axis = 1)

		epochgp = sessdf.groupby(['hemi', 'epoch'])
		epochdf = epochgp.FFTabs.apply(np.mean)
		epochdf = pd.DataFrame(dict(FFTabs = epochdf), index = epochdf.index)

		colors = dict(Pre = 'b', Post = 'r')
		axs = dict(l = 1, r = 2)

		for (k, v) in epochdf.iterrows():
			ax = fig.add_subplot(1, 2, axs[k[0]])
			ax.set_title('%s hemisphere' % k[0])
			ax.plot(fft_freq, np.log10(v['FFTabs']), color = colors[k[1]])

		ax = fig.add_subplot(121)
		fig.suptitle('%s (%s ear lesion)' % (subjID, lesionear[subjID]))
		fig.savefig(os.path.join(studydir, 'Analysis', 'silent_pre_post_%s.png' % subjID))
		fig.suptitle('')
		fig.clf()

def make_silent_sheets(epoch = 'Pre'):

	sesspaths = glob.glob(os.path.join(studydir, 'Sessions', epoch, '*'))
	for sesspath in sesspaths:
		fpaths = glob.glob(os.path.join(sesspath, 'fileconversion', '*.h5'))

		for fpath in fpaths:

			absol, relat = os.path.split(fpath)
			fname, _ = os.path.splitext(relat)
			
			f = pd.HDFStore(fpath, 'r')
			df = f['df']
			f.close()

			df = df[df.rr==-1]
			df = df[df.good]
			gp = df.groupby('hemi')

			nsamp = df.lfp.values[0].size
			Fs = nsamp / trial_duration

			lfp_mean = gp.lfp.apply(np.mean)
			lfp_err = gp.lfp.apply(np.std)

			fft_mean = gp.fft.apply(np.mean)
			fft_err = gp.fft.apply(np.std)

			t = np.linspace(0, trial_duration, nsamp)
			f = np.linspace(0, Fs/2., df.fft.values[0].size)

			fig = plt.figure()
			ax_lfp = []; ax_fft = []
			for i, ((k, lfp_mean_), (k, lfp_err_), (k, fft_mean_), (k, fft_err_)) in enumerate(zip(lfp_mean.iterkv(), lfp_err.iterkv(),\
				fft_mean.iterkv(), fft_err.iterkv())):
				# lfp_err_ = lfp_err[k]
				ax_lfp.append(fig.add_subplot(2, 2, i+1));
				misc.errorfill(t, lfp_mean_, lfp_err_, ax = ax_lfp[-1])
				ax_lfp[-1].set_title(k)

				ax_fft.append(fig.add_subplot(2, 2, i+3));
				misc.errorfill(f, np.log(fft_mean_), np.log(fft_err_), ax = ax_fft[-1])

			misc.sameyaxis(ax_lfp); misc.sameyaxis(ax_fft)
			misc.samexaxis(ax_lfp, [0, stim_duration]); misc.samexaxis(ax_fft, [f.min(), f.max()])

			figpath = os.path.join(studydir, 'Sheets', '%s_silent.png' % fname)
			fig.savefig(figpath)
			plt.close(fig)

def make_noise_rr_sheets():

	fpaths = glob.glob(os.path.join(studydir, 'Sessions', 'Pre', 'both',  '*.h5'))

	for fpath in fpaths:

		absol, relat = os.path.split(fpath)
		fname, _ = os.path.splitext(relat)

		f = pd.HDFStore(fpath, 'r')
		df = f['df']
		f.close()
		df = df[df.rr>0] # subset only non-silent trials
		df = df[df.good] # use only good trials
		urrs = np.unique(df.rr)
		nrrs = urrs.size

		gp = df.groupby(('rr', 'hemi'))

		nsamp = df.lfp.values[0].size
		Fs = nsamp / trial_duration

		lfp_mean = gp.lfp.apply(np.mean) 
		lfp_err = gp.lfp.apply(np.std)

		fft_mean = gp.fft.apply(np.mean)
		fft_err = gp.fft.apply(np.std)

		t = np.linspace(0, trial_duration, nsamp)
		f = np.linspace(0, Fs/2., df.fft.values[0].size)

		fig = plt.figure(figsize = (14, 8.8))
		ax_lfp = []; ax_fft = []
		for i, ((k, lfp_mean_), (k, lfp_err_), (k, fft_mean_), (k, fft_err_)) in enumerate(zip(lfp_mean.iterkv(), lfp_err.iterkv(),\
			fft_mean.iterkv(), fft_err.iterkv())):

			rr, hemi = k
			if hemi=='l':
				j=1
			else: j=2

			ax_lfp.append(fig.add_subplot(nrrs, 4, ((i/2)*4)+j))
			misc.errorfill(t, lfp_mean_, lfp_err_, ax = ax_lfp[-1])
			ax_lfp[-1].set_title(k)

			ax_fft.append(fig.add_subplot(nrrs, 4, ((i/2)*4)+2+j))
			misc.errorfill(f, np.abs(fft_mean_), np.abs(fft_err_), ax = ax_fft[-1])


		misc.sameyaxis(ax_lfp); misc.sameyaxis(ax_fft)
		misc.samexaxis(ax_lfp, [0, stim_duration]); misc.samexaxis(ax_fft, [0, 100.])
		[a.set_yticklabels('') for a in ax_lfp[1::2]]
		fig.savefig(os.path.join(studydir, 'Sheets', '%s_noise.png'%fname))
		plt.close(fig)

def make_noise_singleburst_sheets():

	fpaths = glob.glob(os.path.join(studydir, 'Sessions', 'Pre', 'both',  '*.h5'))

	fig = plt.figure(figsize = (14, 8.8))
	ax_lfp = []

	for j, fpath in enumerate(fpaths):

		absol, relat = os.path.split(fpath)
		fname, _ = os.path.splitext(relat)

		f = pd.HDFStore(fpath, 'r')
		df = f['df']
		f.close()
		df = df[np.logical_and(df.rr>0, df.rr<16, df.good)] # subset only non-silent trials

		gp = df.groupby('hemi')

		nsamp = df.lfp.values[0].size
		Fs = nsamp / trial_duration

		lfp_mean = gp.lfp.apply(np.mean) 
		lfp_err = gp.lfp.apply(np.std)

		fft_mean = gp.fft.apply(np.mean)
		fft_err = gp.fft.apply(np.std)

		t = np.linspace(0, trial_duration, nsamp)
		f = np.linspace(0, Fs/2., df.fft.values[0].size)

		for i, ((k, lfp_mean_), (k, lfp_err_), (k, fft_mean_), (k, fft_err_)) in enumerate(zip(lfp_mean.iterkv(), lfp_err.iterkv(),\
			fft_mean.iterkv(), fft_err.iterkv())):

			ax_lfp.append(fig.add_subplot(len(fpaths), 2, (j*2)+i+1))
			misc.errorfill(t, lfp_mean_, lfp_err_, ax = ax_lfp[-1])
			if i==0:
				ax_lfp[-1].set_title('%s / %s' % (fname, k))
			else:
				ax_lfp[-1].set_title(k)

	misc.sameyaxis(ax_lfp, [-0.0004, 0.0004])
	misc.samexaxis(ax_lfp, [0.05, 0.175])
	[a.set_yticklabels('') for a in ax_lfp[1::2]]
	fig.savefig(os.path.join(studydir, 'Sheets', 'noise_singleburst.png'))
	plt.close(fig)


def add_lfp_filts():

	for epoch in epochs:
		sesspaths = glob.glob(os.path.join(studydir, 'Sessions', epoch, '*'))
		for sesspath in sesspaths:
			fpaths = glob.glob(os.path.join(sesspath, 'fileconversion', '*.h5'))
			for fpath in fpaths:
				f = pd.HDFStore(fpath, 'a')
				if not 'lfp_filt' in f.keys():
					add_lfp_filt(f)
				f.close()

def add_lfp_filts():
	add_field_from_fcn(add_lfp_filt, 'lfp_filt')

def add_ffts():
	add_field_from_fcn(add_fft, 'fft')

def add_field_from_fcn(fcn, fieldname, **kwargs):

	for epoch in epochs:
		print epoch
		sesspaths = glob.glob(os.path.join(studydir, 'Sessions', epoch, '*'))
		for sesspath in sesspaths:
			print '\t%s\n' % os.path.split(sesspath)[1]
			fpaths = glob.glob(os.path.join(sesspath, 'fileconversion', '*.h5'))
			for fpath in fpaths:
				f = h5py.Filep(fpath, 'a')
				if not fieldname in f.keys():
					fcn(f, **kwargs)
				f.close()
			print '\n'
		print '\n'

def add_lfp_filt(f, fs = 384.384384384):

	lfp = f['lfp'].value
	nchan, ntrials, nsamp = lfp.shape

	lfp_filt = np.empty_like(lfp)

	for i in range(nchan):
		for j in range(ntrials):
			lfp_filt[i, j, :] = EEG.remove_60hz(lfp[i, j, :])

	f.create_dataset('lfp_filt', data = lfp_filt, compression = 'gzip')


def add_fft(f, fs = 384.384384384):
	
	lfp = f['lfp'].value
	nchan, ntrials, nsamp = lfp.shape

	# calcualte FFT
	FFT = []
	for i in range(nchan):
		FFT.append(np.fft.rfft(lfp[i, ...]))

	# calculate fft frequencies
	nfreq = FFT[0].shape[1]
	fft_freq = np.fft.fftfreq(nsamp, d = 1./fs)	
	fft_freq = fft_freq[:nfreq]

	FFT = np.asarray(FFT)

	f.create_dataset('fft', data = FFT, compression = 'gzip')
	f.create_dataset('fft_freq', data = fft_freq)

	return f

def plot_signal_quality():
	'''
	Plot "signal quality" measure for each animal.
	This is simply the histogram of maximum trial amplitudes.
	Movement artifacts are usually large amplitude (>0.002) spikes in the signal.
	'''
	fpaths = glob.glob(os.path.join(studydir, 'Sessions', 'Pre', 'both', '*.h5'))

	axs = []
	fig = plt.figure()
	npaths = len(fpaths)
	nrows = np.ceil(np.sqrt(npaths))
	for i, fpath in enumerate(fpaths):
		d = pd.HDFStore(fpath, 'r')
		df = d['df']
		d.close()
		axs.append(fig.add_subplot(nrows, nrows, i+1))
		axs[-1].hist(df.lfp.apply(np.max), bins = np.arange(0, 0.01, 0.0005))
		axs[-1].set_title(os.path.split(fpath)[-1])

	misc.samexaxis(axs)
	plt.show()
	fig.savefig(os.path.join(studydir, 'Analysis', 'signal_quality.png'))


def mark_good_trials(epoch = 'Pre'):

	subjpaths = glob.glob(os.path.join(studydir, 'Combined', '*'))
	for subjpath in subjpaths:
		
		d = pd.HDFStore(subjpath)

		df = d['df']
		df = mark_good(df)
		d['df'] = df

		d.close()

def mark_good(df):
	
	df['lfp_max'] = df.lfp.apply(np.abs)
	df['lfp_max'] = df.lfp_max.apply(np.max)
	df['good'] = df.lfp_max<0.0005

	return df

def plot_trials(df):

	plt.close('all');
	fig1 = plt.figure();
	fig2 = plt.figure();
	ax1 = []; ax2 = []
	hemis = ['r', 'l']
	for j, hemi in enumerate(hemis[0]):
		df_ = df[df.hemi==hemi]
		ntrials = df_.lfp.size
		nsamp = df_.lfp[df_.index[0]].size
		for k, i in enumerate(range(ntrials)[::40]):
			ax1.append(fig1.add_subplot(5, 5, k+1))
			ix = df_.index[i]
			ax1[-1].hist(df_.lfp[ix], bins = 100)
			ax1[-1].set_xticklabels('')
			ax1[-1].set_yticklabels('')
			ax1[-1].set_title(ix)

			ax2.append(fig2.add_subplot(5, 5, k+1))
			ix = df_.index[i]
			ax2[-1].plot(df_.lfp[ix])
			ax2[-1].set_xticklabels('')
			ax2[-1].set_yticklabels('')

	misc.samexaxis(ax1)
	misc.sameyaxis(ax2)

def convert_to_pd():

	'''
	Loop order:
		Subject, stimulus type, (initialize file), epoch, session
	'''

	for subjID in subjIDs:
		for stimtype in stimtypes: # loop through stimtypes
			LFP = []; LFPFILT = []; FFT = []; STIMID = []; HEMI = []; SUBJID = []; EPOCH = []; SESSID = [];
			if stimtype == 'noise':
				ATTEN = [];

			for epoch in epochs:

				sesspaths = glob.glob(os.path.join(studydir, 'Sessions', epoch, '*'))
				
				for sesspath in sesspaths:

					fpath, = glob.glob(os.path.join(sesspath, 'fileconversion', '%s*%s.h5' % (subjID, stimtype)))

					absol, relat = os.path.split(fpath)
					blockID, _ = os.path.splitext(relat)
					sessID = os.path.split(os.path.split(absol)[0])[-1]
					subjID, blocknum, stimtype = blockID.split('_')

					fin = h5py.File(fpath, 'r')
					lfp = fin['lfp'].value
					lfp_filt = fin['lfp_filt'].value
					fft = fin['fft'].value
					fft_freq = fin['fft_freq'].value
					stimparams = fin['stimID'].value
					fin.close()

					nchan, ntrials, nsamp = lfp.shape

					for i in xrange(nchan):

						[LFP.append(lfp[i, j, :]) for j in xrange(ntrials)]
						[LFPFILT.append(lfp_filt[i, j, :]) for j in xrange(ntrials)]
						[FFT.append(fft[i, j, :]) for j in xrange(ntrials)]
						if stimtype == 'noise':
							ATTEN.extend(stimparams[:, 1].tolist())
						HEMI.extend([hemis[i]]*ntrials)
						SUBJID.extend([subjID]*ntrials)
						EPOCH.extend([epoch]*ntrials)
						SESSID.extend([sessID]*ntrials)

			if stimtype == 'noise':
				d = dict(lfp = LFP, lfpfilt = LFPFILT, FFT = FFT, atten = ATTEN, hemi = HEMI, subj = SUBJID, epoch = EPOCH, sess = SESSID)
			elif stimtype == 'silent':
				d = dict(lfp = LFP, lfpfilt = LFPFILT, FFT = FFT, hemi = HEMI, subj = SUBJID, epoch = EPOCH, sess = SESSID)
			df = pd.DataFrame(d)

			outpath = os.path.join(studydir, 'Combined', '%s_%s.h5' % (subjID, stimtype))
			fout = pd.HDFStore(outpath)
			fout['df'] = df
			fout.close()

def compare_pre_post_noise_all():
	pass

def compare_pre_post_noise(animalID):

	for epoch in epochs:
		print epoch
		sesspaths = glob.glob(os.path.join(studydir, 'Sessions', epoch, '*'))
		for sesspath in sesspaths:
			fpath, = glob.glob(os.path.join(sesspath, 'both', '%s*noise.h5' % animalID))
			d = pd.HDFStore(fpath)
			df = d['df']
			d.close()




def summarize_sessions():
	for epoch in epochs:
		print epoch
		sesspaths = glob.glob(os.path.join(studydir, 'Sessions', epoch, '*'))
		sesspaths = np.sort(sesspaths)
		for sesspath in sesspaths:
			summarize_session(sesspath)

def summarize_session(sesspath):

	print '\t%s' % os.path.split(sesspath)[1]
	fpaths = glob.glob(os.path.join(sesspath, 'fileconversion', '*noise*.h5'))
	# tmp = [os.path.split(i)[-1] for i in fpaths]
	# for i in tmp: print i

	f = h5py.File(fpaths[0])

	nlfp = f['lfp'].value.shape[2]
	ustimid = np.unique(f['stimID'].value[:, 1])
	print ustimid

	print '\t\t#LFP samples:\t%u' % nlfp
	print '\t\t'
	f.close()

	


'''
def combine_pre_post():


	for subjID in subjIDs:
		for epoch in epochs:
		print subjID, epoch
		sesspaths = os.path.join(studydir, 'Sessions', epoch)
		for sesspath in sesspaths:'''
			
'''
EEG.fileconvert_all(experiment = 'awakeeeg', epoch = 'Pre')
convert_to_pd(epoch = 'Pre')

'''

