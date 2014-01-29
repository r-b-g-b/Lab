import os, glob, h5py, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import Spikes, RF, RR, misc
reload(RR);
import fileconversion

# studydir = '/Volumes/BOB_SAGET/Fmr1_voc'
studydir = '/Users/robert/Desktop/Fmr1_voc'
sessionpaths = glob.glob(os.path.join(studydir, 'Sessions', '[A-Za-z]*'))

ix2freq = 1000 * 2**(np.arange(0, 64)/10.)
nrrs = 7
t_rrtf = np.arange(0, 6, 0.001)
t_rf = np.arange(0, 0.333, 0.001)


def rr_make_contactsheets():
	'''
	loop through all the sessions and plot the rrtfs
	'''

	fig = plt.figure(figsize = (30, 18));
	txt_suptitle = fig.suptitle('')
	ax_cfrrtf = fig.add_axes((0.76, 0.76, 0.24, 0.23));
	ax_cfvs = ax_cfrrtf.twinx();
	ax_cfcircpsthall = fig.add_axes((0.62, (11/14.)-0.02, 0.1, (1/7.)+0.04), polar = True)
	ax_cfcircpsthall.set_xticklabels(''); ax_cfcircpsthall.set_yticklabels('');
	ax_rf = fig.add_axes((0.67, 0.51, 0.33, 0.23));
	ax_rfrast = fig.add_axes((0.67, 0.25, 0.33, 0.24));
	ax_rfrast.set_xticklabels('');
	ax_rfpsth = fig.add_axes((0.67, 0.01, 0.33, 0.24));

	ax_cfrr = [fig.add_axes((0.03, 1-((i+1)/7.), 0.35, 1/7.)) for i in np.arange(nrrs)]
	ax_cfalignedpsth = [fig.add_axes((0.38, 1-((i+1)/7.), 0.17, 1/7.)) for i in np.arange(nrrs)]
	ax_cfcircpsth = [fig.add_axes((0.53, 1-((i+1)/7.), 0.1, 1/7.), polar = True) for i in np.arange(nrrs)]
	# ax_noiserr = [fig.add_subplot(nrrs, 3, i) for i in np.arange(1, 3*nrrs, 3)]

	for sessionpath in sessionpaths:

		session = os.path.split(sessionpath)[1]
		unitinfos = fileconversion.get_session_unitinfo(sessionpath, onlycomplete = ('RF', 'RR', 'VOC'))
			
		for unitkey in unitinfos.keys():
			
			txt_suptitle.set_text('%s %s' % (session, unitkey))

			unitinfo = unitinfos[unitkey]

			rf_ix = unitinfo['stimtype'].index('RF')
			
			f_rf = h5py.File(unitinfo['fpath'][rf_ix], 'r')
			rf_rast = f_rf['rast'].value
			rf_stimparams = f_rf['stimID'].value
			cf_ix = f_rf['cf'].value
			f_rf.close()
			
			cf = ix2freq[20:][int(cf_ix)]

			''' calculate and plot RF, psth, and sorted raster'''
			rf = RF.calc_rf(rf_rast, rf_stimparams)
			rf_psth = Spikes.calc_psth(rf_rast)
			RF.plot_rf(rf, cf = cf_ix, axes_on = False, ax = ax_rf) # plot RF
			ax_rf.axvline(cf_ix, color = 'r', lw = 1.5)
			Spikes.plot_sorted_raster(rf_rast, rf_stimparams, ax = ax_rfrast) # plot raster
			ax_rfpsth.plot(t_rf, Spikes.exp_smoo(rf_psth, tau = 0.005)) # plot PSTH

			''' calcualte and plot RRTFs for CF and noise stimuli '''
			rr_ix = unitinfo['stimtype'].index('RR')
			
			f_rr = h5py.File(unitinfo['fpath'][rr_ix], 'r')
			rr_rast = f_rr['rast'].value
			rr_stimparams = f_rr['stimID'].value
			f_rr.close()

			# find the played CF
			rr_ufreqs = np.unique(rr_stimparams[:, 0])
			urrs = np.unique(rr_stimparams[:, 1])
			npips = (urrs*4).astype(int)
			rr_freq, rr_ufreq_ix, _ = misc.closest(rr_ufreqs, cf, log = True)

			ax_rf.axvline(RF.calc_freq2ix(rr_freq), color = 'g', lw = 1.5)
			# calculate the PSTHs for each repetition rate
			tmp = Spikes.calc_psth_by_stim(rr_rast, rr_stimparams)
			rr_cfpth = tmp[0][rr_ufreq_ix, :, :]
			# rrtf_noisepsth = tmp[0][0, :, :]

			# plot the aligned psths
			RR.aligned_psth_separate_all(rr_rast, rr_stimparams, rr_freq, npips, axs = ax_cfalignedpsth)
			[a.set_yticklabels('') for a in ax_cfalignedpsth]
			[a.set_xticklabels('') for a in ax_cfalignedpsth[:-1]]

			# plot circular psths
			r, V, theta = RR.circ_psth_all(rr_rast, rr_stimparams, rr_freq, npips, axs = ax_cfcircpsth)
			[a.set_yticklabels('') for a in ax_cfcircpsth]
			[a.set_xticklabels('') for a in ax_cfcircpsth]

			# plot all circular summed vector strengths
			ax_cfcircpsthall.plot(theta, V, '.-')
			[ax_cfcircpsthall.plot([0, th], [0, v], color = 'b', alpha = 1-(i/10.)) for i, (th, v) in enumerate(zip(theta, V))]


			# plot RRTF
			rrtf = RR.calc_rrtf_all(rr_rast, rr_stimparams, rr_freq, urrs, npips)
			ax_cfrrtf.plot(rrtf, '.-', ms = 10)
			ax_cfvs.plot(V*np.cos(theta), 'g.-', ms = 10)
			for tick in ax_cfvs.yaxis.get_major_ticks():
				tick.set_pad(-5)
				tick.label2.set_horizontalalignment('right')

			# plot repetition rate PSTHs
			for i in xrange(nrrs):
				# RR.plot_rrtf(t_rrtf, rrtf_noisepsth[i, :], urrs[i], int(4*urrs[i]), onset = 0.05, duration = 0.025, ax = ax_noiserr[i])
				RR.plot_rrtf(t_rrtf, rr_cfpth[i, :], urrs[i], int(4*urrs[i]), onset = 0.05, duration = 0.025, ax = ax_cfrr[i])

			# ax_noiserr[0].set_title('Noise RRTFs')
			ax_cfrr[0].set_title('CF RRTFs (%.0f kHz)' % (cf/1000))
			# [a.set_xlim(0, 4.5) for a in ax_noiserr]
			[a.set_xlim(0, 4.5) for a in ax_cfrr]
			misc.sameyaxis(ax_cfrr+ax_cfalignedpsth)

			figsavepath = os.path.join(studydir, 'Sheets', 'RRTFs', '%s_%s_RRTF.png' % (session, unitkey))
			print figsavepath
			fig.savefig(figsavepath)
			[a.cla() for a in fig.get_axes()] # clear all axes


def voc_make_contactsheets():

	tmp = np.load(os.path.join(studydir, 'stims', 'voc_43to80kHz.npz'))
	stims = tmp['P']; T = tmp['T']; F = tmp['F']

	nbins, nfreqs, nstims = stims.shape
	plot_chunk_size = nbins/3

	fig = plt.figure()
	ax_stim = [fig.add_subplot(4, 1, i+1) for i in xrange(nstims)]
	for i in xrange(nstims):
		ix = np.arange(i*plot_chunk_size, ((i+1)*plot_chunk_size))
		voc.plot_spec(stims[ix, :, 0].T**0.5, F = F, T = T[ix], ax = ax_stim[i])
	

		
		



		
		