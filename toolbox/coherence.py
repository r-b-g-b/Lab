import os, glob, h5py, re
import numpy as np
import pandas as pd
from fileconversion import get_session_unitinfo
import nitime as nt
import nitime.viz as viz
from nitime.analysis import SNRAnalyzer

studydir = '/Volumes/BOB_SAGET/Fmr1_voc'
# studydir = '/Users/robert/Desktop/Fmr1_voc'
sessionpaths = glob.glob(os.path.join(studydir, 'Sessions', 'voc_*'))
sessions = [os.path.split(s)[1] for s in sessionpaths]

def calc_coherence_all(lb = 0, ub = 100, prefix = 'VOC'):

	# dtype = np.dtype([('gen', 'S2'), ('exp', 'S3'), ('sess', 'S20'), ('unit', 'i4'), ('stim', 'i4'), ('Psignal', '1600f4'), ('Pnoise', '1600f4'), ('snr', '1600f4'), ('coh', '1600f4'), ('F', '1600f4'), ('info_snr', 'f4')])
	# coh_info = np.empty(0, dtype = dtype)
	
	npsth = 16000
	Gens = []
	Exps = []
	Sessions = []
	Unitnums = []
	Stims = []
	Info_snrs = []

	for session in sessions:
		print session
		_, gen, exp, _ = session.split('_')
	
		unitinfos = get_session_unitinfo(session, onlycomplete = ('RF', 'RR', 'VOC'))
		# fpaths = glob.glob(os.path.join(studydir, 'Sessions', session, 'fileconversion', '%s*.h5' % prefix))
		nunits = len(unitinfos)

		# loop through penetrations
		for j, unitkey in enumerate(unitinfos.keys()):

			unitinfo = unitinfos[unitkey]
			fpath = unitinfo['fpath'][unitinfo['stimtype'].index('VOC')]
			unitnum = unitinfo['unitnum']
		
			print '%i\t(%i of %i)' % (unitnum, j+1, nunits)

			f = h5py.File(fpath, 'r')
			rast = f['rast'].value
			stimparams = f['stimID'].value
			stimparams = stimparams[:, 0][..., np.newaxis]
	
			ustim = np.unique(stimparams[:, 0])
			nstim = ustim.size
		
			stim = []; info_snr = []
			# loop through stimuli
			for i in range(nstim):

				ix = stimparams[:, 0]==ustim[i]
				rast_ = rast[ix, :npsth]
				t1 = nt.TimeSeries(rast_, sampling_interval = 0.001)
				snr1 = SNRAnalyzer(t1)
				
				F = snr1.mt_frequencies
				freq_ix = F<ub
				F = F[freq_ix]
				info = snr1.mt_information[freq_ix]
				pnoise = snr1.mt_noise_psd[freq_ix]
				psignal = snr1.mt_signal_psd[freq_ix]
				psnr = snr1.mt_snr[freq_ix]
				coh = snr1.mt_coherence[freq_ix]
				fstep = np.diff(F[:2])[0]
				stim.append(ustim[i])
				info_snr.append((-np.log2(1-coh)).sum()*fstep)

								
			Gens.append(gen); Exps.append(exp); Sessions.append(session); Unitnums.append(unitnum)
			Stims.append(stim); Info_snrs.append(info_snr);


				# coh_info.resize(coh_info.size+1)
				# coh_info[-1] = np.array((gen, exp, session, unitnum, ustim[i], psignal, pnoise, psnr, coh, F, info_snr), dtype = dtype)

			df = pd.DataFrame(dict(gen = Gens, exp = Exps, sess = Sessions, unit = Unitnums, stim = Stims, info_snr = Info_snrs))

			d = pd.HDFStore(os.path.join(studydir, 'Analysis', 'coh_info_%s.h5' % session))
			d['df'] = df
			d.close()

	# np.savez(os.path.join(studydir, 'Analysis', 'tmp_coh_info.npz'))
	# return coh_info
	return df
	

def load_coherence():

	fpath = os.path.join(studydir, 'Analysis', 'coh_data.h5')
	if os.path.exists(fpath):
		print 'Loading %s' % os.path.split(fpath)[1]
		d = pd.HDFStore(fpath, 'r')
		df = d['df'].copy()
		d.close()
		return df
	else:
		print 'Loading individual coh_info files'
		dfs = []
		fpaths = glob.glob(os.path.join(studydir, 'Analysis', 'coh_info_*.h5'))
		for fpath in fpaths:
			d = pd.HDFStore(fpath)
			dfs.append(d['df'].copy())
			d.close()

		return pd.concat(dfs)







