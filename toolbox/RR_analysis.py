import os, glob, h5py
import numpy as np
import pandas as pd
import Spikes, RF, RR
import fileconversion
import misc

studydir = '/Users/robert/Desktop/Fmr1_voc'
sessionpaths = glob.glob(os.path.join(studydir, 'Sessions', '[A-Za-z]*'))
sessions = [os.path.split(i)[1] for i in sessionpaths]

def analyze():

	gens = []; exps = []; sesss = []; unitnums = []; cfs = []; noise_rrtfs = []; cf_rrtfs = []

	for session in sessions:

		_, gen, exp, _ = session.split('_')
		unitinfos = fileconversion.get_session_unitinfo(session, onlycomplete = ('RF', 'RR', 'VOC'))

		for unitkey in unitinfos.iterkeys():

			unitinfo = unitinfos[unitkey]
			
			print unitinfo['session'], unitinfo['unitnum']

			rr_ix = unitinfo['stimtype'].index('RR')
			f = h5py.File(unitinfo['fpath'][rr_ix], 'r')
			rast = f['rast'].value
			stimparams = f['stimID'].value
			cf_ix = f['cf'].value
			f.close()

			ufreqs = np.unique(stimparams[:, 0])
			nfreqs = ufreqs.size

			urrs = np.unique(stimparams[:, 1])

			cf = RF.ix2freq[20:][np.round(cf_ix).astype(int)]

			_, ix, err = misc.closest(ufreqs, cf, log = True)

			if np.abs(err)>0.5:
				print 'Repetition rate stimuli not present for this unit (nearest tone %2.2f octaves away)' % err
				continue

			ufreq_played = ufreqs[ix]

			npips = urrs*4
			noise_rrtf = RR.calc_rrtf_all(rast, stimparams, 0., urrs, npips = npips)
			cf_rrtf = RR.calc_rrtf_all(rast, stimparams, ufreq_played, urrs, npips = npips)

			gens.append(gen)
			exps.append(exp)
			sesss.append(session)
			unitnums.append(unitinfo['unitnum'])
			cfs.append(cf_ix)
			noise_rrtfs.append(noise_rrtf)
			cf_rrtfs.append(cf_rrtf)

	df = pd.DataFrame(dict(gen = gens, exp = exps, sess = sesss, unit = unitnums, cf = cfs, noise_rrtf = noise_rrtfs, cf_rrtf = cf_rrtfs))
	d = pd.HDFStore(os.path.join(studydir, 'Analysis', 'rrtf_data.h5'))
	d['df'] = df
	d.close()
	return df
