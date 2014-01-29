import os, glob, h5py, re
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
from scipy.ndimage import gaussian_filter
import RF; reload(RF)
import RR; reload(RR)
import Spikes; reload(Spikes)
import misc; reload(misc)


verbose = True

basedir = '/Volumes/BOB_SAGET/Fmr1_RR/'
ix2freq = 1000 * 2**(np.arange(0, 64)/10.)
lfp_fs = 1603/4.2
lfp_P = 1./lfp_fs
lfp_t = np.arange(0, 4.2, lfp_P)

def add_field():

	savepath = '/Users/robert/Documents/Work/Bao/Fmr1_RR/Analysis/RRTFs.h5'
	fsave = h5py.File(savepath, 'a')

	sessnames = fsave.keys()

	for sessname in sessnames:
	
		fsess = fsave[sessname]
		rrblocks = fsess.keys()

		for rrblock in rrblocks:
	
			rrpath = os.path.join(basedir, 'Sessions', 'full_window', 'good', sessname, 'fileconversion', rrblock + '.h5')
			p = re.compile('RR')
			rfpath_ = p.sub('RF', rrblock)
			rfpath = os.path.join(basedir, 'Sessions', 'full_window', 'good', sessname, 'fileconversion', rfpath_+'.h5')
		
			p = re.compile('(\d+)')
			unitnum = np.int32(p.findall(rrblock)[0])

			funit = fsess[rrblock]
			# d = funit.require_dataset('unit', (), 'i4')
			# d.write_direct(np.array([unitnum]))

			if verbose:
				print rrblock
				print rfpath
			
			# LOAD RF FILECONVERSION FILE
			frf = h5py.File(rfpath, 'r')
			cf_ix = np.int32(frf['cf'].value)
			frf.close()

			cf = ix2freq[20:][cf_ix]
			
			# LOAD RR FILECONVERSION FILE
			frr = h5py.File(rrpath, 'r')
			rast = frr['rast'].value
			stimparams = frr['stimID'].value
			frr.close()
			
			# CALCULATE FREQ PLAYED TO THIS UNIT'S CF
			freqs = np.unique(stimparams[:, 0])
			rrs = np.unique(stimparams[:, 1])
			freq, _, freq_err = misc.closest(freqs, cf, log = True)
			freq = np.round(freq, 0)
			assert freq_err<0.3
			
			# COMPUTE NEW FIELD
			on, off = Spikes.calc_on_off(funit['ev_psth'].value)
			# vs, vs_p = RR.calc_vs_all(rast, stimparams, [0.], rrs)
			# rrtf = RR.calc_rrtf_all(rast, stimparams, freq, rrs)

			# ADD NEW FIELD
			d = funit.require_dataset('on_halfmax', (1,), float)
			d.write_direct(np.array([on]))
			d = funit.require_dataset('off_halfmax', (1,), float)
			d.write_direct(np.array([off]))
		
			# open the RR file, get rast and stimparams, then close it
			# frf = h5py.File(rfpath, 'r')
			# RF.add_rf_analysis(frf, funit)
			# frf.close()

	fsave.close()




def h5_to_numpy(f, fields, types):
	
	dtype = np.dtype(zip(fields+['gen', 'exp', 'sess', 'unit', 'cf'], types+['S2', 'S3', 'S20', 'S20', 'f4']))
	db = np.empty(0, dtype = dtype)
	for sess in f:
		fsess = f[sess]
		gen, exp, _ = sess.split('_')
		for unit in fsess:
			db_ = np.empty(1, dtype = dtype)
			funit = fsess[unit]
			db_['gen'] = gen
			db_['exp'] = exp
			db_['sess'] = sess
			db_['unit'] = unit
			db_['cf'] = funit['cf_com'].value
			for field in fields:
				db_[field] = funit[field].value
			
			db.resize(db.size+1)
			db[-1] = db_
			
	return db
	
def gen_exp_analysis_2d(f, key = 'rrtf_noise', key_dtype = '7f4'):
	
	db = h5_to_numpy(f, [key, 'rrs'], [key_dtype, '7f4'])

	fig = plt.figure()
	ax = fig.add_subplot(111)	
	gens = ['ko', 'wt']
	exps = ['nai', 'exp']
	pltopts = {'ko' : {'nai' : {'color' : 'r', 'ls' : '-'}, 'exp' : {'color' : 'r', 'ls' : '--'}}, 'wt' : {'nai' : {'color' : 'b', 'ls' : '-'}, 'exp' : {'color' : 'b', 'ls' : '--'}}}
	# t = np.arange(7)
	rrs = db['rrs'][0]
	leg = []
	for i, gen in enumerate(gens):
		for j, exp in enumerate(exps):
			ix = np.vstack((db['gen']==gen, db['exp']==exp)).all(0)
			db_ = db[ix]
			nunits = db_.size
			# y = st.nanmean(db_[key], 0)
			yerr = st.nanstd(db_[key], 0) / np.sqrt(nunits)
			ax.errorbar(rrs, y, yerr = yerr, color = pltopts[gen][exp]['color'], ls = pltopts[gen][exp]['ls'])
			
			leg.append('-'.join((gen, exp)))
			
	ax.legend(leg)
	# ax.set_title('Evoked PSTHs')
	# ax.set_xlabel('Time (ms)')
	# ax.set_ylabel('Firing rate (spks/s)')
	plt.show()




def gen_exp_analysis_1d(f, key = 'on_halfmax', key_dtype = 'f4'):
	db = h5_to_numpy(f, [key, 'rrs'], [key_dtype, '7f4'])
	fig = plt.figure()
	ax = fig.add_subplot(111)	
	gens = ['wt', 'ko']
	exps = ['nai', 'exp']
	pltopts = {'ko' : {'nai' : {'color' : 'r', 'hatch' : None}, 'exp' : {'color' : 'r', 'hatch' : '/'}}, 'wt' : {'nai' : {'color' : 'b', 'hatch' : None}, 'exp' : {'color' : 'b', 'hatch' : '/'}}}
	# t = np.arange(7)
	leg = []
	for i, gen in enumerate(gens):
		for j, exp in enumerate(exps):
			ix = np.vstack((db['gen']==gen, db['exp']==exp)).all(0)
			db_ = db[ix]
			nunits = db_.size
			# y = st.nanmean(data[ix], 0)
			# yerr = st.nanstd(data[ix], 0) / np.sqrt(nunits)
			y = st.nanmean(db_[key], 0)
			yerr = st.nanstd(db_[key], 0) / np.sqrt(nunits)
			ax.bar((i*3)+j, y, yerr = yerr, color = pltopts[gen][exp]['color'], ecolor = 'k', hatch = pltopts[gen][exp]['hatch'])
			
	# 		leg.append('-'.join((gen, exp)))
	# 		
	# ax.legend(leg)
	# # ax.set_title('Evoked PSTHs')
	# # ax.set_xlabel('Time (ms)')
	# # ax.set_ylabel('Firing rate (spks/s)')
	# plt.show()

	
			
	
	