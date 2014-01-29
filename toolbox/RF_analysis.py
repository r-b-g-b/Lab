import numpy as np
import h5py, glob, os, shutil
import Spikes; reload(Spikes)
import RF; reload(RF)
import RR; reload(RR)
import matplotlib.pyplot as plt
import misc; reload(misc);
import scipy.stats as st
import re
import shutil
import itertools

sesss = ['ko_w2_20120802', 'wt_exp_20120314', 'wt_exp_20120320', 'wt_exp_20120508', 
'wt_exp_20120509', 'wt_nai_20120208', 'wt_nai_20120410', 'ko_exp_20120315', 'wt_nai_20120412', 
'ko_exp_20120321', 'wt_w1_20120625', 'ko_exp_20120326', 'wt_w1_20120626', 'ko_exp_20120420', 
'wt_w1_20120628', 'ko_nai_20120202', 'wt_w1_20120629', 'ko_nai_20120206', 'wt_w2_20120619', 
'ko_nai_20120209', 'wt_w2_20120620', 'ko_nai_20120305', 'wt_w2_20120621', 
'ko_nai_20120306', 'wt_w2_20120622', 'ko_w2_20120521', 'wt_w3_20120709', 'ko_w2_20120523']


basedir = '/Volumes/BOB_SAGET/'
nbins = 333
stim_on = 50
ix2freq = 1000 * 2**(np.arange(0, 64)/10.)


# dtype = np.dtype([('group', 'S2'), ('sess', 'S20'), ('unit', 'i8'), ('cf', 'f4'), \
# 	('vs_noise', 'f4'), ('vs_tone', 'f4'), ('p_noise', 'f4'), ('p_tone', 'f4'), \
# 	('dist_il_noise', 'f4'), ('dist_il_tone', 'f4'), \
# 	('dist_in_noise', 'f4'), ('dist_in_tone', 'f4'), \
# 	('power_noise', 'f4'), ('power_tone', 'f4'), \
# 	('rrtf_noise', 'f4'), ('rrtf_tone', 'f4'), \
# 	('rr', 'f4'), ('base_mean', 'f4'), ('evoked_mean', 'f4')])
	
dtype = np.dtype([('gen', 'S2'), ('exp', 'S3'), ('sess', 'S20'), ('unit', 'i8'), \
	('psth', '333f4'), ('ev_psth', '333f4'), \
	('rf', '(8, 43)f4'), ('rf_clean', '(8, 43)f4'), \
	('cf_man', 'f4'), ('com', 'f4'), ('com_tip', 'f4'), \
	('bw', '8f4'), ('bw_lr', '(8, 2)f4'), ('thresh', 'f4'), ('coord', '2f4'), \
	('resp_on', 'f4'), ('resp_off', 'f4'), ('ev_resp_on', 'f4'), ('ev_resp_off', 'f4'), \
	('base_mean', 'f4'), ('ev_mean', 'f4')])
	
p = re.compile('(\d+).h5')

def get_ix2freq():
	return ix2freq

def characterize(sesss = sesss, experiment = 'Fmr1_RR', pplot = True, verbose = False):

	if type(sesss) == str:
		sesss = [sesss]

	# set up figure
	figsize = (12, 12)
	fig = plt.figure(figsize = figsize)

	# loop through sesss
	for sess in sesss:

		DB = np.empty(0, dtype = dtype)
		
		print '%s\n%s\n\n' % (sess, '-'*50)
		
		# build the output directory path
		savedir = os.path.join(basedir, experiment, 'Sessions', sess, 'analysis')
		if not os.path.exists(savedir):
			os.mkdir(savedir)
			
		# WT or KO / CTL or EXP
		gen, exp, date = sess.split('_')

		# find the RF blocks
		pens = glob.glob(os.path.join(basedir, experiment, 'Sessions', sess, 'fileconversion', 'RF*.h5'))

		# load the cfs for this sess
		cfs = np.loadtxt(os.path.join(basedir, experiment, 'Sessions', sess, 'cfs.txt'), ndmin = 1)
		
		# loop through blocks in this sess
		for pen in pens:
			
			absol, relat = os.path.split(pen)
			blockname = os.path.splitext(relat)[0]

			# get unit number from filename
			unitnum = np.int32(p.findall(relat))[0] # unit number
			ix = cfs[:, 0] == unitnum
			if ix.sum() > 0:
				cf_man = cfs[ix, 1][0]
				if verbose:
					print pen
			
				# load the RF block
				f = h5py.File(pen, 'r')
				spktimes = f['spktimes'].value
				# if not np.isnan(spktimes[0]):

				'''--------RF--------'''
				# load the RF block to get the RF
				rf = f['rf'].value;	rast = f['rast'].value; spktimes = f['spktimes'].value; stimparams = f['stimID'].value; spktrials = f['spktrials'].value; coord = f['coord'].value; ntrials = f['rast'].shape[0]; f.close()

				# calculate the psth (spk/s*trial, normalized by number of trials
				psth = Spikes.calc_psth(rast, normed = True) # spk/s
			
				# baseline firing rate
				base_mean = psth[:stim_on].mean() # spk/s
			
				# response onset/offset
				psth_smoo = Spikes.exp_smoo(psth, tau = 0.003) 
				resp_on, resp_off = Spikes.calc_on_off(psth_smoo, stim_on = stim_on)
			
				# rewindowed RF
				rf_rewin = RF.calc_rf(rast, stimparams, resp_on = resp_on + stim_on - 3, resp_off = resp_off + stim_on + 3, normed = True)
			
				# thresholded RF
				rf_thresh = rf_rewin.copy()
				rf_threshold = np.percentile(rf_thresh, 66)# upper quartile
				rf_peak = rf_rewin.max()
				rf_thresh[rf_thresh < rf_threshold] = 0
			
				# find maximum RF cluster
				(rf_clust, clust_sizes) = RF.findmaxcluster(rf_thresh, cf = cf_man, include_diagonal = False)
				# if clust_sizes.max() < 10: # if it's a tiny RF, set it to nans
				# 	rf_clust = np.empty(rf_clust.shape) * np.nan
			
				rf_mask = rf_clust > 0
				# find evoked psth
				ev_psth = RF.calc_evoked_psth(rast, stimparams, rf_mask)
				ev_psth_smoo = Spikes.exp_smoo(ev_psth, tau = 0.003)
				ev_resp_on, ev_resp_off = Spikes.calc_on_off(ev_psth_smoo, stim_on = stim_on)
				ev_mean = ev_psth[ev_resp_on : ev_resp_off].mean()
				
				# bandwidth and threshold
				bw, bw_lr, _, thresh = RF.calc_bw_cf_thresh(rf_mask)

				# center of mass
				com = RF.calc_rf_com(rf_clust)
				tip_top = np.max([thresh-2, 0])
				tip_bottom = thresh
				com_tip = RF.calc_rf_com(rf_clust[tip_top:tip_bottom, :])

				'''PLOT'''
				if pplot:
					
					rf1_ax = fig.add_subplot(221)
					rf2_ax = fig.add_subplot(222)
					psth1_ax = fig.add_subplot(223)
					psth2_ax = fig.add_subplot(224)
					
					RF.plot_RF(rf, bw_lr = bw_lr, thresh = thresh, cf = cf_man, ax = rf1_ax)
					RF.plot_RF(rf_clust, bw_lr = bw_lr, thresh = thresh, cf = cf_man, ax = rf2_ax)
					rf2_ax.axvline(com, color = 'g', ls = '--')
					rf2_ax.axvline(com_tip, color = 'b', ls = '--')
				
					psth1_ax.plot(psth_smoo)
					psth2_ax.plot(ev_psth_smoo)
					psth1_ax.axvline(resp_on+stim_on, color = 'r', ls = '--')
					psth1_ax.axvline(resp_off+stim_on, color = 'r', ls = '--')
				
					figpath = os.path.join(savedir, blockname + '.png')
					fig.savefig(figpath);
					fig.clf()
				'''PLOT'''
			

				DB.resize(DB.size + 1)
				DB[-1] = np.array((gen, exp, sess, unitnum, \
					psth[:333], ev_psth[:333], \
					rf, rf_clust, \
					cf_man, com, com_tip, \
					bw, bw_lr, thresh, coord, \
					resp_on, resp_off, ev_resp_on, ev_resp_off, \
					base_mean, ev_mean), dtype = dtype)
		
			# end unit loop
	
			np.savez(os.path.join(basedir, experiment, 'Sessions', sess, sess + '_RF.npz'), DB = DB)
		if verbose:
			print '\n'*4
		# end sess loop		
	
	
def load_DBs(groups, experiment = 'Fmr1_RR', suffix = 'DB', DB = None, v = False):

	if DB is None:
		DB = np.empty(0, dtype = dtype)

	if groups == 'all': # uses shell-style wildcards
		groups = ['[!_]']

	pat = '_'.join(groups) + '*' + '/'
	dirs = glob.glob(os.path.join(basedir, experiment, 'Sessions', pat))
	for d in dirs:
		_, sess = os.path.split(os.path.split(d)[0])
		try:
			fname = os.path.join(d, sess + '_' + suffix + '.npz')
			db = np.load(fname)['DB']
			for db_ in db:
				DB.resize(DB.size + 1)
				DB[-1] = db_
			nunits = db.size
		except IOError:
			nunits = 0
		if v:
			print '%s\t%i' % (fname, nunits)
	return DB

def export_DB(DB, experiment = 'Fmr1_RR', fname = 'DB'):

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

	np.savetxt(os.path.join(basedir, experiment, 'Analysis', fname+'.txt'), txt, fmt = '%s')