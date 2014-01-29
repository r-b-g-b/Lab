import os, glob, h5py, re
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.stats import linregress
import Spikes, RF, RR
import misc

# studydir = '/Volumes/BOB_SAGET/Fmr1_RR'
# studydir = '/Volumes/BOB_SAGET/Fmr1_voc/'
studydir = '/Users/robert/Desktop/Fmr1_voc/'

ix2freq = 1000 * 2**(np.arange(0, 64)/10.)

cmap = plt.get_cmap('jet')
def cf_to_rgba(i, N = 40):
	if i>N:
		return cmap(1.)
	elif i<0:
		return cmap(0.)
	else:
		return cmap(1.*i/N)
	

def get_ix2freq():
	return ix2freq

def calc_freq2ix(freq, start_ix = 20):

	_, ix, err = misc.closest(ix2freq[start_ix:], freq, log = True)
	if err>1:
		raise 'Large error encountered. (%f)' % err

	return ix

def calc_rf(rast, stimparams, resp_on = 57, resp_off = 80, normed = False, smooth = False):
	'''
	Output:
		rf : spks per second for each freq/atten combo
	'''
	
	
	freqs = stimparams[:, 0]
	attens = stimparams[:, 1]

	ufreqs = np.unique(freqs)
	uattens = np.unique(attens)

	nfreqs = ufreqs.size
	nattens = uattens.size
	
	rf = np.zeros((nattens, nfreqs))
	count = np.zeros((nattens, nfreqs))
	for f in range(nfreqs):
		ix1 = freqs == ufreqs[f] # trial where the frequency was this frequency
		for a in range(nattens):
			ix2 = attens == uattens[a] # trial where the attenuation was this attenuation
			ix = np.logical_and(ix1, ix2) # trial where both were true
			rf[a, f] = rast[ix, resp_on:resp_off].sum()
			count[a, f] = ix.sum()

	if normed:

		duration = resp_off - resp_on
		rf = (1000 / duration) * (rf / count)
	
	if smooth:
		rf = gaussian_filter(rf, 0.5)
	
	return rf

def calc_lfp_rf(lfp, stimparams, onset = 0.055, offset = 0.08):
	
	time_ix = np.linspace(0, 0.333, lfp.shape[1])
	onset_ix = (time_ix < onset).nonzero()[0][-1]
	offset_ix = (time_ix > offset).nonzero()[0][0]
	onset = np.int32(onset*1000)
	offset = np.int32(offset*1000)
	
	freqs = stimparams[:, 0]
	attens = stimparams[:, 1]

	ufreqs = np.unique(freqs)
	uattens = np.unique(attens)

	nfreqs = ufreqs.size
	nattens = uattens.size
	
	rf = np.zeros((nattens, nfreqs))
	for f in range(nfreqs):
		ix1 = freqs == ufreqs[f] # trial where the frequency was this frequency
		for a in range(nattens):
			ix2 = attens == uattens[a] # trial where the attenuation was this attenuation
			ix = np.logical_and(ix1, ix2) # trial where both were true
			rf[a, f] = np.nanmax(lfp[ix, onset_ix:offset_ix]).mean()

	return rf

def add_bf_man(experiment):
	'''
	Runs through all of the units and allows you to click on the CF. Clicking a negative x-coordinate
	results in the unit being discarded (renamed to _RR###.h5)
	For valid units, the CF and threshold are saved in a text file called cfs.txt
	'''

	rf_blocks = glob.glob(os.path.join(studydir, 'Sessions', experiment, 'fileconversion', '*RF*.h5'))
	cfs = []
	badunits = []
	for rf_block in rf_blocks:
				
		fig = plt.figure();
		ax = fig.add_subplot(111);
		
		print rf_block
		f = h5py.File(rf_block, 'r')

		blockname = os.path.splitext(os.path.split(rf_block)[1])[0]
		pennum = np.int32(blockname.split('RF')[1])
	
		rf = f['rf'].value
		f.close()
		
		RF.plot_rf(rf, ax = ax)
		xlim = ax.get_xlim()
		ylim = ax.get_ylim()
		ax.set_xlim([xlim[0]-2, xlim[1]])
		ax.set_ylim([ylim[0], ylim[1]-1])
		ax.set_title(pennum)
		plt.show()
		xy = plt.ginput(1, timeout = -1)[0]
		xy = np.array(xy)
		if np.prod(xy)>0:
			xy = np.array(xy)
			cfs.append(np.hstack((pennum, xy)))
		else:
			badunits.append(pennum)

		plt.close(fig)


		savepath = os.path.join(studydir, 'Sessions', experiment, 'cfs.txt')
		np.savetxt(savepath, cfs)

		# print badunits

	plt.close(fig)

def look_at_map(experiments, onlygood = False, prefix = 'nonA1', ax = None, makechanges = False, verbose = False):
	'''
	Reads the 
	'''
	
	if type(experiments) is str:
		experiments = [experiments]

	# imgpath = os.path.join(studydir, experiment, 'experimentfiles', 'ScopePhoto1.png')
	# img = ndimage.read(imgpath)

	nexperiments = len(experiments)
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.set_axis_bgcolor('0.7')
	for i, experiment in enumerate(experiments):

		ax.set_title(experiment)
		
		if onlygood:
			badunits = []
		else:
			badunits = np.sort(glob.glob(os.path.join(studydir, 'Sessions', experiment, 'fileconversion', '_RF*.h5')))

		goodunits = np.sort(glob.glob(os.path.join(studydir, 'Sessions', experiment, 'fileconversion', 'RF*.h5')))
		units = np.concatenate((goodunits, badunits))
	
		cfs = np.loadtxt(os.path.join(studydir, 'Sessions', experiment, 'cfs.txt'), ndmin = 1)

		
		maxx = -np.inf
		maxy = -np.inf
		minx = +np.inf
		miny = +np.inf

		for unit in units:
			f = h5py.File(unit, 'r')
			if not np.isnan(f['coord'][0]):
				x = f['coord'][0]
				y = f['coord'][1]
			else:
				x = y = 200
			f.close()
		
			unitnum = os.path.splitext(os.path.split(unit)[1])[0]
			p = re.compile('\d+')
			unitnum = np.int32(p.findall(unitnum))

			ix = cfs[:, 0] == unitnum # what is this unit's cf
			if ix.sum() == 0:
				clr = '0.5'
			elif unit in badunits:
				clr = '0.5'
			else:
				cf = cfs[ix, 1]
				clr = cf_to_rgba(cf)[0]
				print clr
				
			unitnum = os.path.splitext(os.path.split(unit)[1])[0]
			p = re.compile('\d+')
			unitnum = p.findall(unitnum)[0]
			unitnum = np.int32(unitnum)
			ax.text(x, y, unitnum, color = clr)
			maxx = np.max((maxx, x))
			minx = np.min((minx, x))
			maxy = np.max((maxy, y))
			miny = np.min((miny, y))

		ax.set_xlim([minx-0.5, maxx+0.5])
		ax.set_ylim([miny-0.5, maxy+0.5])
		# ax.set_xticklabels('')
		# ax.set_yticklabels('')
		fig.savefig(os.path.join(studydir, 'Maps', '%s_map.png' % experiment))
		ax.cla()
		
	# 
	# if verbose:
	# 	print 'Bad units\n...'
	# 	for badunit in badunits:
	# 		print os.path.splitext(os.path.split(badunit)[1])[0]
	# 
	# 	print 'Good units\n...'
	# 	for goodunit in goodunits:
	# 		print os.path.splitext(os.path.split(goodunit)[1])[0]
	# 
	# if makechanges:
	# 	unitnums = raw_input('Enter units you would like to change-->')
	# 
	# 	if len(unitnums) > 0:
	# 	
	# 		unitnums = np.int32(unitnums.split(' '))
	# 		nonA1path = os.path.join(studydir, 'Sessions', experiment, prefix + '.txt')
	# 		nonA1 = np.empty(0)
	# 		if os.path.exists(nonA1path):
	# 			nonA1 = np.loadtxt(nonA1path)
	# 		nonA1 = np.concatenate((nonA1, unitnums))
	# 		np.savetxt(nonA1path, nonA1)
	# 	
	# 	else:
	# 		print 'No changes made.'
	
	return fig, ax

def get_spktimes(spktimes, spktrials, stimparams, param):
	
	ix = make_spk_mask(spktrials, stimparams, param)
	return spktimes[ix], spktrials[ix]
	
	
def get_trials(stimparams, param):
	'''
	Returns the indices of trials in which non-nan entries of param match corresponding entries in stimparams.
	'''
	param = np.asarray(param)
	
	ignoreparams = np.isnan(param)
	param_ = param[~ignoreparams]
	stimparams_ = stimparams[:, ~ignoreparams]
	
	return (stimparams_ == np.tile(param_, (stimparams_.shape[0], 1))).all(1).nonzero()[0]

def make_spk_and_trial_masks(spktrials, stimparams, param):
	ntrials = stimparams.shape[0]
	ix = RF.get_trials(stimparams, param)
	trial_mask = np.zeros(ntrials, dtype = np.bool)
	spk_mask = np.zeros(spktrials.size, dtype = np.bool)
	for ix_ in ix:
		spk_mask[spktrials == ix_] = True
		trial_mask[ix_] = True
		
	return spk_mask, trial_mask
	
def make_spk_mask(spktrials, stimparams, param):
	
	ix = RF.get_trials(stimparams, param)
	spk_mask = np.zeros(spktrials.size, dtype = np.bool)
	for ix_ in ix:
		spk_mask[spktrials == ix_] = True
		
	return spk_mask
	
def make_trial_mask(stimparams, param):
		
	ntrials = stimparams.shape[0]
	ix = RF.get_trials(stimparams, param)
	trial_mask = np.zeros(ntrials, dtype = np.bool)
	for ix_ in ix:
		trial_mask[ix_] = True
		
	return trial_mask
		
		
def remove_no_cf_units(experiments):
	
	if type(experiments) is not list:
		experiments = [experiments]

	for experiment in experiments:
		print experiment

		cfs = np.loadtxt(os.path.join(studydir, 'Sessions', experiment, 'cfs.txt'))
		fpaths = glob.glob(os.path.join(studydir, 'Sessions', experiment, 'fileconversion', 'RF*.h5'))
		p = re.compile('(\d+)\.h5')

		for fpath in fpaths:
			absol, rel = os.path.split(fpath)

			unitnum = np.int32(p.findall(rel))[0]
			cf = cfs[cfs[:, 0] == unitnum, 1]
			associated_fpaths = glob.glob(os.path.join(studydir, 'Sessions', experiment, 'fileconversion', '[A-Za-z]*%3.3u.h5' % unitnum)) # not already underscore leading file names
			if (cf.size) == 0 or np.isnan(cf):
				for apath in associated_fpaths:
					_, rel = os.path.split(apath)
					os.rename(apath, os.path.join(absol, '_'+rel))


def remove_units(sessions, kind = 'nonA1'):
	'''
	Will find the nonA1.txt file and rename all of those RF (and related files)
	to start with an underscore '_'
	'''		
	
	if type(sessions) is str:
		sessions = [sessions]

	for session in sessions:
		print session

		nona1_path = os.path.join(studydir, 'Sessions', session, kind + '.txt')
		if os.path.exists(nona1_path):
			nona1 = np.loadtxt(nona1_path, ndmin = 1)
			absol = os.path.join(studydir, 'Sessions', session, 'fileconversion')
			for nona1_ in nona1:
				files = glob.glob(os.path.join(studydir, 'Sessions', session, 'fileconversion', '[A-Za-z]*%3.3u.h5' % nona1_))
				# rename all h5 files in the fileconversion folder with this unit number
				for f in files:
					absol, rel = os.path.split(f)
					os.rename(f, os.path.join(absol, '_'+rel))
				
				# delete all analysis files in the analysis folder with this unit number
				files = glob.glob(os.path.join(studydir, 'Sessions', session, 'analysis', '*%3.3i*' % nona1_))
				for f in files:
					os.remove(f)
			
						
		else:
			print 'No sites to remove!'		

def threshold_rf(rf, thresh_mag = 0.25):

	rf_thresh = rf.copy()
	rf_thresh[rf < (thresh_mag * rf.max())] = 0

	return rf_thresh

def rf_contact_sheet(experiments, showcleaned = True):

	if type(experiments) == str:
		experiments = [experiments]

	fig = plt.figure(figsize = (16, 12))
	for experiment in experiments:
		cf_path = os.path.join(studydir, 'Sessions', experiment, 'cfs.txt')
		if os.path.exists(cf_path):
			cfs = np.loadtxt(cf_path)
		else:
			cfs = None

		units = glob.glob(os.path.join(studydir, 'Sessions', experiment, 'fileconversion', 'RF*.h5'))
		units = np.sort(units)
		nunits = len(units)
		ncols = np.ceil(np.sqrt(nunits))

		fig.suptitle(experiment)
		for i, unit in enumerate(units):
			absol, rel = os.path.split(unit)
			unitname = re.findall('(\d+).h5', rel)[0]
			f = h5py.File(unit, 'r')
			rf = f['rf'].value
			f.close()
			if i == 0:
				rf_shape = rf.shape
				textx = rf_shape[1]/2.
				texty = rf_shape[0] - 1

			if showcleaned:
				# thresholded rf
				rf_thresh = threshold_rf(rf, 0.25)
				rf_clust, _ = findmaxcluster(rf_thresh, include_diagonal = False)

				RF = np.hstack([rf, np.ones((rf.shape[0], 1))*rf.max(), rf_clust])
			else:
				RF = rf

			ax = fig.add_subplot(ncols, ncols, i+1)
			ax.imshow(RF, interpolation = 'nearest', aspect = 'auto', cmap = 'hot')
			if not cfs is None:
				cf = cfs[cfs[:, 0]==np.int32(unitname), 1]
				ax.axvline(cf, color = 'w', lw = 3)
			ax.text(textx, texty, unitname, color = 'k', bbox = dict(facecolor = 'white', alpha = 0.8))
			if i < (nunits-1):
				ax.set_xticks([]);
			if i > 0:
				ax.set_yticks([]);


		fig.savefig(os.path.join(studydir, 'Sheets', 'RFs', 'rf_' + experiment + '.png'))
		fig.clf()


def findmaxcluster(B, cf = None, include_diagonal = False):
	'''
	Finds clusters defined by direct left-right or up-down pixel contact;
	-1 - cell is member of largest cluster
	0 - no activity in this cell
	1 - activity in this cell
	2 - cell neighbors cluster-member
	3 - cell is member of currently evaluated cluster
	'''

	sizes = np.ndarray(0)
	D = B > 0
	# pad with zeros above and below
	D = np.concatenate((np.zeros((1, D.shape[1])), D, np.zeros((1, D.shape[1]))), axis = 0)
	# pad with zeros right and left
	D = np.concatenate((np.zeros((D.shape[0], 1)), D, np.zeros((D.shape[0], 1))), axis = 1)
	D = D > 0
	D = np.float32(D)
	max_clust_size = 0
	if cf is not None:
		cf = np.int32(cf)
		max_cf_match = 0
	# everything equals 0 or 1
	while np.sum(D[D>=0])>0: # while there are entries in D
		# for each cell (except border)
		for n in range(1, D.shape[0]-1):
			for m in range(1, D.shape[1]-1):
				if D[n,m]>0: # if the center cell is filled
					D[n,m]=2 # set it equal to 2
					while np.sum(D==2)>0: # while there are stil 2s
						# for each cell (except border)
						for ll in range(1, D.shape[0]-1): 
							for k in range(1, D.shape[1]-1):
								if D[ll,k]==2: # if the other cell is also 2
									# numcluster n,m,ll,k
									if D[ll-1,k]==1: # if the cell to the left is filled
										D[ll-1,k]=2 # set it equal to 2
									if include_diagonal:
										if D[ll-1,k-1]==1: # find diagonal neighbors as well
											D[ll-1,k-1]=2;
										if D[ll-1,k+1]==1: # find diagonal neighbors as well
											D[ll-1,k+1]=2;
									if D[ll+1,k]==1:  # if the cell to the right is filled
										D[ll+1,k]=2
									if include_diagonal:
										if D[ll+1,k-1]==1: # find diagonal neighbors as well
											D[ll+1,k-1]=2;
										if D[ll+1,k+1]==1: # find diagonal neighbors as well
											D[ll+1,k+1]=2;
									if D[ll,k-1]==1: # downstairs neighbor
										D[ll,k-1]=2
									if D[ll,k+1]==1: # upstairs neighbor
										D[ll,k+1]=2
									D[ll,k]=3; # include the test cell in this cluster

					if cf is not None:
						cf_match = ((D==3).nonzero()[1] == cf).sum()
						if cf_match > max_cf_match:
							swap_winner = True
							max_cf_match = cf_match
						else:
							swap_winner = False
					elif cf is None:
						clust_size = np.sum(D==3)
						if clust_size > max_clust_size: # if this is the largest cluster
							swap_winner = True
							max_clust_size = clust_size
						else:
							swap_winner = False

					if swap_winner:
						D[D==-1]=0 # set previously largest cluster to 0
						D[D==3]=-1 # set current cluster to largest cluster
					else: # if this is not the largest cluster
						D[D==3]=0; # zero out the test cell
						
	C = -1 * B * D[1:-1, 1:-1]
	
	return (C, max_clust_size)

def calc_rf_psth(rast, stimparams):
	
	ufreqs = np.unique(stimparams[:, 0])
	uattens = np.unique(stimparams[:, 1])
	nfreqs, nattens = ufreqs.size, uattens.size
	rf_psth = np.empty((nattens, nfreqs, rast.shape[1]))
	for i in range(nfreqs):
		for j in range(nattens):
			ix = RF.get_trials(stimparams, np.array([ufreqs[i], uattens[j]]))
			rf_psth[j, i, :] = rast[ix, :].mean(0)
	
	return rf_psth
	
	
def calc_evoked_psth(rast, stimparams, rf_mask):

	# # get the unique frequencies and attenuations played
	# ufreqs = np.unique(stimparams[:, 0])
	# uattens = np.unique(stimparams[:, 1])
	# 
	# nfreqs, nattens = ufreqs.size, uattens.size
	# 
	# # initialize an empty array that contains one PSTH for each FREQ/ATTEN combination
	# rf_psth = np.empty((nattens, nfreqs, rast.shape[1]))
	# # initialize the nfreqs x nattens mask that will indication which FREQ/ATTEN combos are in the RF
	# 
	# # build the RF PSTH
	# for i in range(nfreqs):
	# 	for j in range(nattens):
	# 		ix = RF.get_trials(stimparams, (ufreqs[i], uattens[j]))
	# 		rf_psth_ = rast[ix, :].mean(0)
	# 		rf_psth[j, i, :] = rf_psth_
	
	rf_psth = calc_rf_psth(rast, stimparams)

	ev_psth = 1000 * rf_psth[rf_mask].mean(0)
	
	return ev_psth

def calc_rf_com(rf):
	'''
	Calculates the center of mass of the receptive field, one possible measure of CF
	Input:
		B - some receptive field, processed or raw
	Output:
		COM - the center of mass
	'''
	
	nfreq = rf.shape[1]
	x = np.sum(rf, 0)
	com = np.sum(np.arange(nfreq) * x) / np.float32(np.sum(x))
	
	return com

def calc_com(masses, radii):

	com = (masses*radii).sum() / masses.sum()
	return com

def calc_cf(B):
	'''
	input is smoothed/thresholded/clustered rf
	output is the mean tip location
	'''
	tmp = B.nonzero()
	a = tmp[0]
	f = tmp[1]
	thresh_freqs = f[a==a.max()]
	thresh_attens = a[a==a.max()]
	cf = calc_com(B[thresh_attens, thresh_freqs], thresh_freqs)
	
	return cf

def calc_bw_thresh(B):
	
	nattens = B.shape[0]
	# sum across rows and find the last one that is non-zero
	thresh = (B.sum(1) > 0).nonzero()[0].max()
	
	# calculate the bandwidth at each atten level above thresh
	bw = np.empty(nattens) * np.nan
	bw_lr = np.empty((nattens, 2)) * np.nan
	for i, b in enumerate(B[thresh::-1]):
		resp = b.nonzero()[0]
		if resp.size > 1:
			resp_max = resp.max(); resp_min = resp.min()
			bw_lr[i, :] = np.array([resp_min, resp_max])
			bw[i] = resp_max - resp_min + 1
		elif resp.size == 1:
			bw_lr[i, :] = np.array([resp[0], resp[0]])
			bw[i] = 0.1
		elif resp.size == 0:
			bw_lr[i, :] = np.array([0, 0])
			bw[i] = 0
			
	return bw, bw_lr, thresh
		
def plot_rf(rf, bw_lr = None, cf = None, thresh = None, ax = None, cmap = 'gray', axes_on = True, **kwargs):
	
	if type(rf) is h5py.File:
		rast = rf['rast'].value
		stimparams = rf['stimID'].value
		rf = calc_rf(rast, stimparams)
	
	nattens, nfreqs = rf.shape
	if ax is None:
		fig = plt.figure();
		ax = fig.add_subplot(111);
	else:
		fig = ax.get_figure()	
	
	im = ax.imshow(rf, interpolation = 'nearest', aspect = 'auto', cmap = cmap, origin = 'upper', **kwargs)
	if axes_on:
		cb = fig.colorbar(im, ax = ax)
		cb.set_label('Firing rate (spikes/s)')
	
	if bw_lr is not None:
		j = 0
		for i, (bw_l, bw_r) in enumerate(bw_lr[::-1]):
			if not np.isnan(bw_l):
				ax.plot((bw_l, bw_r), (j, j), color = 'r')
				j += 1
				
	if cf is not None:
		ax.axvline(cf, color = 'r', ls = '--')

	if thresh is not None:
		ax.axhline(thresh, color = 'g', ls = '--')
	
	if axes_on:
		ax.set_ylim([rf.shape[0]-0.5, -0.5])
		ax.set_xlim([-0.5, rf.shape[1]-0.5])
		ax.set_xticks(np.arange(0, nfreqs, 10))
		ax.set_xticklabels([('%0.0f' % (f/1000.)) for f in ix2freq[20:][np.arange(0, nfreqs, 10)]])
		ax.set_xlabel('Frequency (kHz)')
		ax.set_yticks(np.arange(0, nattens))
		ax.set_yticklabels(np.arange((nattens-1)*10, -1, -10))
		ax.set_ylabel('Sound Intensity (db SPL)')
	plt.show()
	
	return ax

def myconv(x, wind = None):
	if wind is None:
		wind = np.hanning(5)
		wind = wind / wind.sum()
	
	return np.convolve(x, wind, 'same')	

def calc_latency_by_stim(rast, stimparams):
	
	stim_psth, _ = Spikes.calc_psth_by_stim(rast, stimparams, bins = np.arange(0, 0.334, 0.001))
	stim_psth_smoo = np.apply_along_axis(myconv, 2, stim_psth)
	stim_peak_times = np.apply_along_axis(np.argmax, 2, stim_psth_smoo)
	
	return stim_peak_times

def plot_latency_by_stim(rast, stimparams):
	
	rf = calc_rf(rast, stimparams)
	stim_peak_times = calc_latency_by_stim(rast, stimparams)
	
	x = stim_peak_times.copy()
	x[x>80] = 80
	x[x<50] = 50

	fig = plt.figure()
	ax1 = fig.add_subplot(121)
	ax2 = fig.add_subplot(122)
	RF.plot_rf(rf, ax = ax1)
	RF.plot_rf(x.T, ax = ax2, cmap = 'jet')
	

def add_coords_to_h5py(sessions):
	if not type(sessions) is list:
		sessions = [sessions]
	
	for session in sessions:
		print session
		p = re.compile('(\d+)\.h5')
		coords = np.loadtxt(os.path.join(studydir, 'Sessions', session, 'experimentfiles', '%s.txt' % session), ndmin = 1)
		fpaths = glob.glob(os.path.join(studydir, 'Sessions', session, 'fileconversion', '*.h5'))
		for fpath in fpaths:
			print fpath
			absol, relat = os.path.split(fpath)
			unitnum = int(p.findall(relat)[0])
			coord = coords[coords[:, 0]==unitnum, 1:3]
			if coord is np.nan: coord = np.array((np.nan, np.nan))
			f = h5py.File(fpath, 'r+')
			# if f['coord'].value is not np.array((np.nan, np.nan)):
			try:
				f['coord'].write_direct(coord)
			except:
				print 'Failed!'
			f.close()
			

def calc_fake_strf(rast, stimparams):
	'''
	2-d matrix no.freqs x no.time bins, values are spontaneous firing rate. basically each row is a psth for a particular frequency
	'''
	
	stim_psth, _ = Spikes.calc_psth_by_stim(rast, stimparams)
	fake_strf = stim_psth.mean(1)
	
	return fake_strf

def plot_fake_strf(rast, stimparams):
	
	fig = plt.figure(figsize = (15, 10))
	ax1 = fig.add_axes([0.05, 0.2, 0.9, 0.75])
	ax2 = fig.add_axes([0.05, 0.05, 0.9, 0.15])
	
	fake_strf = calc_fake_strf(rast, stimparams)
	RF.plot_rf(fake_strf, ax = ax1, axes_on = False)
	psth = fake_strf.mean(0)
	psth_smoo = Spikes.hamming_smoo(psth, windlen = 5)
	ax2.plot(psth_smoo)
	ax2.set_xlim([0, rast.shape[1]])
	plt.draw(); plt.show();

	return fig

def add_rf_analysis(frf, funit, stimparams):

	
	rast = frf['rast'].value
	# stimparams = frf['stimID'].value
	

	# calculate, smooth, threshold, and cluster RF
	rf = calc_rf(rast, stimparams, normed = True, smooth = False)
	rf[:, :5] = 0
	nan_ix = np.isnan(rf)
	if nan_ix.sum() > 0:
		print 'NaNs found in rf! Converting NaNs to zeros...'
		print nan_ix.sum()
		rf[np.isnan(rf)] = 0
	rf_smoo = gaussian_filter(rf, 0.5)
	rf_thresh = rf_smoo.copy()
	rf_mask = rf_thresh<0.25*rf_thresh.max()
	rf_thresh[rf_mask] = 0
	(rf_clust, max_clust_size) = findmaxcluster(rf_thresh)
	if max_clust_size < 10:
		print 'FUCKING SMALL RF dude'
	
	# baseline_mean : mean firing rate 0-50ms for all trials	
	baseline_psth = 1000 * rast[:, :50].mean(0)
	baseline_mean = baseline_psth.mean()
	baseline_std = baseline_psth.std()
	fr_thresh = baseline_mean + 2*baseline_std

	# evoked and spontaneous psths
	ev_psth = calc_evoked_psth(rast, stimparams, rf_clust>0)
	ev_psth_smoo = Spikes.hamming_smoo(ev_psth, windlen = 5)
	spont_psth = calc_evoked_psth(rast, stimparams, rf_clust==0)
	

	# psth properties
	# peak firing rate
	resp_max = ev_psth_smoo.max()
	peak_time_ix = ev_psth_smoo.argmax()
	
	# first millisecond after the final sub-threshold millisecond
	resp_on_ix_ = (ev_psth_smoo[57:peak_time_ix]<fr_thresh).nonzero()[0]
	if not resp_on_ix_.any():
		resp_on_ix = peak_time_ix
	else:
		resp_on_ix = resp_on_ix_.max() + 58
	resp_off_ix = np.empty(3)
	resp_off_ix[0] = (ev_psth_smoo[peak_time_ix:]<resp_max/2).nonzero()[0].min() + peak_time_ix
	resp_off_ix[1] = (ev_psth_smoo[peak_time_ix:]<resp_max/4).nonzero()[0].min() + peak_time_ix
	resp_off_ix[2] = (ev_psth_smoo[peak_time_ix:]<fr_thresh).nonzero()[0].min() + peak_time_ix
	
	peak_time = peak_time_ix / 1000.
	resp_on = resp_on_ix / 1000.
	resp_off = resp_off_ix / 1000.
	
	# calculate a decay constant from the peak to the post-peak sub-threshold
	y = np.log(ev_psth[peak_time_ix:resp_off_ix[2]])
	tmp = linregress(np.arange(y.size), y)
	decay_const = tmp[0]
	
	# firing rate statistics on the rf
	rf_spont_mean = rf[rf_clust==0].mean()
	rf_resp_max = rf[rf_clust>0].max() - rf_spont_mean
	rf_resp_median = np.median(rf[rf_clust>0]) - rf_spont_mean
	rf_resp_mean = rf[rf_clust>0].mean() - rf_spont_mean
	rf_size = (rf_clust>0).sum()

	# automated CF measures (center-of-mass and auto)
	# automated BW measures
	cf_com = calc_rf_com(rf_clust)
	cf_auto = calc_cf(rf_clust)
	(bw, bw_lr, thresh) = calc_bw_thresh(rf_clust)

	# # PSTH for freq/attens in RF
	# fig = plt.figure()
	# ax = fig.add_subplot(111)
	# ax.plot(ev_psth_smoo, 'b.-')
	# ax.axhline(baseline_mean+3*baseline_std, color = 'r', ls = '--')
	# ax.axvline(resp_on_ix, color = 'r', ls = '--')
	# ax.axvline(resp_off_ix[0], color = 'r', ls = '--')
	# ax.axvline(resp_off_ix[1], color = 'r', ls = '--')
	# ax.axvline(resp_off_ix[2], color = 'r', ls = '--')
	
	rf_shape = rf.shape
	d = funit.require_dataset('rf', rf_shape, float)
	d.write_direct(rf)
	d = funit.require_dataset('rf_clust', rf_shape, float)
	d.write_direct(rf_clust)
	d = funit.require_dataset('rf_spont_mean', (), float)
	d.write_direct(np.array([rf_spont_mean]))
	d = funit.require_dataset('rf_resp_median', (), float)
	d.write_direct(np.array([rf_resp_median]))
	d = funit.require_dataset('rf_resp_mean', (), float)
	d.write_direct(np.array([rf_resp_mean]))
	d = funit.require_dataset('rf_resp_max', (), float)
	d.write_direct(np.array([rf_resp_max]))
	d = funit.require_dataset('resp_max', (), float)
	d.write_direct(np.array([resp_max]))
	d = funit.require_dataset('peak_time', (), float)
	d.write_direct(np.array([peak_time]))
	d = funit.require_dataset('resp_on', (), float)
	d.write_direct(np.array([resp_on]))
	d = funit.require_dataset('resp_off', (3,), float)
	d.write_direct(np.array([resp_off]))
	d = funit.require_dataset('decay_const', (), float)
	d.write_direct(np.array([decay_const]))	
	d = funit.require_dataset('rf_size', (), int)
	d.write_direct(np.array([rf_size]))
	d = funit.require_dataset('baseline_mean', (), float)
	d.write_direct(np.array([baseline_mean]))
	d = funit.require_dataset('cf_auto', (), float)
	d.write_direct(np.array([cf_auto]))
	d = funit.require_dataset('cf_com', (), float)
	d.write_direct(np.array([cf_com]))
	d = funit.require_dataset('bw', (8,), float)
	d.write_direct(bw)
	d = funit.require_dataset('bw_lr', (8,2), float)
	d.write_direct(bw_lr)
	d = funit.require_dataset('thresh', (), int)
	d.write_direct(np.array([thresh]))
	d = funit.require_dataset('ev_psth', (333,), float)
	d.write_direct(ev_psth[:333])

	
	