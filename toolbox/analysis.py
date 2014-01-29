import scipy.stats as st
import numpy as np
import matplotlib.pyplot as plt
import glob, os
import misc; reload(misc)
import RF; reload(RF)
import itertools

# basedir = '/Users/robert/Desktop/Fmr1_RR'
basedir = '/Volumes/BOB_SAGET/Fmr1_RR/'

def gap_startle_ratio(data, title = None, ax = None, show_all = True, animals = None):
	if ax is None:
		fig = plt.figure();
		ax = fig.add_subplot(111);
	
	freqs = np.unique(data['freq'])
	nfreqs = freqs.size
	gaps = np.unique(data['gap'])
	ngaps = gaps.size
	if animals is None:
		animals = np.unique(data['animal'])
	nanimals = animals.size
	x = np.arange(ngaps)
	ppi = np.empty((nfreqs, ngaps, nanimals))
	for f, freq in enumerate(freqs):
		for a, animal in enumerate(animals):
			dat = data[np.c_[data['animal']==animal, data['freq']==freq].all(1)]
			basestartle = dat[dat['gap']==0]['maxstartle'].mean()
			for r, gap in enumerate(gaps):
				dat_ = dat[dat['gap']==gap]['maxstartle']
				ppi[f, r, a] = dat_.mean() / basestartle
		
		ax.errorbar(x+ngaps*f, st.nanmean(ppi[f, ...], 1), yerr = st.nanstd(ppi[f, ...], 1), lw = 3)

		if show_all:
			ax.plot(x+ngaps*f, ppi[f, ...], color = '0.7')
			
	ax.set_xticks(np.arange(nfreqs*ngaps))
	ax.set_xticklabels(np.tile(gaps, nfreqs))
	ax.axhline(1, color = 'r', ls = '--')
	ax.set_ylabel('PPI')
	ax.set_xlabel('Gap duration (s)')
	ax.set_title(title)

def calc_rate_startle_ratio(dat):
	
	return st.nanmean(dat[dat['cued']==1]['maxstartle']) / st.nanmean(dat[dat['cued']==0]['maxstartle'])
	
def rate_startle_ratio(data, title = None, ax = None, show_all = True):

	if ax is None:
		fig = plt.figure();
		ax = fig.add_subplot(111);
	
	freqs = np.unique(data['freq'])
	nfreqs = freqs.size
	rates = np.unique(data['rate'])
	nrates = rates.size
	animals = np.unique(data['animal'])
	nanimals = animals.size
	x = np.arange(nrates)
	ppi = np.empty((nfreqs, nrates, nanimals))
	for f, freq in enumerate(freqs):
		for r, rate in enumerate(rates):
			for a, animal in enumerate(animals):
				dat_ = data[np.c_[data['freq']==freq, data['rate']==rate, data['animal']==animal].all(1)]
				ppi[f, r, a] = calc_rate_startle_ratio(dat_)
		
		ax.errorbar(x+nrates*f, st.nanmean(ppi[f, ...], 1), yerr = st.nanstd(ppi[f, ...], 1), lw = 3)

		if show_all:
			ax.plot(x+nrates*f, ppi[f, ...], color = '0.7')
			
	ax.set_xticks(np.arange(nfreqs*nrates))
	ax.set_xticklabels(np.tile(rates, nfreqs))
	ax.axhline(1, color = 'r', ls = '--')
	ax.set_ylabel('PPI')
	ax.set_xlabel('Rate (pps)')
	ax.set_title(title)

def startle_ratio_by_indep(indep, data, ax = None, color = None, show_all = False):

	if ax is None:
		fig = plt.figure();
		ax = fig.add_subplot(111);

	ax.axhline(1, c = 'r', ls = '--')
	animals = np.unique(data['animal'])
	nanimals = len(animals)
	indeps = data[indep]
	uindeps = np.unique(indeps)
	nindeps = uindeps.size
	xax = np.arange(nindeps)
	ratio = np.empty((nanimals, nindeps))
	p = []
	for i, animal in enumerate(animals):
		data_ = data[data['animal'] == animal]
		
		for j, uindep in enumerate(uindeps):
			dat = data_[data_[indep]==uindep]
			ratio[i, j] = calc_startle_ratio(dat)

		if show_all:
			if color is None:
				ax.plot(xax, ratio[i, :], 'o-')
			else:
				ax.plot(xax, ratio[i, :], 'o-', color = color)
	
	if color is None:		
		l = ax.errorbar(misc.jitter(xax, 0.1), st.nanmean(ratio, 0), yerr = st.nanstd(ratio, 0), color = 'k', lw = 3)
	else:
		l = ax.errorbar(misc.jitter(xax, 0.1), st.nanmean(ratio, 0), yerr = st.nanstd(ratio, 0), color = color, lw = 3)
	
	ax.set_ylim([0.5, 1.5])
	ax.set_xlim([-1, nindeps])
	ax.set_xticks(np.arange(nindeps))
	ax.set_xticklabels(uindeps)
	ax.set_xlabel(indep)
	ax.set_ylabel('Max startle amplitude (cued/uncued)')
	# plt.show()
	return ax, l


def startle_ratio_by_indep_by_group(indep, data, group, ax = None):

	if ax is None:
		fig = plt.figure();
		ax = fig.add_subplot(111);
	
	ugroups = np.unique(data[group])
	ngroups = ugroups.size
	clrs = 'bgrymc'
	l = []
	for i in range(ngroups):
		data_ = data[data[group] == ugroups[i]]
		ax, l_ = startle_ratio_by_indep(indep, data_, ax = ax, color = clrs[i])
		l.append(l_)
		
	ax.legend(l, ugroups)

def bar_by_indep_by_group_1d(dep, indep, data, group_keys, bins = None, visible = True, ax = None, color = 'b', show_all = False, use_bar = False, ordered_x = False):

	'''
	For 1D data (a value such as bandwidth), plots 
	'''

	'''
	make ugroup_keys, a list of the group names
	make groups_bin, the group values for each data point (binned for easier indexing)	
	'''
	ngroups = len(group_keys)
	if bins is None:
		bins = [None] * ngroups
	# get all levels for all groups
	groups_bin = []
	ugroup_keys = []
	for gr, bin in zip(group_keys, bins):
		if bin is None:
			ugroup_keys.append(list(np.unique(data[gr])))
			groups_bin.append(data[gr])
		else:
			ugroup_keys.append(bin)
			groups_bin.append(misc.bin(data[gr], bin))

	# loop through all combinations of groups
	ndata = data.size
	for ugroup_key in itertools.product(*ugroup_keys):
		ix = np.empty(ndata) * True
		for group_key, ugroup_key_ in zip(group_keys, ugroup_key):
			ix_ = data[group_key] == ugroup_key_
			ix = np.vstack((ix, ix_)).all(0)
			dat
			
			
		data[group[0]==gr[0]]

	if ax is None:
		fig = plt.figure();
		ax = fig.add_subplot(111);

	clrs = 'bgrym'
	groups = data[group]
	ugroups = np.unique(groups)
	ngroups = ugroups.size
	for i in range(ngroups):
		data_ = data[data[group] == ugroups[i]]
		bar_by_indep_1d(dep, indep, data_, visible = visible, ax = ax, color = clrs[i], show_all = show_all, use_bar = use_bar, ordered_x = ordered_x)

	ax.legend(ugroups)

# def bar_by_indep_by_group_1d(dep, indep, data, group, visible = True, ax = None, color = 'b', show_all = False, use_bar = False, ordered_x = False):
# 	
# 	'''
# 	For 1D data (a value such as bandwidth), plots 
# 	'''
# 	
# 	if ax is None:
# 		fig = plt.figure();
# 		ax = fig.add_subplot(111);
# 	
# 	clrs = 'bgrym'
# 	groups = data[group]
# 	ugroups = np.unique(groups)
# 	ngroups = ugroups.size
# 	for i in range(ngroups):
# 		data_ = data[data[group] == ugroups[i]]
# 		bar_by_indep_1d(dep, indep, data_, visible = visible, ax = ax, color = clrs[i], show_all = show_all, use_bar = use_bar, ordered_x = ordered_x)
# 	
# 	ax.legend(ugroups)

def bar_by_indep_2d(dep_key, indep_key, data, ax = None, bins = None, color = 'b', show_all = False):

	x = np.asarray(data[indep_key])
	y = np.asarray(data[dep_key])

	if bins is None:
		x_bin = x
	else:
		x_bin = misc.bin(x, bins)

	bins = np.unique(x_bin)
	nbins = bins.size

	y_mean = np.empty(nbins)
	y_sem = np.empty(nbins)
	for i in range(nbins):
		y_ = y[x_bin == bins[i]]
		y_mean[i] = st.nanmean(y_)
		y_sem[i] = st.nanstd(y_) / np.sqrt(y_.size)

	if ax is None:
		fig = plt.figure();
		ax = fig.add_subplot(111);

	if show_all:
		ax.scatter(x, y, color = color, alpha = 0.25)
		lw = 2
	else:
		lw = 1
	ax.errorbar(bins, y_mean, yerr = y_sem, color = color, lw = lw)

	ax.set_xlim([bins[0]-1, bins[-1]+1])
	plt.show()

def bar_by_indep(dep_key, indep_key, data, ax = None, bins = None, color = 'b', show_all = False):
	
	'''
	passes along the arguments to the 1D case (i.e. independent is CF, dependent is baseline firing rate)
	or the 2D case (i.e. independent is RR, dependent is RRTF)
	'''
	
	if len(data[dep_key].shape) == 1:
		bar_by_indep_1d(dep_key, indep_key, data, ax = ax, bins = bins, color = color, show_all = show_all)
	elif len(data[dep_key].shape) == 2:
		bar_by_indep_2d(dep_key, indep_key, data, ax = ax, bins = bins, color = color, show_all = show_all)
	
def bar_by_indep_1d(dep_key, indep_key, data, ax = None, bins = None, color = 'b', show_all = False):
	
	'''
	the 1D case (i.e. independent is CF, dependent is baseline firing rate)
	'''
	if type(dep_key) is str:
		y = np.asarray(data[dep_key])
	else:
		y = np.asarray(dep_key)
	if type(indep_key) is str:
		x = np.asarray(data[indep_key])
	else:
		x = np.asarray(x)
		
	if bins is None:
		x_bin = x
	else:
		x_bin = misc.bin(x, bins)

	bins = np.unique(x_bin)
	nbins = bins.size
	
	y_mean = np.empty(nbins)
	y_sem = np.empty(nbins)
	for i in range(nbins):
		y_ = y[x_bin == bins[i]]
		y_mean[i] = st.nanmean(y_)
		y_sem[i] = st.nanstd(y_) / np.sqrt(y_.size)
	
	if ax is None:
		fig = plt.figure();
		ax = fig.add_subplot(111);
		
	if show_all:
		ax.scatter(x, y, color = color, alpha = 0.25)
		lw = 2
	else:
		lw = 1
	ax.errorbar(bins, y_mean, yerr = y_sem, color = color, lw = lw)

	ax.set_xlim([bins[0]-1, bins[-1]+1])
	plt.show()
		
def bar_by_indep_by_group_2d(dep_key, indep_key, group_key, data, bins = None, group_bins = None, visible = True, ax = None, color = 'bgmcry', show_all = False):
	
	groups = data[group_key]
	
	if group_bins is None:
		group_bins = groups
	else:
		group_bins = misc.bin(groups, group_bins)
	
	bins = np.unique(group_bins)
	nbins = bins.size
	
	if ax is None:
		fig = plt.figure();
		ax = fig.add_subplot(111);
	
	for i, bin in enumerate(bins):
		data_ = data[group_bins==bin]
		if data_.size > 0:
			bar_by_indep_2d(dep_key, indep_key, data_, visible = visible, ax = ax, color = color[i], show_all = show_all)

def bar_by_indep_2d(dep_key, indep_key, data, visible = True, ax = None, color = 'b', show_all = False, use_bar = False, **kwargs):
	
	'''
	the 2D case (i.e. independent is RR, dependent is RRTF)
	'''
	if type(indep_key) is str:
		x = data[0][indep_key]
	else:
		x = indep_key
	y = data[dep_key]
	nbins = x.size
	

	y_means = st.nanmean(y, 0)
	y_sems = st.nanstd(y, 0) / np.sqrt(y.shape[0])
	
	if visible:
		if ax is None:
			fig = plt.figure();
			ax = fig.add_subplot(111);

		if show_all:
			ax.plot(x, y.T, 'gray')
		
		if use_bar:
			line, = ax.bar(x, y_means, yerr = y_sems, color = color, **kwargs)
		else:
			line, _, _ = ax.errorbar(x, y_means, yerr = y_sems, lw = 2, color = color, **kwargs)
	
	plt.show()
	return line, ax


def bar_by_indep_by_group_3d(dep_key, indep_key, data, visible = True, ax = None, group_keys = 'gen', bins = None, iscategorical = True, show_all = False, use_bar = False, fig = None, **kwargs):
	'''
	plots a comparison over many groups
	attemps to make multiple line plots containing 2 lines per plot for ease-of-viewing
	gen exp sess
	'''
	
	ngroups = len(group_keys) # number of groups
	# bin data for each group
	if bins is None:
		bins = [None] * ngroups
		
	groups = [] # actual entries for each group
	groups_bins = []
	ugroups = []
	nbins_per_group = []
	for bin, group_key in zip(bins, group_keys):

		groups.append(data[group_key])
		if bin is None:
			groups_bins.append(groups[-1]) # binned labels for each data point
			ugroups.append(np.unique(groups[-1])) # all unique binned labels
		else:
			groups_bins.append(misc.bin(groups[-1], bin))
			ugroups.append(bin[:-1])

		nbins_per_group.append(ugroups[-1].size) # how many unique labels in each group

	''' ix is for all possible group combinations, a boolean vector indicating membership of each data point in the data set
		l is passed to itertools.product to iterate through all combinations of group values
	'''
	ix, l = build_index(groups_bins, ugroups, nbins_per_group)

	if fig is None:
		fig = plt.figure(figsize = (11, 14));
		
	naxs1 = nbins_per_group[0] # exp / row
	naxs2 = nbins_per_group[1] # cf / column		

	cm = plt.cm.gist_ncar
	clrs = [cm(i) for i in np.linspace(0, 0.9, nbins_per_group[-1])]

	nlines = np.zeros((nbins_per_group[0], nbins_per_group[1]))
	for i in itertools.product(*l):
		print i
		axisnum = i[0]*naxs2 + i[1] + 1
		ax = fig.add_subplot(naxs1, naxs2, axisnum)
		if ix[i].sum()>0:
			nlines[i[0], i[1]] += 1
			bar_by_indep_2d(dep_key, indep_key, data[ix[i]], ax = ax, color = clrs[i[-1]])

			if i[1] == 0:
				ax.set_ylabel('%s\n%s' % (str(ugroups[0][i[0]]), dep_key), multialignment = 'center')
			if i[0] == 0:
				ax.set_title('%s = %s' % (str(group_keys[1]), str(ugroups[1][i[1]])))
			else:
				ax.set_title('')
	
	axs = fig.get_axes()
	misc.sameyaxis(axs)
	misc.samexaxis(axs)
	done = False
	i = 0
	while not done:
		i += 1
		if nlines[i % nbins_per_group[0], int(np.floor(i/nbins_per_group[0]))] == nbins_per_group[2]:
			axs[i].legend(ugroups[2])
			done = True
	
	plt.show()
	return fig
	
def build_index(groups_bins, ugroups, nbins_per_group):

	nsamp = len(groups_bins[0])
	ix = np.empty((np.hstack([np.asarray(nbins_per_group), nsamp])), dtype = bool)

	
	l = []
	for i in range(len(nbins_per_group)):
		l.append(range(nbins_per_group[i]))
	
	for i in itertools.product(*l):
		ix_ = np.ones(nsamp, dtype = bool)
		for j in range(len(i)):
			ix_[groups_bins[j]!=ugroups[j][i[j]]] = False
		ix[i] = ix_
	
	return ix, l



def maps_by_group(gens = ['wt', 'ko'], exps = ['nai', 'exp', 'w1', 'w2', 'w3']):
	
	for i, (gen, exp) in enumerate(itertools.product(gens, exps)):
		sesss = glob.glob(os.path.join(basedir, 'Sessions', gen+'_'+exp+'*'))
		sesss = [os.path.basename(s) for s in sesss]
		nsess = len(sesss)
		if nsess > 0:
			fig, ax = RF.look_at_map(sesss)
			fig.savefig(os.path.join(basedir, 'maps', '%s_%s' % (gen, exp)))
			plt.close('all')


	return


def dep_by_indep_for_2_groups(dep_key, indep_key, group1_key, group2_key, data, v = False):

	'''
	takes 2-D data and splits it up by two groups
	it plots group1 groups on separate axes and group2 groups as different colors on the same axis
	'''
	
	group1 = np.unique(data[group1_key]); ngroup1 = len(group1)
	group2 = np.unique(data[group2_key]); ngroup2 = len(group2)

	leg = []
	line = []
	for i in range(ngroup1):
		leg.append([])
		line.append([])
	colors = 'bgmyrc'

	fig = plt.figure()
	fig.tight_layout()
	for i, (g1, g2) in enumerate(itertools.product(range(ngroup1), range(ngroup2))):
		if v:
			print group1[g1], group2[g2]
		data_ = data[np.vstack((data[group1_key]==group1[g1], data[group2_key]==group2[g2])).all(0)]
		ax = fig.add_subplot(1, ngroup1, g1+1)
		ax.set_title(group1[g1])
		if data_.size > 0:
			line_, ax_ = bar_by_indep_2d(dep_key, indep_key, data_, ax = ax, color = colors[g2])
			line[g1].append(line_)
			leg[g1].append(group2[g2])

	for g1 in range(ngroup1):
		ax = fig.add_subplot(1, ngroup1, g1+1)
		ax.legend(line[g1], leg[g1])
	
	misc.samexaxis(fig.get_axes())
	misc.sameyaxis(fig.get_axes())
	
	plt.show()
	
	return fig


	
def analyze(DB):
	
	analysisdir = '/Users/robert/Desktop/Fmr1_RR/Analysis'
	fig = analysis.dep_by_indep_for_2_groups('rrtf_tone_norm', 'rr', 'gen', 'exp', DB, v = True)
	fig.savefig(os.path.join(analysisdir, 'rrtf_gen_and_exp.png'))
	fig = analysis.dep_by_indep_for_2_groups('rrtf_tone_norm', 'rr', 'exp', 'gen', DB, v = True)
	fig.savefig(os.path.join(analysisdir, 'rrtf_exp_and_gen.png'))

	group_keys = np.array(['exp', 'cf', 'gen'])
	bins = np.array([None, np.array([0, 30, 60]), None])
	fig = plt.figure(figsize = (14, 11))
	for i in itertools.permutations([0, 1, 2]):
		print group_keys[np.array(i)], bins[np.array(i)]
		analysis.bar_by_indep_by_group_3d('rrtf_tone_norm', 'rr', DB, group_keys = group_keys[np.array(i)], bins = bins[np.array(i)], fig = fig)
		fig.savefig(os.path.join('/Users/robert/Desktop/Fmr1_RR/analysis', '_'.join(group_keys[np.array(i)]) + '.png'))
		fig.clf();
	