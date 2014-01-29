import glob, os
import numpy as np
import misc; reload(misc)
import wavey; reload(wavey)
from matplotlib.mlab import specgram
import matplotlib.pyplot as plt
import struct
import pdb
import h5py
from scipy.stats import nanmean, nanstd
from scipy.signal import sepfir2d
from scipy.sparse import lil_matrix
# from findclusters import findclusters


# basedir = '/Volumes/BOB_SAGET/Vocalization/'
# basedir = '/Users/robert/Desktop/Vocalization'
basedir = '/Volumes/BOB_SAGET/Vocalization'
Fs = 327868 # sample rate for sound recordings
dtype = np.dtype([('gen', 'S2'), ('cage', 'S8'), ('age', 'i4'), ('chunk', 'i4'), ('ncalls', 'i4')])

	
def plot_spec(Pxx, F = None, T = None, ax = None, F_range = None, **kwargs):
	'''
	Plots and correctly formats the spectrogram (axis tick labels)
	Input:
		Pxx : MxN spectrogram array (M frequencies, N time bins)
		F : vector of length M indicating the frequencies included in the spectrogram
		ax : axis on which to plot
		lowpass : boolean, whether or not to exclude frequencies lower than 2 kHz
		kwargs : properties for plotting
	Output:
		ax : axis containing the spectrogram
	'''

	ax, fig = misc.axis_check(ax)
	
	if T is None:
		extent = [0, Pxx.shape[1]-1, None, None]
	else:
		extent = [T[0], T[-1], None, None]
		
	if F is None:
		extent[2:4] = [0, Pxx.shape[0]-1]
	else:
		extent[2:4] = [F[0], F[-1]]

	ax.imshow(Pxx, aspect = 'auto', origin = 'lower', interpolation = 'nearest', extent = extent, **kwargs)
	
	# ax.set_yticks(np.arange(0, F.size, 20))
	# ax.set_yticklabels(F[np.arange(0, F.size, 20)])

	# fig.tight_layout()

	return ax

def load_spec(path):
	'''
	Loads the spectrogram file that was saved by the function sess_spectrogram
	Input:
		path : either the absolute path for the spectrogram file or a list containing the
				cage ID, block number, and sess
	Output:
		P : MxN spectrogram array (M frequencies, N time bins)
		F : vector of length M indicating the frequencies included in the spectrogram
	'''

	exten = os.path.splitext(path)[-1]
	
	if exten == '.npz':
		tmp = np.load(path)
		P = np.array(tmp['P'].min().todense())
		F = tmp['F']
		T = tmp['T']
		tmp.close()
		
	elif exten == '.h5':
		NFFT = 512
		noverlap = 256

		f = h5py.File(path, 'r+')
		s = f['s'].value
		fs = f['fs'].value
		f.close()

		(P, F, T) = specgram(s, Fs = fs, NFFT = NFFT, noverlap = noverlap)

	return P, F, T


def bin2arr(f, sampwidth = 4, ind = (0, None)):
	'''
	Converts the binary file to a numpy vector
	Input:
		f : an open file containing the binary data
		sampwidth : number of bytes per sample; 4-single, 8-double
		ind : optional tuple containing the range of samples to load from the file
	Output:
		x : numpy vector of data
		fs : sample rate obtained from the file header
	'''
	
	''' Get sampling frequency from header '''
	f.seek(0)
	header = struct.unpack('>5d', f.read(40))
	fs = int(header[2])

	''' Read data '''
	f.seek(40 + ind[0] * sampwidth) # seek to starting index
	if ind[1] is None:
		xs = f.read() # read entire file
	else:
		xs = f.read((ind[1]-ind[0]) * sampwidth) # read part of file

	nsamp = len(xs) / sampwidth
	if sampwidth == 4:
		dtype = 'f' # if single precision float
	elif sampwidth == 8:
		dtype = 'd' # if double precision float
	x = struct.unpack('>%s%s' % (nsamp, dtype), xs)
	del xs
	
	x = np.array(x, dtype = 'float32')
	
	return x, fs
	
def convert_bin_to_wav(f):
	'''
	Wrapper that loads a binary file, converts it to a numpy array, and saves it to a wav file
	Input:
		f : an open file containing the binary data
	'''
	
	x, fs = bin2arr(f)
	wavey.write(x, fs, f.name+'.wav', float)
	del x


def power_search(fname):
	'''
	Search for power peaks in a defined frequency band. Validated on several 60s recordings. KO4_0_7_30_2012 : automated-0 calls, manual-0 calls. KO_3_7_30_2012 : automated-0
	Input:
		Pxx : MxN spectrogram array (M frequencies, N time bins)
		F : vector of length M indicating the frequencies included in the spectrogram
		F_range : tuple containing the range of frequencies in which to search
		path_prefix : string containing the path prefix that will be used to save out the files
			(suffix is the index at which the call is located)
	'''
	P, F, T = load_spec(fname)
	
	P[40:43] = 0
	P[153:156, :] = 0
	
	F_ix = (F>50000)
	P = P[F_ix, :]
	F = F[F_ix]
	
	wind = np.hamming(30)
	wind = wind / wind.sum()
	P_sum = P.sum(0)

	P_smoo = np.convolve(P.sum(0), wind, 'same')
	
	starts_x, heights = misc.find_peaks(P_smoo)
		
	dtype = np.dtype([('starts_x', 'i8'), ('heights', 'f8')])
	xy = np.empty(starts_x.size, dtype = dtype)
	xy['starts_x'] = starts_x
	xy['heights'] = heights
	xy.sort(order = 'heights')
	xy = xy[::-1] # reverse array
	upper_ix = np.int32(xy.size * 0.05)
	xy = xy[:upper_ix]
	# 
	# else:
	# 	xy = np.nan * np.empty((1, 2))
		
	savefname = os.path.split(os.path.splitext(fname)[0])[-1] + '_blobs.txt'
	np.savetxt(os.path.join(basedir, 'analysis', 'power_search_output', savefname), xy)
	
	return xy

def run_power_search():
	
	cagedirs = np.sort(glob.glob(os.path.join(basedir, 'Cages/*')))

	for cagedir in cagedirs:
		cage = os.path.split(cagedir)[1]
		print cage
		sessdirs = np.sort(glob.glob(os.path.join(cagedir, 'P*')))
		for sessdir in sessdirs:
			_, sess = os.path.split(sessdir)
			print sess
			chunks = np.sort(glob.glob(os.path.join(sessdir, '*.npz')))
			for chunk in chunks:
				chunk_ = os.path.splitext(os.path.split(chunk)[-1])[0] + '_blobs.txt'
				if not os.path.exists(os.path.join(basedir, 'analysis', 'power_search_output', chunk_)):
					print os.path.split(chunk)[-1]
					power_search(chunk)
	
def compress_file(fname, v = True):
	f = open(fname, 'rb')
	s, fs = bin2arr(f)
	f.close()
	f = h5py.File(fname + '.h5')
	f.create_dataset('s', data = s, compression = 'gzip')
	f.create_dataset('fs', data = fs)
	f.close()
	
def compress_sess(sessdir, v = True):
	
	fnames = np.sort(glob.glob(os.path.join(sessdir, '*')))
	fnames = [fname for fname in fnames if not len(os.path.splitext(fname)[-1])]
	for fname in fnames:
		absol, ext = os.path.splitext(fname)
		if not os.path.exists(absol + '.h5'):
			if v:
				print fname
			try:
				compress_file(fname)
			except:
				print 'File not saved!'
	
def compress_sess_all():
	
	cagedirs = np.sort(glob.glob(os.path.join(basedir, 'Cages', '*')))
	for cagedir in cagedirs:
		_, cage = os.path.split(cagedir)
		print cagedir
		sessdirs = np.sort(glob.glob(os.path.join(cagedir, '*')))
		for sessdir in sessdirs:
			_, sess = os.path.split(sessdir)
			compress_sess(sessdir)
	

def calc_specgram_all():
	
	cagedirs = glob.glob(os.path.join(basedir, 'Cages', '*'))
	cagedirs.sort()
	
	for cagedir in cagedirs:
		cage = os.path.split(cagedir)[1]
		print 'Cage %s' % cage
		sessdirs = glob.glob(os.path.join(cagedir, 'P*'))
		sessdirs.sort()
		for sessdir in sessdirs:
			_, sess = os.path.split(sessdir)
			print 'Session %s' % sess
			chunks = glob.glob(os.path.join(sessdir, '*.h5'))
			chunks.sort()
			for chunk in chunks:
				absol, relat = os.path.split(chunk)
				if not os.path.exists(os.path.join(absol, os.path.splitext(relat)[0]) + '.npz'):
					print os.path.split(chunk)[-1]
					calc_specgram(chunk)
	
def zscore_spec(P):
	
	P_std = P.std(1)
	P_mean = P.mean(1)
	P_z = (P-np.tile(P_mean, (P.shape[1], 1)).T) / np.tile(P_std, (P.shape[1], 1)).T
	
	return P_z
	
def lowpass_spec(P, F):
	
	# only save power between 25 and 110 kHz (Holy and Guo 2005)
	F_range = (25000, 120000)
	F_ix = np.vstack((F>=F_range[0], F<=F_range[1])).all(0)
	P = P[F_ix, :]
	F = F[F_ix]
	
	return P, F

def calc_specgram(fname):
	
	# load raw spectrogram
	P, F, T = load_spec(fname)
	
	# remove low frequencies
	P, F = lowpass_spec(P, F)
	P = zscore_spec(P)
	
	# threshold
	thresh_ = P.mean(1)	# calculate mean for each frequency
	thresh1 = np.tile(thresh_, (P.shape[1], 1)).T
	# P_thresh = P.copy()
	P[P < thresh1] = 0
	# thresh2 = 4E-11
	# Pxx_thresh[Pxx_thresh < thresh2] = 0
	
	P_lil = lil_matrix(P)
	
	newfname = os.path.splitext(fname)[0] + '.npz'
	np.savez(newfname, P = P_lil, F = F, T = T)

	
def summarize_all_sess():
		cages = glob.glob(os.path.join(basedir, 'Cages/*'))
		cages.sort()
		for cage in cages:
			print '%s' % cage
			sesss = glob.glob(os.path.join(cage, '*'))
			sesss.sort()
			for sess in sesss:
				_, sess_ = os.path.split(sess)
				files = glob.glob(os.path.join(cage + '_*.h5'))
				files.sort()
				print '\t%s\t%i' % (sess, len(files))
			print '\n'


def plot_spec_incr(i, ax):
	ax.set_xlim((i, i+1000))
	plt.show()
	return i+1000

def manual_count(fname, ax, fig):
	
	Pxx, F, T = load_spec(fname)
	
	savefname = os.path.splitext(os.path.split(fname)[1])[0] + '_mancalls.txt'
	
	ax = plot_spec(np.power(Pxx, 0.05), ax = ax)
	fig = ax.get_figure()
	hits = np.empty((0, 2), dtype = 'float32')
	for i in np.arange(0, Pxx.shape[1], 500):
		ax.set_xlim((i, i+1000))
		plt.show()
		xy = fig.ginput(n = 0, timeout = 0)
		if len(xy) > 0:
			hits = np.vstack((hits, xy))
			for xy_ in xy:
				ax.axvline(xy_[0], color = 'r', ls = '--')
		np.savetxt(os.path.join(basedir, 'analysis', 'manual_count_output', savefname), hits)

def run_gui_count():
	
	cagedirs = glob.glob(os.path.join(basedir, 'Cages', '*'))
	cagedirs.sort()
	
	for cagedir in cagedirs:
		cage = os.path.split(cagedir)[1]
		print cage
		sessdirs = glob.glob(os.path.join(cagedir, 'P*'))
		sessdirs.sort()
		for sessdir in sessdirs:
			_, sess = os.path.split(sessdir)
			print sess
			chunks = glob.glob(os.path.join(sessdir, '*.npz'))
			chunks.sort()
			for chunk in chunks[::10]:
				print chunk
				chunk_ = os.path.splitext(os.path.split(chunk)[-1])[0] + '_guicalls.txt'
				if not os.path.exists(os.path.join(basedir, 'analysis', 'gui_count_output', chunk_)):
					print os.path.split(chunk)[-1]
					gui_count(chunk)

def run_manual_count():
	
	cagedirs = glob.glob(os.path.join(basedir, 'Cages', '*3pups*'))
	cagedirs.sort()
	
	fig = plt.figure(figsize = (15, 8))
	ax = fig.add_subplot(111)
	for cagedir in cagedirs:
		cage = os.path.split(cagedir)[1]
		print cage
		sessdirs = glob.glob(os.path.join(cagedir, 'P*'))
		sessdirs.sort()
		for sessdir in sessdirs:
			_, sess = os.path.split(sessdir)
			print sess
			chunks = glob.glob(os.path.join(sessdir, '*.npz'))
			chunks.sort()
			for chunk in chunks[::10]:
				print chunk
				chunk_ = os.path.splitext(os.path.split(chunk)[-1])[0] + '_mancalls.txt'
				if not os.path.exists(os.path.join(basedir, 'analysis', 'manual_count_output', chunk_)):
					print os.path.split(chunk)[-1]
					voc_score(chunk)
					plt.close()
	

def plot_blob(blob, Pxx, ax = None):
	
	ax, fig = misc.axis_check(ax)
	
	plot_spec(np.power(Pxx, 0.5), ax = ax)
	ax.set_xlim((max(blob['ij'][0][1] - 100, 0), blob['ij'][-1][1]+100))
	plt.show()



def guided_count(fname, ax, fig):
	
	Pxx, F, T = load_spec(fname)
	relat = os.path.splitext(os.path.split(fname)[-1])[0]
	blobs = np.loadtxt(os.path.join(basedir, 'analysis', 'power_search_output', relat + '_blobs.txt'), ndmin = 2)
	
	savefname = relat + '_autocalls.txt'
	
	plot_spec(np.power(Pxx, 0.5), ax = ax)
	hits = np.empty((0, 2), dtype = 'float32')
	
	skip_ct = 0
	i=0
	while (skip_ct <= 10) and i < len(blobs):
		blob = blobs[i]
		ax.set_xlim((max(blob[0] - 200, 0), blob[0] + 300))
		plt.show()
		xy = fig.ginput(n = 0, timeout = -1)
		i += 1
		if len(xy) > 0:
			skip_ct = 0
			hits = np.vstack((hits, xy))
			for xy_ in xy:
				ax.axvline(xy_[0], color = 'r', ls = '--')

		else:
			skip_ct += 1

		np.savetxt(os.path.join(basedir, 'analysis', 'guided_count_output', savefname), hits)

def run_guided_count():
	
	cagedirs = glob.glob(os.path.join(basedir, 'Cages/*'))
	cagedirs.sort()

	fig = plt.figure(figsize = (15, 8))
	ax = fig.add_subplot(111)
	for cagedir in cagedirs:
		if cagedir.find('3pups')<0:
			cage = os.path.split(cagedir)[1]
			print cage
			sessdirs = glob.glob(os.path.join(cagedir, 'P*'))
			sessdirs.sort()
			for sessdir in sessdirs:
				absol, sess = os.path.split(sessdir)
				print sess
				# chunks = glob.glob(os.path.join(sessdir, '*.npz'))
				# chunks.sort()
				i = 0
				while i < 10:
					chunk = os.path.join(basedir, 'Cages', cage, sess, '%s_%s_%2.2i.npz' % (cage, sess, i))
					chunk_ = '%s_%s_%2.2i_autocalls.txt' % (cage, sess, i)
					# chunk = chunks[i]
					# chunk_ = os.path.splitext(os.path.split(chunk)[-1])[0] + '_autocalls.txt'
					if not os.path.exists(os.path.join(basedir, 'analysis', 'guided_count_output', chunk_)):
						print os.path.split(chunk)[-1]
						try:
							guided_count(chunk, ax, fig)
							ax.cla()
						except:
							print 'Not found!'
					i += 1	
					
	

def compare_auto_vs_man(prefix, P, pplot = False):
	
	auto = np.loadtxt(os.path.join(basedir, 'analysis', 'guided_count_output', prefix + '_autocalls.txt'))
	man = np.loadtxt(os.path.join(basedir, 'analysis', 'manual_count_output', prefix + '_calls.txt'))
	
	auto_only = auto.copy()
	man_only = np.empty((0, 2), dtype = 'float32')
	for man_ in man:
		closest = np.abs(np.min(auto_only[:, 0] - man_[0]))
		closest_ix = np.argmin(auto_only[:, 0] - man_[0])
		if closest < 20:
			auto_only = np.vstack((auto_only[:closest_ix, :], auto_only[closest_ix+1:, :]))
		else:
			man_only = np.vstack((man_only, man_))
	
	if pplot:
		plot_compare_auto_vs_man(auto, man, auto_only, man_only, P)
		
	return auto_only, man_only

def plot_compare_auto_vs_man(auto, man, auto_only, man_only, P):
	
	fig1 = plt.figure()
	ax = fig1.add_subplot(111)

	ax.plot(auto[:, 0], np.ones(auto.shape[0]), 'g*')
	ax.plot(man[:, 0], np.zeros(man.shape[0]), 'b*')
	ax.set_ylim((-10, 10))
	
	
	# plot spectrograms of the missed calls
	fig2 = plt.figure()
	for i, auto_ in enumerate(auto_only):
		ax = fig2.add_subplot(2, auto_only.shape[0], i+1)
		plot_spec(P[:, auto_[0]-100:auto_[0]+1000], ax = ax)
		ax.axvline(100, color = 'r', ls = '--')
		ax.set_yticklabel('')
	for i, man_ in enumerate(man_only):
		ax = fig2.add_subplot(2, man_only.shape[0], auto_only.shape[0]+i+1)
		plot_spec(P[:, man_[0]-100:man_[0]+1000], ax = ax)
		ax.axvline(100, color = 'r', ls = '--')
		ax.set_yticklabels('')		
	fig2.tight_layout()
	
	print 'Manual:\t%i calls\nAuto:\t%i calls' % (man.shape[0], auto.shape[0])
	
	plt.show()
	

def quantify_all():
	cagedirs = glob.glob(os.path.join(basedir, 'Cages/*'))
	cagedirs.sort()

	wind = np.hamming(15)
	wind = wind / wind.sum()
	means = np.array([])
	stds = np.array([])
	sesss = np.array([])
	fhits = np.empty((0, 217))
	fmeans = np.empty((0, 217))
	hists = np.empty((0, 20))
	for cagedir in cagedirs:
		cage = os.path.split(cagedir)[1]
		print cage
		sessdirs = glob.glob(os.path.join(cagedir, 'P*'))
		sessdirs.sort()
		for sessdir in sessdirs:
			_, sess = os.path.split(sessdir)
			print sess
			chunks = glob.glob(os.path.join(sessdir, '*.npz'))
			chunks.sort()
			ix = np.random.randint(0, len(chunks)-1, 2)
			for i in ix:
				chunk = chunks[i]
				print os.path.split(chunk)[-1]
				P, F, T = load_spec(chunk)
				if P.shape[0] == 217:
					P[:, 153:156] = 0
					
					fhits_ = (P>0).sum(1) / np.float32(P.shape[1])
					fhits = np.vstack((fhits, fhits_))
					P_sum = P.sum(0)
					P_smoo = np.convolve(P_sum, wind, 'same')
					hist_ = plt.hist(P_sum, bins = 20, visible = False)
					hists = np.vstack((hists, hist_[0]))
					means = np.append(means, P_sum.mean())
					stds = np.append(stds, P_sum.std())
					sesss = np.append(sesss, os.path.split(chunk)[-1])
					fmeans = np.vstack((fmeans, P.sum(1)))

	return means, stds, sesss, fmeans

def export_DB(DB, fname = 'DB'):

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

	np.savetxt(os.path.join(basedir, fname+'.txt'), txt, fmt = '%s')




fnames = glob.glob(os.path.join(basedir, 'analysis', 'power_search_output/*'))
for fname in fnames:
	absol, relat = os.path.split(fname)
	r_split = relat.split('_')
	if len(r_split)==4:
		if len(r_split[2]) == 1:
			os.rename(fname, os.path.join(basedir, 'analysis', 'power_search_output', '%s_%s_%2.2i_%s' % (r_split[0], r_split[1], int(r_split[2]), r_split[3])))


# P_ = P[:, :20000]
# m, n = P.shape
# P_rav = P_.ravel(order = 'F')
# X = P_rav.reshape((10*m, 2000)).T
# pca = RandomizedPCA(n_components = 4)
# pca.whiten = True
# pca.fit(X)
# X_tran = pca.transform(X)
# 
# scatter(X_tran[:, 0], X_tran[:, 1])
# 
# kmeans = KMeans(k=4)
# kmeans.fit(X_tran)
# 
# scatter(X_tran[kmeans.labels_==0, 0], X_tran[kmeans.labels_==0, 1], color = 'b')
# scatter(X_tran[kmeans.labels_==1, 0], X_tran[kmeans.labels_==1, 1], color = 'g')
# scatter(X_tran[kmeans.labels_==2, 0], X_tran[kmeans.labels_==2, 1], color = 'r')
# 
# labels = kmeans.labels_
# c_ix = [[]]*4
# for i in range(4):
# 	c_ix[i] = (labels==i).nonzero()[0]
# c_ix = [(labels==0).nonzero()[0], (labels==1).nonzero()[0], (labels==2).nonzero()[0]]
# for c_ix_ in c_ix:
# 	fig = figure()
# 	for i in range(1, min(17, len(c_ix_))):
# 		ax = fig.add_subplot(4, 4, i)
# 		voc.plot_spec(X[c_ix_[i], :].reshape((m, 10)), ax = ax)
	
# def analyze_call_counts():
	
	

# dictmap = {'KO3' : {'7_18_2012' : 'P08', '7_19_2012' : 'P09'}, 'KO4' : {'7_27_2012' : 'P09', '7_30_2012' : 'P13', '7_31_2012' : 'P14', '8_1_2012' : 'P15', '8_2_2012' : 'P16', '8_6_2012' : 'P20'}, 'KO5' : {'7_25_2012' : 'P08', '7_26_2012' : 'P09', '7_27_2012' : 'P10', '7_30_2012' : 'P13', '7_31_2012' : 'P14', '8_2_2012' : 'P16', '8_6_2012' : 'P20', '8_7_2012' : 'P21'}, 'WT5' : {'8_8_2012' : 'P09', '8_9_2012' : 'P10'}, 'WT6' : {'8_13_2012' : 'P09', '8_14_2012' : 'P10', '8_16_2012' : 'P12', '8_17_2012' : 'P13'}}
# 
# fnames = glob.glob(os.path.join(basedir, 'analysis', 'power_search_output', '*'))
# for fname in fnames:
# 	absol, relat = os.path.split(fname)
# 	relat, exten = os.path.splitext(relat)
# 	if relat.find('P')<0:
# 		r_split = relat.split('_')
# 		cage = r_split[0].split('-')[0]
# 		sess = '_'.join(r_split[2:5])
# 		age = dictmap[cage][sess]
# 		num = r_split[1]
# 		newfname = os.path.join(absol, '_'.join((cage, age, num, 'blobs')) + '.txt')
# 		os.rename(fname, newfname)
# 

# '''Rename old naming system files to new one (with P-- in the filename)'''
# cagedirs = glob.glob(os.path.join(basedir, 'Cages/*'))
# 
# for cagedir in cagedirs:
# 	cage = os.path.split(cagedir)[1]
# 	print cage
# 	sessdirs = glob.glob(os.path.join(cagedir, 'P*'))
# 	for sessdir in sessdirs:
# 		_, sess = os.path.split(sessdir)
# 		print sess
# 		chunkpaths = glob.glob(os.path.join(sessdir, '*'))
# 		for chunkpath in chunkpaths:
# 			absol, chunk = os.path.split(chunkpath)
# 			chunk_, exten = os.path.splitext(chunk)
# 			chunk_split = chunk_.split('_')
# 			if len(chunk_split) == 5:
# 				newchunk_ = '_'.join((chunk_split[0], sess, '%2.2i' % int(chunk_split[1])))
# 				newchunk = newchunk_ + exten
# 				newchunkpath = os.path.join(absol, newchunk)
# 				print chunk, newchunk
# 				os.rename(chunkpath, newchunkpath)
# 			else:
# 				print 'Problem with %s' % chunk

# blob_starts = np.array([b['ij'][0][1] for b in blobs])
# plot(blob_starts, 1.1*np.ones(blob_starts.shape[0]), 'g.')
# plot(x[:, 0], np.ones(x.shape[0]), '.')
# 
# 
# cagedirs = glob.glob(os.path.join(basedir, 'Cages/*'))
# 
# for cagedir in cagedirs:
# 	cage = os.path.split(cagedir)[1]
# 	print cage
# 	sessdirs = glob.glob(os.path.join(cagedir, 'P*'))
# 	for sessdir in sessdirs:
# 		_, sess = os.path.split(sessdir)
# 		print sess
# 		chunkpaths = glob.glob(os.path.join(sessdir, '*'))
# 		for chunkpath in chunkpaths:
# 			absol, chunk = os.path.split(chunkpath)
# 			chunk_, exten = os.path.splitext(chunk)
# 			chunk_split = chunk_.split('_')
# 			if len(chunk_split[-1]) == 1:
# 				newchunk_ = '_'.join(chunk_split[:2]) + '_0' + chunk_split[-1] 
# 				newchunk = newchunk_ + exten
# 				newchunkpath = os.path.join(absol, newchunk)
# 				print chunk, newchunk
# 				os.rename(chunkpath, newchunkpath)
# 			else:
# 				print 'Problem with %s' % chunk


		
# fig = plt.figure()
# ax1 = fig.add_subplot(211)
# voc.plot_spec(np.power(P, 0.05), ax = ax1)
# ax2 = fig.add_subplot(212)
# ax2.plot(P_smoo)
# ax2.axhline(thresh, color = 'r', ls = '--')
# for start_x, start_y in zip(starts_x, starts_y):
# 	ax1.set_xlim((start_x-100, start_x+200))
# 	ax2.set_xlim((start_x-100, start_x+200))
# 	ax1.plot(start_x, start_y, 'r*', ms = 20)
# 	plt.show()
# 	fig.ginput(n = 0)

# Pxx, F, T = load_spec(fname)
# Pxx[153:156, :] = 0
# hits = Pxx>0
# # convolve to get rid of speckles
# hr = np.ones(3, dtype = 'float32') / 3.
# hits_filt = sepfir2d(hits, hr, hr)
# hits2 = hits_filt > 4./9
# 
# tmp = hits2.nonzero()
# hits_ij = zip(tmp[0], tmp[1])
# blobs__ = findclusters(hits_ij)
# 
# # eliminate small blobs
# blobs_ = [b for b in blobs__ if b['size']>10]
# 
# # eliminate narrow blobs (usually noise)
# blobs = []
# for blob in blobs_:
# 	y, x = map(list, zip(*blob['ij']))
# 	y_range = max(y) - min(y)
# 	x_range = max(x) - min(x)
# 	if (float(y_range) / x_range) < 2:
# 		blobs.append(blob)
# 
# Pxx_Tfs = Pxx.shape[1] / T[-1].astype('float32')
# Pxx_Ffs = Pxx.shape[0] / (F[-1] - F[0]).astype('float32')
# t_thresh = 0.02 * Pxx_Tfs
# f_thresh = 50000 * Pxx_Ffs
# 
# # combine close blobs
# if len(blobs) > 1:
# 	newblobs = [blobs[0]]
# 	for i in range(1, len(blobs)):
# 		blob0 = newblobs[-1]
# 		blob1 = blobs[i]
# 		if ((blob1['ij'][0][1] - blob0['ij'][-1][1]) < t_thresh) and \
# 			((blob1['ij'][0][0] - blob0['ij'][-1][0]) < f_thresh):
# 			blob0['ij'].extend(blob1['ij'])
# 			blob0['size'] += blob1['size']
# 		else:
# 			newblobs.append(blob1)
# else:
# 	newblobs = blobs
# 	
# # print '%i blobs' % len(newblobs)
# 
# savefname = os.path.split(os.path.splitext(fname)[0])[-1] + '_blobs.pickle'
# f = open(os.path.join(basedir, 'analysis', 'power_search_output', savefname), 'w')
# cPickle.dump(newblobs, f)
# f.close()