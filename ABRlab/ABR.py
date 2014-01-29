import os, re
from glob import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt
import misc

class ABRAnalysis(object):
	'''
	A class dealing with analysis of ABR data.
	
	Initialize using:
	
		abr = ABR.ABRAnalysis(basedir)

	where 'basedir' is the full path to a folder containing these direcories
		'data' : containing raw text file outputs from SigGen
		'fileconversion' : made following ABRAanalysis.convert_all
			containing the text files converted to '.npz' files
		'Analysis' : where all analysis files will be saved
		'Figures' : where all analysis figures will be saved

	'''

	def __init__(self, basedir):

		self.basedir = basedir
		self.savedir= os.path.join(basedir, 'Analysis')
		
		self._checkDirectories()

		# load all of the ABR files in the fileconversion folder
		# self.fpaths = glob(os.path.join(self.basedir, 'fileconversion', '*.npz'))

		# dtype = np.dtype([('gen', 'S2'), ('exp', 'S3'), ('animal', 'i4'), ('freq', 'f4'), ('atten', 'f4'), ('lat', '5f4'), ('ampl', '5f4')])

		plt_opts = {'wt' : {'color' : 'b', 'x_offset' : 0, 'y_offset' : 0.00001}, 'ko' : {'color' : 'r', 'x_offset' : 0.09, 'y_offset' : 0}}

	def manual_peaks_all(self):

		fpaths = glob(os.path.join(self.basedir, 'fileconversion', '*.npz'))
		for fpath in fpaths:
			d = np.load(fpath)['arr_0']
			manual_peaks(d, ax)

	def manual_threshold_all(self):
		
				# if the path already exists, load it
		savepath = os.path.join(self.basedir, 'Analysis', 'thresholds.csv')
		if os.path.exists(os.path.join(self.basedir, 'Analysis', 'thresholds.csv')):
			print 'Thresholds already calculated and saved in %s' % savepath
			print 'Delete this file if you want to relabel the thresholds.'
			return

		fpaths = glob(os.path.join(self.basedir, 'fileconversion', '*.npz'))
		# shuffle the files to prevent experimenter bias
		fpaths = np.array(fpaths)
		fpaths = fpaths[np.random.permutation(fpaths.size)]
		
		for fpath in fpaths:
			d = np.load(fpath)['arr_0']
			self.manual_threshold(d)
	
	def manual_threshold(self, d):
		'''
		Takes a np.array of data for one ABR run (one animal), plots the
		average traces for each attenuation (frequencies separately) in a
		vertical stack for simple discernment of the ABR threshold.

		The user enters the actual attenuation (displayed in the legend) into
		a raw_input prompt in the terminal.

		Output is saved to a comma-separated value file
		'''

		# construct the output path
		savepath = os.path.join(self.basedir, 'Analysis', 'thresholds.csv')
		
		# else initialize threshold to be empty
		df = pd.DataFrame(columns=['gen', 'exp', 'animal', 'freq', 'thresh'])

		ufreqs = np.unique(d['freq'])
		nsamples = 244
		duration = 0.099424
		sample_rate = nsamples / duration
		t = np.arange(0, duration, 1./sample_rate)

		for freq in ufreqs:
			print 'Frequency: %f' % freq
			d_ = d[d['freq']==freq]
			
			nattens = d_.size
			y = self.filter_abr(d_['data'].T)

			paa = plotABRAttens(t, y, d_['atten'])
			thresh = paa.thresh
			# thresh = int(raw_input('Threshold--> '))
			
			df = df.append(dict(gen=d_['gen'][0],
				exp=d_['exp'][0],
				animal=d_['animal'][0],
				freq=freq,
				thresh=thresh), ignore_index=True)

		df.to_csv(savepath)

	def manual_peaks(self, D, ax):

		savepath = os.path.join(self.basedir, 'peaks.npz')
		dtype = np.dtype([('gen', 'S2'), ('exp', 'S3'), ('animal', 'i4'), ('freq', 'f4'), ('attens', '2f4'), ('peaks', '(2,5)f4')])

		if os.path.exists(savepath):
			peaks = np.load(savepath)['arr_0']
		else:
			peaks = np.empty(0, dtype = dtype)

		# ufreqs = np.unique(D['freq'])
		ufreqs = [8000., 16000.]
		
		nsamples = 244
		duration = 0.099424
		sample_rate = nsamples / duration
		t = np.arange(0, duration, 1./sample_rate)
		peak_times = (np.array([0.9, 1.6, 2.3, 3.2, 4.2]) / 1000.)+0.001
		for freq in ufreqs:
			
			d_ = D[D['freq']==freq]
			
			attens_ = d_['atten']
			atten = np.sort(d_['atten'])[np.array([0, 5])]
			d_ = d_[(np.tile(d_['atten'], (2, 1)) == np.tile(atten, (d_['atten'].shape[0], 1)).T).any(0)]
			
			nattens = d_.size

			y = d_['data'].T
			y = y - np.tile(y[0, :], (244, 1))
			y_offset = np.tile(np.linspace(0.00005, 0, nattens), (244, 1))
			
			y = y + y_offset
			
			# ax = fig.add_subplot(111)
			ax.plot(t, y)
			for i in range(len(peak_times)):
				ax.axvline(peak_times[i], color = 'r', ls = '--')
			ax.legend(d_['atten'])
			plt.show()
			
			xy = plt.ginput(n = -1, timeout = -1)
			
			ax.cla()
			
			peak_ = np.array(xy)[:, 0]
			peak_ = np.array(peak_)
			peak_.resize(peak_.size/5, 5)
			peak = np.empty((2, 5)) * np.nan
			peak[:peak_.shape[0], :] = peak_
			
			attens = np.empty(2) * np.nan
			attens[:nattens] = d_['atten']
			print peak
			peaks.resize(peaks.size+1)
			peaks[-1] = np.array((d_['gen'][0], d_['exp'][0], d_['animal'][0], freq, d_['atten'], peak), dtype = dtype)

		np.savez(savepath, peaks)

	def convert_all(self):
		
		fpaths = glob(os.path.join(self.basedir, 'data', '*.txt'))

		for fpath in fpaths:
			
			print fpath
			self.convert(fpath)

	def convert(self, fpath):

		# unpack the file path
		relat, absol = os.path.split(fpath)
		fname_, _ = os.path.splitext(absol)
		# get the animal info from the file path
		gen, exp, animal = fname_.split('_')
		
		# parse the ABR text output
		x = open(fpath)
		x = x.readlines()
		x = x[13:]
		p = re.compile('(\d+)')

		dtype = np.dtype([('gen', 'S2'), ('exp', 'S3'), ('animal', 'i4'), \
			('freq', 'f4'), ('atten', 'f4'), \
			('data', '244f4')])
		
		X = np.array([], dtype = dtype)
		done = False
		while len(x)>0:

			npts = np.int32(p.findall(x[3])[0])
			assert npts == 244
			freq = np.int32(p.findall(x[10])[0])
			atten = np.int32(p.findall(x[11])[0])
			data_ = x[12:12+npts]
			data = [float(d[:-2]) for d in data_]
		
			X.resize(X.size+1)
			X[-1] = np.array((gen, exp, animal, freq, atten, data), dtype = dtype)

			x = x[12+npts+1:]

		savepath = os.path.join(self.basedir, 'fileconversion', '_'.join((gen, exp, str(animal), 'ABR'))+'.npz')
		np.savez(savepath, X)
		
	def filter_abr(self, x, duration = 0.00999424, nsamples = 244, Wn = 0.02, btype = 'high'):

		b, a = butter(10, Wn = Wn, btype = btype)
		if len(x.shape)==1:
			x_filt = filtfilt(b, a, x)
		elif len(x.shape)>1:
			x_filt = np.empty_like(x)
			for i in range(x.shape[1]):
				x_filt[:, i] = filtfilt(b, a, x[:, i])
		return x_filt

	def plot_peak_gen_vs_exp(self, x, measure = 'ampl'):
		'''
		'''
		gens = ['wt', 'ko']
		exps = ['nai', 'exp']
		freqs = [8000., 16000.]

		plt_opts = {'wt' : {'color' : 'b', 'x_offset' : 0}, 'ko' : {'color' : 'r', 'x_offset' : 1}}

		fig = plt.figure(figsize = (14, 8))
		ax = []
		for i in range(10):
			ax.append(fig.add_subplot(2, 5, i+1))
			ax[-1].set_title('Peak %i' % ((i%5)+1))

		a = np.arange(5)
		for f, freq in enumerate(freqs):
			x2 = x[x['freq']==freq]
			attens = np.unique(x2['atten'])
			for atten in attens[:1]:
				x3 = x2[x2['atten']==atten]
			
				for gen in gens:
					ampl = np.empty((len(exps), 5))
					ampl_err = np.empty((len(exps), 5))
					for j, exp in enumerate(exps):
						x4 = x3[np.vstack((x3['gen']==gen, x3['exp']==exp)).all(0)]
						ampl[j, :] = (x4[measure]*10).mean(0)
						ampl_err[j, :] = (x4[measure]*10).std(0) / np.sqrt(x4.size)
					for i in range(ampl.shape[1]):
						ampl_ = ampl[:, i]
						ampl_err_ = ampl_err[:, i]
						ax[(f*5)+i].errorbar(np.arange(2), ampl_, yerr = ampl_err_, color = plt_opts[gen]['color'])

		misc.sameyaxis(ax)
		misc.samexaxis(ax, [-1, 2])
		for i in range(10):
			ax[i].set_xticks([0, 1])
			ax[i].set_xticklabels(exps)
		for i in range(1, 10):
			ax[i].set_yticklabels([])

		plt.show()	

	def plot_threshold_gen_vs_exp(self, x):
		'''
		'''
		gens = ['wt', 'ko']
		exps = ['nai', 'exp']
		freqs = np.unique(x['freq'])[:-1]

		plt_opts = {'wt' : {'color' : 'b', 'x_offset' : 0}, 'ko' : {'color' : 'r', 'x_offset' : 1}}

		fig = plt.figure(figsize = (14, 4))
		ax = []
		for i, freq in enumerate(freqs):
			ax.append(fig.add_subplot(1, 4, i+1))
			ax[-1].set_title('%u kHz' % (freq/1000))
			for j, gen in enumerate(gens):
				Y = np.empty(2)
				Y_err = np.empty(2)
				for k, exp in enumerate(exps):
					x_ = x[np.vstack((x['gen']==gen, x['exp']==exp, x['freq']==freq)).all(0)]
					y = x_['thresh_calib'].mean()
					y_err = x_['thresh_calib'].std() / np.sqrt(x_.size)
					
					Y[k] = y
					Y_err[k] = y_err
			
				ax[-1].errorbar([1, 2], Y, yerr = Y_err, color = plt_opts[gen]['color'])
		
		misc.sameyaxis(ax)
		misc.samexaxis(ax, [0, 3])
		for i in range(4):
			ax[i].set_xticks([1, 2])
			ax[i].set_xticklabels(exps)
		for i in range(1, 4):
			ax[i].set_yticklabels([])

		plt.show()	

	def calc_latency_delay(self):

		dtype = dtype = np.dtype([('gen', 'S2'), ('exp', 'S3'), ('animal', 'i4'), ('freq', 'f4'), ('atten', 'f4'), ('lat', '5f4'), ('lat_del', '5f4'), ('ampl', '5f4')])
		
		x2 = np.empty(0, dtype = dtype)
		for x_ in x:
			lat_del = np.diff(np.hstack((0, x_['lat'])))
			x2.resize(x2.size+1)
			x2[-1] = np.array((x_['gen'], x_['exp'], x_['animal'], x_['freq'], x_['atten'], x_['lat'], lat_del, x_['ampl']), dtype = dtype)

		np.savez(os.path.join(self.basedir, 'peaks.npz'), x2)

	def _checkDirectories(self):
		'''
		Makes sure the correct directory tree is in place within the base directory.
		'''
		subdirs = ['data', 'fileconversion', 'Analysis', 'Figures']
		for subdir in subdirs:
			directory = os.path.join(self.basedir, subdir)
			if not os.path.exists(directory):
				os.mkdir(directory)
				print 'Folder %s did not exist...' % subdir
				print 'Creating %s' % directory

class plotABRAttens(object):
	'''
	A simple interactive class for manually marking the
	ABR threshold. A figure showing the average ABR trace
	will appear. The user is expected to click on the red
	dot corresponding to the threshold dB SPL. If no
	threshold is present, the user can click the red dot
	at the bottom of the figure.

	Initialize using:
		t : time points of ABR samples
		abr : array of ABR traces (ntraces x nsamples)
		attens : the attenuation for each of the ntraces
		y_offset_mag : amplitude of offset for displaying
			multiple ABR traces on one plot

	'''
	def __init__(self, t, abr, attens, y_offset_mag=0.0003):


		plt.ion()
		# initialize figure
		self.fig, self.ax = plt.subplots(1, figsize=(12, 14))

		# incorporate attens
		self.attens = attens
		
		# offset abr by y_offset_mag for plotting
		nattens = abr.shape[1]
		nsamples = t.size
		y_offset = np.tile(np.linspace(y_offset_mag, 0, nattens), (nsamples, 1))
		abr = abr + y_offset

		# plot offset abr
		self.ax.plot(t, abr)
		# plot big red dots for UI
		self.ax.plot(np.zeros(nattens+1), np.hstack((abr[0, :], -y_offset_mag/nattens)), 'ro', ms=10, picker=5)
		self.ax.legend(self.attens)

		# draw
		plt.draw()
		plt.show()

		# connect to event handler
		self.fig.canvas.mpl_connect('pick_event', self.picker_event)
		# pause program execution for user input
		self.fig.canvas.start_event_loop(timeout=-1)

	def picker_event(self, event):
		if len(event.ind)>1:
			print 'More than one trace selected!\n'
			print

		try:
			print 'Atten: %s\n' % self.attens[event.ind[0]]
			self.thresh = self.attens[event.ind[0]]
		except IndexError:
			print np.nan
			self.thresh = np.nan

		self.fig.canvas.stop_event_loop()
		plt.close(self.fig)



