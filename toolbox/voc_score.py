from matplotlib.widgets import Button
import os, glob
import numpy as np
import matplotlib.pyplot as plt
from voc import plot_spec, load_spec, zscore_spec, lowpass_spec

plt.ion()
plt.rcParams['keymap.save'] = ''
plt.rcParams['keymap.yscale'] = ''
plt.rcParams['keymap.grid'] = ''
basedir = '/Volumes/BOB_SAGET/Vocalizations'
cages = None

class Score():

	def __init__(self, basedir, cages = None):
	
		self.hits = np.empty((0, 2), dtype = 'float32')
		self.counter = -1

		if cages is None:
			cages_ = glob.glob(os.path.join(basedir, 'Cages', '*'))
			cages = [os.path.split(c)[-1] for c in cages_]
		
		self.initiate_figure()
		
		fnames = []
		savepaths = []
		for cage in cages:
			print 'Cage: %s' % cage
			sessdirs = glob.glob(os.path.join(basedir, 'Cages', cage, 'P*'))
			print 'Sessions:'
			print sessdirs
			for sessdir in sessdirs:
				sess = os.path.split(sessdir)[-1]
				for i in np.arange(0, 60, 10): # set which files to do
					relat = '%s_%s_%2.2d' % (cage, sess, i)
					new_fname = os.path.join(basedir, 'Cages', cage, sess, relat+'.npz')
					new_savepath = os.path.join(basedir, 'analysis', 'gui_count_output', relat+'.txt')
					if not os.path.exists(new_savepath) and os.path.exists(new_fname):
						fnames.append(new_fname)
						savepaths.append(new_savepath)
					
		fnames.sort()
		savepaths.sort()
		for (fname, savepath) in zip(fnames, savepaths):
			print fname, savepath
		
		self.fnames = fnames
		self.savepaths = savepaths
		
	def keypress(self, event):
		
		old_xlim = self.ax1.get_xlim()
		new_xlim = old_xlim
		xrange = old_xlim[1] - old_xlim[0]
		xcenter = np.mean(old_xlim)
		if event.key == 'alt+d': #small forward
			self.ax1.set_xlim([old_xlim[0]+(xrange/4.), old_xlim[1]+(xrange/4.)])
		elif event.key == 'alt+s': #small back
			self.ax1.set_xlim([old_xlim[0]-(xrange/4.), old_xlim[1]-(xrange/4.)])
		elif event.key == 'alt+f': #large forward
			self.ax1.set_xlim([old_xlim[1], old_xlim[1]+xrange])
		elif event.key == 'alt+a': # large back
			self.ax1.set_xlim([old_xlim[0]-xrange, old_xlim[0]])

		elif event.key == 'alt+e': #zoom in
			self.ax1.set_xlim([xcenter-0.25*xrange, xcenter+0.25*xrange])
		elif event.key == 'alt+w': #zoom out
			self.ax1.set_xlim([xcenter-1.5*xrange, xcenter+1.5*xrange])

		elif event.key == 'alt+g':
			xy = plt.ginput(n = 0, timeout = -1)
			self.hits = np.vstack((self.hits, xy))
			for xy_ in xy:
				self.ax1.axvline(xy_[0], color = 'r', ls = '--')
				self.ax2.axvline(xy_[0], color = 'r', ls = '--')
			plt.show()
			np.savetxt(self.savepaths[self.counter], self.hits)
			
		elif event.key == 'alt+q':
			
			self.ax1.cla()
			self.ax2.cla()
			
			self.counter += 1
			
			print self.fnames[self.counter]
			print self.savepaths[self.counter]

			self.hits = np.empty((0, 2), dtype = 'float32')
			np.savetxt(self.savepaths[self.counter], self.hits)

			P, F, T = load_spec(self.fnames[self.counter])
			P, F = lowpass_spec(P, F)
			P = zscore_spec(P)
			P[P<0] = 0
			P = P**0.5

			wind = np.hamming(10)
			wind = wind / wind.sum()
			P_sum = (np.convolve(P.sum(0), wind, 'same'))

			nfreqs, ntimes = P.shape
			Tmax = T.max()
			print P.shape
			plot_spec(P[:, ::10], F, T[::10], ax = self.ax1)
			self.ax1.set_title(self.fnames[self.counter])
			self.ax1.set_xlim([0, 1.5])

			self.ax2.plot(T, P_sum)
			self.ax2.set_xlim([0, Tmax])
			
			plt.show()

		plt.draw()
		plt.show()
	
	def initiate_figure(self):
		print 'initiate_figure'
		
		fig = plt.figure(figsize = [15, 8])
		ax1 = fig.add_axes([0, 0.2, 1, 0.75])
		ax2 = fig.add_subplot(616)

		fig.subplots_adjust(bottom = 0.05)
		fig.subplots_adjust(left = 0)
		fig.subplots_adjust(right = 1)
		self.fig = fig
		self.ax1 = ax1
		self.ax2 = ax2
		self.fig.canvas.mpl_connect('key_press_event', self.keypress)

	if __name__=='__main__':
		Score()