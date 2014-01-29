#!/usr/bin/env

import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.widgets import Button
from scipy.ndimage import imread
from tkFileDialog import askopenfilename, askdirectory
from tkSimpleDialog import askstring

cmap = plt.get_cmap('jet')
def cf_to_rgba(i, N = 40):
	if i>N:
		return cmap(1.)
	elif i<0:
		return cmap(0.)
	else:
		return cmap(1.*i/N)

class Penetrations():
	
	def __init__(self, ax):
		
		self.pens = []
		self.cfs = []
		self.text = []
		self.npens = 0


		'''Get file path and display scopephoto'''
		impath = askopenfilename()

		img = imread(impath) # load scopephoto
		ax.imshow(img[::-1, :, :], origin = 'bottom') # show scopephoto
		
		# Set output directory for this experiment as the relative path of the scopephoto
		outputpath = os.path.split(impath)[0]

		'''Place the outputpath and note text on the axis'''
		xlim = np.array(ax.get_xlim())
		ylim = np.array(ax.get_ylim())
		
		colwidth = (xlim[1]-xlim[0]) / 20.
		rowwidth = (ylim[0]-ylim[1]) / 20.
		
		xloc = xlim[0] + np.abs(np.diff(xlim))/2.
		yloc = ylim[0] - rowwidth
		self.outputpath = ax.text(xloc, yloc, outputpath)
		
		yloc = ylim[0] - 2*rowwidth
		self.note = ax.text(xloc, yloc, 'Note')
		
		ax.figure.canvas.mpl_connect('button_press_event', self.button_press_callback)
		ax.figure.canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)
			
		ax.callbacks.connect('xlim_changed', self.axis_update)
		ax.callbacks.connect('ylim_changed', self.axis_update)
		
		self.ax = ax
		self.selected = None
		self.textselected = None
		
	def button_press_callback(self, event):
		# if we are trying to move a penetration
		if event.key is None:
			if self.selected is None:
				selected = [p.contains(event)[0] for p in self.pens]				
				if sum(selected)==1:
					self.selected = np.int32(np.array(selected).nonzero()[0])
					return
				else: return
			else:
				self.selected = None
			
		# if we are trying to relabel a penetration
		elif event.key=='z':
			selected = [p.contains(event)[0] for p in self.pens]
			if sum(selected)==1:
				self.textselected = np.int32(np.array(selected).nonzero()[0])
				self.edit_cf()
				return
			if self.note.contains(event)[0]:
				self.edit_note()
				return
			if self.outputpath.contains(event)[0]:
				self.edit_outputpath()
				return

		return
		
	def button_release_callback(self, event):
		self.selected = None
		return
		
	def motion_notify_callback(self, event):
		if event.inaxes!=self.ax: return
		if self.selected is None: return
		
		self.pens[self.selected].set_position((event.xdata, event.ydata))
		plt.draw()

		return

	def edit_cf(self):
		if self.textselected is None: return
		value = askstring('Assign CF', 'New CF:')
		try:
			freq = np.int32(value)
		except ValueError: return
		
		if freq==0:
			clr = (0, 0, 0)
		else:
			clr = cf_to_rgba(freq)
		self.cfs[self.textselected] = np.int32(freq)
		
		self.pens[self.textselected].set_color(clr)
		plt.draw()
		self.textselected = None

	def edit_outputpath(self):
		
		outputpath = askdirectory()
		print outputpath
		if os.path.exists(outputpath):
			self.outputpath.set_text(outputpath)
		else: return
		
		plt.draw()

	def edit_note(self):
		
		note = askstring('Note', '')
		self.note.set_text(note)
		
		plt.draw()


	def add_pen(self, event):

		# add number to the ScopePhoto
		# get current view
		xlim = self.ax.get_xlim()
		ylim = self.ax.get_ylim()
		colwidth = (xlim[1]-xlim[0]) / 20.
		rowwidth = (ylim[0]-ylim[1]) / 20.
		self.npens += 1
		x = xlim[0] - (rowwidth * ((self.npens%20)))
		y = ylim[1] - (colwidth * (np.floor(self.npens/20.)+1))
		self.pens.append(self.ax.text(x, y, np.str(self.npens)))
		self.cfs.append(0)
				
		plt.draw()
		return


	def rem_pen(self, event):
		# Remove CF label box
		if self.npens==0: return
		# Remove movable text on ScopePhoto
		del self.pens[-1]
		del self.cfs[-1]
		self.npens -= 1
		plt.show()
		return
		
	def axis_update(self, event):

		xlim = np.array(self.ax.get_xlim())
		ylim = np.array(self.ax.get_ylim())
		
		colwidth = (xlim[1]-xlim[0]) / 20.
		rowwidth = (ylim[0]-ylim[1]) / 20.
		
		xloc = xlim[0] + np.abs(np.diff(xlim))/2.
		yloc = ylim[0] - rowwidth
		self.outputpath.set_position([xloc, yloc])
		
		yloc = ylim[0] - 2*rowwidth
		self.note.set_position([xloc, yloc])

	def swap_e1(self, event):

		for i in xrange(0, self.npens-4, 4):
			colors = [p.get_color() for p in self.pens[i:i+4]]
			texts = [p.get_text() for p in self.pens[i:i+4]]

			self.pens[i].set_color(colors[1])
			self.pens[i+1].set_color(colors[0])


			self.pens[i].set_text(texts[1])
			self.pens[i+1].set_text(texts[0])

		plt.show()
		plt.draw()

	def swap_e2(self, event):

		for i in xrange(0, self.npens-4, 4):

			colors = [p.get_color() for p in self.pens[i:i+4]]
			texts = [p.get_text() for p in self.pens[i:i+4]]

			self.pens[i+2].set_color(colors[3])
			self.pens[i+3].set_color(colors[2])

			self.pens[i+2].set_text(texts[3])
			self.pens[i+3].set_text(texts[2])

		plt.draw()

	def save(self, event):
		
		if self.npens==0: return
		outputpath = self.outputpath.get_text()
		try:
			xy = np.empty((self.npens, 4), dtype = np.float32)
			for i in range(self.npens):
				xy_ = self.pens[i].get_position()
				pennum = int(self.pens[i].get_text())
				xy[i, :] = np.array([pennum, xy_[0], xy_[1], self.cfs[i]])
			np.savetxt(os.path.join(outputpath, 'markpens.txt'), xy)
			ax.get_figure().savefig(os.path.join(outputpath, 'scopephotograph.png'))
			print xy

		except IOError:
			print '%s does not exist!' % outputpath
			
	def load(self, event):
		loadpath = askopenfilename()
		tmp = np.loadtxt(loadpath)
		
		self.npens = tmp.shape[0]
		self.pens = []
		self.cfs = list(tmp[:, 3].astype(int))
		for (u, x, y, cf) in tmp:
			# add number to the ScopePhoto
			if cf==0:
				clr = (0, 0, 0)
			else:
				clr = cf_to_rgba(cf)
			self.pens.append(self.ax.text(x, y, str(int(u)), color = clr))
				
		plt.draw()
			
			

if __name__ == '__main__':
	
	# make figure
	fig = plt.figure(figsize = [7.175, 5.4])
	fig.subplots_adjust(bottom = 0)
	fig.subplots_adjust(left = 0)
	fig.subplots_adjust(right = 1)
	ax = fig.add_axes([0, 0.1, 1, 0.9])
	ax.set_xticks([])
	ax.set_yticks([])

	# make buttons
	ax_plus = fig.add_axes([0.09, 0.02, 0.06, 0.04])
	ax_minus = fig.add_axes([0.02, 0.02, 0.06, 0.04])
	ax_save = fig.add_axes([0.85, 0.02, 0.06, 0.04])
	ax_load = fig.add_axes([0.92, 0.02, 0.06, 0.04])
	ax_swap_e1 = fig.add_axes([0.3, 0.02, 0.06, 0.04])
	ax_swap_e2 = fig.add_axes([0.37, 0.02, 0.06, 0.04])

	# make program object
	pens = Penetrations(ax)

	# link buttons to callbacks
	b_plus = Button(ax = ax_plus, label = '+')
	b_plus.on_clicked(pens.add_pen)
	b_minus = Button(ax = ax_minus, label = '-')
	b_minus.on_clicked(pens.rem_pen)
	b_save = Button(ax = ax_save, label = 'Save')
	b_save.on_clicked(pens.save)
	b_load = Button(ax = ax_load, label = 'Load')
	b_load.on_clicked(pens.load)
	b_swap_e1 = Button(ax = ax_swap_e1, label = 'swap e1')
	b_swap_e1.on_clicked(pens.swap_e1)
	b_swap_e2 = Button(ax = ax_swap_e2, label = 'swap e2')
	b_swap_e2.on_clicked(pens.swap_e2)

	plt.show()