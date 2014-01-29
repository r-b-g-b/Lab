from ggplot import *
import inspect
import os, shutil
from glob import glob
import numpy as np
from scipy import stats as st
import pandas as pd
import Gap; reload(Gap)
import matplotlib as mpl
from matplotlib import pyplot as plt
import misc; reload(misc)
import pdb

color_cycle = mpl.rcParams['axes.color_cycle']
colors = dict(control='k', salicylate='r', prenihl='k', postnihl='g',
	postnihl5d='y', thalid='r', vehicle='b', preinjection='k',
	pre='k', tnfa='r', thalidwashout='y', vehiclewashout='m', ee='r', post='g')
lss = dict(control='-', salicylate='--', prenihl='-', postnihl='--',
	postnihl5d='--', thalid='--', vehicle='--', preinjection='-',
	pre='-', tnfa='--', thalidwashout='--', vehiclewashout='--', ee='--', post='--')

def agg_nanmean(ser):
	return np.mean(ser[~pd.isnull(ser)])

def agg_nansem(ser, axis=0):
	return st.sem(ser[~pd.isnull(ser)])

class GapAnalysis(object):

	def __init__(self, studyID=None, cageID=None, kind='gapdetection'):

		self.basedir = '/Volumes/BOB_SAGET/TNFalpha/tinnitus/behavior'
		self.kind = kind
		if (studyID is not None) and (cageID is not None):
			self._select_study(studyID, cageID)

		if not os.path.exists(os.path.join(self.studydir, 'Analysis')):
			os.mkdir(os.path.join(self.studydir, 'Analysis'))

		self.figdir = os.path.join(self.studydir, 'Analysis', kind)
		if not os.path.exists(self.figdir):
			os.mkdir(self.figdir)

		self.postdate_bins = [-100, 0, 1, 2, 4, 6, 8, 10, 12]
		self._load_experiment()

		# if a date, postdate1 is present in the dataframe, calculate the number of days
		# between each session and postdate1
		if hasattr(self.df, 'postdate1'):
			self.df = self._add_condition_postdate(self.df)
		
		notepath = os.path.join(self.studydir, 'notes.txt')
		if os.path.exists(notepath):
			with open(notepath) as f:
				self._notes = f.read()
		else:
			self._notes = None

		self.figformat = 'png'

	def pairwise_compare_conditions_by_day_by_freq(self, control='pre', bins=None):
		'''
		Plots pairwise differences (each animal's change in gapratio performance
		from their control data) by day as a function of frequency
		'''

		df = self.df

		if not hasattr(df, 'condition_postdate'):
			df = self._add_condition_postdate(df)

		# what post days are present
		udays = np.unique(df.postdate_bin)
		ndays = len(udays)-1 # control should have its own date bin

		'''
		Compute the control mean gap ratio for each animal/frequency
		'''
		# DataFrame with only the control data (all animals, freqs)
		df_control = df[df.condition==control]
		# group by animal and frequency
		animalgp_control = df_control.groupby(('animalID', 'freq'))
		# control mean gap ratio for each animal/frequency
		animalmeans_control = animalgp_control.agg(dict(gapratio=np.mean))

		'''
		Compute the manipulation mean gap ratio for each animal/postdate_bin/frequency 
		'''
		df_manip = df[df.postdate1>0]
		postdategp = df_manip.groupby(('animalID', 'freq', 'postdate_bin'))
		postdate_manip = postdategp.agg(dict(gapratio=np.mean)).unstack('postdate_bin')
		
		# make a dataframe, rows are animal/freq, columns are postdate_bin
		# loop through postdates, add a column with change in gap ratio for each one
		gapdiff = pd.DataFrame(index=animalmeans_control.index)
		for (key, df_) in postdate_manip.iteritems():
			gapdiff[key[1]] = df_ - animalmeans_control.gapratio

		condition_by_animal = df_manip.groupby('animalID').condition.first()

		cond = []
		for (key, value) in gapdiff.iterrows():
			cond.append(condition_by_animal[key[0]])
		gapdiff['condition'] = cond

		# add frequency column (from index)
		gapdiff['freq'] = gapdiff.index.get_level_values('freq')

		condgp = gapdiff.groupby(('condition', 'freq'))
		condmean = condgp.agg(agg_nanmean)
		condsem = condgp.agg(agg_nansem)

		ufreqs = np.unique(df.freq)
		x = range(len(ufreqs))

		fig, axs = plt.subplots(1, ndays, figsize=(12, 8))
		for i, ((ix, y_), (ix2, yerr_)) in enumerate(zip(condmean.iteritems(), condsem.iteritems())):
			assert ix==ix2
			for (cond, cond_y), (cond, cond_yerr) in zip(y_.unstack('condition').iteritems(), yerr_.unstack('condition').iteritems()):
				axs[i].errorbar(x, cond_y, yerr=cond_yerr, c=colors[cond], label=cond)
			axs[i].set_title('Postdate %s' % ix[1])

		[a.set_xticks(x) for a in axs]
		[a.set_xlim([-0.5, len(x)-0.5]) for a in axs]
		[a.set_ylim([-0.7, 0.7]) for a in axs]
		[a.set_xticklabels(np.int32(ufreqs/1000)) for a in axs]
		axs[0].set_ylabel('$\Delta$ Gap ratio')
		[a.set_xlabel('Frequency (kHz)') for a in axs]
		[a.axhline(0., c='k', ls='--') for a in axs]
		fig.suptitle('%s %s' % (self.studyID, os.path.split(self.studydir)[1]))
		axs[-1].legend(loc='best')

		fig.savefig(os.path.join(self.figdir, 'pairwise_compare_conditions_by_day_by_freq.%s' % self.figformat))

		plt.close(fig)
		return gapdiff


	def pairwise_compare_conditions_by_day(self, control = 'pre', bins=None):
		'''
		Plots change in gap ratio pairwise for each animal, pre- and post-manipulation
		
		The gap ratio is calculated as the average for each this animal, each frequency
		and then averaged for all frequencies.
		Change in gap ratio is the difference between this value for the pre- and post-
		manipulation days.
		
		bins : Post-manipulation data is averaged in bins of 2 days, by default, but arbitrary bin
			sizes can be used
		'''

		df = self.df
		df_control = df[df.condition==control]
		animalgp_control = df_control.groupby(('animalID', 'freq'))
		animalmeans_control = animalgp_control.agg(dict(gapratio=np.mean))

		df_manip = df[df.postdate1>0]
		# if bins is None:
			# bins = np.arange(df_manip.postdate1.min()-1, df_manip.postdate1.max()+1, 2)
		# df_manip['postdate_bin'] = pd.cut(df_manip.postdate1, bins=bins)
		postdategp = df_manip.groupby(('animalID', 'freq', 'postdate_bin'))
		postdate_manip = postdategp.agg(dict(gapratio=np.mean)).unstack('postdate_bin')
		
		gapdiff = pd.DataFrame(index=animalmeans_control.index)
		for (key, df_) in postdate_manip.iteritems():
			gapdiff[key[1]] = df_ - animalmeans_control.gapratio

		condition_by_animal = df_manip.groupby('animalID').condition.first()

		cond = []
		for (key, value) in gapdiff.iterrows():
			cond.append(condition_by_animal[key[0]])
		gapdiff['condition'] = cond

		condgp = gapdiff.groupby('condition')
		condmean = condgp.agg(agg_nanmean)
		condsem = condgp.agg(agg_nansem)

		styles = dict(tnfa=dict(color='r', hatch='///'), vehicle=dict(color='b', hatch=None))
		x = 2*np.arange(len(condmean.columns))
		width=0.8
		fig, ax = plt.subplots()
		for i,  ((key, y), (_, yerr)) in enumerate(zip(condmean.iterrows(), condsem.iterrows())):
			ax.bar(x-width+(i*width), y, yerr=yerr, width=width, fill=False, hatch=styles[key]['hatch'], ecolor='k', label=key)

		ax.legend(loc = 'best')
		ax.set_xticks(x)
		ax.set_xticklabels(condmean.columns)
		ax.set_xlabel('Days post injection')
		ax.set_ylabel('Change in gap detection ratio')
		ax.axhline(0, color = 'k')
		ax.set_title('%s %s' % (self.studyID, os.path.split(self.studydir)[1]))

		fig.savefig(os.path.join(self.figdir, 'pairwise_compare_condition_by_day.%s' % self.figformat))

		plt.close(fig)
		return gapdiff

	def pairwise_compare_nihl_effect_by_freq(self, control='prenihl'):
		
		df = self.df

		conditions = np.unique(df.condition).tolist()
		conditions.remove(control)

		# make animal/condition/freq groups
		animalgp = df.groupby(('animalID', 'condition', 'freq'))
		# take means of gap performance for those groups
		animalmeans = animalgp.agg(dict(gapratio=np.mean))

		# put the control and manipulation columns side-by-side
		condstack = animalmeans.gapratio.unstack('condition')

		manipdiff = pd.DataFrame(index=condstack.index)
		for cond in conditions:
			manipdiff[cond+'diff'] = condstack[cond] - condstack[control]

		manipdiffgp = manipdiff.groupby(level='freq')

		condmeans = manipdiffgp.mean()
		condsems = manipdiffgp.agg(st.sem)

		return condmeans, condsems

	def pairwise_compare_nihl_effect(self, control='prenihl'):
		
		df = self.df

		conditions = np.unique(df.condition).tolist()
		conditions.remove(control)

		# make animal/condition/freq groups
		animalgp = df.groupby(('animalID', 'condition', 'freq'))
		# take means of gap performance for those groups
		animalmeans = animalgp.agg(dict(gapratio=np.mean))

		# put the control and manipulation columns side-by-side
		condstack = animalmeans.gapratio.unstack('condition')

		manipdiff = pd.DataFrame(index=condstack.index)
		for cond in conditions:
			manipdiff[cond+'diff'] = condstack[cond] - condstack[control]

		manipdiff = manipdiff.stack(0)
		manipdiffgp = manipdiff.groupby(level=2)

		condmeans = manipdiffgp.mean()
		condsems = manipdiffgp.agg(st.sem)

		return condmeans, condsems


	def pairwise_compare_conditions_by_genotype_by_freq(self, control='pre'):
		
		df = self.df

		conditions = np.unique(df.condition).tolist()
		conditions.remove(control)

		# make animal/condition/freq groups
		animalgp = df.groupby(('animalID', 'condition', 'freq'))
		# take means of gap performance for those groups
		animalmeans = animalgp.agg(dict(gapratio=np.mean,
										gen=lambda x: x.iloc[0]))

		# put the control and manipulation columns side-by-side
		condstack = animalmeans.gapratio.unstack('condition')
		genstack = animalmeans.gen.unstack('condition')[control]

		manipdiff = pd.DataFrame(index=condstack.index)
		for cond in conditions:
			manipdiff[cond+'diff'] = condstack[cond] - condstack[control]

		manipdiff['gen'] = genstack
		manipdiff.reset_index(inplace=True)
		manipdiffgp = manipdiff.groupby(('gen', 'freq'))

		condmeans = manipdiffgp.mean().unstack('gen')
		condsems = manipdiffgp.agg(st.sem).unstack('gen')

		return condmeans, condsems

	def pairwise_compare_conditions_by_freq(self, control = 'prenihl', conditions = ['tnfa', 'vehicle']):

		df = self.df

		ucond = np.unique(df.condition)
		ncond = len(ucond)

		# make animal/condition/freq groups
		animalgp = df.groupby(('animalID', 'condition', 'freq'))
		# take means of gap performance for those groups
		animalmeans = animalgp.agg(dict(gapratio=np.mean))

		# put the control and manipulation columns side-by-side
		condstack = animalmeans.unstack('condition')
		
		# add a column indicating which condition this animal is
		animalcondition = []
		for key, value in condstack.iterrows():
			if not pd.isnull(value['gapratio']['tnfa']):
				animalcondition.append('tnfa')
			else:
				animalcondition.append('vehicle')
		condstack['animalcondition'] = animalcondition

		# make a column that is the result after manipulation
		# (combines multiple columns for each manipulation)
		# subtract the control from the manipulation values (for each frequency)
		manip = condstack.gapratio.vehicle.copy()
		manip.update(condstack.gapratio.tnfa)
		condstack['manip'] = manip
		condstack['manipdiff'] = condstack.manip - condstack.gapratio[control]

		# add a frequency column (turns frequency index into column)
		condstack['freq'] = condstack.index.get_level_values('freq')

		# group changes in gap detection by condition and frequency
		condgp = condstack.groupby(('animalcondition', 'freq'))
		condmeans = condgp.manipdiff.apply(np.mean).unstack('animalcondition')
		condsems = condgp.manipdiff.apply(st.sem).unstack('animalcondition')

		# line plot (change in gap detection VS frequency)
		ufreqs = np.unique(condstack.freq)
		x = range(len(ufreqs))
		fig, ax = plt.subplots()
		styles = dict(tnfa=dict(color='r', hatch='///'), vehicle=dict(color='b', hatch=None))
		for (key, y), (_, yerr) in zip(condmeans.iteritems(), condsems.iteritems()):
			misc.errorfill(x, y, yerr, ax=ax, label=key, color=styles[key]['color'])
		ax.set_xticks(x)
		ax.set_xticklabels((ufreqs/1000.).astype(int))
		ax.legend()

		fig.savefig(os.path.join(self.figdir, 'pairwise_compare_conditions_by_freq.%s' % self.figformat))

		plt.close(fig)

	def pairwise_compare_conditions(self, control = 'pre'):
		'''
		Performs a pairwise comparison between a control condition and several "manipulated"
		conditions where the manipulations were performed on DIFFERENT ANIMALS.
		
		This script plots the paired pre-post difference for each animal, and the manipulation
		group is indicated by different marker type/line styles
		'''

		df = self.df

		ucond = np.unique(df.condition)
		ncond = len(ucond)

		# make animal/condition/freq groups
		animalgp = df.groupby(('animalID', 'condition', 'freq'))
		# take means of gap performance for those groups
		animalmeans = animalgp.agg(dict(gapratio = np.mean))

		# put the control and manipulation columns side-by-side
		condstack = animalmeans.unstack('condition')
		
		# add a column indicating which condition this animal is
		animalcondition = []
		for key, value in condstack.iterrows():
			if not pd.isnull(value['gapratio']['tnfa']):
				animalcondition.append('tnfa')
			else:
				animalcondition.append('vehicle')
		condstack['animalcondition'] = animalcondition

		# make a column that is the result after manipulation
		# (combines multiple columns for each manipulation)
		# subtract the control from the manipulation values (for each frequency)
		manip = condstack.gapratio.vehicle.copy()
		manip.update(condstack.gapratio.tnfa)
		condstack['manip'] = manip
		condstack['manipdiff'] = condstack.manip - condstack.gapratio[control]

		# add a frequency column (turns frequency index into column)
		condstack['freq'] = condstack.index.get_level_values('freq')

		# calculate change collapsing over freq
		condgp = condstack.groupby('animalcondition')
		condmeans = condgp.manipdiff.apply(np.mean)
		condsems = condgp.manipdiff.apply(st.sem)

		# bar plot (change in gap detection, collapsing frequency)
		fig, ax = plt.subplots()
		for i, ((key, y), (_, yerr)) in enumerate(zip(condmeans.iteritems(), condsems.iteritems())):
			ax.bar(i+0.6, y, yerr=yerr, width=0.8, hatch=styles[key]['hatch'], fill=False, ecolor='k', label=key)
		ax.axhline(0, color='k')
		ax.set_title('%s %s' % (self.studyID, os.path.split(self.studydir)[1]))
		ax.legend(loc='upper center')
		ax.set_xlim([0, 3])

		fig.savefig(os.path.join(self.figdir, 'pairwise_compare_conditions.%s' % self.figformat))

		plt.close(fig)

	def group_compare_conditions_by_day(self, control='pre', close=False):
		'''
		Plot group averages for condition for each day bin
		'''

		df = self.df

		uconds = np.unique(df.condition)
		manip = [i for i in uconds if i!=control]

		if not hasattr(df, 'condition_postdate'):
			df = self._add_condition_postdate(df)

		udays = np.unique(df.postdate_bin)
		ndays = len(udays)-1 # control should have its own date bin

		animalgp = df.groupby(('animalID', 'condition', 'postdate_bin', 'freq'))
		animaldf = animalgp.agg(dict(gapratio = np.mean))
		condgp = animaldf.groupby(level = ('condition', 'postdate_bin', 'freq'))

		ufreqs = np.unique(df.freq)
		x = np.arange(ufreqs.size)
		y = condgp.agg(np.mean)
		yerr = condgp.agg(st.sem)
		ymean_control = y.ix[control]
		yerr_control = yerr.ix[control]
		ymean_manip = y.ix[manip].unstack('postdate_bin')
		yerr_manip = yerr.ix[manip].unstack('postdate_bin')

		fig, axs = plt.subplots(1, ndays, figsize=(12, 8))
		for i, ((ix, y_), (ix2, yerr_)) in enumerate(zip(ymean_manip.iteritems(), yerr_manip.iteritems())):
			assert ix==ix2
			for (cond, cond_y), (cond, cond_yerr) in zip(y_.unstack('condition').iteritems(), yerr_.unstack('condition').iteritems()):
				axs[i].errorbar(x, cond_y, yerr=cond_yerr, c=colors[cond], label=cond)
			axs[i].errorbar(x, ymean_control.values.flatten(), yerr=yerr_control.values.flatten(), c=colors[control], marker='.', label=control)
			axs[i].set_title('Postdate %s' % ix[1])

		[a.set_xticks(x) for a in axs]
		[a.set_xlim([-0.5, len(x)-0.5]) for a in axs]
		[a.set_ylim([0, 1.5]) for a in axs]
		[a.set_xticklabels(np.int32(ufreqs/1000)) for a in axs]
		axs[0].set_ylabel('Gap ratio')
		[a.set_xlabel('Frequency (kHz)') for a in axs]
		[a.axhline(1., c='r', ls='--') for a in axs]
		fig.suptitle('%s %s' % (self.studyID, os.path.split(self.studydir)[1]))
		axs[-1].legend(loc='best')

		figpath = os.path.join(self.figdir, 'group_compare_conditions_by_day.%s' % self.figformat)
		fig.savefig(figpath)

		if close:
			plt.close(fig)

		return fig, axs

	def group_compare_conditions_by_freq(self, conditions = None, ax = None, close=False):

		df = self.df

		animalgp = df.groupby(('animalID', 'condition', 'freq'))
		animaldf = animalgp.agg(dict(gapratio = np.mean))
		condgp = animaldf.groupby(level = ('condition', 'freq'))

		ufreqs = np.unique(df.freq)
		x = np.arange(ufreqs.size)
		y = condgp.agg(np.mean)
		yerr = condgp.agg(st.sem)

		fig = plt.figure(figsize = (10, 8))
		ax = fig.add_subplot(111)
		for j, ((i, y_), (i, yerr_)) in enumerate(zip(y.unstack().iterrows(), yerr.unstack().iterrows())):
			ax.errorbar(x, y_, yerr=yerr_, label = i, marker = '.', color='brc'[j])
		
		ax.set_xticks(x)
		ax.set_xlim([0, len(x)-1])
		ax.set_xticklabels(np.int32(ufreqs/1000))
		self._format_axis(ax)
		ax.set_title('%s %s' % (self.studyID, os.path.split(self.studydir)[1]))
		ax.legend(loc='best')

		figpath = os.path.join(self.figdir, 'group_compare_conditions.%s' % self.figformat)
		fig.savefig(figpath)

		if close:
			plt.close(fig)
		return fig, ax


	def group_compare_conditions(self, conditions = None, ax = None, close=False):

		df = self.df

		animalgp = df.groupby(('animalID', 'condition', 'freq'))
		animaldf = animalgp.agg(dict(gapratio = np.mean))
		condgp = animaldf.groupby(level = ('condition', 'freq'))

		ufreqs = np.unique(df.freq)
		x = np.arange(ufreqs.size)
		y = condgp.agg(np.mean)
		yerr = condgp.agg(st.sem)

		fig = plt.figure(figsize = (10, 8))
		ax = fig.add_subplot(111)
		for ((i, y_), (i, yerr_)) in zip(y.unstack().iterrows(), yerr.unstack().iterrows()):
			ax.errorbar(x, y_, yerr=yerr_, label = i, marker = '.')
		
		ax.set_xticks(x)
		ax.set_xlim([0, len(x)-1])
		ax.set_xticklabels(np.int32(ufreqs/1000))
		self._format_axis(ax)
		ax.set_title('%s %s' % (self.studyID, os.path.split(self.studydir)[1]))
		ax.legend(loc='best')

		figpath = os.path.join(self.figdir, 'group_compare_conditions.%s' % self.figformat)
		fig.savefig(figpath)

		if close:
			plt.close(fig)
		return fig, ax


	def single_subject_compare_conditions_by_day(self):
		'''
		For each subject, plot the mean gap ratio for each condition
		binned by date. (same as single_subject_compare_conditions, but
		uses binned date as an additional condition)
		'''
		fig, ax = plt.subplots(figsize = (10, 8))

		df = self.df

		for animalID in self.animalIDs:

			animaldf = df[df.animalID==animalID]
			condgp = animaldf.groupby(('condition_postdate', 'freq'))

			ufreqs = np.unique(animaldf.freq)
			x = np.arange(ufreqs.size)
			y = condgp.agg(dict(gapratio = np.mean))
			yerr = condgp.agg(dict(gapratio = st.sem))

			for ((i, y_), (i, yerr_)) in zip(y.unstack().iterrows(), yerr.unstack().iterrows()):
				misc.errorfill(x, y_, yerr_, label = i, ax = ax, marker = 'o')

			ax.set_xticks(x)
			ax.set_xticklabels(np.int32(ufreqs/1000))
			self._format_axis(ax)
			ax.set_title('%s %s' % (self.studyID, os.path.split(self.studydir)[1]))
			ax.legend()

			figpath = os.path.join(self.figdir, 'single_subject_compare_conditions_by_day_%s.%s' % (animalID, self.figformat))
			fig.savefig(figpath)
			ax.cla();

		plt.close(fig)		

	def single_subject_compare_conditions(self):
		'''
		For each subject, plot the condition means for each condition
		on the same plot for comparison.
		'''
		fig, ax = plt.subplots(figsize = (10, 8))

		for animalID in self.animalIDs:

			animaldf = self.df[self.df.animalID==animalID]
			condgp = animaldf.groupby(('condition', 'freq'))

			ufreqs = np.unique(animaldf.freq)
			x = np.arange(ufreqs.size)
			y = condgp.agg(dict(gapratio = np.mean))
			yerr = condgp.agg(dict(gapratio = st.sem))

			for ((i, y_), (i, yerr_)) in zip(y.unstack().iterrows(), yerr.unstack().iterrows()):
				# ax.errorbar(x, y_, yerr_, label = i, marker = 'o', color = colors[i], ls = lss[i])
				misc.errorfill(x, y_, yerr_, label = i, ax = ax, marker = 'o', color = colors[i], ls = lss[i])

			ax.set_xticks(x)
			ax.set_xticklabels(np.int32(ufreqs/1000))
			self._format_axis(ax)
			ax.set_title('%s %s' % (self.studyID, os.path.split(self.studydir)[1]))
			ax.legend()

			figpath = os.path.join(self.figdir, 'single_subject_compareconditions_%s.%s' % (animalID, self.figformat))
			fig.savefig(figpath)
			ax.cla();

		plt.close(fig)

	def single_subject_daily_results_by_condition(self):
		'''
		For each condition, plots the daily results and a condition mean
		'''

		df = self.df
		
		animals = np.unique(df.animalID)
		ufreqs = np.unique(df.freq)
		x = range(len(ufreqs))

		fig, ax = plt.subplots()
		for animal in animals:
		
			animaldf = df[df.animalID==animal]
			uconds = np.unique(animaldf.condition)

			for cond in uconds:
				conddf = animaldf[animaldf.condition==cond]
				df1 = conddf.pivot(index='freq', columns='sess', values='gapratio')
				# plot individual days
				for i, (sess, gapratio) in enumerate(df1.iteritems()):
					if i>6: ls='--'
					else: ls='-'
					ax.plot(x, gapratio, label=sess, ls=ls)
				
				# plot mean
				ax.plot(x, df1.mean(1), c='k', lw=2)
				ax.set_title('%s: %s' % (animal, cond))
				ax.set_xlim([0, len(x)-1]); ax.set_ylim([0., 1.4])
				ax.axhline(1., c='r', ls='--')
				ax.set_xticks(x)
				ax.set_xticklabels((ufreqs/1000.).astype(int))
				ax.set_xlabel('Frequency (kHz)'); ax.set_ylabel('Gap ratio')
				ax.legend(loc='best', fontsize='small')

				# save out figure
				fig.suptitle('%s %s' % (self.studyID, os.path.split(self.studydir)[1]))
				figname = 'single_subject_dailyresults_%s_%s.%s' % (animal, cond, self.figformat)
				fig.savefig(os.path.join(self.figdir, figname))

				ax.cla()

		plt.close(fig)

	def single_subject_startleampl(self):

		fig = plt.figure()
		animalIDs = pd.read_csv(os.path.join(self.studydir, 'dobs.csv'))['animalID'].values
		for animalID in animalIDs:
			fnames = glob(os.path.join(self.studydir, 'data', 'Gap', Gap.reverse_animalID(animalID), '*.txt'))
			nfnames = len(fnames)
			nrows, ncols = misc.get_subplot_grid(nfnames)
			for i, fname in enumerate(fnames):
				ax = fig.add_subplot(nrows, ncols, i+1)
				df = Gap.txt2pd(fname)
				Gap.plot_startleampl(df, ax)
				ax.set_ylim([0, 8])
				ax.set_title('%s %s' % (self.studyID, os.path.split(self.studydir)[1]))
			figpath = os.path.join(self.figdir, 'startleampl_%s.%s' % (animalID, self.figformat))
			fig.savefig(figpath)
			fig.clf();

		plt.close(fig)

	def test_pairwise():

		means_post, sems_post = _calc_pairwise_gapratio(df, ['postnihl'])
		means_tha, sems_tha = _calc_pairwise_gapratio(df, ['thalid1', 'thalidwashout1'])
		means_veh, sems_veh = _calc_pairwise_gapratio(df, ['vehicle1', 'vehiclewashout1'])

		means = means_post.join(means_tha).join(means_veh)
		sems = sems_post.join(sems_tha).join(sems_veh)
		col_order = ['diff_postnihl', 'diff_thalid1', 'diff_vehicle1', 'diff_thalidwashout1', 'diff_vehiclewashout1']
		mean = means[col_order]
		sems = sems[col_order]

		styles = dict(postnihl=dict(color='white', hatch=None),
			thalid1=dict(color='gray', hatch='\\\\'),
			thalidwashout1=dict(color='white', hatch='\\\\'),
			vehicle1=dict(color='gray', hatch='---'),
			vehiclewashout1=dict(color='white', hatch='---'))
		fig, ax = plt.subplots(figsize=(12, 4))
		x = np.arange(len(means))
		for i, ((key, y), (_, yerr)) in enumerate(zip(means.iteritems(), sems.iteritems())):
			label = key[5:]
			ax.bar(x+i/6., y, yerr=yerr, width = 1/6., ecolor='k', label=label, **styles[label])
		ax.axhline(0, color='k')
		ax.legend(loc='best', fontsize='small')


		def _calc_pairwise_gapratio(df, group, control = 'prenihl'):

			group1 = list([control])
			group1.extend(group)
			ix = [cond in group1 for cond in df.condition]
			df1 = df.ix[ix]
			means1 = _group_animal_condition_freq(df1)

			control_means = means1.gapratio.prenihl
			manip_means = means1.gapratio.filter(group)

			for k, v in manip_means.iteritems():
				manip_means['diff_%s' % k] = v-control_means

			cond_group = manip_means.filter(regex='diff_').groupby(level='freq')
			cond_means = cond_group.agg(np.mean)
			cond_sems = cond_group.agg(st.sem)

			return cond_means, cond_sems


		def group_animal_condition_freq(df):
			gp = df.groupby(('animalID', 'condition', 'freq'))
			return gp.agg(dict(gapratio=np.mean)).unstack('condition').dropna()




	def _format_axis(self, ax):
		
		ax.set_xlabel('Frequency (kHz)')
		ax.set_ylabel('Gap ratio')
		ax.set_ylim([0, 1.5])
		ax.axhline([1.0], color='r', ls='--')
		# ax.legend()

	def _add_condition_postdate(self, df, bins=None):

		if bins is None: bins = self.postdate_bins

		df['postdate_bin'] = pd.cut(df.postdate1, bins=bins)
		df['condition_postdate'] = df.condition+df.postdate_bin
		return df

	def _load_experiment(self, conditions = None, onlygood = True):
		'''

		'''
		if conditions is None:
			cond_patts = ['']
		else:
			if type(conditions) is str:
				conditions = [conditions]

			cond_patts = []
			for condition in conditions:	
				cond_patts.append('_'+condition+'_')

		fnames = []
		for cond_patt in cond_patts:
			if onlygood:
				patt = os.path.join(self.studydir, self.kind, '[A-Za-z]*%s*' % cond_patt)
			else:
				patt = os.path.join(self.studydir, self.kind, '*%s*' % cond_patt)
			fnames.extend(glob(patt))

		dfs = []
		if len(fnames)==0:
			print 'No files found with pattern\n%s' % patt
		else:
			for fname in fnames:
				df_ = pd.read_csv(fname)
				dfs.append(df_)
			df = pd.concat(dfs)

			df = df[np.vstack((5000<df.freq, df.freq<40000)).all(0)]
			if hasattr(df, 'postdate'):
				df = df[df.postdate1<10]
			df.index = range(len(df))
			self.df = df

	def _select_study(self, studyID, cageID):

		studydir = os.path.join(self.basedir, studyID, cageID)
		animalIDs = pd.read_csv(os.path.join(studydir, 'dobs.csv'))['animalID'] 
		animalIDs = [i for i in animalIDs if i[0]!='_']
		
		self.animalIDs = animalIDs
		self.studydir = studydir
		self.studyID = studyID
		self.cageID = cageID

	def _compute_gapdiff(self, control='pre'):

		df = self.df

		if not hasattr(df, 'condition_postdate'):
			df = self._add_condition_postdate(df)

		'''
		Compute the control mean gap ratio for each animal/frequency
		'''
		# DataFrame with only the control data (all animals, freqs)
		df_control = df[df.condition==control]
		# group by animal and frequency
		animalgp_control = df_control.groupby(('animalID', 'freq'))
		# control mean gap ratio for each animal/frequency
		animalmeans_control = animalgp_control.agg(dict(gapratio=np.mean))

		'''
		Compute the manipulation mean gap ratio for each animal/postdate_bin/frequency 
		'''
		df_manip = df[df.postdate1>0]
		postdategp = df_manip.groupby(('animalID', 'freq', 'postdate_bin'))
		postdate_manip = postdategp.agg(dict(gapratio=np.mean)).unstack('postdate_bin')

		condition_by_animal = df_manip.groupby('animalID').condition.first()
		
		# make a dataframe, rows are animal/freq, columns are postdate_bin
		# loop through postdates, add a column with change in gap ratio for each one
		gapdiff = pd.DataFrame(index=animalmeans_control.index)
		for (key, df_) in postdate_manip.iteritems():
			gapdiff[key[1]] = df_ - animalmeans_control.gapratio

		gapdiff.columns.name = 'postdate_bin'
		gapdiff = gapdiff.stack(level='postdate_bin')
		gapdiff = pd.DataFrame({'gapdiff': gapdiff})
		gapdiff.reset_index(inplace=True)

		# add condition by animal
		gapdiff['condition'] = gapdiff.animalID.map(condition_by_animal)

		self.gapdiff = gapdiff

	def describe_age(self):
		print 'Age---\n\tMean: %.2f\n\tMin: %i\n\tMax: %i\n' % (self.df.age.mean(), self.df.age.min(), self.df.age.max())

	def analyze(self):

		methods = inspect.getmembers(self, predicate=inspect.ismethod)
		methods = [i for i in methods if not i[0].startswith('_') and not i[0]=='analyze']
		for mname, meth in methods:
			try:
				meth()
				print mname
			except:
				print 'Could not run %s' % mname

	def notes(self):
		print self._notes
	
def gather_studies(basedir):
	
	IDs = []
	studypaths = glob(os.path.join(basedir, '[A-Za-z]*'))
	for studypath in studypaths:
		_, studyID = os.path.split(studypath)
		cagepaths = glob(os.path.join(studypath, '*'))
		for cagepath in cagepaths:
			_, cageID = os.path.split(cagepath)
			IDs.append((studyID, cageID))

	return IDs

def describe_ages(basedir):
	IDs = gather_studies(basedir)
	for studyID, cageID in IDs:
		gap = GapAnalysis(studyID, cageID)
		print studyID, cageID
		try:
			gap.describe_age()
		except:
			pass



