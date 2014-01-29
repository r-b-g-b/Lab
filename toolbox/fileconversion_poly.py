# fileconversion.py
# this should be a script that takes the blocks from
# fileconversion and makes them into python structures of a similar nature
# the blocks should be named with 
# prefix -- specifies the stimulus set used (e.g. b for RF or v for vocalizations)
# location ID -- a number indicating which location this was, continuous starting with 1
#		e.g. if you recorded two blocks in the same location but with different stimulus
#		sets, they should have different prefixes, but the same number
# .mat -- the file extension
# note: remove all other characters from the file name, e.g. b15.mat is OK, b15_R3.mat is NOT

#you should have the file0 in a folder and load in the matlab version of the file0 using pylab's
#matlab loader

import pylab
import scipy.io
import numpy as np
import glob
import os
import h5py
import Spikes; reload(Spikes)
import RF; reload(RF)
from pytables_classes import *
import re

sitemapdict = {'fourbyfour' : np.array([[1, 7, 13, 14], [3, 4, 10, 16], [2, 8, 12, 11], [6, 5, 9, 15]])}
prefix = {'b' : 'RF', 'r' : 'RR', 'q' : 'QU', 'v' : 'VOC', 'i' : 'IS'}
basedir = '/Volumes/BOB_SAGET/Fmr1_KO_ising/Sessions/good/'

nchanperblockdict = {'b' : 4, 'i' : 16}
# basedir = '/Users/robert/Desktop/tetrode_test'
# basedir = '/Volumes/BOB_SAGET/Fmr1_RR/Sessions'
# basedir = '/Volumes/BOB_SAGET/for_shaowen/'

def fileconversion(experiments, v = True):
	
	if type(experiments) is str:
		experiments = [experiments]
	
	# compile regexp patterns
	p_blocknum = re.compile('(?<=[a-zA-Z])[\d]+')
	
	for experiment in experiments:
		print experiment

		#get [b|r|etc]##.mat file names
		files = glob.glob(os.path.join(basedir, experiment, 'data', '[tbsyrpvi]*.mat'))
		nfiles = len(files)
		
		# make the destination directory if it doesn't exist (~/fileconversion/)
		if not os.path.exists(os.path.join(basedir, experiment, 'fileconversion')):
			os.mkdir(os.path.join(basedir, experiment, 'fileconversion'))

		blockinfo = []

		# parse block file name into
		# blockID - the block name
		# blocknum - the block number
		# blockrep - if multiple blocks were recorded at one site, which one was this
		for ff, file0 in enumerate(files):
			file0 = os.path.basename(file0)
			blockID = os.path.splitext(file0)[0] #removes extension
			tmp = blockID.split('_')
			blocknum = tmp[0]
			if len(tmp)==1:
				blockrep = 0
			else:
				blockrep = int(tmp[1])
			blocknum = int(p_blocknum.findall(blocknum)[0])
			blockinfo.append((blocknum, blockrep, blockID))
		blockinfo.sort()
	
		# load cfs and coords
		cfs = load_cfs(experiment, v = v)
		coords = load_coords(experiment, v = v)
		
		#loop through penetrations
		for blocknum, blockID, blockrep in blockinfo:

			fileconversion_block(experiment, blocknum, blockID, blockrep, coords = coords, cfs = cfs, v = v, nchanperblock = nchanperblockdict[blockID[0]], sitemap = sitemapdict['fourbyfour'])

		# end block loop

	# end experiment loop

def fileconversion_block(experiment, blocknum, blockID, blockrep, coords = None, cfs = None, v = True, nchanperblock = 4, sitemap = None):
	
	if cfs is None:
		cfs = load_cfs(experiment)
	if coords is None:
		coords = load_coords(experiment)
	
	tempfile = scipy.io.loadmat(os.path.join(basedir, experiment, 'data', blockID + '.mat')) # load the .mat file0
	Data0 = tempfile['Data'] # save the Data file0 to a temporary variable

	nchan = Data0['trial'][0][0][0][0]['CH'][0].size

	print '\t%s' % blockID

	blockstr = '%s%3.3i' % (prefix[blockID[0]], blocknum)
	print '\t%s' % blockstr

	savepath = os.path.join(basedir, experiment, 'fileconversion', blockstr + '.h5')
	b_ = h5py.File(savepath, 'w')

	for cc in range(nchan):
			
		# u_ = h5py.File(savepath, 'w')	
		u_ = b_.create_group
		unit(u_, Data0, cc, blockID)
		# try:
		# 	add_coords(u_, coords, unitnum)
		# 	add_cfs(u_, cfs, unitnum)
		# except:
		# 	pass
		u_.close()


def load_coords(experiment, v = True):
	
	try:
		coords = np.loadtxt(os.path.join(basedir, experiment, 'experimentfiles', experiment + '.txt'), ndmin = 1)
		if v:
			print 'Found coordinates at %s' % os.path.join(basedir, experiment, 'experimentfiles', experiment + '.txt')
	except:
		coords = np.nan
		if v:
			print 'Coordinates not found'
			
	return coords
	
def load_cfs(experiment, v = True):

	try:
		cfs = np.loadtxt(os.path.join(basedir, experiment, 'cfs.txt'), 'float32', ndmin = 1)
		if v:
			print 'Found CFs at %s' % os.path.join(basedir, experiment, 'cfs.txt')
	except:
		cfs = np.nan
		if v:
			print 'CFs not found'
	
	return cfs

def unit(u_, Data0, cc, blockID):
	
	# number of time bins to include in the LFP array
	nlfpsamp = 0
	for tt, trial in enumerate(Data0['trial'][0][0][0]):

		thislfpsamp = trial['LFP'].shape[1]
		if thislfpsamp>nlfpsamp:
			nlfpsamp = thislfpsamp
	
	
	ntrials = Data0['trial'][0][0][0].size # find number of trials
	nstimID = Data0['trial'][0][0][0][0]['Epoch_Value'][0].size
	
	# initialze LFP, spike times, spike trials, spike waveform
	lfp = np.ndarray((0, nlfpsamp), dtype = 'float32')
	spktimes = np.ndarray(0)
	spktrials = np.ndarray(0)
	spkwaveform = np.ndarray((0, 22))

	# initialize frequency and attenuation IDs
	stimID = np.ndarray((0, nstimID), dtype = 'float32')

	ttt = 0 # valid trial counter
	for tt in range(ntrials):
		trial = Data0['trial'][0][0][0][tt]

		thisstimID = np.float32(trial['Epoch_Value'][0])

		# if not ((blockID.startswith('b')) and (thisstimID[0] < 2)):

		# get the LFP for this trial and pad it with nans so it can fit in a matrix (since
		# some of the trials have +/-1 data point for LFP)
		lfpchannel = trial['LFP'][cc]
		lfpchannel = np.concatenate((lfpchannel, np.zeros(nlfpsamp - len(lfpchannel)) * np.nan))
		lfp = np.vstack((lfp, lfpchannel))

		spktime = trial['CH'][0][cc]['latency']
		if np.prod(spktime.shape) > 0:
			spktimes = np.append(spktimes, spktime)
			spktrials = np.append(spktrials, np.ones(spktime.size) * ttt)
			spkwaveform = np.concatenate((spkwaveform, trial['CH'][0][cc]['spkwaveform'].T), 0)
			
		# add to Epoch_Value
		stimID = np.vstack((stimID, thisstimID))

		ttt += 1 # increment valid trial counter

			# end if valid ID
	
	# end trial loop
	
	if spktimes.size == 0: # if no spikes
		print 'No spikes detected for this unit.'
		spktimes = np.array([np.nan])
		spktrials = np.array([np.nan])
		spkwaveform = np.array([np.nan])
		rast = np.array([np.nan])
		
	else:
		# filter out unwanted trials
		remID = np.array([np.nan, np.nan])
		if blockID.startswith('b'):
			remID = np.array([1., 70.])
		elif blockID.startswith('r'):
			remID = np.array([0., 0.])
	
		spktimes, spktrials, spkwaveform, lfp, stimID = \
			remove_trials(spktimes, spktrials, spkwaveform, lfp, stimID, remID)
	
		# create raster
		ntrials = stimID.shape[0]
		nbins = np.ceil(1000 * spktimes.max())+1
		rast = Spikes.calc_rast(spktimes, spktrials, ntrials, nbins)

	
	# save out to file
	u_.create_dataset('chan', data = cc)
	u_.create_dataset('blockID', data = blockID)
	# add stimulus ID datasets to this stimset on this unit
	u_.create_dataset('stimID', data = stimID)
	u_.create_dataset('lfp', data = lfp, compression = 'gzip')
	u_.create_dataset('spktimes', data = spktimes, compression = 'gzip')
	u_.create_dataset('spktrials', data = spktrials, compression = 'gzip')
	u_.create_dataset('spkwaveform', data = spkwaveform, compression = 'gzip')
	u_.create_dataset('rast', data = rast, compression = 'gzip')
		
	if blockID.startswith('b'):
		rf = RF.calc_rf(rast, stimID)
		u_.create_dataset('rf', data = rf, compression = 'gzip')
	
	
def add_coords(u_, coords, unitnum):
	
	try:
		coord = coords[coords[:, 0] == unitnum, 1:3][0]
	except:
		coord = np.nan
	u_.create_dataset('coord', data = coord)

	
def add_cfs(u_, cfs, unitnum):
	
	try:
		ix = cfs[:, 0] == unitnum
		cf = cfs[ix, 1][0]
	except:
		cf = np.nan
	u_.create_dataset('cf', data = cf)

def remove_trials(spktimes, spktrials, spkwaveform, lfp, stimID, remID):
	
	spk_mask, trial_mask = RF.make_spk_and_trial_masks(spktrials, stimID, remID)
	
	if spk_mask.sum()>0:
		spktimes = spktimes[~spk_mask]
		spktrials = spktrials[~spk_mask]
		spkwaveform = spkwaveform[~spk_mask, :]
		lfp = lfp[~trial_mask, :]
		stimID = stimID[~trial_mask, :]
	
	return spktimes, spktrials, spkwaveform, lfp, stimID
	
def sort_redos(d, v = False):
	'''
	Input:
		d - session main directory, e.g. /Users/robert/Desktop/Fmr1_RR/Sessions/ko_nai_20120306
	'''
	# get a list of the files
	os.chdir(os.path.join(basedir, d, 'data'))
	files = glob.glob('[b|r]*.mat')
	files.append('^^^') # special treatment for the last file
	p = re.compile('_R[\d]*')
	for f in range(len(files)-1):
		blockid1 = files[f][:3]
		blockid2 = files[f+1][:3]
		if blockid1 == blockid2:
			if v:
				print 'Renamed %s to %s' % (files[f], '_' + files[f])
			os.rename(files[f], '_' + files[f])
		else:
			if files[f].find('_R') > -1:
				if v:
					print 'Renamed %s to %s' % (files[f], p.sub('', files[f]))
				os.rename(files[f], p.sub('', files[f]))