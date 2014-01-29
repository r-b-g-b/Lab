import numpy as np
import glob, os, re
import matplotlib as mpl
import matplotlib.pyplot as plt
import csv
import misc; reload(misc);
import cPickle as pickle

# basedir = '/Users/robert/Dropbox/behavior/Auditory maze'
basedir = '/Volumes/BOB_SAGET/Auditory maze'

def plot_heatmap(data, epoch = 'pre', ax = None, gridsize = 30, rel_to_target = False):
	'''
	plots a heatmap from x and y locations in one of the trial epochs.
	uses hexagonal bins
	'''
	xy = np.empty((0, 2), dtype = np.int32)
	x_key = epoch + '_x'
	y_key = epoch + '_y'
	if rel_to_target is True:
		for dat in data:
			xy = np.vstack((xy, np.vstack((dat[x_key] - dat['target_x'], dat[y_key] - dat['target_y'])).T))

	elif rel_to_target is False:
		for dat in data:
			xy = np.vstack((xy, np.vstack((dat[x_key], dat[y_key])).T))
	
	ax, _ = misc.axis_check(ax)
	
	ax.set_aspect('equal')
	ax.hexbin(xy[:, 0], xy[:, 1], gridsize = gridsize)
	
	if rel_to_target is False:
		misc.target(data[0]['spout_x'], data[0]['spout_y'], ax)
	if rel_to_target is True:
		misc.target(0, 0, ax)

	plt.show();
	
def plot_trial(dat, info, ax = None):
	
	ax, _ = misc.axis_check(ax)
	
	# plot the pre-track
	plot_track(dat['pre_x'], dat['pre_y'], ax = ax, color = 'b')
	ax.plot(dat['pre_x'][0], dat['pre_y'][0], marker = '*', ms = 20, color = 'b')

	# if there is a post-track, plot it
	if len(dat['post_x']) > 0:
		plot_track(dat['post_x'], dat['post_y'], ax = ax, color = 'g')
		ax.plot(dat['post_x'][0], dat['post_y'][0], marker = '*', ms = 20, color = 'g')

	# draw the target radius
	cir = mpl.patches.Circle((dat['target_x'], dat['target_y']), radius = info['target_radius'], color = 'r', fill = False, linestyle = 'dotted')
	ax.add_patch(cir)
	
	# draw the arena radius
	cir = mpl.patches.Circle((info['arena_x'], info['arena_y']), radius = info['arena_radius'], color = 'k', fill  = False)
	ax.add_patch(cir)
	misc.target(dat['target_x'], dat['target_y'], ax, color = 'r')
	misc.target(dat['spout_x'], dat['spout_y'], ax, color = 'g')
	
	# set the axis limits
	ax.set_xlim([info['arena_x']-info['arena_radius']-20, info['arena_x']+info['arena_radius']+20])
	ax.set_ylim([info['arena_y']-info['arena_radius']-20, info['arena_y']+info['arena_radius']+20])
	
	plt.show()

def plot_trials(data, info):

	ax = []
	fig = plt.figure(figsize = (11, 10));
	for i, dat in enumerate(data):
		subpltnum = np.mod(i, 25) + 1
		ax_ = fig.add_subplot(5, 5, subpltnum)
		ax_.set_xticklabels('')
		ax_.set_yticklabels('')
		plot_trial(dat, info, ax_)
		ax.append(ax_)
		if subpltnum == 25:
			fig.tight_layout();
			fig.savefig(os.path.join(basedir, 'analysis', 'paths', info['sess'] + '_trials_%2.2i.png' % np.int32(np.floor(i/25))))
			fig.savefig(os.path.join(basedir, 'analysis', 'paths', info['sess'] + '_trials_%2.2i.eps' % np.int32(np.floor(i/25))))
			fig = plt.figure(figsize = (11, 10));


		
def plot_track(x, y, ax = None, **kwargs):
	
	'''
	shows the path for one trial
	'''
	
	ax, _ = misc.axis_check(ax)
	ax.set_aspect('equal')
	
	ax.plot(x, y, '.-', **kwargs)
	plt.show()
	
def pre_vs_post_heatmap(data):
	'''
	heatmap of locations visited while searching for the target and after the target was located
	'''
	fig = plt.figure();
	ax1 = fig.add_subplot(211);
	ax2 = fig.add_subplot(212);
	plot_heatmap(data, epoch = 'pre', ax = ax1)
	plot_heatmap(data, epoch = 'post', ax = ax2)
	
def corr_vs_incorr_targets(data):
	
	'''
	shows the locations of targets that were correctly found versus targets that were missed
	'''
	
	target_corr = np.array([dat for dat in data if dat['target'] == 1])
	target_incorr = np.array([dat for dat in data if dat['target'] == 0])
	target_loc = np.array([dat for dat in data])
	
	fig = plt.figure();
	
	# correct trial target locations
	ax1 = fig.add_subplot(131);
	plot_heatmap(target_corr, epoch = 'target', ax = ax1, gridsize = 20)

	# incorrect trial target locations
	ax2 = fig.add_subplot(132);
	plot_heatmap(target_incorr, epoch = 'target', ax = ax2, gridsize = 20)
	
	# all trial target locations
	ax3 = fig.add_subplot(133);
	plot_heatmap(target_loc, epoch = 'target', ax = ax3, gridsize = 20)
	misc.sameyaxis([ax1, ax2, ax3])
	misc.samexaxis([ax1, ax2, ax3])
	
def calc_speed(x, y, t):
	
	return np.sqrt((np.diff(x)**2) + (np.diff(y)**2)) / (t[-1] - t[0])
	
def session_speed(data, epoch = 'pre'):
	'''
	returns the average speed for each trial
	'''
	x_key = epoch + '_x'
	y_key = epoch + '_y'
	t_key = epoch + '_t'
	speed = np.empty(len(data))
	for i, dat in enumerate(data):
		speed[i] = calc_speed(dat[x_key], dat[y_key], dat[t_key]).mean()
		
	return speed
	
def calc_pct_correct(data):
	
	ntrials = len(data)
	ncorrect = np.sum([dat['target'] for dat in data])
	nreward = np.sum([dat['reward'] for dat in data])
	
	pct_correct = 100*ncorrect / ntrials
	pct_reward = 100*nreward / ntrials
	pct_reward_of_correct = 100*nreward / ncorrect
	
	print '%30s\t%i' % ('no. trials:', ntrials)
	print '%30s\t%3.3f%%' % ('targets found:', pct_correct)
	print '%30s\t%3.3f%%' % ('reward found:', pct_reward)
	print '%30s\t%3.3f%%' % ('rewards per target found:', pct_reward_of_correct)
	print
	return pct_correct, pct_reward, pct_reward_of_correct
	
def save_session(sessname):

	tmp = sessname.split('_')
	subj = tmp[0]
	date = '_'.join(tmp[1:])
	data = []
	info = {}
	trials = []
	f = csv.reader(open(os.path.join(basedir, sessname + '.csv'), 'rU'), delimiter = ',')
	for l, line in enumerate(f):
		line = np.array(line, dtype = np.float32)
		if l == 0:
			info['sess'] = sessname
			info['arena_x'] = line[0]
			info['arena_y'] = line[1]
			info['arena_radius'] = line[2]
			info['trial_time'] = line[3]
			info['trial_hold'] = line[4]
			info['target_radius'] = line[5]
			info['reward_time'] = line[6]
			info['reward_hold'] = line[7]
			info['spout_radius'] = line[8]
			info['pip_rate'] = line[9]
			info['pip_freq'] = line[10]
			info['difficulty'] = line[11]
		else:

			data.append({})

			# performance (found target? got reward?)
			data[-1]['target'] = line[0]
			data[-1]['reward'] = line[1]
			# target and spout location and arena radius
			data[-1]['target_x'] = line[2]
			data[-1]['target_y'] = line[3]
			data[-1]['spout_x'] = line[4]
			data[-1]['spout_y'] = line[5]

			
			# table cell that separates pre from post epochs
			splitpoint_ = (line[7:] == -1).nonzero()[0]
			if splitpoint_.size > 0:
				splitpoint = splitpoint_[0]+7

				# how many points pre?
				pre_xy = line[7:splitpoint]
				pre_bins = pre_xy.size / 3.
				# how many points post?
				post_xy = line[(splitpoint+1):]
				post_bins = post_xy.size / 3.

				# save location and time pre and post
				data[-1]['pre_x'] = pre_xy[:pre_bins]
				data[-1]['pre_y'] = pre_xy[pre_bins:2*pre_bins]
				data[-1]['pre_t'] = pre_xy[2*pre_bins:3*pre_bins]
				data[-1]['post_x'] = post_xy[:post_bins]
				data[-1]['post_y'] = post_xy[post_bins : 2*post_bins]
				data[-1]['post_t'] = post_xy[2*post_bins : 3*post_bins]

			else:
				# how many points pre?
				pre_xy = line[7:]
				pre_bins = pre_xy.size / 3.
				# save location and time pre and post
				data[-1]['pre_x'] = pre_xy[:pre_bins]
				data[-1]['pre_y'] = pre_xy[pre_bins:2*pre_bins]
				data[-1]['pre_t'] = pre_xy[2*pre_bins:3*pre_bins]
				data[-1]['post_x'] = np.empty(0)
				data[-1]['post_y'] = np.empty(0)
				data[-1]['post_t'] = np.empty(0)


	outfname = os.path.join(basedir, os.path.splitext(os.path.split(sessname)[-1])[0] + '.p')
	
	db = {'data' : data, 'info' : info}
	pickle.dump(db, open(outfname, 'w'))

	return data, info


def load_session(sessname):
	
	db = pickle.load(open(os.path.join(basedir, sessname + '.p'), 'r'))
	return db['data'], db['info']


def remove_rest_periods(data, speed_thresh = 3, trial_thresh = 5):
	
	speeds = session_speed(data)
	
	# 1 : speed over threshold -- 0 : speed under threshold
	speeds_thresh = np.array(speeds > speed_thresh, dtype = np.float32)
	# more than 5 trials in a row with subthreshold speeds
	speeds_thresh_2 = np.array(np.convolve(speeds_thresh, np.ones(trial_thresh), 'same') > 0, dtype = np.float32)
	data_ = [dat for (dat, spee) in zip(data, speeds_thresh_2) if spee > 0]
	
	return data_


'''
Things to calculate

	1. Percent targets hit, rewards gotten
	
	1. speed as a function of trial (did the animal rest for any long periods of time?)
	2. average speed (do animals differ in this?)
	3. correct as a function of speed (do faster animals get more right?)
	
	Heatmaps
	1. correct target locations versus incorrect target locations

'''