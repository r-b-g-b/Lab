import glob, os

subjs = ['Blue', 'Blue2', 'Orange', 'Orange2', 'Black', 'Black2']
basedir = '/Volumes/BOB_SAGET/Auditory maze'
# load and compare percentages for each session
fig = plt.figure();
axs = []
for i, subj in enumerate(subjs):
	print subj
	pcts = []
	fnames = glob.glob(os.path.join(basedir, subj + '_*.csv'))
	fnames.sort()
	for fname in fnames:
		sessname = os.path.splitext(os.path.split(fname)[-1])[0]
		print sessname
		data, info = Maze.load_session(sessname)
		Maze.plot_trials(data, info)
# 		data_ = Maze.remove_rest_periods(data)
# 		pct_data = Maze.calc_pct_correct(data)
# 		pct_data_ = Maze.calc_pct_correct(data_)
# 		pcts.append(np.hstack((pct_data, pct_data_)))
# 
# 	pcts = np.array(pcts)
# 	ax = fig.add_subplot(2, 3, i+1)
# 	ax.plot(pcts[:, 0])
# 	ax.plot(pcts[:, 3])
# 	axs.append(ax)
# 
# misc.sameyaxis(axs)

# process and save each session
for subj in subjs:
	fnames = glob.glob(os.path.join(basedir, subj + '_*.csv'))
	for fname in fnames:
		sessname = os.path.splitext(os.path.split(fname)[-1])[0]
		data, info = Maze.save_session(sessname)
