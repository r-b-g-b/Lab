import h5py, os, glob, re
import RF, RR, Spikes, misc
from fileconversion import load_cfs

basedir = '/Volumes/BOB_SAGET/Fmr1_voc/voc_ko_nai_20130116'
experiment = 'voc_ko_nai_20130116'

rf_paths = glob.glob(os.path.join(basedir, 'fileconversion', 'RF*.h5'))

rr_paths = glob.glob(os.path.join(basedir, 'fileconversion', 'RR*.h5'))
voc_paths = glob.glob(os.path.join(basedir, 'fileconversion', 'VOC*.h5'))
ix2freq = RF.get_ix2freq()
p = re.compile('RF(\d+).h5')

for rf_path in rf_paths:
	penno = p.findall(rf_path)[0]
	rr_path = [f for f in rr_paths if penno in f][0]
	voc_path = [f for f in voc_paths if penno in f][0]
	
	rf_file = h5py.File(rf_path, 'r')
	rf_rast = rf_file['rast'].value
	rf_stimparams = rf_file['stimID'].value
	rf_file.close()

	cfs = load_cfs(experiment)
	cf_ix = np.int32(np.round(RF.find_cf(cfs, np.int32(penno))))
	cf = ix2freq[20:][cf_ix]

	# perform analysis 
	if len(rr_path) > 0:
		
		rr_file = h5py.File(rr_path, 'r')
		rr_rast = rr_file['rast'].value
		rr_stimparams = rr_file['stimID'].value
		rr_file.close()
		
		# rr_rast = rr_rast[1:, :]
		# rr_stimparams = rr_stimparams[:-1, :]
		ufreqs = np.unique(rr_stimparams[:, 0])
		urrs = np.unique(rr_stimparams[:, 1])
		freq_played, freq_ix_played, _ = misc.closest(ufreqs, cf, log = True)
		rr_bins = np.arange(0, 6, 0.05)
		rr_psth_stim, usp = Spikes.calc_psth_by_stim(rr_rast, rr_stimparams, bins = rr_bins)
		cf_psth = rr_psth_stim[freq_ix_played, :, :]
		noise_psth = rr_psth_stim[0, :, :]
		
		fig = plt.figure()
		nrrs = cf_psth.shape[0]
		
		ax = []
		for i in range(nrrs):
			ax.append(fig.add_subplot(nrrs, 1, i+1))
			ax[-1].plot(rr_bins[:-1], cf_psth[i, :])
			RR.plot_tone_pips(urrs[i], 6, 0.05, 0.025, ax = ax[-1], color = 'r')
			if i>0:
				ax[-1].set_xticklabels('')
			
		misc.sameyaxis(ax)
	
	if len(voc_path) > 0:
		
		voc_file = h5py.File(voc_path, 'r')
		voc_rast = voc_file['rast'].value
		voc_stimparams = voc_file['stimID'].value
		voc_file.close()
		
		# voc_rast = voc_rast[1:, :]
		# voc_stimparams = voc_stimparams[:-1, :]
		voc_stimparams = voc_stimparams[:, 0]
		voc_stimparams = np.hstack((voc_stimparams[..., np.newaxis], np.zeros((voc_stimparams.size, 1))))
		uvocs = np.unique(voc_stimparams[:, 0])
		nvocs = uvocs.size
		
		voc_bins = np.arange(0, 20, 0.05)
		voc_psth_stim, usp = Spikes.calc_psth_by_stim(voc_rast, voc_stimparams, bins = voc_bins)
		
		fig = plt.figure()
		ax = []
		for i in range(nvocs):
			ax.append(fig.add_subplot(nvocs, 1, i+1))
			ax[-1].plot(voc_bins[:-1], voc_psth_stim[i, 0, :])
			if i>0:
				ax[-1].set_xticklabels('')
			
		misc.sameyaxis(ax)
		
