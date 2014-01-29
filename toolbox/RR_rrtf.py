import os, glob, h5py, re
import numpy as np
import RF; reload(RF)
import RR; reload(RR)
import Spikes; reload(Spikes)
import misc; reload(misc)

verbose = True

ix2freq = 1000 * 2**(np.arange(0, 64)/10.)

basedir = '/Volumes/BOB_SAGET/Fmr1_voc/'
savepath = os.path.join(basedir, 'Analysis', 'RRTFs.npz')

dtype = np.dtype([('gen', 'S2'), ('exp', 'S3'), ('sess', 'S20'), ('unit', 'i4'), ('cf', 'f4'), ('rrtf', '7f4'), ('rrs', '7f4')])
db = np.empty(0, dtype = dtype)

sesspaths = glob.glob(os.path.join(basedir, 'Sessions', 'voc_*'))

for sesspath in sesspaths:
	sess = os.path.split(sesspath)[-1]
	_, gen, exp, sess = sess.split('_')
	rrpaths = glob.glob(os.path.join(basedir, 'Sessions', sesspath, 'fileconversion', 'RR*.h5'))
	cfs = np.loadtxt(os.path.join(basedir, 'Sessions', sesspath, 'cfs.txt'))

	for rrpath in rrpaths:
	
		absol, relat = os.path.split(rrpath)
		p = re.compile('RR')
		rfpath_ = p.sub('RF', relat)
		rfpath = os.path.join(absol, rfpath_)
		blockname = os.path.splitext(relat)[0]
				
		p = re.compile('(\d+).h5')
		unitnum = np.int32(p.findall(relat))[0] # unit number
		if verbose:
			print rrpath
			print rfpath

		# open the RR file, get rast and stimparams, then close it
		frr = h5py.File(rrpath, 'r')
		rast = frr['rast'].value
		stimparams = frr['stimID'].value
		frr.close()

		# if not np.isnan(spktimes[0]):
		cf = cfs[cfs[:, 0]==unitnum, 1][0]
		cf_hz = ix2freq[20:][int(cf)]
		freqs = stimparams[:, 0]
		rrs = stimparams[:, 1]
		ufreqs = np.unique(freqs)
		urrs = np.unique(rrs)
		nrrs = urrs.size
	
		# now we determine which of the frequencies we played is closest to this neuron's CF
		thisfreq, thisfreq_ix, thisfreq_err = misc.closest(ufreqs, cf_hz, log = True)
		if np.abs(thisfreq_err) > 0.2:
			print 'No close frequency found!'
		thisfreq = ufreqs[thisfreq_ix]
	
		# isolate the parts of the raster for this frequency and build a psth for each RR
		ix = RF.get_trials(stimparams, np.array([thisfreq, np.nan]))
		thisrast = rast[ix, :1050]
		thisstims = stimparams[ix, :]
		psths, ustims = Spikes.calc_psth_by_stim(thisrast, thisstims)
	
		rrtf = RR.calc_rrtf_all(thisrast, thisstims, thisfreq, urrs)
	
		db.resize(db.size+1)
		db[-1] = np.array((gen, exp, sess, unitnum, cf_hz, rrtf, urrs), dtype = dtype)

np.savez(savepath, db)