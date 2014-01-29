from ggplot import *
from scipy.stats import sem
import pandas as pd

cutgp = diamonds.groupby('cut')
fit = cutgp.price.mean()
fit_se = cutgp.price.apply(sem)
model = pd.DataFrame({'fit': fit, 'fit_se': fit_se, \
	'fit_se_max': fit+fit_se, 'fit_se_min': fit-fit_se})

# reorder levels
model = model.ix[['Fair', 'Good', 'Very Good', 'Premium', 'Ideal']]
model.index.name = 'cut'

model.reset_index(inplace=True)
model['x'] = range(len(model))

# plot it
p = ggplot(aes(x='x', y='fit', ymin = 'fit_se_min', \
	ymax='fit_se_max', color = 'cut'), data=model)
print p + geom_linerange()


from pandas import DataFrame
from numpy.random import randn
from ggplot import *

letters = ['a','a','a','b','b','b','c','c','c','d','d','d']*2
colors = ['r','g','b','r','g','b','r','g','b','r','g','b']*2
values1 = randn(24)
values2 = randn(24)
ymin = values1 - rand(24)
ymax = values1 + rand(24)
df = pd.DataFrame({'letter': letters, 'color': colors, \
    'ix': np.arange(24), \
    'value1': values1, \
    'ymin': ymin, \
    'ymax': ymax, \
    'value2': values2})


from pandas import DataFrame, series

d = {'x': np.arange(5), 'y': np.arange(5, 0, -1)}
df = DataFrame(d)

np.random.seed(0)
d_err = {'x': rand(5), 'y': rand(5)}
df_err= DataFrame(d_err)

# normal
df.plot()

df.plot(yerr=df_err)

df.x.plot(yerr=df_err)
df.y.plot(yerr=df_err)


from ggplot import *
import analyze_gap

gap = analyze_gap.GapAnalysis('thalidomide', 'WT')

import os
from glob import glob
studydir = '/Volumes/BOB_SAGET/TNFalpha/tinnitus/behavior/thalidomide/WT'
animalID = '112Black'
fpaths = glob(os.path.join(studydir, 'data', 'Gap', animalID, '*.txt'))
def plot_session_summary(fname):

for fpath in fpaths:
	df = Gap.txt2pd(fpath)
	df = df[df.freq>5000]
	df = df[df.freq<40000]
	df['trial'] = np.arange(len(df))
	df = Gap.add_startleampl(df)
	gg = ggplot(aes(x='trial', y='startleampl', color='cued'), data=df) \
	+ geom_point() \
	+ facet_wrap('freq', scales='fixed') \
	+ stat_smooth()
	gg.draw(); plt.draw();
	absol, relat = os.path.split(fpath)
	fname, _ = os.path.splitext(relat)
	figpath = os.path.join(studydir, 'Analysis', 'Sheets', fname+'.png')
	ggsave(figpath)
	plt.close('all')


# RAT VOCALIZATIONS

basedir = '/Volumes/BOB_SAGET/for_shaowen'
stimnames = np.loadtxt('/Volumes/BOB_SAGET/for_shaowen/stimulus/stimnames.txt', 'S')
for i, stimname in enumerate(stimnames):
	stimnames[i] = os.path.splitext(stimname.split('\\')[-1])[0]
	
nstims = stimnames.size


# make giant contact sheets of psths and rasters
plt.close('all')
voc_fnames = glob.glob('VOC*.h5')
onset_offset = np.loadtxt(os.path.join(basedir, 'sounds', 'onset_offset.txt'))
fig = plt.figure(figsize = (16, 20))
for voc_fname in voc_fnames:
	
	f = h5py.File(voc_fname, 'r')
	rast = f['rast'].value
	stimparams = f['stimID'].value
	f.close()

	rast_by_stim, ustimparams = Spikes.calc_rast_by_stim(rast, stimparams)
	psth, ustimparams = Spikes.calc_psth_by_stim(rast, stimparams)

	ax = []
	for i in range(nstims):
		
		onset = 50 + 1000*onset_offset[i, 0]
		offset = 50 + 1000*onset_offset[i, 1]
		ax_ = fig.add_subplot(10, 4, i+1)
		Spikes.plot_raster(rast_by_stim[:, i, 0, :], ax_)
		ax_.plot(5*Spikes.hamming_smoo(psth[i, 0, :])-5)
		ax_.set_title(stimnames[i])
		ax_.set_xlim([onset-50, offset+50])
		ax_.axvline(50+1000*onset_offset[i, 0], color = 'r', ls = '--')
		ax_.axvline(50+1000*onset_offset[i, 1], color = 'r', ls = '--')
		ax.append(ax_)
		if i>0:
			ax_.set_xticklabels('')
			ax_.set_yticklabels('')

	misc.sameyaxis(ax)
	fig.suptitle(voc_fname)
	# fig.tight_layout()
	plt.show()
	fig.savefig(os.path.join('/Volumes/BOB_SAGET/for_shaowen/analysis', '%s_psth_rast.png' % os.path.splitext(voc_fname)[0]))
	fig.clf();


stim_fnames = glob.glob('/Volumes/BOB_SAGET/for_shaowen/sounds/*.wavtxt')
for stim_fname in stim_fnames:
	x = np.loadtxt(stim_fname)
	h5_fname = os.path.splitext(stim_fname)[0] + '.h5'
	f = h5py.File(h5_fname, 'w')
	f.create_dataset('s', data = x, compression = 'gzip')
	f.create_dataset('fs', data = 200000)
	f.close()
	


# add spectrogram to sound databases
stim_fnames = glob.glob('/Volumes/BOB_SAGET/for_shaowen/sounds/*.h5')
NFFT = 512
noverlap = 256
for stim_fname in stim_fnames:
	f = h5py.File(stim_fname, 'r+')
	P, F, T = specgram(f['s'].value, Fs = f['fs'].value, NFFT = NFFT, noverlap = noverlap)
	f.create_dataset('P', data = P, compression = 'gzip')
	f.create_dataset('F', data = F, compression = 'gzip')
	f.create_dataset('T', data = T, compression = 'gzip')
	f.close()
	
# mark sound onset times
ax, fig = misc.axis_check(None)
# onset_offset = np.zeros((40, 2), np.float32)

for i, stimname in enumerate(stimnames):
	f = h5py.File(os.path.join(basedir, 'sounds', stimname + '.h5'), 'r')
	P = f['P'].value
	F = f['F'].value
	T = f['T'].value
	f.close()
	
	voc.plot_spec(P**0.3, F, T, ax = ax)
	ax.set_title(stimname)
	fig.savefig(os.path.join(basedir, 'sounds', stimname+'.png'))
	ax.cla()
	# x = plt.ginput(2, timeout = -1)
	# onset_offset[i, :] = np.array([x[0][0], x[1][0]])
	

straight_series = np.arange(29, 40)
c_series = np.arange(2, 11)
def calc_tuning_curve(psth, ixs):
	resp = np.zeros(ixs.size)
	for i, ix in enumerate(ixs):
		onset = 58 + 1000*onset_offset[ix, 0]
		# offset = 68 + 1000*onset_offset[ix, 1]
		offset = onset + 30
		resp[i] = psth[ix, 0, onset:offset].mean()
	
	return resp


# plot tuning curves for various stimulus series
ixs = c_series
voc_fnames = glob.glob('VOC*.h5')
# voc_fnames = ['VOC017.h5', 'VOC029.h5']
onset_offset = np.loadtxt(os.path.join(basedir, 'sounds', 'onset_offset.txt'))
plt.close('all')
fig = plt.figure()
ax = []
for i, voc_fname in enumerate(voc_fnames):
	
	print voc_fname
	f = h5py.File(voc_fname, 'r')
	rast = f['rast'].value
	stimparams = f['stimID'].value
	f.close()

	rast_by_stim, ustimparams = Spikes.calc_rast_by_stim(rast, stimparams)
	psth, ustimparams = Spikes.calc_psth_by_stim(rast, stimparams)

	ax_ = fig.add_subplot(4, 2, i+1)
	ax_.plot(calc_tuning_curve(psth, ixs))
	ax_.set_title(voc_fname)
	ax.append(ax_)
	
fig.suptitle('synthesized morph series')
fig.savefig(os.path.join(basedir, 'analysis', 'synth_morph_series.png'))

# make tuning sheet for one neuron
psth_all = np.zeros(0)
for i, ix in enumerate(ixs):
	onset = 58 + 1000*onset_offset[ix, 0]
	# offset = 68 + 1000*onset_offset[ix, 1]
	offset = onset + 30

	psth_ = psth[ix, 0, onset:offset]
	psth_all = np.hstack((psth_all, psth_, np.zeros(100-psth_.size)))


tuning_curve = calc_tuning_curve(psth, ixs)
x = np.arange(15, 15+9*100, 100)

fig = plt.figure(figsize = (10, 8))
ax1 = fig.add_axes([0.1, 0.5, 0.8, 0.4])
ax1.set_title(os.path.splitext(voc_fname)[0])
ax1.plot(x, tuning_curve, 'ro-')
ax1.set_ylabel('Mean firing rate (spikes / second)')
ax1.set_xticks(x)
ax1.set_xticklabels('')
ax1.set_xlim([x[0]-30, x[-1]+30])
ax2 = fig.add_axes([0.1, 0.1, 0.8, 0.4])
ax2.plot(Spikes.hamming_smoo(psth_all, windlen = 10))
ax2.set_xlabel('Stimulus ID')
ax2.set_ylabel('Firing rate (spikes / second)')
ax2.set_ylim([0, 1])
ax2.set_xticks(x)
ax2.set_xticklabels(stimnames[ixs])
ax2.set_xlim([x[0]-30, x[-1]+30])


# plot RFs
plt.close('all');
fig = plt.figure()
rf_fnames = ['RF017.h5', 'RF018.h5', 'RF019.h5', 'RF020.h5']
for i, rf_fname in enumerate(rf_fnames):
	f = h5py.File(rf_fname, 'r')
	rast = f['rast'].value
	stimparams = f['stimID'].value
	f.close()
	rf = RF.calc_rf(rast, stimparams)
	ax = fig.add_subplot(4, 2, i+1)
	RF.plot_RF(rf, ax = ax)
	ax.set_title(rf_fname)



dtype = np.dtype([('gen', 'S2'),  ('exp', 'S3'), ('sess', 'S20'), ('cf', 'f4'), ('rrtf', '7f4')])
db = np.empty(0, dtype = dtype)
for sess in f.iterkeys():
	gen, exp, _ = sess.split('_')
	for unit in f[sess].iterkeys():
		print f[sess][unit]['rrs']
		rrtf = f[sess][unit]['rrtf'].value
		if rrtf.size != 7:
			rrtf = np.append(rrtf, np.nan)
		
		db.resize(db.size+1)
		db[-1] = np.array((gen, exp, sess, f[sess][unit]['cf'].value, rrtf), dtype = dtype)

clrs = {'wt' : {'nai' : 'b-', 'exp' : 'b--'}, 'ko' : {'nai' : 'r-', 'exp' : 'r--'}}
gens = ['wt', 'ko']
exps = ['nai', 'exp']
rrs = np.array([1.5, 2, 3, 4, 6, 8, 12])
fig = plt.figure()
ax = fig.add_subplot(111)
leg = []
for gen, exp in itertools.product(gens, exps):
	clrs_ = clrs[gen][exp]
	data = db[np.vstack((db['gen']==gen, db['exp']==exp)).all(0)]['rrtf']
	ax.errorbar(rrs, data.mean(0), yerr = st.sem(data), color = clrs_[0], linestyle = clrs_[1:])
	leg.append(' '.join((gen, exp)))

ax.set_title('Repetition rate transfer function')
ax.set_xlabel('Repetition rate (pips/sec)')
ax.set_ylabel('RRTF')
ax.legend(leg)
plt.show()

ax = fig.add_subplot(122)
leg = []
for gen, exp in itertools.product(gens, exps):
	clrs_ = clrs[gen][exp]
	data = db[np.vstack((db['gen']==gen, db['exp']==exp, db['cf']>30)).all(0)]['rrtf']
	ax.errorbar(rrs, data.mean(0), yerr = sem(data), color = clrs_[0], linestyle = clrs_[1:])
	leg.append(' '.join((gen, exp)))

ax.set_title('RRTF for units outside of exposure range (4-32 kHz)')
ax.set_xlabel('Repetition rate (pips/sec)')
ax.set_ylabel('RRTF')
ax.legend(leg)
plt.show()

lfp = f['lfp'].value
rast = f['rast'].value
stimparams = f['stimID'].value
f.close()

lfp_stim = calc
fs_lfp = lfp.shape[1] / 4.2
P_lfp = 1./fs_lfp
lfp_t = np.arange(0, 4.2, P_lfp)

fig = plt.figure()
for i in range(psth_stim.shape[0]):
	for j in range(7):
		ax = fig.add_subplot(4, 7, (i-1)*7 + j)
		ax.plot(lfp_t, lfp_stim[i, j, :])
		ax.set_ylabel('%f %f' % (ustimparams[0][i], ustimparams[1][j]))

plt.show()
	
# do statistics by animal
sesss = np.unique(db['sess'])
rrtfs = {'wt_nai' : np.empty((0, 7)), 'wt_exp' : np.empty((0, 7)), \
'ko_nai' : np.empty((0, 7)), 'ko_exp' : np.empty((0, 7)),}
for sess in sesss:
	db_ = db[db['sess'] == sess]
	gen_ = db_[0]['gen']
	exp_ = db_[0]['exp']
	rrtf_ = db_['rrtf'].mean(0)
	key = gen_+'_'+exp_
	rrtfs[key] = np.vstack((rrtfs[key], rrtf_[np.newaxis, ...]))

dtype2 = np.dtype([('gen', 'S2'), ('exp', 'S3'), ('rr', 'f4'), ('rrtf', 'f4')])
db2 = np.empty(0, dtype = dtype2)
rrs = np.array([1.5, 2, 3, 4, 6, 8, 12])
for key in rrtfs:
	gen, exp = key.split('_')
	for rrtf in rrtfs[key]:
		for rr, rrtf_ in zip(rrs, rrtf):
			db2.resize(db2.size+1)
			db2[-1] = np.array((gen, exp, rr, rrtf_), dtype = dtype2)
	

clrs = {'wt' : {'nai' : 'b-', 'exp' : 'b--'}, 'ko' : {'nai' : 'r-', 'exp' : 'r--'}}
gens = ['wt', 'ko']
exps = ['nai', 'exp']

fig = plt.figure()
ax = fig.add_subplot(111)
leg = []
for gen, exp in itertools.product(gens, exps):
	clrs_ = clrs[gen][exp]
	data = rrtfs[gen+'_'+exp]
	ax.errorbar(rrs, data.mean(0), yerr = st.sem(data), color = clrs_[0], linestyle = clrs_[1:])
	leg.append(' '.join((gen, exp)))

ax.set_title('Repetition rate transfer functions (RRTFs)')
ax.set_xlabel('Repetition rate (pips/sec)')
ax.set_ylabel('RRTF')
ax.legend(leg)
path = misc.get_path()
fig.savefig(path, format = 'eps')

plt.show()


clrs = {'wt' : {'nai' : 'b-', 'exp' : 'b--'}, 'ko' : {'nai' : 'r-', 'exp' : 'r--'}}
exp = 'nai'
leg = []
fig = plt.figure()
ax = fig.add_subplot(111)

for gen in gens:
	clrs_ = clrs[gen][exp]
	data = db[db['gen']==gen]['rrtf']
	ax.errorbar(rrs, data.mean(0), yerr = st.sem(data), color = clrs_[0], linestyle = clrs_[1:])
	leg.append(gen)

ax.set_title('Repetition rate transfer functions (RRTFs)')
ax.set_xlabel('Repetition rate (pips/sec)')
ax.set_ylabel('RRTF')
ax.legend(leg)
path = misc.get_path()
fig.savefig(path, format = 'eps')

plt.show()

b1 = 2.
b2 = 3.
b12 = 4.
b3 = 5.

x0 = 10.
x1 = np.random.rand(100)
x2 = np.random.rand(100)
c = np.random.rand(100)

y = x0 + b1*x1 + b2*x2 + b12*x1*x2 + b3*c + np.random.rand(100)
f = open('/Users/robert/Desktop/test.csv', 'w')
f.write(','.join(('y', 'x1', 'x2', 'c')))
f.write('\n')
for i in range(100):
	f.write('%4.4f,%4.4f,%4.4f,%4.4f\n' % (y[i], x1[i], x2[i], c[i]))
f.close()
	

clrs = {'wt' : {'nai' : 'b-', 'exp' : 'b--'}, 'ko' : {'nai' : 'r-', 'exp' : 'r--'}}
gens = ['wt', 'ko']
exps = ['nai', 'exp']
leg = []
t = np.arange(-50, 333-50)
fig = plt.figure()
ax = fig.add_subplot(111)
for gen in gens:
	data = db_[db_['gen']==gen]['ev_psth']
	ax.errorbar(t, data.mean(0), yerr = st.sem(data, 0), color = clrs[gen]['nai'][0], ls = clrs[gen]['nai'][1:])
	leg.append(gen+'-nai')
ax.set_xlim((5, 35))
plt.show()


stimparams = f['stimID'].value
for i in range(16):
	rast = f['site%2.2i' % i]['rast'].value
	fig = plt.figure()
	ax = fig.add_subplot(111)
	Spikes.plot_sorted_raster(rast, stimparams, ax = ax)


dtype = np.dtype([('gen', 'S2'), ('exp', 'S3'), ('animal', 'i4'), ('freq', 'f4'), ('atten', 'f4'), ('peaks', '5f4')])

d2 = np.empty(0, dtype = dtype)
for d_ in d:
	gen = d_['gen']
	exp = d_['exp']
	animal = d_['animal']
	freq = d_['freq']
	for a in range(d_['attens'].size):
		atten = d_['attens'][a]
		peaks = d_['peaks'][a]
		d2.resize(d2.size+1)
		d2[-1] = np.array((gen, exp, animal, freq, atten, peaks), dtype = dtype)


savepath = os.path.join(basedir, fname+'.txt')

assert not os.path.exists(savepath)

f = open(savepath, 'w')
line = ''
dat = data[0]
for key_ in dat.dtype.descr:
	key = key_[0]
	newfield = dat[key]
	if newfield.size==1:
		newfield = [newfield]
	for i, newfield_ in enumerate(newfield):
		header = key
		if len(newfield)>1:
			header = '%s%2.2i' % (header, i+1)
		line += '%s\t' % header
line += '\n'
f.write(line)

for dat in data:
	line = ''
	for key_ in dat.dtype.descr:
		key = key_[0]
		newfield = dat[key]
		if newfield.size==1:
			newfield = [newfield]
		for newfield_ in newfield:
			if not type(newfield_) in [str, np.string_]:
				newfield_ = '%.5e' % newfield_
			line += '%s\t' % newfield_
	line += '\n'
	f.write(line)

f.close()

plot_opts = {'wt' : {'color' : 'b'}, 'ko' : {'color' : 'r'}, 'nai' : {'ls' : '-'}, 'exp' : {'ls' : '--'}}

peak_no = np.arange(1, 6)

gens = np.unique(d['gen'])
exps = np.unique(d['exp'])

fig = plt.figure()
ax = [[]]*2
ax[0] = fig.add_subplot(211)
ax[1] = fig.add_subplot(212)
for g, e in itertools.product(gens, exps):
	d_ = d[np.vstack((d['gen']==g, d['exp']==e)).all(0)]
	for a in range(2):
		dat = d_['peaks'][:, i, :]
		y = dat.mean(0)
		yerr = st.sem(dat, 0)
		ax[a].errorbar(peak_no, y, yerr = yerr, color = plot_opts[g]['color'], ls = plot_opts[e]['ls'])

plt.show()


ons = []
for db_ in db:
	on, off = Spikes.calc_on_off(db_['ev_psth'])
	ons.append(on)

basedir = '/Volumes/BOB_SAGET/fmr1_forRobert/naive/'
sess_fnames = glob.glob(os.path.join(basedir, '*'))
for sess_fname in sess_fnames:
	vor_fname = glob.glob(os.path.join(sess_fname, 'VOR2*.txt'))[0]
	vor = np.loadtxt(vor_fname)
	vor = vor[vor[:,0]>0]
	cf = vor[:, :3]
	np.savetxt(os.path.join(sess_fname, 'cf.txt'), cf)


f = h5py.File('/Users/robert/Documents/Work/Bao/Fmr1_Heesoo/RFs.h5')
basedir = '/Volumes/BOB_SAGET/fmr1_forRobert/naive/'
sess_fnames = glob.glob(os.path.join(basedir, '*'))
for sess_fname in sess_fnames:
	print sess_fname
	absol, sess = os.path.split(sess_fname)
	s = f.create_group(sess)
	pen_fnames = glob.glob(os.path.join(sess_fname, 'fileconversion', 'RF*.h5'))
	for pen_fname in pen_fnames:
		print pen_fname
		absol, pen = os.path.split(pen_fname)
		pen = os.path.splitext(pen)[0]
		frf_ = h5py.File(pen_fname, 'r')
		stimparams = frf_['stimID']
		# print stimparams.shape
		frf = frf_['site00']
		funit = s.create_group(pen)
		RF.add_rf_analysis(frf, funit, stimparams)
		frf_.close()

f.close()


fig = plt.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
basedir = '/Volumes/BOB_SAGET/fmr1_forRobert/naive/'
for s_ in f:
	s = f[s_]
	for u_ in s:
		u = s[u_]
		RF.plot_RF(u['rf'], ax = ax1)
		RF.plot_RF(u['rf_clust'], ax = ax2)
		title = '_'.join((s_, u_))
		ax1.set_title(title)
		fig.savefig(os.path.join(basedir, 'analysis', title + '.png'), format = 'png', dpi = fig.dpi)
		ax1.cla()
		ax2.cla()
				

fig = plt.figure(figsize = (16, 4))
ax = fig.add_subplot(111)
t = np.arange(333)
dat = db[db['gen']=='wt']['ev_psth']
misc.errorfill(dat, color = 'b', ax = ax)
dat = db[db['gen']=='ko']['ev_psth']
misc.errorfill(dat, color = 'r', ax = ax)

y = db['ev_psth'][:, 60:100]
y = y / np.tile(y.max(1), (40, 1)).T

pca = PCA(n_components = 2)
y_ = pca.fit_transform(y)

km = KMeans(k = 2)
km.fit(y_)

fig = plt.figure()
ax = fig.add_subplot(111)
t = np.arange(333)
dat = db[km.labels_==0]['ev_psth']
misc.errorfill(dat, color = 'b', ax = ax)
dat = db[km.labels_==1]['ev_psth']
misc.errorfill(dat, color = 'r', ax = ax)

peak_time = np.argmax(db['ev_psth'], 1)

fig = plt.figure()
ax = fig.add_subplot(111)
ix = np.argsort(db['resp_on'])
for i in ix[:20]:
	if db[i]['gen']=='wt':
		color = 'b'
	if db[i]['gen']=='ko':
		color = 'r'
	y = db[i]['ev_psth']
	y = y / y.max()
	wind = np.hanning(10)
	wind = wind / wind.sum()
	y = np.convolve(y, wind, 'same')
	ax.plot(y, color = color)


for i in ix[:3]:
	y = db_ko[i]['ev_psth']
	y = y / y.max()
	wind = np.hanning(10)
	wind = wind / wind.sum()
	y = np.convolve(y, wind, 'same')
	plot(y)



fig = plt.figure()
ax = fig.add_subplot(111)
for db_ in db:
	RF.plot_RF(db_['rf_clust'], bw_lr = db_['bw_lr'], ax = ax)
	title = '_'.join((db_['sess'], db_['unit']))
	ax.set_title(title)
	fig.savefig(os.path.join('/Users/robert/Documents/Work/Bao/Fmr1_RR/analysis/RFs', title+'.png'))
	ax.cla()
	




fig = plt.figure()
ax = fig.add_subplot(111)
for s_ in f:
	s = f[s_]
	for u_ in s:
		u = s[u_]
		RF.plot_RF(u['rf_clust'].value, bw_lr = u['bw_lr'].value, ax = ax)
		ax.plot(u['cf'].value, u['thresh'].value-2, '*', color = 'r', ms = 10)
		xy = plt.ginput()
		bw.append((s_, u_, xy))


bw = []
fig = plt.figure()
ax = fig.add_subplot(111)
for (s_, u_) in su:
		u = f[s_][u_]
		RF.plot_RF(u['rf_clust'].value, bw_lr = u['bw_lr'].value, ax = ax)
		ax.plot(u['cf'].value, u['thresh'].value-2, '*', color = 'r', ms = 10)
		plt.show()
		xy_ = plt.ginput(2)
		xy = (xy_[0][0], xy_[1][0])
		bw.append((s_, u_, xy))
		ax.cla()


for (s_, u_, bw_lr) in bw:
	u = f[s_][u_]
	u['bw'][2] = np.abs(bw_lr[1]-bw_lr[0])


basedir = '/Volumes/BOB_SAGET/Fmr1_voc/voc_ko_nai_20130116'
files = glob.glob(os.path.join(basedir, 'fileconversion', 'RR*.h5'))
for file0 in files:
	f = h5py.File(file0, 'r')
	stimparams = f['stimID'].value
	for siteno in range(4):
		sitekey = 'site%2.2u' % siteno
		rast = f[sitekey]['rast'].value
		rast = rast[:-1, :]
		stimparams = stimparams[2:, :]
		psth, usp = Spikes.calc_psth_by_stim(rast, stimparams, bins = np.arange(0, 6, 0.050))

		fig = plt.figure()
		for i in range(3):
			for j in range(7):
				ax = fig.add_subplot(7, 3, i+(j*3))
				ax.plot(psth[i, j, :])
		
	f.close()




voc_path = voc_paths[2]
for shift in range(-3, 3):	
	voc_file = h5py.File(voc_path, 'r')
	voc_rast = voc_file['rast'].value
	voc_stimparams = voc_file['stimID'].value
	voc_stimparams = voc_stimparams[:, 0]
	voc_stimparams = np.hstack((voc_stimparams[..., np.newaxis], np.zeros((voc_stimparams.size, 1))))
	voc_file.close()
	
	if shift>0:
		voc_rast = voc_rast[shift:, :]
		voc_stimparams = voc_stimparams[0:-shift, :]
	elif shift<0:
		voc_rast = voc_rast[:shift, :]
		voc_stimparams = voc_stimparams[(-shift):, :]
	elif shift==0:
		pass
	ufreqs = np.unique(voc_stimparams[:, 0])
	urrs = np.unique(voc_stimparams[:, 1])
	freq_played, freq_ix_played, _ = misc.closest(ufreqs, cf, log = True)
	voc_bins = np.arange(0, 20, 0.05)
	voc_psth, usp = Spikes.calc_psth_by_stim(voc_rast, voc_stimparams, bins = voc_bins)
	# cf_psth = voc_psth_stim[freq_ix_played, :, :]
	# noise_psth = voc_psth_stim[0, :, :]
	
	fig = plt.figure()
	nrrs = voc_psth.shape[0]
	
	ax = []
	for i in range(nrrs):
		ax.append(fig.add_subplot(nrrs, 1, i+1))
		ax[-1].plot(voc_bins[:-1], voc_psth[i, 0, :])
		# RR.plot_tone_pips(urrs[i], np.int32(urrs[i]*6), 0.05, 0.025, ax = ax[-1], color = 'r')
		if i<2:
			ax[-1].set_xticklabels('')
		
	misc.sameyaxis(ax)

ix2freq = 1000 * 2**(np.arange(0, 64)/10.)

sesspaths = glob.glob(os.path.join(basedir, 'Sessions', 'full_window', 'good', '*'))
for sesspath in sesspaths:
	rrpaths = glob.glob(os.path.join(sesspath, 'fileconversion', 'RF*.h5'))
	for rrpath in rrpaths:
		f = h5py.File(rrpath, 'r')
		rast = f['rast'].value
		stimparams = f['stimID'].value
		cf = f['cf'].value
		f.close()
		rf = RF.calc_rf(rast, stimparams)


dtype = np.dtype([('gen', 'S2'), ('exp', 'S3'), ('sess', 'S20'), ('cf', 'f4'), ('resp_on', 'f4'), ('resp_off', 'f4')])
db = np.empty(0, dtype = dtype)
for s_ in f:
	s = f[s_]
	gen, exp, _ = s_.split('_')
	for u_ in s:
		u = s[u_]
		db.resize(db.size+1)
		db[-1] = np.array((gen, exp, s_, u['cf'].value, u['resp_on'].value, u['resp_off'].value[0]), dtype = dtype)
		
d = dict(gen = db['gen'], exp = db['exp'], resp_on = db['resp_on'], resp_off = db['resp_off'])

df = pd.DataFrame(d)



ev_psth_norm = np.apply_along_axis(zscore, 1, db['ev_psth'])	


gens = ['wt', 'ko']
exps = ['nai', 'exp']
data_mean = np.empty((2, 2), dtype = float)
data_sem = np.empty((2, 2), dtype = float)
for g, gen in enumerate(gens):
	for e, exp in enumerate(exps):
		db_ = db[np.vstack((db['gen']==gen, db['exp']==exp)).all(0)]
		sesss = np.unique(db_['sess'])
		dat = []
		for sess in sesss:
			db__ = db_[db_['sess']==sess]
			dat.append(db__['resp_on'].mean())
		data_mean[g, e] = np.mean(dat)
		data_sem[g, e] = np.std(dat) / np.sqrt(len(dat))
			
			
fig = plt.figure()
rrs = np.array([1.5, 2, 3, 4, 6, 8, 12])
for i in range(7):
	ax = fig.add_subplot(7, 1, i+1)
	ax.set_xlim([0, 4])
	RR.plot_tone_pips(rrs[i], 6, 0.05, 0.025, ax = ax)
	ax.plot([0, 4], [0.05, 0.05], 'k')
	ax.set_ylim([0, 0.1])
	ax.set_yticks([])
	ax.set_xticks([])
	



sesspaths = glob.glob(os.path.join(basedir, 'voc_ko_*'))
for sesspath in sesspaths:
	rf_paths = glob.glob(os.path.join(sesspath, 'fileconversion', 'RF*.h5'))
	fig = plt.figure()
	fig.suptitle(os.path.split(sesspath)[-1])
	for i, rf_path in enumerate(rf_paths[:4]):
		ax = fig.add_subplot(2, 2, i+1)
		f = h5py.File(rf_path, 'r')
		RF.plot_rf(f, ax = ax)
		f.close()






# plot VOC psth contact sheets
basedir = '/Volumes/BOB_SAGET/Fmr1_voc'
sesspaths = glob.glob(os.path.join(basedir, 'voc_ko_nai*'))

bins = np.arange(0, 20.001, 0.001)
fig = plt.figure(figsize = (10.0125, 11.6875))
for sesspath in sesspaths:
	_, sess = os.path.split(sesspath)
	print sess
	voc_paths = glob.glob(os.path.join(sesspath, 'fileconversion', 'VOC*.h5'))
	rf_paths = glob.glob(os.path.join(sesspath, 'fileconversion', 'RF*.h5'))
	voc_paths.sort()
	rf_paths.sort()
	fpaths = zip(voc_paths, rf_paths)
	for i, (voc_path, rf_path) in enumerate(fpaths):
		figno = np.floor(i/5.)
		rowno = i%5

		f = h5py.File(voc_path, 'r')
		voc_rast = f['rast'].value
		voc_stimparams = f['stimID'].value
		f.close()
		voc_stimparams = voc_stimparams[:, 0][..., np.newaxis]
		voc_psths, usp = Spikes.calc_psth_by_stim(voc_rast, voc_stimparams, bins = bins)
		
		f = h5py.File(rf_path, 'r')
		rf_rast = f['rast'].value
		rf_stimparams = f['stimID'].value
		f.close()
		rf = RF.calc_rf(rf_rast, rf_stimparams)
		
		ax = fig.add_subplot(5, 4, (rowno*4)+1)
		RF.plot_RF(rf, ax = ax)
		ax.set_title(os.path.split(rf_path)[-1])
		
		for j, voc_psth in enumerate(voc_psths):
			colno = j+2
			ax = fig.add_subplot(5, 4, (rowno*4)+colno)
			voc_psth_smoo = Spikes.hamming_smoo(voc_psth, windlen = 20)
			assert voc_psth_smoo.max()<10
			ax.plot(bins[:-1], voc_psth_smoo)
			ax.set_ylim([0,5])
		if rowno==4 or i==len(fpaths):
			fig.savefig(os.path.join(basedir, 'Analysis', '%s_voc_psth_%i.png' % (sess, figno)))
			fig.clf()

fig = plt.figure()
ax = fig.add_subplot(432)
RF.plot_RF(rf, ax = ax)		
for j, voc_psth in enumerate(voc_psths):
	ax = fig.add_subplot(4, 1, j+2)	
	voc_psth_smoo = Spikes.hamming_smoo(voc_psth, windlen = 20)
	ax.plot(bins[:-1], voc_psth_smoo)
	ax.set_ylim([0, 5])
			
bins = np.arange(0, 4.501, 0.001)
fig = plt.figure()
# fig.tight_layout()
fig.set_facecolor('w')
ax = []
for i in range(7):
	ax.append(fig.add_axes([0, 1-(i/7.), 1, (1./7)]))
	# ax.append(fig.add_subplot(7, 1, i+1))
	psth_smoo = Spikes.hamming_smoo(psths[1, i, :], windlen = 20)
	ax[-1].plot(bins[:-1], psth_smoo)
	ax[-1].set_ylabel('%i' % rrs[i])
	ax[-1].set_xticklabels([])
	ax[-1].set_yticklabels([])
	ax[-1].set_axis_off()
misc.sameyaxis(ax)
misc.samexaxis(ax, [0, 4])
for i in range(7):
	ax_ = ax[i]
	RR.plot_tone_pips(rrs[i], 6, 0.050, 0.025, ax = ax_, color = 'r')
plt.draw()
plt.show()




# plot abr amplitudes

gens = ['wt', 'ko']
exps = ['nai', 'exp']

fig = plt.figure()
ax = fig.add_subplot(111)
for gen in gens:
	for exp in exps:

		data = x[np.vstack((x['gen']==gen, x['exp']==exp, x['freq']==8000, x['atten']==38)).all(0)]['ampl']
		ampl = data.mean(0)
		ampl_sem = data.std(0) / np.sqrt(data.shape[0])
			
		ax.errorbar(np.arange(1, 6), ampl, yerr = ampl_sem)
			
			
for path in fpaths:
	absol, relat = os.path.split(path)
	if relat.find(' ')>-1:
		shutil.move(path, os.path.join(absol, relat.replace(' ', '0')))



# Gens = np.empty(0, 'S2')
# Exps = np.empty(0, 'S3')
# Sesss = np.empty(0, 'S20')
# Unitnums = np.empty(0, np.int32)
# Cfs = np.empty(0, np.float32)
# Peak_times = np.empty(0, np.float32)
# On_times = np.empty(0, np.float32)
# Xloc = np.empty(0, np.float32)
# Yloc = np.empty(0, np.float32)

p = re.compile('(\d+)')
onset_time_ix = 50

dtype = np.dtype([('gen', 'S2'), ('exp', 'S3'), ('sess', 'S20'), ('unit', 'i4'), ('cf', 'f4'), ('peak_time', 'f4'), ('on_time', 'f4'), ('cf_on_time', 'f4'), ('x', 'f4'), ('y', 'f4'), ('psth', '333f4'), ('cf_psth', '333f4'), ('rf', '(8,43)f4')])

db = np.empty(0, dtype = dtype)


for experiment in experiments:
	rfpaths = glob.glob(os.path.join(basedir, experiment, 'fileconversion', 'RF*.h5'))
	cfs = np.loadtxt(os.path.join(basedir, experiment, 'cfs.txt'))
	# xys = np.loadtxt(os.path.join(basedir, experiment, 'xys.txt'))
	_, gen, exp, sess = experiment.split('_')
	for rfpath in rfpaths:
		absol, relat = os.path.split(rfpath)
		fname, ext = os.path.splitext(relat)
		unitnum = np.int32(p.findall(fname)[0])
		cf = cfs[cfs[:, 0]==unitnum, 1][0]
		# xy = xys[xys[:, 0]==unitnum, 1:][0]
		xy = [None, None]
		f = h5py.File(rfpath, 'r')
		rast = f['rast'].value
		stimparams = f['stimID'].value
		f.close()
		
		# ix = stimparams[:, 1]>10
		psth = rast[:, :333].sum(0)
		psth_smoo = Spikes.hamming_smoo(psth, windlen = 3)
		
		peak_time = (psth_smoo[55:100].argmax() + 5) / 1000.
		on_time, _ = Spikes.calc_on_off(psth_smoo)
		on_time = on_time / 1000.

		stim_psth, _ = Spikes.calc_psth_by_stim(rast[:, :333], stimparams)
		cf_int = np.int32(cf)
		
		cf_psth = stim_psth[cf_int-3:cf_int+2, :, :].sum(0).sum(0)
		cf_psth = cf_psth / cf_psth.max()
		cf_psth_smoo = Spikes.hamming_smoo(psth, windlen = 3)
		cf_on_time, _ = Spikes.calc_on_off(cf_psth_smoo)
		cf_on_time = cf_on_time / 1000.
		
		rf = RF.calc_rf(rast, stimparams)
		
		db.resize(db.size+1)
		db[-1] = np.array((gen, exp, sess, unitnum, cf, peak_time, on_time, cf_on_time, xy[0], xy[1], psth, cf_psth, rf), dtype = dtype)
		
		# Gens = np.hstack((Gens, gen))
		# Exps = np.hstack((Exps, exp))
		# Sesss = np.hstack((Sesss, sess))
		# Unitnums = np.hstack((Unitnums, unitnum))
		# Cfs = np.hstack((Cfs, cf))
		# Peak_times = np.hstack((Peak_times, peak_time))
		# On_times = np.hstack((On_times, on_time))
		# Xloc = np.hstack((Xloc, xy[0]))
		# Yloc = np.hstack((Yloc, xy[1]))
		
# df = pd.DataFrame(dict(gen = Gens, exp = Exps, sess = Sesss, unit = Unitnums, cf = Cfs, peak_time = Peak_times, on_time = On_times, x = Xloc, y = Yloc))
		


grouped = df.groupby('exp')
means = grouped.mean()
stds = grouped.std()
for m, s in zip(means.peak_time, stds.peak_time):
	print '%.5f (%.5f)' % (m, s)
		
		
for experiment in experiments:
	markpenspath = os.path.join(basedir, experiment, 'experimentfiles', experiment+'.txt')
	recsites = np.loadtxt(markpenspath)
	cfs = np.loadtxt(os.path.join(basedir, experiment, 'cfs.txt'))
 	goodunits = cfs[:, 0]
	recsites = recsites[np.array([r in goodunits for r in recsites[:, 0]]).nonzero()[0]]
	unitnums = recsites[:, 0]
	xloc = recsites[:, 1]
	yloc = recsites[:, 2]
	xmin = xloc.min()
	xmax = xloc.max()-xmin
	ymin = yloc.min()
	ymax = yloc.max()-ymin
	
	xloc = (xloc - xmin) / xmax
	yloc = (yloc - ymin) / ymax
	
	assert (xloc.min()==0) and (xloc.max()==1) and (yloc.min()==0) and (yloc.max()==1)
	
	xy = np.hstack((unitnums[..., np.newaxis], xloc[..., np.newaxis], yloc[..., np.newaxis]))

	np.savetxt(os.path.join(basedir, experiment, 'xy.txt'), xy)
	

basedir = '/Volumes/BOB_SAGET/Fmr1_Heesoo/A1/'
condpaths = glob.glob(os.path.join(basedir, '*'))
for condpath in condpaths:
	absol, cond = os.path.split(condpath)
	sesspaths = glob.glob(os.path.join(basedir, cond, '*'))
	for sesspath in sesspaths:
		_, sess = os.path.split(sesspath)
		shutil.move(sesspath, os.path.join(basedir, '_'.join((cond, sess))))
	
for i in range(158, 178):
	
	df_ = df.ix[i]
	rast, stimparams = misc.get_rast_stimparams_from_df(df_, unitprefix = 'RF')
	fig = RF.plot_fake_strf(rast, stimparams)
	fig.savefig(os.path.join(basedir, 'Analysis', 'fake_strf_%3.3i.png' % i))
	plt.close(fig)
	

def mysem(x):
	return x.std() / np.sqrt(x.size)
	

x = np.load('/Volumes/BOB_SAGET/Fmr1_RR/latency.npz')['arr_0']
df = misc.np_to_pd(x)
df['cf_bin'] = pd.cut(df.cf, [0, 20, 30, 40])

ganimal = df.groupby(['gen', 'exp', 'sess', 'cf_bin'])
animal_mean = ganimal.agg({'peak_time' : np.mean, 'on_time' : np.mean})
gcond = animal_mean.groupby(level = ['cf_bin', 'exp', 'gen'])
cond_mean = gcond.agg(np.mean)
cond_sem = gcond.agg(mysem)

cond_mean.plot(kind = 'bar')



'''


'''
for experiment in experiments:
	coh_info = coherence.calc_coherence_all(experiment)
	np.savez(os.path.join(basedir, 'Analysis', '%s_coh_info2.npz' % experiment), coh_info)
		

sessions = glob.glob(os.path.join(studydir, 'Sessions', 'voc_*'))
for session in sessions:
	print session
	fpaths = glob.glob(os.path.join(studydir, 'Sessions', session, 'fileconversion', '*.h5'))
	for fpath in fpaths:
		f = h5py.File(fpath, 'r')
		absol, relat = os.path.split(fpath)
		if len(f['coord'].shape)==0:
			print fpath
			newfpath = os.path.join(absol, 'new_%s' % relat)
			newf = h5py.File(newfpath)

			for key in f.iterkeys():
				if not key=='coord':
					f.copy(key, newf)
				else:
					newf.create_dataset('coord', data = np.array((np.nan, np.nan)))

			newf.close()
			f.close()
			shutil.move(fpath, os.path.join(absol, 'bk_%s' % relat))
			shutil.move(newfpath, os.path.join(absol, relat))				

		else:
			f.close()



	

		
sessions = ['voc_ko_exp_20130128', 'voc_ko_exp_20130130', 'voc_ko_exp_20130201', 'voc_ko_nai_20130116', 'voc_ko_nai_20130117', 
'voc_ko_nai_20130118', 'voc_ko_nai_20130122']

sessions = ['voc_wt_nai_20130305','voc_wt_nai_20130307',
'voc_wt_nai_20130308','voc_wt_nai_20130315',
'voc_wt_nai_20130321','voc_wt_nai_20130325']


for session in sessions:
	print session
	x = np.loadtxt(os.path.join(studydir, 'Sessions', session, 'experimentfiles', '%s.txt' % session))
	print x
	x[:, 2] = 480 - x[:, 2]
	np.savetxt(os.path.join(studydir, 'Sessions', session, 'experimentfiles', 'new_%s.txt' % session), x)
		
		
		
class C(object):
    def __init__(self):
        self._x = 'sweet'

    @property
    def x(self):
    	return self._x

    def __len__(self):
    	return 4
    # def getx(self):
    #     return self._x
    # def setx(self, value):
    #     self._x = value
    # def delx(self):
    #     del self._x
    # x = property(getx, setx, delx, "I'm the 'x' property.")
		
		
class C(object):
    def __init__(self):
        self._x = None

    def getx(self):
    	print 'getting x you bitch'
        return self._x
    def setx(self, value):
        self._x = value
    def delx(self):
        del self._x
    x = property(getx, setx, delx, "I'm the 'x' property.")		
		
'''Load in translation between attenuation and db SPL'''		
tmp = np.loadtxt('/Users/robert/Dropbox/FMRP/JN_3/ABR_stuff/threshold_key.csv', 'str', delimiter = ',')
colnames = tmp[0, 1:]
rownames = tmp[1:, 0]
x = tmp[1:, 1:]
d = dict()
for i in xrange(len(colnames)):
	d[colnames[i]] = np.int32(np.vstack((rownames, x[:, i])))

'''Translate atten to db SPL'''
x = pd.read_csv('/Users/robert/Dropbox/FMRP/JN_3/ABR_stuff/thresholds2.csv')
threshs = []
for ind, x_ in x.iterrows():
	# if x_['thresh']==0:
	# 	thresh = 100
	# else:
	freqkey = str(x_['freq'])
	try:
		ix = (d[freqkey][1, :]==x_['thresh']).nonzero()[0][0]
		thresh = float32(d[freqkey][0, ix])
	except:
		thresh = 75
	threshs.append(thresh)		
		

xls = pd.ExcelFile('/Users/robert/Dropbox/FMRP/JN_3/thresholds.xls')
x = xls.parse('thresholds')

animalnum = x.animal.astype(np.str)
x['animalID'] = x.gen+x.exp+animalnum
		
gengp = x.groupby(('gen', 'exp', 'freq'))
genmeans = gengp.agg(dict(threshdbspl = np.mean))
generrs = gengp.agg(dict(threshdbspl = misc.pd_sem))

means = genmeans.unstack().values
stds = generrs.unstack().values
labels = genmeans.unstack().index
fig = plt.figure()
ax = fig.add_subplot(111)
freqs = ['4k', '8k', '16k', '32k', '45k']
for i in xrange(4):
	ax.errorbar(range(len(freqs)), means[i, :], yerr = stds[i, :], label = labels[i], mew = 0, marker = 'o', ms = 5, capsize = 0)

ax.set_xticks(range(len(freqs)))
ax.set_xticklabels(freqs)
ax.legend()
plt.show()
fig.savefig('tmp.eps')

animalpaths = glob.glob('*')
for animalpath in animalpaths:
	os.chdir(animalpath)
	fnames = glob.glob('*')
	for fname in fnames:
		d = pd.HDFStore(fname)
		df = d['df']
		df = Gap.add_gapratio(df)
		d['df'] = df
		d.close()

	os.chdir('..')


dir1 = '/Volumes/BOB_SAGET/TNFalpha/tinnitus/behavior/incoming/xxx'
dir2 = '/Volumes/BOB_SAGET/TNFalpha/tinnitus/behavior/thalidomide/Gap/data'
fnames1 = glob.glob(os.path.join(dir1, '*.txt'))
fnames2 = glob.glob(os.path.join(dir2, '*.txt'))
set1 = []
set2 = []
for fname in fnames1:
	absol, relat = os.path.split(fname)
	animalID, _, _, mo, da, yr, _ = relat.split('_')
	set1.append('_'.join((mo, da, yr, animalID)))

for fname in fnames2:
	absol, relat = os.path.split(fname)
	animalID, _, _, mo, da, yr, _ = relat.split('_')
	set2.append('_'.join((mo, da, yr, animalID)))

set1 = set(set1); set2 = set(set2);
diffset = list(set1-set2)
diffset.sort()



'''combine two parts of the same gap detect session'''
# dir1 = '/Volumes/BOB_SAGET/TNFalpha/tinnitus/behavior/tmp'
dir1 = '/Users/robert/Desktop/tmp'
os.chdir(dir1)
fnames1 = glob('*3_0.txt')
fnames2 = glob('*3_1.txt')
fnames1.sort()
fnames2.sort()
for fname1, fname2 in zip(fnames1, fnames2):
	fname, ext = os.path.splitext(fname1)
	outfname = fname[:-4]+ext
	header = np.loadtxt(fname1, 'S')[0]
	x1 = np.loadtxt(fname1, 'S', skiprows = 1)
	x2 = np.loadtxt(fname2, 'S', skiprows = 1)
	x = np.vstack((x1, x2))
	np.savetxt(outfname, x.astype(str), '%s', header = '\t'.join(header), delimiter = '\t', comments = '')



fpaths = glob.glob('*')
fpaths.sort()
fig = plt.figure();
for i, fpath in enumerate(fpaths):
	df = Gap.txt2pd(fpath)
	ax = fig.add_subplot(3, 3, i+1)
	Gap.plot_startleampl(df, ax)
	ax.set_ylim([0, 8])
	ax.set_title(os.path.split(fpath)[1])
	ax.set_xticklabels('')


plt.close('all')
ylim = (37.65, 37.85)
xlim = (-122.6, -122.3)
map_ = basemap.Basemap(resolution = 'h', projection='merc', \
	llcrnrlat = ylim[0], urcrnrlat = ylim[1], \
	llcrnrlon = xlim[0], urcrnrlon = xlim[1])
map_.drawcoastlines()
x, y = map_.makegrid(10, 10)
map_.drawmeridians(x, labels = [0, 0, 0, 1])
map_.drawparallels(y, labels = [1, 0, 0, 0])
map_.fillcontinents(color = 'coral')
# map_.drawcountries()
plt.show()

#Let's create a basemap of Europe
	x1 = -20.
	x2 = 40.
	y1 = 32.
	y2 = 64.
	m = basemap.Basemap(resolution='h',projection='merc', llcrnrlat=y1,urcrnrlat=y2,llcrnrlon=x1,urcrnrlon=x2,lat_ts=(x1+x2)/2)
	m.drawcountries(linewidth=0.5)
	m.drawcoastlines(linewidth=0.5)
m.drawparallels(np.arange(y1,y2,2.),labels=[1,0,0,0],color='black',dashes=[1,0],labelstyle='+/-',linewidth=0.2) # draw parallels
m.drawmeridians(np.arange(x1,x2,2.),labels=[0,0,0,1],color='black',dashes=[1,0],labelstyle='+/-',linewidth=0.2) # draw meridians

class Coordinate(object):
	def __init__(self, x, y):
		self.x = x
		self.y = y
	def __repr__(self):
		return 'Coord: ' +str(self.__dict__)


ix = np.vstack((df.condition=='postnihl', df.condition=='prenihl',
	df.condition=='thalid1', df.condition=='thalidwashout1',
	df.condition=='vehicle1', df.condition=='vehiclewashout1')).any(0)
df = df[ix]


im_pre = imread('pre-lesion.png')
im_post = imread('post-lesion.png')

fig, ax = plt.subplots()
ax.imshow(im_post, extent=(0, 1, 0.4, 1.1))
plt.show()
xy = plt.ginput(7)

df = pd.DataFrame(columns=('group', 'condition', 'freq', 'gapratio', 'yerr'))

freq = [5, 7, 10, 14, 20, 28, 40]
condition = 'post'
group = 'lesion'
x, gapratio = zip(*post_lesion_y)
x, y = zip(*post_lesion_yerr)
yerr = np.abs(np.array(gapratio)-np.array(y))
d = dict(freq=freq, gapratio=gapratio, yerr=yerr)
df_ = pd.DataFrame(d)
df_['condition'] = condition
df_['group'] = group
df = df.append(df_, ignore_index=True)

conditions = ['pre', 'post']
groups = ['control', 'lesion']

for condition in conditions:
	df1 = df[df.condition==condition]
	fig, ax = plt.subplots()
	ax.set_title(condition)
	for group in groups:
		df2 = df1[df1.group==group]
		ax.errorbar(range(7), df2.gapratio, yerr=df2.yerr, label=group)
	ax.set_ylim([0, 1.5])
	ax.legend()

ix = ['2' not in i for i in df.condition]
df = df.ix[ix]


df = df[df.condition]

ucond = np.unique(df.condition).values
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

# group changes in gap detection by condition and frequency
condgp = condstack.groupby(('animalcondition', 'freq'))
condmeans = condgp.manipdiff.apply(np.mean).unstack('animalcondition')
condsems = condgp.manipdiff.apply(st.sem).unstack('animalcondition')

# line plot (change in gap detection VS frequency)
ufreqs = np.unique(condstack.freq).values
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


		
# calculate difference caused by a manipulation
# input is a df with gap ratio at each frequency for each session
#     freq  animalID  gapratio condition gen  age        sess  postdatenihl  \
# 1   7071  Black112  0.765494   prenihl  WT   76  2013-07-01            -2   
# 2  10000  Black112  0.646148   prenihl  WT   76  2013-07-01            -2   
# 3  14142  Black112  0.783128   prenihl  WT   76  2013-07-01            -2   
# 4  20000  Black112  0.617598   prenihl  WT   76  2013-07-01            -2   
# 5  28284  Black112  0.974314   prenihl  WT   76  2013-07-01            -2   

#           DOB  datenihl  
# 1  2013-04-16  7/3/2013  
# 2  2013-04-16  7/3/2013  
# 3  2013-04-16  7/3/2013  
# 4  2013-04-16  7/3/2013  
# 5  2013-04-16  7/3/2013  
# with one "control" condition and one manipulation condition
# should include multiple animals
# output will be the difference in gap detection at each freq for individual animals






		






