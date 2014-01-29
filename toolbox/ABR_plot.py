import os, glob, re
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import ABR
import misc


basedir = '/Users/robert/Documents/Work/Bao/Fmr1_Heesoo/ABR'
duration = 0.099424
nsamples = 244

'''
redo peak amplitudes using peak-trough height
'''
t = np.linspace(0, duration, nsamples)
x = np.load(os.path.join(basedir, 'peaks.npz'))['arr_0']
for i in range(x.size):
	x_ = x[i]
	gen = x_['gen']
	exp = x_['exp']
	animal = x_['animal']
	freq = x_['freq']
	atten = x_['atten']
	print x_['ampl']
	print x_['lat']
	abr = np.load(os.path.join(basedir, 'fileconversion', '%s_%s_%1.1i_ABR.npz' % (gen, exp, animal)))['arr_0']
	
	abr_ = abr[np.vstack((abr['gen']==gen, abr['exp']==exp, abr['animal']==animal, abr['freq']==freq, abr['atten']==atten)).all(0)]['data'][0, :]
	abr_filt = ABR.filter_abr(abr_)
	
	lats = np.hstack((x_['lat'], x_['lat'][-1]+0.0015))*10 # *10 because i got the times wrong
	for j in range(lats.size-1):
		[_, lat_ix_pe, _] = misc.closest(t, lats[j]) # peak ix
		[_, lat_ix_tr, _] = misc.closest(t, lats[j+1]) # trough ix
		ampl_peak = abr_filt[(lat_ix_pe-2):(lat_ix_pe+2)].max()
		ampl_trough = abr_filt[lat_ix_pe:lat_ix_tr].min()
		ampl_ = ampl_peak - ampl_trough
		
		x[i]['ampl'][j] = ampl_
		x[i]['lat'][j] = lats[j]


	print x[i]['ampl']
	print x[i]['lat']
	# fig = plt.figure()
	# ax = fig.add_subplot(111)
	# ax.plot(abr_filt)
	# ax.axvline(lat_ix_pe, color = 'r', ls = '--')
	# ax.axvline(lat_ix_tr, color = 'r', ls = '--')
	# ax.axhline(ampl_peak, color = 'r', ls = '--')
	# ax.axhline(ampl_trough, color = 'r', ls = '--')
	


'''redo peak amplitudes using new filtering...'''
t = np.linspace(0, duration, nsamples)
for i in range(x.size):
	print x[i]['ampl']
	x_ = x[i]
	gen = x_['gen']
	exp = x_['exp']
	animal = x_['animal']
	freq = x_['freq']
	atten = x_['atten']

	abr = np.load('%s_%s_%1.1i_ABR.npz' % (gen, exp, animal))['arr_0']
	
	abr_ = abr[np.vstack((abr['gen']==gen, abr['exp']==exp, abr['animal']==animal, abr['freq']==freq, abr['atten']==atten)).all(0)]['data'][0, :]
	
	ampl = np.empty(5)
	for j, lat in enumerate(x_['lat']):
		lat = lat*10
		[_, lat_ix, _] = misc.closest(t, lat)
		ampl_ = abr_[(lat_ix-2):(lat_ix+2)].max()
		x[i]['ampl'][j] = ampl_

	print x[i]['ampl']

# basedir = '/Volumes/BOB_SAGET/Fmr1_Heesoo/ABR'
basedir = '/Users/robert/Documents/Work/Bao/Fmr1_Heesoo/ABR'
dtype = np.dtype([('gen', 'S2'), ('exp', 'S3'), ('animal', 'i4'), ('freq', 'f4'), ('atten', 'f4'), ('lat', '5f4'), ('ampl', '5f4')])

plt_opts = {'wt' : {'color' : 'b', 'x_offset' : 0, 'y_offset' : 0.00001}, 'ko' : {'color' : 'r', 'x_offset' : 0.09, 'y_offset' : 0}}


'''

plot raw ABR traces in a nice array

'''
fig = plt.figure(figsize = [6., 8.7625])

gens = ['wt', 'ko']
exps = ['nai', 'exp']
duration = 0.00999424
nsamples = 244

t = np.linspace(0, duration, nsamples)
freq = 16000.
# t = t[20:-20]

for e, exp in enumerate(exps):
	ax = fig.add_subplot(1, 2, e+1)
	for gen in gens:
	for animal in range(1, 5):
		# animal = np.random.randint(1, 5)
		fig = plt.figure()
		ax = fig.add_subplot(111)
		abr = np.load('%s_%s_%1.1i_ABR.npz' % (gen, exp, animal))['arr_0']

		abr_ = abr[abr['freq']==freq]

		y = abr_['data'][:10:2, :].T
		y_filt = ABR.filter_abr(y)
		y_filt = y_filt - np.tile(y_filt.mean(0), (244, 1))

		ABR.plot_abr_attens(t, y_filt, abr_['atten'][:10:2], ax = ax, color = plt_opts[gen]['color'], y_offset_mag = 0.0005)
		
		ax.set_xlim([0.001, 0.007])

ax.plot([0.0015, 0.0015, 0.0025], [-0.000002, -0.000052, -0.000052], ls = '-', color = 'k')


plt.show()


'''

plot raw ABR traces and vlines at peaks

'''

fig = plt.figure()
ax = fig.add_subplot(111)	


gens = ['wt', 'ko']
duration = 0.099424
nsamples = 244
t = np.linspace(0, duration, nsamples)
# t = t[20:-20]

for gen in gens:
	exp = 'nai'
	animal = 1
	abr = np.load('%s_%s_%1.1i_ABR.npz' % (gen, exp, animal))['arr_0']

	freq = 16000.
	x_ = x[np.vstack((x['gen']==gen, x['exp']==exp, x['animal']==animal, x['freq']==freq)).all(0)]
	
	atten = x_['atten'][0]
	
	abr_ = abr[np.vstack((abr['freq']==freq, abr['atten']==atten)).all(0)]

	# y = abr_['data'].T
	y = ABR.filter_abr(abr_['data'].T)
	ax.plot(t, y, color = plt_opts[gen]['color'])
	for lat in x_['lat'][0]:
		ax.axvline(lat*10, color = 'r', ls = '--')

# ax.set_xlim([0.01, 0.06])
plt.show()




'''

make ABR table for paper

'''
measures = ['lat_del', 'ampl']
scales = [1000., 1000.]
freqs = [8000., 16000.]
gens = ['wt', 'ko']
exps = ['nai', 'exp']

calib = np.loadtxt('/Users/robert/Documents/Work/Bao/abr_calib.txt', 'S')
calib_header = calib[0, :]
calib = np.int32(calib[1:, :])

f = open('/Users/robert/Desktop/table2.txt', 'w')
for measure in measures:
	for freq in freqs:
		calib_f_ix = '%uk' % (freq/1000) == calib_header
		attens = [calib[0, calib_f_ix][0], calib[5, calib_f_ix][0]]
		for atten in attens:
			for gen in gens:
				for exp in exps:
					x_ = x[np.vstack((x['gen']==gen, x['exp']==exp, x['freq']==freq, x['atten']==atten)).all(0)]
					x_mean = (1000.*x_[measure]).mean(0)
					x_sem = (1000.*x_[measure]).std(0) / np.sqrt(x_.size)
					f.write( '%s\t%u\t%u\t%s\t%s\t%.3f\t(%.3f)\t%.3f\t(%.3f)\t%.3f\t(%.3f)\t%.3f\t(%.3f)\t%.3f\t(%.3f)\t\n' % (measure, freq, atten, gen, exp, x_mean[0], x_sem[0], x_mean[1], x_sem[1], x_mean[2], x_sem[2], x_mean[3], x_sem[3], x_mean[4], x_sem[4]))
					
f.close()



'''



'''


dtype = x1.dtype
x = np.empty(0, dtype = x1.dtype)
for x1_ in x1:
	ix = np.vstack((x2['gen']==x1_['gen'], x2['exp']==x1_['exp'], x2['animal']==x1_['animal'], x2['freq']==x1_['freq'])).all(0)
	x.resize(x.size+1)
	x[-1] = np.array((x1_['gen'], x1_['exp'], x1_['animal'], x1_['freq'], (x1_['thresh']+x2[ix]['thresh'])/2.), dtype = dtype)

dtype1 = x.dtype
dtype2 = np.dtype([('gen', 'S2'), ('exp', 'S3'), ('animal', 'i4'), ('freq', 'f4'), ('thresh', 'f4'), ('thresh_calib', 'f4')])

dtype_calib = np.dtype([('db', 'f4'), ('4k', 'f4'), ('8k', 'f4'), ('16k', 'f4'), ('32k', 'f4'), ('45k', 'f4')])
header = np.loadtxt('/Users/robert/Documents/Work/Bao/abr_calib.txt', delimiter = '\t', skiprows = 1, dtype = dtype_calib)

x2 = np.empty(0, dtype = dtype2)
for x_ in x:
	
	freq = x_['freq']
	freq_str = '%uk' % (freq/1000)
	thresh = x_['thresh']

	thresh_calib = header['db'][header[freq_str]==thresh][0]
	
	x2.resize(x2.size+1)
	x2[-1] = np.array((x_['gen'], x_['exp'], x_['animal'], x_['freq'], thresh, thresh_calib), dtype = dtype2)
	