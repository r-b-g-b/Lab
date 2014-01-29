import numpy as np
from struct import pack, unpack

valid_ckID = ['RIFF', 'WAVE', 'fmt ', 'fact', 'data']

def read_header(f):
	
	'''WAVE File Format'''
	
	f.seek(0)
	ckID, cksize = read_ckinfo(f)
	if ckID != 'RIFF':
		print 'File format is not RIFF'
	
	ckID = unpack('<4s', f.read(4))[0]
	
	if ckID != 'WAVE':
		print 'File format is not WAVE'
	
	ckID, cksize = read_ckinfo(f)
	if ckID != 'fmt ':
		print 'File format not found'
		
	(wFormatTag, nChannels, nSamplesPerSec, nAvgBytesPerSec, nBlockAlign, wBitsPerSample) = \
		unpack('<hhiihh', f.read(16))
	
	if (cksize == 18) or (cksize == 40):
		cbSize = unpack('<h', f.read(2))
	if cksize == 40:
		(wValidBitsPerSample, dwChannelMask, SubFormat) = unpack('<hi16s', f.read(22))
	
	return wFormatTag, nChannels, nSamplesPerSec, nAvgBytesPerSec, nBlockAlign, wBitsPerSample

# def find_ck(f, target_ckID = 'data'):
# 
# 	f.seek(0)
# 	
# 	while ckID != target_ckID:
# 		ckID, cksize = read_ckinfo(f)

def find_data_chunk(f):

	data = 'xxxx'
	while data!='data':
		data = data[1:] + f.read(1)

	return f.tell()-4

def read_ckinfo(f):
	
	ckID, cksize = unpack('<4si', f.read(8))
	return ckID, cksize
	
def read(fname):
	'''
	
	Read in a wave file
	Input:
		fname : path of wave file
	Output:
		x : Python list of data points
		nSamplesPerSec : sampling frequency of data
	
	'''
	
	f = open(fname, 'rb')
	
	wFormatTag, nChannels, nSamplesPerSec, nAvgBytesPerSec, nBlockAlign, wBitsPerSample \
		= read_header(f)
	
	if wFormatTag == 1:
		sampwidth = 4
		fmt = 'i'
	elif wFormatTag == 3:
		sampwidth = 4
		fmt = 'f'
	else:
		print 'File format %s not supported' % wFormatTag
	
	i = find_data_chunk(f)
	f.seek(i)
 	ckID, cksize = read_ckinfo(f)
	if ckID != 'data':
		print 'Data chunk not found'
		wef
	else:
		nsamp = cksize / sampwidth
		xs = f.read(cksize)[:(nsamp*sampwidth)]
		x = np.array(unpack('<%i%s' % (nsamp, fmt), xs), dtype = 'float32')
		
	f.close()

	return x, nSamplesPerSec

def write(data, fs, filepath, dtype):

	if data.ndim > 1:
		nchan = data.shape[1]
	else:
		nchan = 1
	
	nsamp = data.shape[0]

	if dtype is float:
		sampwidth = 4 # for float
		wFormatTag = 3
	elif dtype is int:
		sampwidth = 4
		wFormatTag = 1
	
	bytespersec = sampwidth * fs
	bitspersample = sampwidth * 8
	nbytes = nsamp * sampwidth

	f = open(filepath, 'wb')

	'''WAVE File Format'''
	f.write('RIFF') #ckID
	f.write(pack('i', 4+nsamp)) #cksize
	f.write('WAVE') #WAVEID

	'''Format chunk'''
	f.write('fmt ') #ckID
	f.write(pack('ihhiihh', 16, wFormatTag, nchan, fs, bytespersec, sampwidth, bitspersample)) # cksize

	'''Data chunk'''
	f.write('data') #ckID
	f.write(pack('i', nbytes)) # cksize
	f.write(pack('%if' % nsamp, *data))
	# ix = np.hstack((np.arange(0, nsamp, 1E6), nsamp))

	# for i in range(len(ix)):
	# 	chunksamp = ix[i+1] - ix[i]
	# 	# f.write(pack('%if' % nbytes, data.astype(dtype).tostring()))
	# 	f.write(pack('%if' % chunksamp, *data[ix[i]:ix[i+1]]) #data.astype(dtype).tostring())

	f.close()
