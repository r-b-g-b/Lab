from sklearn.svm import SVC
from sklearn.decomposition import PCA
from voc import load_spec, plot_spec, lowpass_spec
from scipy.stats import zscore

hits_fname = '/Volumes/BOB_SAGET/Vocalization/analysis/gui_count_output/KO4_P09_20.txt'
spec_fname = '/Volumes/BOB_SAGET/Vocalization/Cages/KO4/P09/KO4_P09_20.h5'

P, F, T = load_spec(spec_fname)
P, F = lowpass_spec(P, F)

nfreqs, ntimes = P.shape

hits = np.loadtxt(hits_fname)
times = hits[:, 0]
ix = time2ix(times, T.max(), ntimes)

X = np.zeros((500, 42300))
Y = np.zeros(500)
for i in range(500):
	X[i, :] = P[:, i*150:i*150+300].ravel()
	Y[i] = np.vstack(((ix-(i*150))>0, (ix-(i*150))<150)).all(0).sum()>0
	
X_ = zscore(X, axis = 0)


plot_spec(P[:, ix[0]-100 : ix[0]+200])
svc = SVC()
svc.fit(X, Y)


svc.predict(X)



def time2ix(t, Tmax, Imax):
	return np.int32((t/Tmax)*Imax)