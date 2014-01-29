import os, glob, re, shutil
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import datetime
import misc; reload(misc);

freqs = np.array([5, 7.1, 10, 14.1, 20, 28.3, 40])

def fileconvert_all(studydir):
    '''
    Converts all of the text file outputs in */Gap/data and converts them to pandas formatted hdf5 files.
    Stores the results in */Gap/fileconversion
    '''

    # dobs = pd.read_csv(os.path.join(studydir, 'dobs.csv'))
    pairpulsedir = os.path.join(studydir, 'pairpulse')
    if not os.path.exists(pairpulsedir):
        os.mkdir(pairpulsedir)

    animalpaths = glob.glob(os.path.join(studydir, 'data', 'PPI', '[0-9]*'))
    for animalpath in animalpaths:
        fpaths = glob.glob(os.path.join(animalpath, '[A-Za-z]*.txt'))
        for fpath in fpaths:

            absol, relat = os.path.split(fpath)
            animalID, gen, condition, mo, da, yr, _ = relat.split('_')

            animalinfo = get_animalinfo(animalID, studydir)
            dob_str = animalinfo.DOB.values[0]

            dob = misc.str2date(dob_str, delimiter = '/', format = 'MMDDYYYY')

            if hasattr(animalinfo, 'date1'):
                date1_str = animalinfo.date1.values[0]
                date1 = misc.str2date(date1_str, delimiter = '/', format = 'MMDDYYYY')

            mo = '%2.2u' % int(mo)
            da = '%2.2u' % int(da)
            yr = '%4.4u' % int(yr)

            sess_str = '_'.join((yr, mo, da))
            sess_date = misc.str2date(sess_str, delimiter = '_', format = 'YYYYMMDD')
            age = (sess_date - dob).days
            if hasattr(animalinfo, 'date1'):
                postdate1 = (sess_date - date1).days
            outpath = os.path.join(pairpulsedir, '%s.csv' % '_'.join((animalID, gen, condition, yr, mo, da)))

            if not os.path.exists(outpath):
                gapratio = fileconvert(fpath)
                if hasattr(animalinfo, 'date1'):
                    df = pd.DataFrame(dict(gapratio = gapratio, animalID = animalID, gen = gen, condition = condition, sess = sess_str, age = age, postdate1 = postdate1))
                else:
                    df = pd.DataFrame(dict(gapratio = gapratio, animalID = animalID, gen = gen, condition = condition, sess = sess_str, age = age))
                df.to_csv(outpath)

def fileconvert(fpath):
    '''
    Takes the text file output for one animal and converts it to a pandas formatted hdf5 file
    '''

    df = txt2pd(fpath)

    df = add_startleampl(df)

    uncued_ampl = df[df.cued=='uncued'].startleampl.mean()
    
    freqgp = df.groupby(('freq', 'cued'))
    startlemean = freqgp.agg(dict(startleampl = np.mean))
    startlemean = startlemean.unstack('cued')
    gapratio = startlemean['startleampl']['cued'] / uncued_ampl

    return gapratio[1:]

def txt2pd(fpath):

    print fpath

    cueddict = {0:'uncued', 1:'cued'}

    header = np.loadtxt(fpath, 'S', delimiter = '\t')[0]
    startletrace_ix = len(header)
    x = np.loadtxt(fpath, 'f', skiprows = 1)
    freq = x[:, 0]
    cued = x[:, 1]
    cuedstr = map(lambda x: cueddict[x], cued)
    # ampl = x[:, 2]
    # holdtime = x[:, 4]
    resp = [i for i in x[:, 5:]]
    df = pd.DataFrame(dict(freq = freq, cued = cuedstr, resp = resp))

    return df

def calc_gap_detection(df):
    df = df[df.cued=='cued']
    gp = df.groupby('freq')
    gapdetection = gp.agg(dict(gapratio = np.mean))

    return gapdetection

def plot_startleampl(df, ax = None):
    '''
    plots the startle amplitude as a function of trial number
    uncued: red
    cued: blue
    '''
    if not hasattr(df, 'startleampl'):
        df = add_startleampl(df)

    ax, fig = misc.axis_check(ax)

    cuegp = df.groupby('cued')
    colors = dict(cued ='b', uncued='r')
    for i, j in cuegp:
        ax.plot(j.index, j.startleampl, '.', color = colors[i])

def plot_startleamp_per_freq(df):

    gp = df.groupby(('freq', 'cued'))
    for i, j in cuegp:
        freqgp.groupby()
        print i

def plot_cued_vs_uncued(df):
    cuegp = df.groupby('cued')
    respmean = cuegp.resp.apply(np.mean)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i, j in respmean.iteritems():
        ax.plot(j, label = i)

    ax.legend()

def plot_group_results(df):

    df = df[df.cued==1]

    # group by session to get session mean gap detection
    gp_sess = df.groupby(['sess', 'freq'])
    gapratio_sess = gp_sess.agg(dict(gapratio = np.mean, freq = get_first, gen = get_first, exp = get_first, animalID = get_first))

    # combine sessions for each animal to give animal performance by condition
    gp_animal = gapratio_sess.groupby(['animalID', 'exp', 'freq'])
    gapratio_animal = gp_animal.agg(dict(gapratio = np.mean, freq = get_first, gen = get_first, exp = get_first))
    gp_exp = gapratio_animal.groupby(['freq', 'exp'])
    gapratio_mean = gp_exp.gapratio.apply(np.mean).unstack()
    gapratio_sem = gp_exp.gapratio.apply(apply_sem).unstack()
    fig = plt.figure();
    ax = fig.add_subplot(111);
    for (k, v), (kerr, verr) in zip(gapratio_mean.iteritems(), gapratio_sem.iteritems()):
        misc.errorfill(np.arange(freqs.size), v, verr, ax = ax, color = colors[k], marker = '.', ls = '-', label = k)

    ax.legend()

    maxy = ax.get_ylim()[1]
    ax.set_ylim([0, maxy])
    ax.set_ylabel('Gap ratio')
    ax.set_xlabel('Frequency (kHz)')
    ax.set_xticks(range(freqs.size))
    ax.set_xticklabels(freqs)

    return ax

def plot_gap_ratio():
    '''
    Somehow plots the gap ratio for all animals averaged over multiple sessions
    '''
    df = df[df.cued==1]
    
    sess_group = df.groupby(['sess', 'freq'])
    sess_gapratio = sess_group.agg(dict(gapratio = np.mean, exp = get_first))
    exp_group = sess_gapratio.groupby(['exp', 'freq'])
    gap_mean = exp_group.gapratio.apply(np.mean)
    plt.close('all')
    gap_mean.plot()

    gp = df.groupby(['animalID', 'freq', 'exp'])
    
    results_mean = gp.gapratio.apply(np.mean)
    results_err = gp.gapratio.apply(np.std)
    
    misc.pd_errorbar(results_mean, results_err)
    fig = plt.figure();
    ax = fig.add_axes([0.1, 0.3, 0.8, 0.6])
    
    results.plot(kind = 'bar', ax = ax)
    return
    
def add_startleampl(df):
    '''
    Finds the startle amplitude for each trial in a dataframe
    '''
    startleampl = []
    for i, d in df.iterrows():
        if d['cued']=='uncued':
            startleampl_ = d['resp'][100:200]
        elif d['cued']=='cued':
            startleampl_ = d['resp'][200:300]
        startleampl.append(startleampl_.max()-startleampl_.min())

    df['startleampl'] = startleampl
    return df


def get_animalinfo(animalID, studydir):
    fname = os.path.join(studydir, 'dobs.csv')
    df = pd.read_csv(fname)
    ncols = len(df.columns)
    ix = df.animalID==animalID    
    info = df[ix]
    # try:
    #     dob = misc.str2date(dob.values[0], delimiter = '/')
    # except IndexError:
    #     print '%s not found!' % animalID
    #     return
    return info

def get_age(animalID, date = 'today'):
    if date=='today':
        date = datetime.datetime.now().date()
    dob = get_animalinfo(animalID)['DOB']
    return (date-dob).days

def reverse_animalID(animalID):
    p = re.compile('[0-9]')
    ix = p.search(animalID).start()
    cagenum = animalID[ix:]
    animal = animalID[:ix]
    return cagenum+animal

def move_to_animal_folder():
    '''Move individual files to folders named for the animal'''
    fnames = glob.glob('*.txt')
    for fname in fnames:
        animalID = fname.split('_')[0]
        outdir = reverse_animalID(animalID)
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        shutil.move(fname, os.path.join(outdir, fname))
        