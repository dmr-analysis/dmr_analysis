import scipy.stats as sc
import itertools as it
import numpy as np
import pandas as pd
import os
import tempfile
try:
    import matlab.engine
except ImportError:
    pass

def ranksum(x, y, technique=None, eng=None):
    if list(set(x)) == list(set(y)):
        return 1
    nx = len(x)
    ny = len(y)
    ns = min(nx, ny)

    #jbw
    if technique == 'Pranksum':
        if ns < 10:
            technique = 'full_enumeration'
        else:
            technique = 'mannwhitneyu'

    if technique == 'full_enumeration':
        if nx <= ny:
            smsample = list(x)
            lgsample = list(y)
            same_order = True
        else:
            smsample = list(y)
            lgsample = list(x)
            same_order = False
        # slow for big sample sizes
        ranks = sc.rankdata(lgsample + smsample)
        srank = ranks[:ns]
        w = np.sum(srank)
        df = pd.DataFrame(list(it.combinations(ranks, ns))).sum(axis=1)
        len_sr = float(len(df))

        plo = np.sum(df <= w)/len_sr
        phi = np.sum(df >= w)/len_sr
        p_tail = min(plo, phi)
        p = min(2*p_tail, 1)
    elif technique == 'mannwhitneyu':
        p = sc.mannwhitneyu(x, y, use_continuity=True, alternative='two-sided').pvalue
    elif technique == 'R_wilcox':
        # rpy2
        #    import rpy2.robjects as robjects
        #    wt_r = robjects.FloatVector(x)
        #    ko_r = robjects.FloatVector(y)
        #    wilcox_result_new = robjects.r['wilcox.test'](wt_r, ko_r)
        #    p = wilcox_result_new[2][0]
        p = exact_wilcox_test(x, y, 'two.sided')
    elif technique == 'M_ranksum':
        x = list(x)
        y = list(y)
        x = matlab.double([x])
        y = matlab.double([y])
        p = eng.ranksum(x,y)

    return p

def exact_wilcox_test(x, y, side):
    tf = tempfile.NamedTemporaryFile(delete=False)
    tf.close()

    with open(tf.name, 'w') as f:
        f.writelines("%s\n" % l for l in np.concatenate((x, y)))

    command = "Rscript --default-packages=stats,utils -e \"options(warn=-1)\" -e \"data=read.table(\\\"" + tf.name + \
        """\\\"); wilcox.test(data[1:""" + str(len(x)) + ",1], data[" + str(len(x) + 1) + ":" + str(len(x) + len(y)) + \
        """,1], correct = TRUE, alternative='""" + side + "')\$p.value \" "

    wilcox_p = os.popen(command).read()

    os.remove(tf.name)
    return float(wilcox_p.split(" ")[-1])
