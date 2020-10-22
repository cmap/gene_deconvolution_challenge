import pickle
import os
import sys
import numpy as np
import itertools

import matplotlib.pyplot as plt

from copy import deepcopy
from math import log, exp, pi, sqrt

from gene import em, gio, conf
from gene.deconv import *

from scipy.stats import norm
from scipy.stats.mstats import spearmanr



def main(experiment, *, num=None):
    ds = gio.load_prepared_dataset(conf.RES, experiment)

    with open(os.path.join(conf.PWD, 'model.bin'), 'rb') as f:
        bins, rf = pickle.load(f);

    keys = list(k for k in ds.input.keys() if k[0] in conf.CODE2GENE)

    featuresA = []
    featuresB = []

    for k in keys:
        f,_ = np.histogram(ds.input[k], bins=bins)

        featuresA.append(np.hstack(([1], f)))
        featuresB.append(np.hstack(([0], f)))

    featuresA = np.vstack(featuresA)
    featuresB = np.vstack(featuresB)

    estA = rf.predict(featuresA)
    estB = rf.predict(featuresB)


    result = {}

    for i, k in enumerate(keys):
        code, well = k
        geneA, geneB = conf.CODE2GENE[code]

        result[(geneA, well)] = estA[i]
        result[(geneB, well)] = estB[i]

    gio.save_gct(os.path.join(conf.RES, 'output2', experiment + '.gct'), result)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        conf.DEBUG = True
        no = 0
        main(conf.EXPERIMENTS[no], num=no)
    else:
        conf.DEBUG = False
        no = int(sys.argv[1])
        main(conf.EXPERIMENTS[no], num=no)
