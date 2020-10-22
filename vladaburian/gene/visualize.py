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

from sklearn.decomposition import PCA

from gene import speedup


SKIP_LIST = [
    {13, 16, 25, 29, 30, 36, 39, 43, 44, 49, 53, 56, 65, 68, 69, 71, 81, 83, 88, 91, 93, 97, 98, 99, 101, 104, 109, 110, 123, 124, 126, 127, 128, 130, 131, 132, 133, 134, 136, 137, 138, 140, 142, 147, 153, 160, 163, 164, 165, 166, 172, 175, 176, 178, 186, 187, 191, 193, 200, 201, 202, 203, 205, 208, 210, 212, 215, 216, 221, 222, 224, 229, 232, 235, 237, 239, 241, 245, 247, 248, 250, 251, 255, 256, 260, 262, 263, 266, 275, 277, 281, 282, 284, 286, 287, 293, 294, 297, 301, 303, 304, 305, 310, 311, 319, 320, 322, 324, 327, 328, 329, 331, 332, 334, 336, 340, 342, 344, 349, 350, 359, 360, 361, 362, 363, 364, 365, 368, 369, 371, 373, 374, 375, 377, 378, 383, 385, 389, 392, 402, 405, 407, 411, 412, 418, 424, 426, 431, 432, 435, 437, 440, 443, 444, 447, 450, 453, 454, 458, 459, 461, 463, 464, 465, 466, 470, 471, 477, 478, 479, 482, 486, 487, 490, 491, 495},
    {13, 16, 23, 25, 29, 30, 31, 33, 36, 37, 43, 44, 45, 47, 48, 49, 51, 52, 53, 56, 60, 65, 68, 69, 71, 72, 74, 77, 81, 83, 84, 85, 86, 87, 88, 91, 93, 95, 97, 99, 101, 104, 108, 109, 110, 113, 117, 118, 123, 124, 126, 127, 128, 129, 130, 131, 132, 133, 134, 136, 137, 138, 140, 141, 142, 145, 147, 149, 153, 154, 156, 160, 162, 163, 164, 165, 166, 170, 172, 173, 175, 176, 178, 179, 186, 187, 190, 191, 193, 198, 200, 201, 202, 203, 205, 206, 210, 211, 212, 215, 216, 217, 219, 221, 222, 224, 225, 229, 232, 235, 237, 239, 240, 241, 245, 247, 248, 250, 251, 255, 256, 259, 260, 262, 263, 266, 268, 269, 272, 273, 275, 277, 278, 280, 281, 282, 284, 286, 287, 290, 293, 294, 297, 299, 301, 303, 304, 305, 307, 310, 311, 312, 313, 315, 316, 317, 319, 320, 322, 324, 328, 329, 331, 332, 333, 334, 335, 336, 340, 342, 344, 346, 349, 350, 352, 355, 358, 359, 360, 361, 362, 363, 364, 365, 367, 368, 369, 371, 373, 374, 375, 377, 378, 379, 383, 385, 386, 388, 389, 390, 392, 393, 395, 399, 402, 403, 404, 405, 406, 407, 411, 412, 413, 418, 422, 424, 426, 427, 428, 431, 432, 435, 436, 437, 443, 444, 447, 450, 453, 455, 458, 459, 461, 463, 464, 465, 466, 470, 471, 472, 473, 474, 477, 478, 479, 480, 482, 485, 486, 487, 490, 491, 492, 493, 494, 495, 497, 500},
]

INCLUDE_LIST = [
    {13, 16, 25, 29, 30, 36, 37, 39, 43, 44, 45, 49, 52, 53, 56, 65, 68, 69, 71, 83},
    {}
]



def plt_input(keys='an'):
    key = [None]

    def press(evt, key=key):
        key[0] = evt.key

    fig = plt.gcf()
    cid = fig.canvas.mpl_connect('key_press_event', press)

    while True:
        if not plt.waitforbuttonpress():
            continue

        if key[0] and key[0] in keys:
            break

    fig.canvas.mpl_disconnect(cid)

    return key[0]

def plt_gmm(k, w, mu, sigma):
    from scipy.stats import norm

    minX = min(mu - 3 * sigma)
    maxX = max(mu + 3 * sigma)

    x = np.linspace(minX, maxX, 1000)
    ysum = np.zeros((len(x),))

    for j in range(k):
        y = w[j] * norm.pdf(x, mu[j], sigma[j])
        ysum += y
        plt.plot(x, y, '--')

    plt.plot(x, ysum, '-')


def bdr(w, mu, sigma):
    if mu[0] > mu[1]:
        w = [w[1], w[0]]
        mu = [mu[1], mu[0]]
        sigma = [sigma[1], sigma[0]]

    a0 = (sigma[0]**2) * (2*log(w[0]) - log(2*pi * (sigma[0]**2)))
    a1 = (sigma[1]**2) * (2*log(w[1]) - log(2*pi * (sigma[1]**2)))

    r = (mu[1] + mu[0]) / 2 + (a1 - a0) / (2 * (mu[1] - mu[0]))

    err0 = w[0] * (1 - norm.cdf(r, mu[0], sigma[0]))
    err1 = w[1] * norm.cdf(r, mu[1], sigma[1])

    return r, err0 + err1



def bdr2(mu, sigma):
    if len(mu) == 1:
        return None, None

    if mu[0] > mu[1]:
        mu = [mu[1], mu[0]]
        sigma = [sigma[1], sigma[0]]


    x = np.linspace(mu[0], mu[1], 10000)
    y = speedup.pdf(x, 2, mu, sigma)
    y = (y[:,0] > y[:,1])

    for i, yi in enumerate(y):
        if not yi:
            r = x[i]
            break

    err0 = (1 - norm.cdf(r, mu[0], sigma[0]))
    err1 = norm.cdf(r, mu[1], sigma[1])

    return r, err0 + err1





def ex(items, key):
    return [item[key] for item in items]


def tform(X):
    return np.sqrt(X)



def main(experiment, *, num=None):
    #conf.DEBUG = True

    SKIP = SKIP_LIST[num]

    ds = gio.load_prepared_dataset(conf.RES, experiment)
    r, codes, wells, code2gene, gt, cmap = ds.input, ds.codes, ds.wells, ds.code2gene, ds.gt, ds.ref

    result = {}

    for code, (geneA, geneB) in conf.CODE2GENE.items():
        result.update({(geneA, well): ds.ref[(code, well)][0] for well in ds.wells})
        result.update({(geneB, well): ds.ref[(code, well)][1] for well in ds.wells})

    #codes = [16]
    #codes = codes[:20]

    inp = ''
    cors = []

    for code, (geneA, geneB) in conf.CODE2GENE.items():

        values_gt = [v for k,v in ds.gt.items() if k[0] == code]
        values_ref = [v for k,v in ds.ref.items() if k[0] == code]

        cntGeneALess = sum(v[0] < v[1] for v in values_gt)
        isGeneALess = cntGeneALess > (0.5 * len(values_gt))
        isGeneALess = conf.IS_A_LESS[code]

        print("GT:", isGeneALess, cntGeneALess, len(values_gt))

        values_gt = list(itertools.chain(values_gt))
        values_ref = list(itertools.chain(values_ref))

        cor, _ = spearmanr([x[0] for x in values_gt], [x[0] for x in values_ref])
        cors.append(cor)
        print("COR = ", cor)

        cor, _ = spearmanr([x[1] for x in values_gt], [x[1] for x in values_ref])
        cors.append(cor)
        print("COR = ", cor)

        if (code in SKIP) and 0:
            #pass
            continue

        if (code not in INCLUDE_LIST[num]) and 1:
            pass
            #continue

        if code > 83:
            pass
            #continue

        if 0 and code not in (41,):
            continue


        errors = []
        errorsRef = []

        print("code", code)

        presult = []
        values = []

        for well in wells:
            values = np.hstack((values, r[(code, well)]))

        values = np.sort(values)
        # values = np.array(values)
        #values = np.log(values[values > 0])
        values = np.sqrt(values)
        #values = np.log(values + 100)

        allvalues = values

        _, _, gest = estimate3(values, noise=.05)

        if not (0 < sum(gest['mu']) < 1e5):
            assert(False)
            print("cant estimate peaks", gest)
            print("skipping code", code)
            continue

        if (abs(gest['mu'][0] - gest['mu'][1])) < 2*(sum(gest['sigma'])):
            #print("very close peaks", gest)
            #print("skipping code", code)
            #continue
            pass

        if conf.DEBUG:
            #_, _, _, _, gest = estimate3(values)

            print("EST:", gest)

            plt.clf()
            plt.hist(values, bins=200, density=True)
            plt.grid()
            plt.title("{} (press 'a' to inspect)".format(code))
            plt_gmm(2, gest['w'], gest['mu'], gest['sigma'])

            thold, err = bdr2(gest['mu'], gest['sigma'])
            print("BDR:", thold, err)

            plt.axvline(thold, color='red')

            inp = plt_input()

        for well in wells:

            if conf.DEBUG and 0:
                print('code', code, 'well', well)

            values = r[(code, well)]
            #values = np.log(values[values > 0])
            values = np.sqrt(values)
            #values = np.log(values + 100)

            value1A, value1B = gt[(code, well)]
            value2A, value2B = cmap[(code, well)]

            #estA, estB, est = estimate1(values)
            estA, estB, est = estimate3(values, .05)
            #estA, estB, _ = estimate4(values, (gest['mu'][0] < gest['mu'][1]))
            #estA, estB, _ = estimate4(values, isGeneALess)

            pA, pBx = norm.pdf([estA, estB], gest['mu'][0], gest['sigma'][0])
            pAx, pB = norm.pdf([estA, estB], gest['mu'][1], gest['sigma'][1])

            #print("G PROB:", pA, pB, pAx, pBx, "mul:", pA*pB, pAx*pBx)


            misclassified = False

            if (value1A < value1B) != (estA < estB):
                misclassified = True
                #print("MISCLASSIFIED !")

                #estA, estB = estB, estA

            exchange = False
            #exchange = pA*pB < pAx*pBx
            exchange = (pB < pBx) and (pA < pAx)
            #exchange = isGeneALess != (estA < estB)
            #exchange = (gest['mu'][0] < gest['mu'][1]) != (estA < estB)
            #exchange = misclassified

            #if misclassified != exchange:
            #    print('false positive')

            if exchange and 1:
                estA, estB = estB, estA

            if exchange != misclassified:
                if exchange:
                    print(well, "FALSE NEGATIVE")
                else:
                    print(well, "FALSE POSITIVE")

            else:
                if exchange:
                    print(well, "CORRECTED")

            if (pB > pBx) != (pA > pAx):
                print(well, "same category")


            if not (0 < (estA + estB) < 1e5):
                print("BAD ESTIMATION", estA, estB)
                estA = np.median(values)
                estB = estA


            if value2A == value2B:
                pass
                #xestA = np.median(values)
                #estB = estA
            #    estA, estB = value2A, value2B

            #estA, estB = value2A, value2B

            #if estA - value1A > 1500:
            #    estA += 1 * (estA - value1A)
            #    print('mod ', code, well)

            #estA = estA ** 2 - ((estA / 30.0) * 250.0)
            #estA = sqrt(estA)

            #presult.append(dict(A=estA**2, B=estB**2, Ap=est['w'][0], Bp=est['w'][1]))

            #result[(geneA, well)] = exp(estA) - 100
            #result[(geneB, well)] = exp(estB) - 100

            #result[(geneA, well)] = estA
            #result[(geneB, well)] = estB

            estA = estA**2
            estB = estB**2


            #estA = round(estA / acc) * acc
            #estB = round(estB / acc) * acc

            #estA = 16000 ** round(min(1.0, log(estA+1) / log(16000)), 2)
            #estB = 16000 ** round(min(1.0, log(estB+1) / log(16000)), 2)

            result[(geneA, well)] = estA
            result[(geneB, well)] = estB

            if conf.DEBUG:
                errors.append((estA - value1A))
                errors.append((estB - value1B))

            #errorsRef.append((value2A - value1A) ** 2)
            #errorsRef.append((value2B - value1B) ** 2)

            #print(well, code)
            #continue
            if inp != 'a':
                continue

            plt.clf()
            plt.hist(values, bins=50, color='skyblue')
            
            plt.axvline(value1A, color='red')
            plt.axvline(value1B, color='red', linestyle='--')
            plt.axvline(value2A, color='green')
            plt.axvline(value2B, color='green', linestyle='--')
            
            plt.axvline(estA, color='magenta')
            plt.axvline(estB, color='magenta', linestyle='--')

            plt.title("{} {} (genes {}, {}) {}".format(code, well, geneA, geneB, ('MISCLASSIFIED' if misclassified else '')))

            plt.grid()
            inp = plt_input()

            pass

        #cor1, _ = spearmanr(
        #    [v[0] for k, v in ds.gt.items() if k[0] == code],
        #    [v[1] for k, v in ds.gt.items() if k[0] == code])
        #
        #cor2, _ = spearmanr(ex(presult, 'A'), ex(presult, 'B'))
        #
        #print("COR BETWEEN GENES:", cor1, cor2)

        #prA = np.array(ex(presult, 'A'))
        #prB = np.array(ex(presult, 'B'))
        #prAB = (prA + prB) / 2
        #prA /= prAB
        #prB /= prAB

        #prAerr = np.array([v[0] for k,v in ds.gt.items() if k[0] == code])
        #prAerr -= np.power(prA, 2)

        #prA = np.sort(prA)
        #prB = np.sort(prB)

        #isAless = sum(prA < 0.999) > sum(prA > 1.001)

        for i, well in enumerate(wells):
            if 0 and isAless != (prA[i] < prB[i]):
                pass
                #result[(geneA, well)] = prB[i] ** 2
                #result[(geneB, well)] = prA[i] ** 2

                #result[(geneA, well)] = np.median(prB) ** 2
                #result[(geneB, well)] = np.median(prA) ** 2

        if conf.DEBUG and 0:
            plt.cla()
            corA, _ = spearmanr([x[0] for x in values_gt], prA)
            corB, _ = spearmanr([x[1] for x in values_gt], prB)

            print("MY COR:", corA, corB)

            #plt.scatter(prA, np.linspace(0, 1, len(prA)), c='red', alpha=.5)
            #plt.scatter(prB, np.linspace(0, 1, len(prB)), c='green', alpha=.5)

            plt.scatter(prA, prAerr, c='red', alpha=.5)

            #plt.scatter(prA, prA + prB, c='red', alpha=.5)
            #plt.scatter(prB, prA + prB, c='green', alpha=.5)
            #plt.scatter(pca_input[0,:], pca_input[1,:])
            plt.grid()
            plt.title("{} (genes {}, {})".format(code, geneA, geneB))
            plt_input()


        if conf.DEBUG:
            """
            plt.cla()
            plt.hist(allvalues, 100)
            plt.hist(ex(presult, 'A'), 100)
            plt.hist(ex(presult, 'B'), 100)
            plt.grid()
            plt.title("{} (genes {}, {})".format(code, geneA, geneB))
            plt_input()
            """

            plt.clf()
            ax1 = plt.gca()
            #fig, ax1 = plt.subplots()
            ax1.hist(allvalues * allvalues, 100)

            ax2 = ax1.twinx()
            #ax2.hist([ex(presult, 'A'), ex(presult, 'B')], 20, color=['red', 'green'], stacked=True, alpha=.8)
            ax2.hist(ex(presult, 'A'), 40, color='red', alpha=.7)
            ax2.hist(ex(presult, 'B'), 40, color='green', alpha=.7)

            #fig.tight_layout()
            #plt.show()
            plt_input()

        if conf.DEBUG and 0:
            plt.clf()
            plt.plot(sorted(errors))
            plt.grid()
            plt_input()


        #errorMed = np.median(errors)
        #errorRefMed = np.median(errorsRef)

        #if errorMed > errorRefMed:
        #    print("BAD RESULT", errorMed, errorRefMed)

    #errors = sorted(errors)
    #print("ERROR MEDIAN:", np.median(errors))
    #print("ERROR MEAN:", np.mean(errors))

    print("median cor:", np.median(cors))

    gio.save_gct(os.path.join(conf.RES, 'output2', experiment + '.gct'), result)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        conf.DEBUG = True
        no = 1
        main(conf.EXPERIMENTS[no], num=no)
    else:
        conf.DEBUG = False
        no = int(sys.argv[1])
        main(conf.EXPERIMENTS[no], num=no)
