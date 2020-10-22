
import argparse
import logging
import multiprocessing
import os
import sys
import numpy as np

from scipy.stats import norm

from gene import conf, speedup
from gene.gio import load_measurements
from gene.deconv import estimate3


logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', stream=sys.stdout, level=logging.INFO)


def process_barcode(measures, code, geneA, geneB, params):
    X = np.hstack(tuple(measures.values()))
    X = np.sort(X)

    X = np.sqrt(X)

    _, _, gest = estimate3(X, noise=.05)

    if not (0 < sum(gest['mu']) < 1e5):
        gest = None

    resultA = {}
    resultB = {}

    for well, values in measures.items():
        values = np.sqrt(values)

        estA, estB, _ = estimate3(values, noise=params['noise'], single_thold=params['single_thold'])

        if gest:
            pA, pBx = norm.pdf([estA, estB], gest['mu'][0], gest['sigma'][0])
            pAx, pB = norm.pdf([estA, estB], gest['mu'][1], gest['sigma'][1])

            exchange = (pB < pBx) and (pA < pAx)

            if exchange:
                estA, estB = estB, estA

        estA = estA**2
        estB = estB**2

        resultA[well] = estA
        resultB[well] = estB

    #wells = list(resultA.keys())
    #single = set(w for w,x in resultA.items() if x == resultB[w])

    #medianA = np.median(list(x for well,x in resultA.items() if not well in single))
    #medianB = np.median(list(x for well,x in resultB.items() if not well in single))

    #for well in single:
    #    resultA[well] = medianA
    #    resultB[well] = medianB

    return geneA, resultA, geneB, resultB


def process(dspath, out, plate, params):
    logging.info("Reading input data")
    ds = load_measurements(dspath)
    logging.info("Processing")

    pool = multiprocessing.Pool()
    results = []

    for code, (geneA, geneB) in sorted(conf.CODE2GENE.items()):
        measures = ds[code]
        results.append(pool.apply_async(process_barcode, (measures, code, geneA, geneB, params)))

    result = {}

    for part in results:
        geneA, resultA, geneB, resultB = part.get()
        result[geneA] = resultA
        result[geneB] = resultB

    pool.close()
    pool.join()

    logging.info("Writing output")
    speedup.save_gct(os.path.join(out, plate + '.gct'), result)

    logging.info("DONE")



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dspath', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--plate', required=True)

    args, _ = parser.parse_known_args()

    params = {
        'noise': 0.02,
        'single_thold': 1.6,
    }

    process(args.dspath, args.out, args.plate, params)


if __name__ == '__main__':
    sys.exit(main())
