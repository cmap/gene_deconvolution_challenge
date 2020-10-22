
import logging
import numpy as np
import os
import re
import subprocess
import sys

from copy import deepcopy

from gene import conf, cmd_main


logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', stream=sys.stdout, level=logging.INFO)


def run_test(params):
    logging.info("Processing with params: %s", params)

    for plate in conf.EXPERIMENTS:
        cmd_main.process(os.path.join(conf.RES, 'input', plate), os.path.join(conf.RES, 'output-tune'), plate, params)

    logging.info("Evaluating results")

    r = subprocess.run(
        [
            'docker', 'run', '--rm',
            '-v', os.path.join(conf.RES, 'output-tune') + ':/workdir',
            '-v', os.path.join(conf.RES, 'ground-truth') + ':/ground-truth',
            'cmap/scorer', conf.EXPERIMENTS[0], conf.EXPERIMENTS[1], '383', '344'
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )

    m = re.search(b'OVERALL SCORE = (\d+)', r.stderr)
    score = int(m[1])

    logging.info("Score is %s", score)

    return score


def main():
    params = {
        'noise': 0.02,
        'single_thold': 1.6,
    }

    X = [.0, .01, .02, .03, .04, .05, .06, .07, .08, .09, .1, .12, .15, .2, .25]
    X = [.012, .014, .016, .018, .022, .024, .026, .028]
    r = []

    for x in X:
        p = deepcopy(params)
        p['noise'] = x

        r.append(run_test(p))

    logging.info("Results: %s", r)


if __name__ == '__main__':
    sys.exit(main())
