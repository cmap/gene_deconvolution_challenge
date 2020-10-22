import numpy as np
import os
import pickle
import re
import types
import multiprocessing

from collections import defaultdict

from gene import speedup


def unique(*args):
    return sorted(set(*args))


def save_gct(filename, dataset):
    genes = unique(k[0] for k in dataset.keys())
    wells = unique(k[1] for k in dataset.keys())

    with open(filename, 'w') as f:
        f.write('#1.3\n')
        f.write('{:d}\t{:d}\t0\t0\n'.format(len(genes), len(wells)))
        f.write('id\t' + '\t'.join(wells) + '\n')

        for gene in genes:
            line = str(gene) + '\t' + '\t'.join('{:.1f}'.format(dataset[(gene, well)]) for well in wells) + '\n'
            f.write(line)


def load_measurement_single(filename):
    with open(filename) as f:
        lines = iter(f)
        next(f)

        result = defaultdict(list)

        for code, fi in (line.split() for line in lines):
            result[int(code)].append(float(fi))

        result = {k:np.array(sorted(v)) for k,v in result.items()}

        return result


def load_measurements(root):
    plate = os.path.split(root)[-1]
    return speedup.load_lbx(root, plate);

    _,_,filenames = next(os.walk(root))

    result = {}
    result_async = []

    pool = multiprocessing.Pool()

    for filename in filenames:
        well = filename[-7:-4]
        assert(re.match('[A-Z][0-9][0-9]', well))

        result_async.append((well, pool.apply_async(load_measurement_single, (os.path.join(root, filename),))))

    for well, r in result_async:
        result.update(((code, well), values) for code, values in r.get().items())

    return result


def load_code2gene(filename):
    with open(filename) as f:
        lines = iter(f)
        next(lines)

        r = defaultdict(lambda: [None, None])
        
        for line in lines:
            code, gene, high = line.split()
            r[int(code)][0 if high == '0' else 1] = int(gene)

        return dict(r.items())


def load_gct(filename):
    with open(filename) as f:
        lines = f.readlines()

    result = {}

    wells = lines[2].split()[1:]

    for line in lines[3:]:
        cols = line.split()
        gene = int(cols[0])
        
        for i, well in enumerate(wells):
            result[(gene, well)] = float(cols[i + 1])

    return result


def uni2duo(code2gene, dataset):
    result = {}
    
    wells = unique(k[1] for k in dataset.keys())
    
    for code in code2gene.keys():
        for well in wells:
            geneA, geneB = code2gene[code]
            
            result[(code, well)] = (
                dataset.get((geneA, well), None),
                dataset.get((geneB, well), None)
            )

    return result


def load_dataset(root, plate, *, with_gt=False, with_ref=False):
    ds = types.SimpleNamespace()

    #ds.input = load_measurements(os.path.join(root, 'input', experiment))
    ds.input = load_measurements(os.path.join(root))
    ds.codes = unique(k[0] for k in ds.input.keys())
    ds.wells = unique(k[1] for k in ds.input.keys())

    #ds.code2gene = load_code2gene(os.path.join(root, 'input', 'barcode_to_gene_map.txt'))

    if with_gt:
        ds.gt = load_gct(os.path.join(root, 'ground-truth', plate + '_DECONV_UNI.gct'))
        ds.gt = uni2duo(ds.code2gene, ds.gt)

    if with_ref:
        ds.ref = load_gct(os.path.join(root, 'output', plate + '.gct'))
        ds.ref = uni2duo(ds.code2gene, ds.ref)

    return ds


def cmd_prepare_dataset(root, experiment):
    ds = load_dataset(root, experiment, with_gt=True, with_ref=True)

    with open(os.path.join(root, 'input', experiment + '.pickle'), 'wb') as f:
        pickle.dump(ds, f, pickle.HIGHEST_PROTOCOL)


def load_prepared_dataset(root, experiment):
    with open(os.path.join(root, 'input', experiment + '.pickle'), 'rb') as f:
        return pickle.load(f)
