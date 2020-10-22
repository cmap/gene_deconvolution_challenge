import numpy as np
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed
from . import params

histColumns = ['hist_' + str(i) for i in range(params.HIST_BINS)]

# def _normalizedHist(x):
#     h, _ = np.histogram(x, bins=params.HIST_BINS)
#     h = np.array(h) / np.max(h)
#     return pd.DataFrame(h.reshape((1, -1)), columns=histColumns)

def _generateHist(experiment, features, normalizeBySum=False):
    bins = params.HIST_BINS

    experiment = pd.merge(experiment, features, on='barcode_id')
    experiment = experiment[(experiment['FI'] >= experiment['min_value']) & (experiment['FI'] <= experiment['max_value'])]
    experiment['FI_norm'] = (experiment['FI'] - experiment['min_value']) / (experiment['max_value'] - experiment['min_value'])
    experiment['hist_idx'] = np.minimum(np.floor(experiment['FI_norm'] * bins).astype(np.int32),
                                        bins - 1)

    h = experiment.groupby(['barcode_id', 'hist_idx'])['hist_idx'].agg(['count'])
    h = h.unstack(fill_value=0)
    if normalizeBySum:
        h = h.div(h.sum(axis=1), axis=0)
    else:
        h = h.div(h.max(axis=1), axis=0)
    h.columns = h.columns.droplevel()
    for i in range(params.HIST_BINS):
        if i not in h.columns:
            h.insert(i, i, 0.0)
    #h['barcode_id'] = h.index.values

    renameCols = {}
    for i in range(bins):
        renameCols[i] = 'hist_' + str(i)
    h.rename(columns=renameCols, inplace=True)

    return h

def _findGlobalMinMax(experiments):
    all_experiments = pd.concat(experiments)
    # df = df[df['FI']>0]
    df = all_experiments.groupby('barcode_id')
    df = df['FI'].agg(['min', 'max'])
    df.rename(columns = {'min': 'min_value', 'max': 'max_value'}, inplace=True)
    # df['barcode_id'] = df.index.values

    hist = _generateHist(all_experiments, df, normalizeBySum=True)
    hist = pd.merge(df, hist, on='barcode_id')


    hist[histColumns] = hist[histColumns] > params.CUT_THRESHOLD
    def apply_func(series):
        nz = np.nonzero(series[histColumns])[0]
        min_v = series['min_value']
        max_v = series['max_value']
        min_nz = nz[0]*(max_v-min_v)/params.HIST_BINS + min_v
        max_nz = (nz[-1]+1)*(max_v-min_v)/params.HIST_BINS + min_v
        return pd.Series([min_nz, max_nz], index=['min_value', 'max_value'])


    # def apply_func(series):
    #     v = series[histColumns].values
    #     left_s = np.argmax(np.cumsum(v)>params.CUT_THRESHOLD_SIDE)
    #     right_s = np.argmax(np.cumsum(v[::-1])>params.CUT_THRESHOLD_SIDE)
    #     min_v = series['min_value']
    #     max_v = series['max_value']
    #     min_nz = left_s*(max_v-min_v)/params.HIST_BINS + min_v
    #     max_nz = (len(v)-right_s)*(max_v-min_v)/params.HIST_BINS + min_v
    #     return pd.Series([min_nz, max_nz, series['barcode_id']], index=['min_value', 'max_value', 'barcode_id'])

    # hist['barcode_id'] = hist.index.values
    df = hist.apply( apply_func, axis=1)


    return df


def _preprocessOneExperiment(c, experiment, barcodeToGene, globalFeatures=None, gt = None):
    groups = experiment.groupby('barcode_id')
    # h = groups['FI'].apply(_normalizedHist)

    if globalFeatures is None:
        features = groups['FI'].agg({'min_value': np.min,
                                     'max_value': np.max})
    else:
        features = globalFeatures

    h = _generateHist(experiment, features)

    features = pd.merge(features, h, on='barcode_id')
    features = pd.merge(features, barcodeToGene, on='barcode_id')
    features['col'] = c

    if gt is not None:
        gt.rename(index=str, columns={c: "gene_pos"}, inplace=True)
        features = pd.merge(features, gt, left_on='gene_low', right_on='id')
        features.rename(index=str, columns={'gene_pos': 'gene_pos_low'}, inplace=True)
        features.drop(columns=['id'], inplace=True)

        features = pd.merge(features, gt, left_on='gene_high', right_on='id')
        features.rename(index=str, columns={'gene_pos': 'gene_pos_high'}, inplace=True)
        features.drop(columns=['id'], inplace=True)

        features['gene_pos_low_norm'] = (features['gene_pos_low'] - features['min_value']) / (
                features['max_value'] - features['min_value'])
        features['gene_pos_high_norm'] = (features['gene_pos_high'] - features['min_value']) / (
                features['max_value'] - features['min_value'])

    return features


def preprocessTrainData(experiments, columns, groundTruth, barcodeToGene, globalScaling):
    features = None
    if globalScaling:
        features = _findGlobalMinMax(experiments)
    allFeatures = Parallel(n_jobs=-1, verbose=1)(delayed(_preprocessOneExperiment)
                                                 (c, experiments[idx], barcodeToGene, features, groundTruth[['id', c]]) for idx, c in enumerate(columns))

    return pd.concat(allFeatures)


def preprocessTestData(experiments, columns, barcodeToGene, globalScaling):
    features = None
    if globalScaling:
        features = _findGlobalMinMax(experiments)
    allFeatures = Parallel(n_jobs=-1, verbose=1)(delayed(_preprocessOneExperiment)
                                                 (c, experiments[idx], barcodeToGene, features) for idx, c in enumerate(columns))

    return pd.concat(allFeatures)
