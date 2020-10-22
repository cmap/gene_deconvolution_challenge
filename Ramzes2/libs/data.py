import numpy as np
import pandas as pd

from tqdm import tqdm
from subprocess import call

from . import params

import os

INPUT_DIR = 'data/input/'
INPUT_DIR_SAMPLE1 = INPUT_DIR + 'DPK.CP001_A549_24H_X1_B42/'
INPUT_DIR_SAMPLE2 = INPUT_DIR + 'LITMUS.KD017_A549_96H_X1_B42/'
INPUT_BARCODE_TO_GENE = INPUT_DIR + 'barcode_to_gene_map.txt'

GROUND_TRUTH_DIR = 'data/ground-truth/'
GROUND_TRUTH_SAMPLE1 = GROUND_TRUTH_DIR + 'DPK.CP001_A549_24H_X1_B42_DECONV_UNI.gct'
GROUND_TRUTH_SAMPLE2 = GROUND_TRUTH_DIR + 'LITMUS.KD017_A549_96H_X1_B42_DECONV_UNI.gct'

INPUT_SAMPLES = [INPUT_DIR_SAMPLE1, INPUT_DIR_SAMPLE2]
GROUND_TRUTH_SAMPLES = [GROUND_TRUTH_SAMPLE1, GROUND_TRUTH_SAMPLE2]


def getInputFileNames(inputPath):
    fileNames = []

    for fileName in sorted(os.listdir(inputPath)):
        if fileName.endswith('txt'):
            fileNames.append(fileName)
    return fileNames


def getColumnNames(fileNames):
    cols = []
    for fileName in fileNames:
        col = fileName[:-4].split('_')[-1]
        cols.append(col)
    return cols


def getInputDataFromPath(inputPath):
    fileNames = getInputFileNames(inputPath)
    cols = getColumnNames(fileNames)

    basePath = inputPath
    experiments = []

    for fileName in tqdm(fileNames):
        a = pd.read_csv(basePath + fileName, sep='\t', engine='c', na_filter=False, keep_default_na=False, low_memory=False,
                        dtype={'barcode_id': np.int32, 'FI': np.float32})
        experiments.append(a)

    return experiments, cols


def getInputData(sampleId=0):
    basePath = INPUT_SAMPLES[sampleId]
    return getInputDataFromPath(basePath)



def getBarcodeToGeneIds():
    barcodeToGene = pd.read_csv(INPUT_BARCODE_TO_GENE, sep='\t',
                                dtype={'barcode_id': np.int32, 'gene_id': np.int32, 'high_prop': np.int32})

    lowPart = barcodeToGene.loc[barcodeToGene['high_prop']==0][['barcode_id', 'gene_id']]
    lowPart.rename(index=str, columns={"gene_id": "gene_low"}, inplace=True)

    highPart = barcodeToGene.loc[barcodeToGene['high_prop'] == 1][['barcode_id', 'gene_id']]
    highPart.rename(index=str, columns={"gene_id": "gene_high"}, inplace=True)

    merged = pd.merge(lowPart, highPart, on='barcode_id')

    return merged


def getGroundTruthData(sampleId=0):
    basePath = GROUND_TRUTH_SAMPLES[sampleId]

    res = pd.read_csv(basePath, sep='\t', comment='#', skiprows=2)
    return res

def getInputDataCppParser(inputPath, testName):
    preprocessedFileName = '/tmp/'+testName+'.csv'

    # if inputPath[0] != '/':
    #     inputPath = '../../' + inputPath

    retCode = call(['cpp_parser/bin/cpp_parser', INPUT_BARCODE_TO_GENE, inputPath, preprocessedFileName, '32', '1', str(params.CUT_THRESHOLD)])
    print(retCode)

    # 'barcode_id', 'col', 'min_value', 'max_value', 'gene_low', 'gene_high',
    # 'hist_0', 'hist_1', 'hist_2', 'hist_3', 'hist_4', 'hist_5', 'hist_6',
    # 'hist_7', 'hist_8', 'hist_9', 'hist_10', 'hist_11', 'hist_12',
    # 'hist_13', 'hist_14', 'hist_15', 'hist_16', 'hist_17', 'hist_18',
    # 'hist_19', 'hist_20', 'hist_21', 'hist_22', 'hist_23', 'hist_24',
    # 'hist_25', 'hist_26', 'hist_27', 'hist_28', 'hist_29', 'hist_30',
    # 'hist_31'],
    # dtype = 'object'

    dtype = {'barcode_id': np.int32, 'col': str, 'min_value': np.float32, 'max_value': np.float32, 'gene_low': np.int32, 'gene_high':np.int32}
    for i in range(params.HIST_BINS):
        dtype['hist_' + str(i)] = np.float32

    res = pd.read_csv(preprocessedFileName, engine='c', sep=';', na_filter=False, keep_default_na=False, low_memory=False, dtype=dtype)
    print(res.columns)

    return res