from libs import data, preprocessing, generators, postprocessing, time_log, params
from libs import seg_models

import argparse

time_log.log('Start')

parser = argparse.ArgumentParser(description='predict cmap')
parser.add_argument('--dspath', dest='inputPath', default=data.INPUT_DIR_SAMPLE1)
parser.add_argument('--out', dest='outputPath', default="output")
parser.add_argument('--plate', dest='testName', default="test")
parser.add_argument('--create_subdir', dest='create_subdir', default="true")
parser.add_argument('--python-cache', dest='python_cache', default="false")

args = parser.parse_args()

if args.python_cache=='true':
    exit(0)

inputPath = args.inputPath
if inputPath[-1] != '/':
    inputPath += "/"

outputPath = args.outputPath
if outputPath[-1] != '/':
    outputPath += '/'

testName = args.testName

time_log.log('Data reading')

barcodeToGene = data.getBarcodeToGeneIds()
cols = data.getColumnNames(data.getInputFileNames(inputPath))

if params.USE_CPP_PARSER:
    experiments = data.getInputDataCppParser(inputPath, testName)
else:
    experiments, _ = data.getInputDataFromPath(inputPath)

time_log.log('Preprocessing')

if params.USE_CPP_PARSER:
    preprocessedData = experiments
else:
    preprocessedData = preprocessing.preprocessTestData(experiments, cols, barcodeToGene, True)
    preprocessedData.reset_index(inplace=True)

time_log.log('Prediction')

model = seg_models.getModel(predict_only=True)
model.load_weights('model.hdf5')
batchSize = 1024

valSequence = generators.DataSequence(preprocessedData, batch_size=batchSize, shuffle=False, generateY=False)
model.load_weights('model.hdf5')

predict = model.predict_generator(valSequence, verbose=1)

time_log.log('Postprocessing')
res = postprocessing.getResultTable(preprocessedData, predict)

time_log.log('Save results')
postprocessing.saveToFile(outputPath, testName, res)

time_log.log('Finish', True)
