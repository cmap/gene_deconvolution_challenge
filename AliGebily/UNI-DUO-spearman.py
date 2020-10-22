import time
start_time = time.time()


# %% import libraries and set input parameters
from os import listdir
from os.path import isfile, join, dirname, realpath
from sys import argv

import numpy as np
from scipy import stats
from os import listdir
from os.path import isfile, join
from numpy import percentile




script_path = dirname(realpath(__file__))
print('args', argv)

arguments_names = []
arguments_values = []

for i in range(1, len(argv)):
    if(argv[i].startswith('--')):
        arguments_names.append(argv[i].lower())
    else:
        arguments_values.append(argv[i])
if(len(arguments_names) != len(arguments_values)):
    print('Number of input paramters is not matched with number of provided values')
    exit()

arguments = {}
for i in range(len(arguments_names)):
    arguments[arguments_names[i]] = arguments_values[i]

plate = arguments['--plate']
predictionspath = arguments['--predictionspath']
groundpath = arguments['--groundpath']

print('predictionspath, groundpath, plate: ', predictionspath, groundpath,  plate)

ground_truth_file = join(groundpath, plate+'_DECONV_UNI.gct')
predictions_file = join(predictionspath, plate+'.gct')

experiment_count = 10050  # 376
data_ground_truth = {}
ground_truth_stream = open(ground_truth_file, "r")
print(ground_truth_stream)
ground_truth_headers = []
index = -1
for line in ground_truth_stream:
    index = index+1
    if(index < 3):  # version, info, headers
        if(index == 2):  # headers
            ground_truth_headers = line[:-1].split("	")
    else:
        parts = line[:-1].split("	")
        gene_id = parts[0]
        data_ground_truth[gene_id] = [
            float(value) for value in parts[1:experiment_count]]
ground_truth_stream.close()
print(len(data_ground_truth))


data_predictions = {}
predictions_stream = open(predictions_file, "r")
print(predictions_stream)
predictions_headers = []
index = -1
for line in predictions_stream:
    index = index+1
    if(index < 3):  # version, info, headers
        if(index == 2):  # headers
            predictions_headers = line[:-1].split("	")
    else:
        parts = line[:-1].split("	")
        gene_id = parts[0]
        data_predictions[gene_id] = [float(value)
                                     for value in parts[1:experiment_count]]
predictions_stream.close()
print(len(data_predictions))

# merge ground with predictions
ground_truth_predictions = {}
for gene_di, values in data_ground_truth.items():
    ground_truth_predictions[gene_di] = [
        values[:len(data_predictions[gene_di])], data_predictions[gene_di]]

import scipy.stats as stats
import statistics as statistics
means = []
for gene, values in ground_truth_predictions.items():
    means.append(statistics.mean(values[0]))

from numpy import percentile
means_medians = percentile(means, [33, 66])
# print(means_medians)
all_correlations = []
high_correlations = []
medium_correlations = []
low_correlations = []
index = 0
for gene, values in ground_truth_predictions.items():
    c = stats.spearmanr(
        values[0], values[1]).correlation
    all_correlations.append(c)
    m=statistics.mean(values[0])
    if(m > means_medians[1]):
        high_correlations.append(c)
        # print("high", m)
    elif(m < means_medians[0]):
        medium_correlations.append(c)
        # print("low", m)
    else:
        low_correlations.append(c)
        # print("medium", m)

print(statistics.median(all_correlations))
# print(statistics.median(high_correlations))
# print(statistics.median(medium_correlations))
# print(statistics.median(low_correlations))
