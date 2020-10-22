

# %% import libraries and set input parameters
import numpy as np
import os
from os import listdir
from os.path import isfile, join, dirname, realpath
from sys import argv
from numpy import percentile
import time
from statistics import mean
import multiprocessing
from multiprocessing import Pool

start_time = time.time()


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

inputpath = arguments['--inputpath']
csvspath = arguments['--csvspath']
plates = arguments['--plates'].split(',')
groundpath = arguments['--groundpath']
print('inputpath, groundpath, csvspath, plates: ', inputpath, groundpath,  csvspath, plates)

# ------- define input paramters
resources_directory = join(script_path, 'resources')
print(resources_directory)

# # Load barcode_to_gene_map
data_barcode_to_gene_map = {}
barcode_to_gene_map_path = join(resources_directory, 'barcode_to_gene_map.txt')
barcode_to_gene_map_stream = open(barcode_to_gene_map_path, 'r')

print(barcode_to_gene_map_stream)
index = -1
for line in barcode_to_gene_map_stream:
    index = index+1
    if(index == 0):  # file headers
        continue
    parts = line[:-1].split("	")
    barcode_id = parts[0]
    gene_id = parts[1]
    # indicates which gene of the pair is mixed in higher (=1) or lower (=0) proportion.
    high_prop = int(parts[2])
    temp_barcode_to_gene_map = None
    if(barcode_id not in data_barcode_to_gene_map):
        temp_barcode_to_gene_map = {}
        data_barcode_to_gene_map[barcode_id] = temp_barcode_to_gene_map
    else:
        temp_barcode_to_gene_map = data_barcode_to_gene_map[barcode_id]
    if high_prop == 0:
        temp_barcode_to_gene_map["g0"] = gene_id
    else:
        temp_barcode_to_gene_map["g1"] = gene_id
barcode_to_gene_map_stream.close()
print(len(data_barcode_to_gene_map))


def generate_csv(plate):
    ground_truth_file = join(groundpath, plate+'_DECONV_UNI.gct')
    percentile_size = 1  # represents %
    # # load experiments measurements from files
    platePath = join(inputpath, plate)
    onlyfiles = [experimentFile for experimentFile in listdir(platePath)]
    data_all_experiments = {}
    experimentId_start_index = len(plate)+1
    experimentFileIndex = -1
    for experimentFile in onlyfiles:
        experimentFileIndex = experimentFileIndex+1
        experimentMeasurementFileStream = open(
            join(platePath, experimentFile), "r")
        # -4 represents ".txt"
        experimentId = experimentFile[experimentId_start_index:-4]
        currentExperimentMeasurements = {}
        data_all_experiments[experimentId] = currentExperimentMeasurements
        experimentIndex = -1
        for line in experimentMeasurementFileStream:
            experimentIndex = experimentIndex + 1
            if(experimentIndex == 0):  # data headers
                continue
            parts = line[:-1].split("	")
            barcode_id = parts[0]  # measurement id
            FI = round(100*float(parts[1]))  # observation
            if(barcode_id not in currentExperimentMeasurements):
                currentExperimentMeasurements[barcode_id] = []
            currentExperimentMeasurements[barcode_id].append(FI)
        experimentMeasurementFileStream.close()
        for barcode_id, observations in currentExperimentMeasurements.items():
            observations = currentExperimentMeasurements[barcode_id]
            observations.sort()
            # assume 50%  percent increase is outlier
            while(observations[0] < 1):
                # print(observations)
                # print("lower", len(observations),
                #         observations[0], observations[1])
                del observations[0]
                # print("lower2", len(observations),
                #         observations[0], observations[1])
            while(observations[-1] > (1.5 * observations[-2])):
                # print(observations)
                # print("upper", len(observations),
                #         observations[-2], observations[-1])
                del observations[-1]
                # print("upper2", len(observations),
                #         observations[-2], observations[-1])
            currentExperimentMeasurements[barcode_id] = observations
        print("experimentId:", experimentId, "total observations:", experimentIndex +
              1, "total measurements", len(data_all_experiments[experimentId]))  # 490
        # if(experimentFileIndex == 2):
        #     break
    print("data_all_experiments count:", len(data_all_experiments))
    # # Load ground truth
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
            data_ground_truth[gene_id] = parts
    ground_truth_stream.close()
    print(len(data_ground_truth))

    # merge data int tabular format
    data_table_columns = ['count', 'mean', 'variance']
    data_table = [data_table_columns]

    percentiles_arr = []
    percentiles_count = int(100/percentile_size + 1)  # add boundaries
    for i in range(percentiles_count):
        percentiles_arr.append(i*percentile_size)
        data_table_columns.append('per'+str(i*percentile_size))

    data_table_columns.append('y0')
    data_table_columns.append('y1')

    # 0 index in headers is the id, so skip it
    for headerIndex in range(1, len(ground_truth_headers)):
        experimentId = ground_truth_headers[headerIndex]
        if(experimentId in data_all_experiments):
            experimentMeasurements = data_all_experiments[experimentId]
            for barcode_id, observations in experimentMeasurements.items():
                measurment = data_barcode_to_gene_map[barcode_id]
                if('g0' in measurment and 'g1' in measurment):
                    g0 = measurment["g0"]
                    g1 = measurment["g1"]

                    nobs = len(observations)
                    mu = round(mean(observations))
                    variance = 0
                    for ob in observations:
                        variance += ((ob-mu)**2)
                    variance = round(variance/nobs)
                    data_row = [nobs, mu, variance]

                    observations_percentiles_values = percentile(
                        observations, percentiles_arr)
                    for i in range(percentiles_count):
                        observations_percentiles_values[i] = int(
                            observations_percentiles_values[i])
                        data_row.append(observations_percentiles_values[i])
                    y0 = data_ground_truth[g0][headerIndex]
                    y1 = data_ground_truth[g1][headerIndex]
                    data_row.append(y0)
                    data_row.append(y1)
                    data_table.append(data_row)
        print(experimentId)

    # write data to file
    import csv
    print("Writing to csv")
    with open(join(csvspath, plate+".csv"), "w+") as my_csv:
        csvWriter = csv.writer(my_csv, delimiter=',', lineterminator='\n')
        csvWriter.writerows(data_table)


if __name__ == '__main__':
    # threads_count = len(plates)  # multiprocessing.cpu_count()
    # po = Pool(threads_count)
    # results = po.map_async(generate_csv,
    #                        ((plate) for plate in plates)).get()  # get will start the processes and execute them
    # po.terminate()  # kill the spawned processes
    for plate in plates:
        generate_csv(plate)

    print("Finished training")
    print("--- %s seconds ---" % ((time.time() - start_time)))
