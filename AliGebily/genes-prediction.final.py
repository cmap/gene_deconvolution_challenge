# %% import libraries and set input parameters
from time import time
start_time = time()

from os import listdir
from os.path import isfile, join, dirname, realpath
from pandas import DataFrame
import xgboost as xgb
from pickle import load
from numpy import percentile, nan, append, exp
from sys import argv
from multiprocessing import cpu_count, Pool
from statistics import mean


class gene_map(object):
    __slots__ = ('g0', 'g1')
    def __init__(self):
        self.g0 = None
        self.g1 = None


class gene_info(object):
    __slots__ = ('prop', 'data')
    def __init__(self, prop, data):
        self.prop = prop
        self.data = data


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

dspath = arguments['--dspath']
out = arguments['--out']
plate = arguments['--plate']
print('dspath, out, plate: ', dspath, out, plate)

dspath_parent_directory = dirname(dspath)

# ------- define input paramters
resources_directory = join(script_path, 'resources')
print(resources_directory)
separator = '	'

percentile_size = 1  # represents %, and it should be the same used while training
percentiles_arr = []
percentiles_count = int(100/percentile_size + 1)  # add boundaries
data_columns = ['count', 'mean', 'variance']
for i in range(percentiles_count):
    percentiles_arr.append(i*percentile_size)
    col = 'per'+str(i*percentile_size)
    data_columns.append(col)
data_columns.append('dpk_litmus')
data_columns.append('y1_larger_than_y0_flag1')
data_columns.append('y1_larger_than_y0_flag2')

dpk_litmus = nan
if(plate.lower().startswith('dpk')):
    dpk_litmus = 0
    print('dpk_litmus=0')
elif(plate.lower().startswith('litmus')):
    dpk_litmus = 1
    print('dpk_litmus=1')


# ------------------------------------------------------------------------------------

# ------- Load barcode_to_gene_map
data_barcode_to_gene_map = {}
barcode_to_gene_map_path = join(
    dspath_parent_directory, 'barcode_to_gene_map.txt')
if(isfile(barcode_to_gene_map_path) == False):  # if not found, find it in resources folders
    barcode_to_gene_map_path = join(
        resources_directory, 'barcode_to_gene_map.txt')

barcode_to_gene_map_stream = open(barcode_to_gene_map_path, 'r')
print(barcode_to_gene_map_stream)
index = -1
for line in barcode_to_gene_map_stream:
    index = index+1
    if(index == 0):  # file headers
        continue
    parts = line[:-1].split(separator)
    barcode_id = parts[0]
    gene_id = parts[1]
    # indicates which gene of the pair is mixed in higher (=1) or lower (=0) proportion.
    high_prop = int(parts[2])
    temp_barcode_to_gene_map = None
    if(barcode_id not in data_barcode_to_gene_map):
        temp_barcode_to_gene_map = gene_map()
        data_barcode_to_gene_map[barcode_id] = temp_barcode_to_gene_map
    else:
        temp_barcode_to_gene_map = data_barcode_to_gene_map[barcode_id]
    if high_prop == 0:
        temp_barcode_to_gene_map.g0 = gene_id
    else:
        temp_barcode_to_gene_map.g1 = gene_id
barcode_to_gene_map_stream.close()
print(len(data_barcode_to_gene_map))
# ------------------------------------------------------------------------------------
modelY0File = join(resources_directory, r'xgboost-percentiles-y0.dat')
loaded_model_y0 = load(open(modelY0File, 'rb'))
modelY1File = join(resources_directory, r'xgboost-percentiles-y1.dat')
loaded_model_y1 = load(open(modelY1File, 'rb'))
experimentId_start_index = len(plate)+1


def predict(experimentFiles):
    experimentIds = []
    data_all_experiments = {}
    genes_values = {}  # row per gene
    experimentFileIndex = -1
    for experimentFile in experimentFiles:
        experimentFileIndex += 1
        experimentMeasurementFileStream = open(join(dspath, experimentFile), 'r')
        # -4 represents '.txt'
        experimentId = experimentFile[experimentId_start_index:-4]
        experimentIds.append(experimentId)
        currentExperimentMeasurements = {}
        data_all_experiments[experimentId] = currentExperimentMeasurements
        measurement_index = -1
        for line in experimentMeasurementFileStream:
            measurement_index = measurement_index + 1
            if(measurement_index == 0):  # data headers
                continue
            parts = line[:-1].split(separator)
            barcode_id = parts[0]  # measurement id
            FI = round(100*float(parts[1]))  # observation, scale to 10
            if(barcode_id not in currentExperimentMeasurements):
                currentExperimentMeasurements[barcode_id] = []
            currentExperimentMeasurements[barcode_id].append(FI)
        # reformat data to be predicted
        for barcode_id in currentExperimentMeasurements.keys():
            if (barcode_id not in data_barcode_to_gene_map):
                print('barcode_id: ', barcode_id, 'in experiment ', experimentFile,
                      'is skipped because it has no mapping data with genes in barcode_to_gene_map.txt')
                continue

            genes = data_barcode_to_gene_map[barcode_id]
            # skip barcodes that don't have the two genes
            if(genes.g0 == None or genes.g1 == None):
                # print('barcode_id: ', barcode_id, 'in experiment ', experimentFile,
                #       'is skipped because it deosn't have the two genes')
                continue

            observations = currentExperimentMeasurements[barcode_id]
            observations.sort()
            # assume 50%  percent increase is outlier
            while(observations[0] < 1):
                del observations[0]
            while(observations[-1] > (1.5 * observations[-2])):
                del observations[-1]

            currentExperimentMeasurements[barcode_id] = observations
            observations.sort()

            nobs = len(observations)
            mu = round(mean(observations))
            variance = 0
            for ob in observations:
                variance += ((ob-mu)**2)
            variance = round(variance/nobs)
            data_row = [nobs, mu, variance]
            currentExperimentMeasurements[barcode_id] = percentile(
                observations, percentiles_arr)
            observations_percentiles_values = currentExperimentMeasurements[barcode_id]
            for i in range(percentiles_count):
                observations_percentiles_values[i] = int(observations_percentiles_values[i])
                data_row.append(observations_percentiles_values[i])
            data_row.append(dpk_litmus)

            a16 = (
                observations_percentiles_values[16]+observations_percentiles_values[17])/2
            a50 = (
                observations_percentiles_values[49]+observations_percentiles_values[50]+observations_percentiles_values[51])/3
            a83 = (
                observations_percentiles_values[83]+observations_percentiles_values[84])/2
            y1_larger_than_y0 = ((a50-a16)/a16) / ((a83-a50)/a83)
            data_row.append(y1_larger_than_y0 if y1_larger_than_y0 >
                            1 else (-1/y1_larger_than_y0))
            data_row.append((mu-a50)/a50)
            g0 = genes.g0
            g1 = genes.g1
            if(g0 not in genes_values):
                genes_values[g0] = []
            if(g1 not in genes_values):
                genes_values[g1] = []

            genes_values[g0].append(gene_info(0, data_row))
            genes_values[g1].append(gene_info(1, data_row))

            # print(experimentId, ':', g0, '=', y0, ', ', g1, '=', y1)
        experimentMeasurementFileStream.close()
        print('experimentId:', experimentId, 'total observations:', measurement_index +
              1, 'total measurements', len(data_all_experiments[experimentId]))  # 490
    # print('data_all_experiments count:', len(data_all_experiments))
    # ------------------------------------------------------------------------------------

    print('................Calculating values................')
    # create dataframe to be predicted
    test0_rows = []
    test1_rows = []
    for gene_id, values in genes_values.items():
        for value in values:
            if(value.prop == 0):
                test0_rows.append(value.data)
            else:
                test1_rows.append(value.data)

    # ------------- Load trained models, and predict
    print('low-prop calculations')
    test0_dataset = DataFrame(test0_rows, columns=data_columns)
    y0_pred = loaded_model_y0.predict(test0_dataset)
    y0_pred = [round(exp(value)) for value in y0_pred]

    print('high-prop calculations')

    test1_dataset = DataFrame(test1_rows, columns=data_columns)
    y1_pred = loaded_model_y1.predict(test1_dataset)
    y1_pred = [round(exp(value)) for value in y1_pred]
    print('preparing data for saving')
    y0_pred_index = 0
    y1_pred_index = 0
    for gene_id, values in genes_values.items():
        new_values = []
        for value in values:
            if(value.prop == 0):
                new_values.append(y0_pred[y0_pred_index])
                y0_pred_index += 1
            else:
                new_values.append(y1_pred[y1_pred_index])
                y1_pred_index += 1
        genes_values[gene_id] = new_values
    return experimentIds, genes_values


if __name__ == '__main__':
    # ------- load experiments measurements from files
    all_experiment_files = listdir(dspath)
    all_experiment_files.sort(reverse=False)
    # all_experiment_files = all_experiment_files[:74]
    # print(all_experiment_files)
    threads_count = cpu_count()
    po = Pool(threads_count)
    all_files_parts = []
    step = int(len(all_experiment_files)/threads_count)
    if(len(all_experiment_files) % threads_count > 0):  # step contains decimal point
        step += 1

    for i in range(threads_count):
        if(i < threads_count-1):
            all_files_parts.append(all_experiment_files[i*step: (i+1)*step])
        else:
            all_files_parts.append(
                all_experiment_files[i*step: len(all_experiment_files)])

    results = po.map_async(predict,
                           ((item)
                            for item in all_files_parts)
                           # this will zip in one iterable object
                           ).get()  # get will start the processes and execute them
    po.terminate()  # kill the spawned processes

    all_genes_values = {}
    all_experimentIds = []
    for thread_experimentIds, thread_genes_values in results:
        for experimentId in thread_experimentIds:
            all_experimentIds.append(experimentId)
        for gene_id, values in thread_genes_values.items():
            if(gene_id not in all_genes_values):
                all_genes_values[gene_id] = []
            gene_values = all_genes_values[gene_id]
            for value in values:
                gene_values.append(value)

    print('loaded experiments', all_experimentIds)
    # ------------------------------------------------------------------------------------

    # Write data to gct
    # order of data: any order is fine, the scorer code anyway re-sorts your solution by gene IDs in the first column to match the order given in the ground truth.
    print('................Saving values to gct................')
    duoFile = join(out, plate+'.gct')
    f = open(duoFile, 'w+')
    f.writelines('#1.3\n')  # version

    gct_rows_count = len(all_genes_values)
    print('rows count', gct_rows_count)

    columns_rows_count = len(all_experimentIds)
    print('columns_rows_count', columns_rows_count)

    # number of rows  ,number of columns, number of row metadata fields, the number of column metadata fields.
    f.write(separator.join(
        [str(gct_rows_count), str(columns_rows_count), '0', '0\n']))

    f.write(separator.join(append(['id'], all_experimentIds)))
    f.write('\n')
    for gene_id, values in all_genes_values.items():
        f.write(separator.join(append([gene_id], values)))
        f.write('\n')
    f.close()
    # ------------------------------------------------------------------------------------

    seconds = ((time() - start_time))
    print('seconds', seconds, 'minutes', round(seconds/60, 3))
    print('Finished')
