
# %%
from os import listdir
from os.path import isfile, join, dirname, realpath
from sys import argv
import pandas as pd
import xgboost as xgb
import pickle
import time
import multiprocessing
from multiprocessing import Pool
import numpy as np


def create_xgb_model(params):
    X, y, saveTo = params
    model = xgb.XGBRegressor(max_depth=10, n_estimators=100,
                             learning_rate=.05, silent=True, n_jobs=2)
    print(model)
    model.fit(X, y)
    print('saving model......')
    pickle.dump(model, open(saveTo, "wb"))
    return model


def skip_rows(index):
    return False# train with full data without skipping any rows
    # return index % 11 > 0# train using only 10% of data for faster development. In production, this line should be commented and previous line(return False) should uncommented


if __name__ == '__main__':
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

    csvspath = arguments['--csvspath']
    plates = arguments['--plates'].split(',')
    print('csvspath, plates: ', csvspath, plates)

    models_directory = join(script_path, 'resources')

    print('loading data......')
    dataset = None
    for plate in plates:
        plate_dataset = pd.read_csv(join(csvspath, plate+'.csv'),
                                    skiprows=skip_rows)
        if(plate.lower().startswith('dpk')):
            plate_dataset['dpk_litmus'] = 0
        elif(plate.lower().startswith('litmus')):
            plate_dataset['dpk_litmus'] = 1
        else:
            plate_dataset['dpk_litmus'] = np.nan
        print(plate)
        print(plate_dataset.head(2))
        if(dataset is None):
            dataset = plate_dataset
        else:
            dataset = dataset.append(plate_dataset)
    
    print(dataset.shape)

    a16 = (dataset['per16']+dataset['per17'])/2
    a50 = (dataset['per49']+dataset['per50']+dataset['per51'])/3
    a83 = (dataset['per83']+dataset['per84'])/2
    dataset['y1_larger_than_y0_flag1'] = ((a50-a16)/a16) / ((a83-a50)/a83)
    dataset['y1_larger_than_y0_flag2'] = (dataset['mean'] - a50)/a50

    dataset['y1_larger_than_y0_flag1'] = [
        x if x > 1 else (-1/x) for x in dataset['y1_larger_than_y0_flag1']]
    
    target0 = np.log(dataset['y0'])
    target1 = np.log(dataset['y1'])
    # train = dataset.drop(['y0', 'y1','per0','per1','per2','per98','per99','per100'], axis=1)
    train = dataset.drop(['y0', 'y1'], axis=1)
    # exit()
    # fit model with training data
    print('fitting models......')
    modelY0File = join(models_directory, 'xgboost-percentiles-y0.dat')
    modelY1File = join(models_directory, 'xgboost-percentiles-y1.dat')
    threads_count = 2  # multiprocessing.cpu_count()
    po = Pool(threads_count)
    results = po.map_async(create_xgb_model,
                           ((args[0], args[1], args[2]) for args in [
                               [train, target0, modelY0File], [train, target1, modelY1File]])).get()  # get will start the processes and execute them
    po.terminate()  # kill the spawned processes

    print("Finished training")
    print("--- %s seconds ---" % ((time.time() - start_time)))
