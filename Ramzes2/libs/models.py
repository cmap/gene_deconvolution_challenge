import numpy as np
import pandas as pd

from keras.models import Model
from keras.layers import Input, Conv1D, Dropout, MaxPool1D, Flatten, Dense

from . import params

def getModel():
    dropoutValue = 0.1
    filters = 16
    inputs = Input((params.HIST_BINS, 2))
    conv1 = Conv1D(filters, 3, activation='relu', padding='same')(inputs)
    conv1 = Conv1D(filters, 3, activation='relu', padding='same')(conv1)
    pool1 = MaxPool1D(pool_size=2)(conv1) # 16
    pool1 = Dropout(dropoutValue)(pool1)

    conv2 = Conv1D(filters*2, 3, activation='relu', padding='same')(pool1)
    conv2 = Conv1D(filters*2, 3, activation='relu', padding='same')(conv2)
    pool2 = MaxPool1D(pool_size=2)(conv2)  # 8
    pool2 = Dropout(dropoutValue)(pool2)

    conv3 = Conv1D(filters*4, 3, activation='relu', padding='same')(pool2)
    conv3 = Conv1D(filters*4, 3, activation='relu', padding='same')(conv3)
    pool3 = MaxPool1D(pool_size=2)(conv3)  # 4
    pool3 = Dropout(dropoutValue)(pool3)

    conv4 = Conv1D(filters*8, 3, activation='relu', padding='same')(pool3)
    conv4 = Conv1D(filters*8, 3, activation='relu', padding='same')(conv4)
    pool4 = MaxPool1D(pool_size=2)(conv4)  # 2
    # pool4 = Dropout(dropoutValue)(pool4)

    conv5 = Conv1D(filters*16, 3, activation='relu', padding='same')(pool4)
    conv5 = Conv1D(filters*16, 3, activation='relu', padding='same')(conv5)

    flatten = Flatten()(conv5)
    output = Dense(2, activation='relu')(flatten)

    model = Model(inputs, output)

    return model
