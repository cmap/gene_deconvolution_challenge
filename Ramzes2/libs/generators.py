import numpy as np
import pandas as pd
from tqdm import tqdm

from keras.utils import Sequence
from . import preprocessing, params

class DataSequence(Sequence):
    def __init__(self, data, batch_size, shuffle=True, generateY=True, fakeX=None, fakeY=None, ae_version=False,
                 use_weights=False, add_scale=False, add_seg_target=False):
        self.use_weights = use_weights
        self.add_scale = add_scale
        self.add_seg_target = add_seg_target

        self.process_data(data, generateY)
        if fakeX is not None:
            self.X = np.concatenate((self.X, fakeX), axis=0)
            if generateY:
                self.Y = np.concatenate((self.Y, fakeY), axis=0)

        self.batch_size = batch_size
        self.shuffle = shuffle
        self.generateY = generateY
        self.ids = list(range(len(self.X)))
        self.ae_version = ae_version

        if shuffle:
            np.random.shuffle(self.ids)

    def process_seg_y(self):
        def seg_target(target):
            bin = int(np.minimum(np.floor(target*params.HIST_BINS), params.HIST_BINS-1))
            shift = target*params.HIST_BINS - bin

            res = np.zeros((params.HIST_BINS, 2))
            res[bin, 0] = 1
            res[:bin, 1] = 1
            res[bin, 1] = shift

            return res

        self.segY = [ [seg_target(pos) for pos in y] for y in tqdm(self.Y, desc="segmentation target")]
        self.segY = np.array(self.segY)

    def process_data(self, data, generateY):
        self.X = data[preprocessing.histColumns].values.astype(np.float32)
        self.weights = data['max_value'].values - data['min_value'].values
        # self.weights /= np.max(self.weights)
        if not self.use_weights:
            self.weights = np.ones_like(self.weights)
        if generateY:
            self.Y = data[['gene_pos_low_norm', 'gene_pos_high_norm']].values.astype(np.float32)
            if self.add_seg_target:
                self.process_seg_y()
            if self.add_scale:
                distance = data['max_value'].values - data['min_value'].values
                self.Y = np.hstack((self.Y, distance.reshape(-1, 1)))

    def __len__(self):
        return int(np.ceil(len(self.ids) / float(self.batch_size)))

    def __getitem__(self, idx):
        batch_ids = self.ids[idx * self.batch_size: min((idx + 1) * self.batch_size, len(self.ids))]

        batchX = self.X[batch_ids]

        pos = np.linspace(0.0, 1.0, batchX.shape[1]+1, dtype=np.float32)[:-1]
        pos = np.repeat(pos.reshape(1, -1), len(batch_ids), axis=0)

        batchX = np.stack([batchX, pos], axis=-1)

        if self.generateY:
            batchY = [self.Y[batch_ids]]
            if self.add_seg_target:
                batchY += [self.segY[batch_ids, 0], self.segY[batch_ids, 1]]
            if self.ae_version:
                batchY += self.X[batch_ids]
            # return batchX, batchY, self.weights[batch_ids]
            return batchX, batchY
        else:
            return batchX

    def on_epoch_end(self):
        if self.shuffle:
            np.random.shuffle(self.ids)
