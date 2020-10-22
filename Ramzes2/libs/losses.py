import keras.backend as K
import tensorflow as tf
from keras import losses


def mse(y_true, y_pred):
    return losses.mean_squared_error(y_true, y_pred)


def mse_mae(y_true, y_pred):
    return losses.mean_squared_error(y_true, y_pred) + losses.mean_absolute_error(y_true, y_pred)


def pearson(y_true, y_pred):
    x = y_true
    y = y_pred
    mx = K.mean(x)
    my = K.mean(y)
    xm, ym = x-mx, y-my
    r_num = K.sum(tf.multiply(xm,ym))
    r_den = K.sqrt(tf.multiply(K.sum(K.square(xm)), K.sum(K.square(ym))))
    r = r_num / r_den

    r = K.maximum(K.minimum(r, 1.0), -1.0)
    return 1 - K.square(r)


def pearson_mse(y_true, y_pred):
    return pearson(y_true, y_pred) + mse(y_true, y_pred)


def kl(y_true, y_pred):
    y_true = y_true / K.reshape(K.sum(y_true, axis=-1), (-1, 1))
    y_pred = y_pred / K.reshape(K.sum(y_pred, axis=-1), (-1, 1))
    y_true = K.clip(y_true, K.epsilon(), 1)
    y_pred = K.clip(y_pred, K.epsilon(), 1)
    return K.sum(y_true * K.log(y_true / y_pred), axis=-1)

def kl_revert(y_true, y_pred):
    y_true = y_true / K.reshape(K.sum(y_true, axis=-1), (-1, 1))
    y_pred = y_pred / K.reshape(K.sum(y_pred, axis=-1), (-1, 1))
    y_true = K.clip(y_true, K.epsilon(), 1)
    y_pred = K.clip(y_pred, K.epsilon(), 1)
    return K.sum(y_pred * K.log(y_pred / y_true), axis=-1)


def scaled_mse(y_true, y_pred):
    return losses.mean_squared_error(K.reshape(y_true[:,-1], (-1, 1))*y_true[:, :-1],
                                     K.reshape(y_true[:, -1], (-1, 1))*y_pred)/100000

def mse_no_scale(y_true, y_pred):
    return losses.mean_squared_error(y_true[:, :-1], y_pred)

def get_weighted_unscaled_scaled_mse(w_unscaled, w_scaled):
    def loss(y_true, y_pred):
        return mse_no_scale(y_true, y_pred)*w_unscaled + scaled_mse(y_true, y_pred)*w_scaled

    return loss
