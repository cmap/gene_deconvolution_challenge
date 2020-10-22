from keras import layers, models, initializers
import keras.backend as K
import tensorflow as tf
import math


def getModel():

    input_img = layers.Input((32, 2))

    dropoutValue = 0.1
    filters = 16
    conv1 = layers.Conv1D(filters, 3, activation='relu', padding='same')(input_img)
    conv1 = layers.Conv1D(filters, 3, activation='relu', padding='same')(conv1)
    pool1 = layers.MaxPool1D(pool_size=2)(conv1)  # 16
    pool1 = layers.Dropout(dropoutValue)(pool1)

    conv2 = layers.Conv1D(filters * 2, 3, activation='relu', padding='same')(pool1)
    conv2 = layers.Conv1D(filters * 2, 3, activation='relu', padding='same')(conv2)
    pool2 = layers.MaxPool1D(pool_size=2)(conv2)  # 8
    pool2 = layers.Dropout(dropoutValue)(pool2)

    conv3 = layers.Conv1D(filters * 4, 3, activation='relu', padding='same')(pool2)
    conv3 = layers.Conv1D(filters * 4, 3, activation='relu', padding='same')(conv3)
    pool3 = layers.MaxPool1D(pool_size=2)(conv3)  # 4
    pool3 = layers.Dropout(dropoutValue)(pool3)

    conv4 = layers.Conv1D(filters * 8, 3, activation='relu', padding='same')(pool3)
    conv4 = layers.Conv1D(filters * 8, 3, activation='relu', padding='same')(conv4)
    pool4 = layers.MaxPool1D(pool_size=2)(conv4)  # 2
    # pool4 = Dropout(dropoutValue)(pool4)

    conv5 = layers.Conv1D(filters * 16, 3, activation='relu', padding='same')(pool4)
    conv5 = layers.Conv1D(filters * 16, 3, activation='relu', padding='same')(conv5)

    flatten = layers.Flatten()(conv5)
    output = layers.Dense(2, activation='relu', name="mus")(flatten)

    # vs = layers.Dense(2, activation='relu')(flatten)
    # ws = layers.Dense(2, activation='relu')(flatten)
    #
    # def get_normal_sum(ws, mus, vs):
    #     vs = tf.maximum(vs, K.epsilon())
    #     def normal_sum(x):
    #         return tf.reduce_sum(ws * tf.exp(-0.5 * (tf.tile([[x]], tf.shape(mus)) - mus) ** 2 / vs ** 2) / tf.sqrt(2 * math.pi) / vs, axis=-1)
    #
    #     return normal_sum
    #
    # probs = layers.Lambda(lambda l: tf.transpose(tf.map_fn(get_normal_sum(l[0], l[1], l[2]), tf.range(0.5/32, 1, 1.0/32)), [1, 0]),
    #                       output_shape=(32,), name="kl")([ws, output, vs])

    bconv4 = layers.Conv1D(filters*16, 3, activation='relu', padding='same')(pool4)
    bconv4 = layers.Conv1D(filters*16, 3, activation='relu', padding='same')(bconv4)
    bup4 = layers.UpSampling1D(2)(bconv4)

    bconv3 = layers.Conv1D(filters*8, 3, activation='relu', padding='same')(bup4)
    bconv3 = layers.Conv1D(filters*8, 3, activation='relu', padding='same')(bconv3)
    bup3 = layers.UpSampling1D(2)(bconv3)

    bconv2 = layers.Conv1D(filters*4, 3, activation='relu', padding='same')(bup3)
    bconv2 = layers.Conv1D(filters*4, 3, activation='relu', padding='same')(bconv2)
    bup2 = layers.UpSampling1D(2)(bconv2)

    bconv1 = layers.Conv1D(filters*2, 3, activation='relu', padding='same')(bup2)
    bconv1 = layers.Conv1D(filters*2, 3, activation='relu', padding='same')(bconv1)
    bup1 = layers.UpSampling1D(2)(bconv1)

    probs = layers.Conv1D(filters, 3, activation='relu', padding='same')(bup1)
    probs = layers.Conv1D(1, 1, activation='relu', padding='same')(probs)
    probs = layers.Flatten(name='recover')(probs)



    model = models.Model(input_img, [output, probs])

    return model