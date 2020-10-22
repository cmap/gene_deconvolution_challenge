from keras import layers, models, activations, losses
import keras.backend as K
# import tensorflow as tf

from . import params


def _getDirectSegmentation(input):
    dropoutValue = 0.1
    filters = 16
    conv1 = layers.Conv1D(filters, 3, activation='relu', padding='same')(input)
    conv1 = layers.Conv1D(filters, 3, activation='relu', padding='same')(conv1)
    conv1 = layers.Dropout(dropoutValue)(conv1)

    conv2 = layers.Conv1D(filters, 3, activation='relu', padding='same')(conv1)
    conv2 = layers.Conv1D(filters, 3, activation='relu', padding='same')(conv2)
    conv2 = layers.Dropout(dropoutValue)(conv2)

    conv3 = layers.Conv1D(filters, 3, activation='relu', padding='same')(conv2)
    conv3 = layers.Conv1D(filters, 3, activation='relu', padding='same')(conv3)
    conv3 = layers.Dropout(dropoutValue)(conv3)

    conv4 = layers.Conv1D(filters, 3, activation='relu', padding='same')(conv3)
    conv4 = layers.Conv1D(filters, 3, activation='relu', padding='same')(conv4)

    conv5 = layers.Conv1D(filters, 3, activation='relu', padding='same')(conv4)
    conv5 = layers.Conv1D(filters, 3, activation='relu', padding='same')(conv5)

    final = conv5

    return final


def _getDenseDirectSegmentation(input):
    filters = 16

    conv1 = layers.Conv1D(filters, 3, activation='relu', padding="same")(input)
    conv1 = layers.Conv1D(filters, 3, activation='relu', padding="same")(conv1)
    conv1 = layers.Conv1D(filters, 3, activation='relu', padding="same")(conv1)
    conv1 = layers.Conv1D(filters, 3, activation='relu', padding="same")(conv1)

    conv2 = layers.Conv1D(filters, 3, activation='relu', padding="same")(conv1)
    conv2 = layers.Conv1D(filters, 3, activation='relu', padding="same")(conv2)
    conv2 = layers.Conv1D(filters, 3, activation='relu', padding="same")(conv2)
    conv2 = layers.Conv1D(filters, 3, activation='relu', padding="same")(conv2)
    conv2 = layers.concatenate([conv1, conv2])

    conv3 = layers.Conv1D(filters, 3, activation='relu', padding="same")(conv2)
    conv3 = layers.Conv1D(filters, 3, activation='relu', padding="same")(conv3)
    conv3 = layers.Conv1D(filters, 3, activation='relu', padding="same")(conv3)
    conv3 = layers.Conv1D(filters, 3, activation='relu', padding="same")(conv3)
    conv3 = layers.concatenate([conv2, conv3])

    conv4 = layers.Conv1D(filters, 3, activation='relu', padding="same")(conv3)
    conv4 = layers.Conv1D(filters, 3, activation='relu', padding="same")(conv4)
    conv4 = layers.Conv1D(filters, 3, activation='relu', padding="same")(conv4)
    conv4 = layers.Conv1D(filters, 3, activation='relu', padding="same")(conv4)
    conv4 = layers.concatenate([conv3, conv4])

    conv5 = layers.Conv1D(filters, 3, activation='relu', padding="same")(conv4)
    conv5 = layers.Conv1D(filters, 3, activation='relu', padding="same")(conv5)
    conv5 = layers.Conv1D(filters, 3, activation='relu', padding="same")(conv5)
    conv5 = layers.Conv1D(filters, 3, activation='relu', padding="same")(conv5)
    conv5 = layers.concatenate([conv4, conv5])

    full_conv = layers.Conv1D(filters*8, 32, activation='relu', padding='same')(conv5)

    return full_conv


def _getUnetSegmentation(input):
    dropoutValue = 0.1
    filters = 16
    conv1 = layers.Conv1D(filters, 3, activation='relu', padding='same')(input)
    conv1 = layers.Conv1D(filters, 3, activation='relu', padding='same')(conv1)
    pool1 = layers.MaxPool1D()(conv1)
    pool1 = layers.Dropout(dropoutValue)(pool1)

    conv2 = layers.Conv1D(filters*2, 3, activation='relu', padding='same')(pool1)
    conv2 = layers.Conv1D(filters*2, 3, activation='relu', padding='same')(conv2)
    pool2 = layers.MaxPool1D()(conv2)
    pool2 = layers.Dropout(dropoutValue)(pool2)

    conv3 = layers.Conv1D(filters*4, 3, activation='relu', padding='same')(pool2)
    conv3 = layers.Conv1D(filters*4, 3, activation='relu', padding='same')(conv3)
    pool3 = layers.MaxPool1D()(conv3)
    pool3 = layers.Dropout(dropoutValue)(pool3)

    conv4 = layers.Conv1D(filters*8, 3, activation='relu', padding='same')(pool3)
    conv4 = layers.Conv1D(filters*8, 3, activation='relu', padding='same')(conv4)
    pool4 = layers.MaxPool1D()(conv4)

    conv5 = layers.Conv1D(filters*16, 3, activation='relu', padding='same')(pool4)
    conv5 = layers.Conv1D(filters*16, 3, activation='relu', padding='same')(conv5)
    
    up4 = layers.UpSampling1D()(conv5)
    up4 = layers.concatenate([up4, conv4])
    up4 = layers.Conv1D(filters * 8, 3, activation='relu', padding='same')(up4)
    up4 = layers.Conv1D(filters * 8, 3, activation='relu', padding='same')(up4)

    up3 = layers.UpSampling1D()(up4)
    up3 = layers.concatenate([up3, conv3])
    up3 = layers.Conv1D(filters * 4, 3, activation='relu', padding='same')(up3)
    up3 = layers.Conv1D(filters * 4, 3, activation='relu', padding='same')(up3)

    up2 = layers.UpSampling1D()(up3)
    up2 = layers.concatenate([up2, conv2])
    up2 = layers.Conv1D(filters * 2, 3, activation='relu', padding='same')(up2)
    up2 = layers.Conv1D(filters * 2, 3, activation='relu', padding='same')(up2)

    up1 = layers.UpSampling1D()(up2)
    up1 = layers.concatenate([up1, conv1])
    up1 = layers.Conv1D(filters, 3, activation='relu', padding='same')(up1)
    up1 = layers.Conv1D(filters, 3, activation='relu', padding='same')(up1)

    final = up1
    return final


def _getPUnetSegmentation(input, levels=5, filters = 16, filters_multiplier = 2, dropoutValue=0.1):

    convs = []

    inp = input
    for i in range(levels):
        f = int(filters*(filters_multiplier**i))
        conv = layers.Conv1D(f, 3, activation='relu', padding='same')(inp)
        conv = layers.Conv1D(f, 3, activation='relu', padding='same')(conv)
        inp = conv
        if i<levels-1:
            pool = layers.MaxPool1D()(conv)
            pool = layers.Dropout(dropoutValue)(pool)
            inp = pool

        convs.append(conv)

    for i in range(levels-2, -1, -1):
        f = int(filters * (filters_multiplier ** i))
        up = layers.UpSampling1D()(inp)
        up = layers.concatenate([up, convs[i]])
        up = layers.Conv1D(f, 3, activation='relu', padding='same')(up)
        up = layers.Conv1D(f, 3, activation='relu', padding='same')(up)

        inp = up

    return inp


def _getPredict(input, name="predict", filters=128):
    precalc = layers.Conv1D(filters, 1, activation='relu', padding='same')(input)
    segment = layers.Conv1D(1, 1, padding='same')(precalc)
    segment = layers.Lambda(lambda x: activations.softmax(x, axis=-2))(segment)
    # segment = layers.Conv1D(1, 1, activation='sigmoid', padding='same')(precalc)
    shifts = layers.Conv1D(1, 1, padding='same')(precalc)
    shifts = layers.Activation(lambda x : activations.relu(x, max_value=1))(shifts)

    res = layers.concatenate([segment, shifts], name=name)
    return res


class PredictValue(layers.Layer):
    def __init__(self, bins_count, **kwargs):
        self.bins_count = float(bins_count)
        super().__init__(**kwargs)

    def call(self, x):
        # argmax = tf.argmax(x[:, :, 0], axis=1)
        # # shift = tf.gather(x[:, :, 1], argmax, axis=-1)
        # arange = K.arange(tf.shape(x)[0], dtype=argmax.dtype)
        # arange = tf.Print(arange, [arange, tf.shape(arange)], message="arange")
        # argmax = tf.Print(argmax, [argmax, tf.shape(argmax), tf.shape(arange)], message="argmax")
        # ids = tf.stack([arange, argmax])
        # ids = tf.Print(ids, [ids, tf.shape(ids)], message="ids")
        # shift = tf.gather_nd(x[:, :, 1], ids)
        # shift = tf.Print(shift, [shift, tf.shape(shift)], message="shift")
        # return K.reshape((tf.cast(argmax, dtype=shift.dtype)+shift)/self.bins_count, (tf.shape(x)[0], 1))

        # seg = x[:, :, 0]
        # shift = x[:, :, 1]
        # mask = tf.cast(tf.equal(seg, tf.reduce_max(seg, axis=1, keepdims=True)), dtype=x.dtype)
        # bin = tf.reduce_sum(K.arange(self.bins_count) * mask, axis=1)
        # shift = tf.reduce_sum(shift * mask, axis=1)
        # return K.reshape((bin+shift)/self.bins_count, (-1, 1))

        seg = x[:, :, 0]
        shift = x[:, :, 1]
        mask = K.cast(K.equal(seg, K.max(seg, axis=1, keepdims=True)), dtype=x.dtype)
        bin = K.sum(K.arange(self.bins_count) * mask, axis=1)
        shift = K.sum(shift * mask, axis=1)
        return K.reshape((bin + shift) / self.bins_count, (-1, 1))



    def compute_output_shape(self, input_shape):
        return (input_shape[0], 1)



def getModel(predict_only=False):
    inp = layers.Input((params.HIST_BINS,2))

    # segmentation = _getDirectSegmentation(inp)
    # segmentation = _getUnetSegmentation(inp)
    # segmentation = _getDenseDirectSegmentation(inp)
    segmentation = _getPUnetSegmentation(inp, levels=3, filters=32, filters_multiplier=1.5)

    low_segment = _getPredict(segmentation, name="lo", filters=64)
    hi_segment = _getPredict(segmentation, name="hi", filters=64)

    low_value = PredictValue(params.HIST_BINS)(low_segment)
    hi_value = PredictValue(params.HIST_BINS)(hi_segment)

    value = layers.concatenate([low_value, hi_value], name="v")

    if predict_only:
        model = models.Model(inp, value)
    else:
        model = models.Model(inp, [value, low_segment, hi_segment])

    return model


def seg_loss_cce(y_true, y_pred):
    return losses.categorical_crossentropy(y_true[:, :, 0], y_pred[:, :, 0])


def seg_loss_shift_loss(y_true, y_pred):
    y_true = y_true[:, :, 1]
    y_pred = y_pred[:, :, 1]

    return losses.mean_squared_error(y_true, y_pred)


def seg_loss_full(y_true, y_pred):
    return seg_loss_cce(y_true, y_pred) + seg_loss_shift_loss(y_true, y_pred)


def seg_loss_none(y_true, y_pred):
    return K.constant(0.0)
