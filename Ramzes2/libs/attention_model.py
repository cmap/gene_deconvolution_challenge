from keras import layers, backend, models, initializers


def _resnext_block(input, filters, count, filters_small):
    res = []
    for i in range(count):
        c = layers.Conv1D(filters_small, 1, padding='same')(input)
        c = layers.Conv1D(filters_small, 3, padding='same')(c)
        res.append(c)
    res = layers.concatenate(res)
    res = layers.Conv1D(filters, 1, padding='same', activation='relu')(res)
    res = layers.add([res, input])
    return res


def _getConvPart(input):
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

    # filters = 32
    # count = 4
    # filters_small = 4
    #
    # conv1 = layers.Conv1D(filters, 3, activation='relu', padding="same")(input)
    # conv2 = _resnext_block(conv1, filters, count, filters_small)
    # conv3 = _resnext_block(conv2, filters, count, filters_small)
    #
    # ful_conv = layers.Conv1D(filters, 32, activation='relu', padding='same')(conv3)
    # return ful_conv

def _getCoords_small(conv, coords, input_count=32):
    coord_conv = layers.Conv1D(1, 1, padding='same')(conv)
    coord_conv = layers.Flatten()(coord_conv)
    coord_conv = layers.Activation('softmax')(coord_conv)
    coord_head = layers.multiply([coords, coord_conv])

    # coord_conv_2 = layers.Conv1D(1, 1, activation='tanh', padding='same')(conv)
    # coord_conv_2 = layers.GlobalAveragePooling1D()(coord_conv_2)
    # coord_head = layers.add([coord_head, coord_conv_2])

    coord_head = layers.Dense(1, activation='relu', use_bias=True, kernel_initializer=initializers.ones())(coord_head)

    return coord_head

def _getCoordsPart(conv, coords, input_count=32):

    coord_head = _getCoords_small(conv, coords, input_count)

    right = layers.Lambda(lambda x: x+1.0/input_count)(coords)
    right_coord_head = _getCoords_small(conv, right, input_count)

    coord_head = layers.concatenate([coord_head, right_coord_head])
    coord_head = layers.Dense(1, use_bias=False, kernel_initializer=initializers.constant(0.5))(coord_head)

    return coord_head


def getModel():

    img_input = layers.Input(shape=(32, 2))

    hist = layers.Lambda(lambda x: x[:, :, 0:1])(img_input)
    coords = layers.Lambda(lambda x: x[:, :, 1])(img_input)

    conv = _getConvPart(hist)

    hi_head = _getCoordsPart(conv, coords)
    low_head = _getCoordsPart(conv, coords)

    output = layers.concatenate([low_head, hi_head])

    model = models.Model(img_input, output)

    return model
