from keras import layers, backend, models, initializers

def identity_block(input_tensor, kernel_size, filters, stage, block):
    """The identity block is the block that has no conv layer at shortcut.

    # Arguments
        input_tensor: input tensor
        kernel_size: default 3, the kernel size of
            middle conv layer at main path
        filters: list of integers, the filters of 3 conv layer at main path
        stage: integer, current stage label, used for generating layer names
        block: 'a','b'..., current block label, used for generating layer names

    # Returns
        Output tensor for the block.
    """
    filters1, filters2, filters3 = filters
    bn_axis = -1
    conv_name_base = 'res' + str(stage) + block + '_branch'
    bn_name_base = 'bn' + str(stage) + block + '_branch'

    x = layers.Conv1D(filters1, (1, ),
                      kernel_initializer='he_normal',
                      name=conv_name_base + '2a')(input_tensor)
    x = layers.BatchNormalization(axis=bn_axis, name=bn_name_base + '2a')(x)
    x = layers.Activation('relu')(x)

    x = layers.Conv1D(filters2, kernel_size,
                      padding='same',
                      kernel_initializer='he_normal',
                      name=conv_name_base + '2b')(x)
    x = layers.BatchNormalization(axis=bn_axis, name=bn_name_base + '2b')(x)
    x = layers.Activation('relu')(x)

    x = layers.Conv1D(filters3, (1, ),
                      kernel_initializer='he_normal',
                      name=conv_name_base + '2c')(x)
    x = layers.BatchNormalization(axis=bn_axis, name=bn_name_base + '2c')(x)

    x = layers.add([x, input_tensor])
    x = layers.Activation('relu')(x)
    return x


def conv_block(input_tensor,
               kernel_size,
               filters,
               stage,
               block,
               strides=(2, )):
    """A block that has a conv layer at shortcut.

    # Arguments
        input_tensor: input tensor
        kernel_size: default 3, the kernel size of
            middle conv layer at main path
        filters: list of integers, the filters of 3 conv layer at main path
        stage: integer, current stage label, used for generating layer names
        block: 'a','b'..., current block label, used for generating layer names
        strides: Strides for the first conv layer in the block.

    # Returns
        Output tensor for the block.

    Note that from stage 3,
    the first conv layer at main path is with strides=(2, 2)
    And the shortcut should have strides=(2, 2) as well
    """
    filters1, filters2, filters3 = filters
    bn_axis = -1
    conv_name_base = 'res' + str(stage) + block + '_branch'
    bn_name_base = 'bn' + str(stage) + block + '_branch'

    x = layers.Conv1D(filters1, (1, ), strides=strides,
                      kernel_initializer='he_normal',
                      name=conv_name_base + '2a')(input_tensor)
    x = layers.BatchNormalization(axis=bn_axis, name=bn_name_base + '2a')(x)
    x = layers.Activation('relu')(x)

    x = layers.Conv1D(filters2, kernel_size, padding='same',
                      kernel_initializer='he_normal',
                      name=conv_name_base + '2b')(x)
    x = layers.BatchNormalization(axis=bn_axis, name=bn_name_base + '2b')(x)
    x = layers.Activation('relu')(x)

    x = layers.Conv1D(filters3, (1, ),
                      kernel_initializer='he_normal',
                      name=conv_name_base + '2c')(x)
    x = layers.BatchNormalization(axis=bn_axis, name=bn_name_base + '2c')(x)

    shortcut = layers.Conv1D(filters3, (1, ), strides=strides,
                             kernel_initializer='he_normal',
                             name=conv_name_base + '1')(input_tensor)
    shortcut = layers.BatchNormalization(
        axis=bn_axis, name=bn_name_base + '1')(shortcut)

    x = layers.add([x, shortcut])
    x = layers.Activation('relu')(x)
    return x

def SEBlock(input, rate=2, trainable=True):
    spatial = layers.Conv1D(1, (1, ), trainable=trainable, kernel_initializer=initializers.ones())(input)
    spatial = layers.Activation("sigmoid")(spatial)
    spatial = layers.multiply([input, spatial])

    filters = input._keras_shape[-1]
    channels = layers.GlobalAveragePooling1D()(input)
    channels = layers.Reshape((1, filters))(channels)
    channels = layers.Conv1D(filters // rate, (1, ), trainable=trainable, kernel_initializer=initializers.ones())(channels)
    channels = layers.Activation("relu")(channels)
    channels = layers.Conv1D(filters, (1, ), trainable=trainable, kernel_initializer=initializers.ones())(channels)
    channels = layers.Activation("sigmoid")(channels)
    channels = layers.multiply([input, channels])

    res = layers.add([spatial, channels])

    return res


def ResNet50(input_shape,
             pooling=None,
             base_filter_count=64):

    f = base_filter_count

    img_input = layers.Input(shape=input_shape)

    x = img_input

    x = conv_block(x, 3, [f, f, f*4], stage=2, block='a', strides=(1, ))
    x = identity_block(x, 3, [f, f, f*4], stage=2, block='b')
    x = identity_block(x, 3, [f, f, f*4], stage=2, block='c')

    x = conv_block(x, 3, [f*2, f*2, f*8], stage=3, block='a')
    x = identity_block(x, 3, [f*2, f*2, f*8], stage=3, block='b')
    x = identity_block(x, 3, [f*2, f*2, f*8], stage=3, block='c')
    x = identity_block(x, 3, [f*2, f*2, f*8], stage=3, block='d')

    x = conv_block(x, 3, [f*4, f*4, f*16], stage=4, block='a')
    x = identity_block(x, 3, [f*4, f*4, f*16], stage=4, block='b')
    x = identity_block(x, 3, [f*4, f*4, f*16], stage=4, block='c')
    x = identity_block(x, 3, [f*4, f*4, f*16], stage=4, block='d')
    x = identity_block(x, 3, [f*4, f*4, f*16], stage=4, block='e')
    x = identity_block(x, 3, [f*4, f*4, f*16], stage=4, block='f')

    x = conv_block(x, 3, [f*8, f*8, f*32], stage=5, block='a')
    x = identity_block(x, 3, [f*8, f*8, f*32], stage=5, block='b')
    x = identity_block(x, 3, [f*8, f*8, f*32], stage=5, block='c')

    if pooling == 'avg':
        x = layers.GlobalAveragePooling1D()(x)
    elif pooling == 'flatten':
        x = layers.Flatten()(x)
    else:
        x = layers.GlobalMaxPooling1D()(x)

    model = models.Model(img_input, x, name='resnet50')

    return model


def small_resnet(input_shape,
             pooling=None,
             base_filter_count=64,
             useSE=False):

    f = base_filter_count

    img_input = layers.Input(shape=input_shape)

    x = img_input

    x = conv_block(x, 3, [f, f, f*4], stage=2, block='a', strides=(1, ))
    x = identity_block(x, 3, [f, f, f*4], stage=2, block='b')
    x = identity_block(x, 3, [f, f, f*4], stage=2, block='c')

    if useSE:
        x = SEBlock(x)

    x = conv_block(x, 3, [f*2, f*2, f*8], stage=3, block='a')
    x = identity_block(x, 3, [f*2, f*2, f*8], stage=3, block='b')
    x = identity_block(x, 3, [f*2, f*2, f*8], stage=3, block='c')

    if useSE:
        x = SEBlock(x)

    x = conv_block(x, 3, [f*4, f*4, f*16], stage=4, block='a')
    x = identity_block(x, 3, [f*4, f*4, f*16], stage=4, block='b')
    x = identity_block(x, 3, [f*4, f*4, f*16], stage=4, block='c')

    if useSE:
        x = SEBlock(x)

    x = conv_block(x, 3, [f*8, f*8, f*32], stage=5, block='a')
    x = identity_block(x, 3, [f*8, f*8, f*32], stage=5, block='b')
    x = identity_block(x, 3, [f*8, f*8, f*32], stage=5, block='c')

    if pooling == 'avg':
        x = layers.GlobalAveragePooling1D()(x)
    else:
        x = layers.GlobalMaxPooling1D()(x)

    model = models.Model(img_input, x, name='resnet50')

    return model


def getModel():
    # resnet = ResNet50((32,2), base_filter_count=8) #loss: 0.0105 - val_loss: 0.0113 (manual stop)
    # resnet = small_resnet((32,2), base_filter_count=4) #loss: 0.0100 - val_loss: 0.0111
    # resnet = small_resnet((32, 2), base_filter_count=4, pooling='avg') # loss: 0.0102 - val_loss: 0.0109 (cosine)
    # resnet = small_resnet((32, 2), base_filter_count=4, pooling='avg', useSE=True) # loss: 0.0102 - val_loss: 0.0109
    # resnet = small_resnet((32, 2), base_filter_count=4, pooling='flatten') # loss: 0.0100 - val_loss: 0.0109
    resnet = small_resnet((32, 2), base_filter_count=16, pooling='avg') # loss: 0.0100 - val_loss: 0.0109

    x = layers.Dense(8*4)(resnet.output)
    x = layers.BatchNormalization()(x)
    x = layers.Activation('relu')(x)
    output = layers.Dense(2, activation='relu')(x)

    model = models.Model(resnet.input, output)

    return model
