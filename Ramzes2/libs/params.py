HIST_BINS = 32

#normalized by sum
# CUT_THRESHOLD = 0.0 # ? (0.8502)
# CUT_THRESHOLD = 0.001 # ? (0.6178)
# CUT_THRESHOLD = 0.005 # ? (0.6548)
CUT_THRESHOLD = 0.0005 # ? (0.6166, 0.6226, 0.6087 => ) (weighted 1:0.01 - 0.6077, 06030 => )
# CUT_THRESHOLD = 0.003 # ? (0.6502)
# CUT_THRESHOLD = 0.0008 # ? (0.6199)
# CUT_THRESHOLD = 0.0003 # ? (0.6149)

#normalized by max
# CUT_THRESHOLD = 0.005 # 0.60644 => 26.60046 lb  (0.6257 as metric only => 28.9 lb)

# CUT_THRESHOLD_SIDE = 0.025

#no global - 0.65416

USE_CPP_PARSER = True