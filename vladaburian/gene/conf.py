
import os

DEBUG = False

PWD = os.path.dirname(__file__)
RES = os.path.join(PWD, '../../competitor_pack_v2/')

EXPERIMENTS = (
    'DPK.CP001_A549_24H_X1_B42',
    'LITMUS.KD017_A549_96H_X1_B42'
)


CODE2GENE = {
    12: [7416, 55847],
    13: [10174, 25803],
    14: [6676, 466],
    15: [1870, 6009],
    16: [8884, 3480],
    17: [5997, 2309],
    18: [57147, 2946],
    19: [6909, 387],
    20: [9093, 3553],
    21: [23142, 427],
    22: [8444, 5898],
    23: [11232, 23365],
    24: [27095, 6657],
    25: [10732, 5054],
    26: [1398, 3108],
    27: [23029, 1950],
    28: [9270, 9517],
    29: [47, 351],
    30: [5889, 25987],
    31: [8243, 23670],
    32: [10730, 4846],
    33: [30, 1452],
    34: [3312, 4776],
    35: [891, 6908],
    36: [7466, 11182],
    37: [89910, 2146],
    38: [56654, 3385],
    39: [501, 672],
    40: [4775, 5710],
    41: [5111, 5261],
    42: [5438, 2115],
    43: [5696, 10765],
    44: [26993, 25839],
    45: [55111, 4172],
    46: [8520, 7015],
    47: [10146, 23161],
    48: [7165, 10298],
    49: [55011, 1802],
    50: [5588, 58533],
    51: [6461, 8726],
    52: [10813, 147179],
    53: [80746, 29937],
    54: [6839, 6184],
    55: [10434, 9854],
    56: [329, 2553],
    57: [5607, 5440],
    58: [22827, 2185],
    59: [51635, 54623],
    60: [84159, 9533],
    61: [5058, 8607],
    62: [998, 207],
    63: [3508, 26054],
    64: [84890, 3315],
    65: [79600, 10682],
    66: [8574, 8573],
    67: [5529, 23378],
    68: [3909, 5289],
    69: [9134, 8553],
    70: [9637, 622],
    71: [8518, 9467],
    72: [6443, 5525],
    73: [58478, 256364],
    74: [4850, 332],
    75: [6616, 2058],
    76: [23210, 7849],
    77: [54733, 868],
    78: [55129, 23244],
    79: [10273, 1677],
    80: [5223, 7043],
    81: [572, 7494],
    82: [5716, 4836],
    83: [22841, 23647],
    84: [2356, 10617],
    85: [6050, 5613],
    86: [64422, 2184],
    87: [4690, 5355],
    88: [8837, 6659],
    89: [55127, 7016],
    90: [23335, 90861],
    91: [10270, 10670],
    92: [6509, 79071],
    93: [83743, 64746],
    94: [6772, 5366],
    95: [22887, 11344],
    96: [51001, 10013],
    97: [26511, 595],
    98: [7874, 5829],
    99: [51070, 5257],
    100: [1738, 4793],
    101: [3416, 22794],
    102: [958, 9761],
    103: [9128, 2736],
    104: [9833, 23326],
    105: [1616, 10972],
    106: [4088, 10523],
    107: [5018, 5290],
    108: [51375, 79090],
    109: [55038, 27032],
    110: [23061, 9097],
    111: [991, 10513],
    112: [64428, 6194],
    113: [11325, 1891],
    114: [9641, 3066],
    115: [6347, 1676],
    116: [23038, 9801],
    117: [57048, 23512],
    118: [94239, 1459],
    119: [51203, 148022],
    120: [388650, 9703],
    121: [6714, 10606],
    122: [652, 1017],
    123: [23039, 1906],
    124: [55179, 9813],
    125: [5827, 8878],
    126: [10245, 3091],
    127: [6813, 7485],
    128: [4331, 2109],
    129: [3251, 23271],
    130: [23212, 9915],
    131: [23536, 9552],
    132: [51282, 93487],
    133: [10227, 2956],
    134: [50813, 11065],
    135: [51024, 9653],
    136: [53343, 200081],
    137: [1465, 3337],
    138: [79006, 11073],
    139: [6988, 23],
    140: [11011, 5641],
    141: [23325, 5743],
    142: [10049, 5321],
    143: [7099, 4144],
    144: [3157, 8720],
    145: [200734, 10782],
    146: [9053, 25932],
    147: [79094, 1981],
    148: [2778, 695],
    149: [54881, 10398],
    150: [54807, 5708],
    151: [10320, 1213],
    152: [323, 3028],
    153: [26001, 6804],
    154: [2042, 226],
    155: [10190, 6777],
    156: [9686, 890],
    157: [843, 6390],
    158: [54915, 2113],
    159: [23658, 4303],
    160: [7750, 10953],
    161: [79643, 22889],
    162: [840, 3098],
    163: [9488, 11044],
    164: [11200, 8985],
    165: [65057, 66008],
    166: [6832, 9217],
    167: [23131, 3329],
    168: [7048, 701],
    169: [4851, 6182],
    170: [6464, 4860],
    171: [10776, 1050],
    172: [6253, 9491],
    173: [79716, 7398],
    174: [23410, 6117],
    175: [29978, 23588],
    176: [9988, 949],
    177: [10845, 7168],
    178: [54957, 9375],
    179: [79170, 22926],
    180: [10610, 4927],
    181: [128, 2523],
    182: [51026, 58497],
    183: [5498, 6856],
    184: [124583, 26036],
    185: [9455, 5899],
    186: [1212, 10206],
    187: [5255, 25793],
    188: [54205, 5467],
    189: [9448, 55556],
    190: [9924, 6774],
    191: [5480, 3978],
    192: [1983, 10681],
    193: [2770, 3638],
    194: [56889, 1788],
    195: [79080, 8349],
    196: [9016, 6275],
    197: [6915, 3895],
    198: [80204, 392],
    199: [9709, 9170],
    200: [116832, 80347],
    201: [670, 5921],
    202: [7867, 23300],
    203: [79143, 836],
    204: [11284, 5721],
    205: [23443, 11188],
    206: [51599, 10898],
    207: [4208, 9276],
    208: [5927, 1647],
    209: [55793, 10285],
    210: [5867, 23386],
    211: [7866, 2263],
    212: [25825, 3156],
    213: [4067, 10237],
    214: [51015, 7077],
    215: [1845, 11157],
    216: [5891, 808],
    217: [23368, 22908],
    218: [3300, 10915],
    219: [64781, 9897],
    220: [51056, 1759],
    221: [8550, 54541],
    222: [50810, 211],
    223: [5583, 5048],
    224: [23011, 84722],
    225: [11168, 30001],
    226: [1385, 348],
    227: [10150, 29890],
    228: [29763, 9868],
    229: [9112, 873],
    230: [23499, 5580],
    231: [8900, 30849],
    232: [85236, 1454],
    233: [5096, 3725],
    234: [55620, 23585],
    235: [5373, 51097],
    236: [51170, 1831],
    237: [54850, 7020],
    238: [80349, 55893],
    239: [355, 6944],
    240: [1062, 22883],
    241: [8804, 5601],
    242: [80758, 6500],
    243: [375346, 10775],
    244: [51385, 4780],
    245: [581, 10904],
    246: [55748, 9917],
    247: [10318, 5627],
    248: [54681, 1399],
    249: [7905, 24149],
    250: [28969, 2222],
    251: [10668, 10131],
    252: [51005, 2958],
    253: [8800, 5715],
    254: [54438, 9650],
    255: [55825, 7319],
    256: [8835, 4817],
    257: [64429, 9143],
    258: [3682, 4791],
    259: [10329, 10494],
    260: [27109, 55324],
    261: [51335, 23149],
    262: [3202, 9943],
    263: [2896, 6599],
    264: [79902, 56924],
    265: [29083, 2690],
    266: [5747, 51031],
    267: [22905, 1019],
    268: [63933, 3566],
    269: [9847, 51422],
    270: [51021, 3815],
    271: [6193, 5788],
    272: [8996, 8624],
    273: [7376, 27346],
    274: [5357, 10489],
    275: [56940, 26128],
    276: [6195, 79921],
    277: [1633, 55746],
    278: [23321, 23636],
    279: [6622, 6709],
    280: [51160, 10892],
    281: [55604, 10652],
    282: [11261, 11004],
    283: [2063, 23530],
    284: [54512, 1277],
    285: [9124, 993],
    286: [79174, 7157],
    287: [23223, 8826],
    288: [39, 1958],
    289: [57215, 4582],
    290: [1052, 9212],
    291: [23139, 3308],
    292: [79947, 7153],
    293: [4616, 8440],
    294: [9710, 6790],
    295: [6990, 11041],
    296: [29082, 8851],
    297: [9670, 6342],
    298: [81533, 5883],
    299: [6696, 10589],
    300: [10525, 5110],
    301: [2353, 2534],
    302: [57192, 4609],
    303: [3398, 5423],
    304: [3775, 1070],
    305: [8312, 10051],
    306: [26064, 6603],
    307: [29103, 9019],
    308: [3206, 4893],
    309: [835, 10969],
    310: [1500, 9797],
    311: [7690, 4043],
    312: [29911, 665],
    313: [2064, 9519],
    314: [65123, 8480],
    315: [5986, 983],
    316: [85377, 10559],
    317: [9903, 8508],
    318: [7088, 9275],
    319: [8324, 8914],
    320: [22796, 54505],
    321: [51116, 4313],
    322: [54442, 6919],
    323: [10695, 2274],
    324: [102, 2288],
    325: [965, 7994],
    326: [3800, 1022],
    327: [10362, 4216],
    328: [5300, 7538],
    329: [8270, 3628],
    330: [27242, 5880],
    331: [9928, 3162],
    332: [1848, 1282],
    333: [1534, 2037],
    334: [5654, 2771],
    335: [10810, 7264],
    336: [7105, 1786],
    337: [23014, 5873],
    338: [55012, 780],
    339: [154, 9926],
    340: [2887, 8446],
    341: [4312, 7159],
    342: [9183, 1026],
    343: [11151, 7082],
    344: [7511, 10058],
    345: [10046, 10797],
    346: [1829, 11230],
    347: [11000, 30836],
    348: [976, 27336],
    349: [3597, 9918],
    350: [178, 79073],
    351: [994, 3383],
    352: [2625, 9688],
    353: [6119, 23659],
    354: [26136, 8396],
    355: [6697, 23077],
    356: [5971, 10857],
    357: [2048, 16],
    358: [22934, 2817],
    359: [9289, 60528],
    360: [4931, 5359],
    361: [1111, 1514],
    362: [4200, 1509],
    363: [5154, 9181],
    364: [93594, 3611],
    365: [5982, 4783],
    366: [26520, 10123],
    367: [4998, 8895],
    368: [10038, 727],
    369: [60493, 9805],
    370: [51495, 4125],
    371: [4482, 4792],
    372: [23224, 1635],
    373: [9690, 4154],
    374: [80212, 26292],
    375: [2920, 5770],
    376: [29928, 25874],
    377: [22809, 11319],
    378: [50814, 1153],
    379: [55256, 4794],
    380: [56997, 11098],
    381: [81544, 902],
    382: [6894, 8974],
    383: [64080, 533],
    384: [55818, 823],
    385: [637, 6812],
    386: [22823, 3693],
    387: [3122, 4891],
    388: [5287, 10954],
    389: [91137, 5909],
    390: [23097, 9221],
    391: [9851, 847],
    392: [23522, 11072],
    393: [1123, 8091],
    394: [5603, 3486],
    395: [51569, 10221],
    396: [2769, 7158],
    397: [55699, 10818],
    398: [2624, 8678],
    399: [5900, 11014],
    400: [8204, 1994],
    401: [9133, 3454],
    402: [10921, 4282],
    403: [6304, 1605],
    404: [51742, 10644],
    405: [2954, 3482],
    406: [27244, 25976],
    407: [23013, 10491],
    408: [5211, 51053],
    409: [8050, 9267],
    410: [79961, 5547],
    411: [2017, 6597],
    412: [55148, 5836],
    413: [2886, 6118],
    414: [960, 2065],
    415: [10451, 2961],
    416: [57406, 310],
    417: [9261, 7296],
    418: [51719, 8321],
    419: [7074, 1429],
    420: [10059, 2548],
    421: [51071, 6499],
    422: [5796, 642],
    423: [55008, 23076],
    424: [8569, 10099],
    425: [23597, 8503],
    426: [8869, 9817],
    427: [8727, 7027],
    428: [63874, 9126],
    429: [8821, 3930],
    430: [57178, 10450],
    431: [5468, 54386],
    432: [10557, 664],
    433: [831, 596],
    434: [9531, 10180],
    435: [2767, 1445],
    436: [1861, 230],
    437: [10165, 9246],
    438: [7106, 6850],
    439: [2908, 5092],
    440: [57761, 23463],
    441: [23338, 1978],
    442: [7852, 10007],
    443: [55033, 896],
    444: [57019, 4925],
    445: [91949, 6251],
    446: [84617, 5777],
    447: [25, 1846],
    448: [899, 5782],
    449: [3303, 5106],
    450: [79850, 5925],
    451: [51466, 2673],
    452: [2852, 644],
    453: [55837, 6284],
    454: [2131, 4605],
    455: [5427, 6810],
    456: [5331, 54499],
    457: [25966, 142],
    458: [55608, 11007],
    459: [1662, 1643],
    460: [9702, 25805],
    461: [11031, 3033],
    462: [10057, 2745],
    463: [5985, 9961],
    464: [3280, 3925],
    465: [55958, 2582],
    466: [51465, 5792],
    467: [2542, 5236],
    468: [5347, 4016],
    469: [10641, 4651],
    470: [11137, 3988],
    471: [9842, 26227],
    472: [10493, 7982],
    473: [10112, 3551],
    474: [1001, 5108],
    475: [23635, 8870],
    476: [29916, 2195],
    477: [2810, 8731],
    478: [5993, 1956],
    479: [11142, 10276],
    480: [79187, 1666],
    481: [51382, 1021],
    482: [8318, 4864],
    483: [9695, 4232],
    484: [4638, 5050],
    485: [9712, 813],
    486: [64943, 3964],
    487: [9738, 26020],
    488: [6793, 481],
    489: [7358, 5019],
    490: [58472, 291],
    491: [50865, 5699],
    492: [2264, 10973],
    493: [1029, 9697],
    494: [23200, 5831],
    495: [51293, 8202],
    496: [1027, 57804],
    497: [10962, 8061],
    498: [10153, 5566],
    500: [874, 57149]
}

IS_A_LESS = {12: False, 13: False, 14: True, 15: True, 16: True, 17: False, 18: False, 19: True, 20: False, 21: True, 22: True, 23: False, 24: False, 25: True, 26: False, 27: False, 28: False, 29: True, 30: True, 31: False, 32: False, 33: True, 34: False, 35: False, 36: True, 37: True, 38: False, 39: False, 40: True, 41: False, 42: False, 43: True, 44: True, 45: True, 46: False, 47: False, 48: True, 49: True, 50: True, 51: True, 52: False, 53: True, 54: True, 55: False, 56: False, 57: True, 58: False, 59: False, 60: True, 61: True, 62: False, 63: True, 64: True, 65: True, 66: False, 67: False, 68: True, 69: True, 70: False, 71: False, 72: True, 73: False, 74: True, 75: True, 76: False, 77: False, 78: True, 79: False, 80: False, 81: True, 82: False, 83: True, 84: True, 85: False, 86: False, 87: False, 88: False, 89: False, 90: True, 91: True, 92: True, 93: True, 94: False, 95: False, 96: False, 97: True, 98: True, 99: True, 100: False, 101: True, 102: True, 103: False, 104: True, 105: True, 106: False, 107: False, 108: False, 109: False, 110: True, 111: False, 112: True, 113: True, 114: True, 115: False, 116: True, 117: True, 118: False, 119: False, 120: False, 121: True, 122: True, 123: True, 124: True, 125: True, 126: True, 127: True, 128: True, 129: False, 130: False, 131: True, 132: True, 133: True, 134: True, 135: False, 136: True, 137: True, 138: True, 139: False, 140: True, 141: False, 142: False, 143: True, 144: False, 145: True, 146: True, 147: True, 148: False, 149: False, 150: True, 151: True, 152: True, 153: False, 154: True, 155: False, 156: True, 157: True, 158: False, 159: False, 160: True, 161: True, 162: False, 163: True, 164: True, 165: False, 166: False, 167: True, 168: True, 169: True, 170: True, 171: False, 172: True, 173: True, 174: True, 175: True, 176: False, 177: True, 178: True, 179: True, 180: True, 181: False, 182: False, 183: True, 184: False, 185: True, 186: True, 187: True, 188: False, 189: False, 190: False, 191: False, 192: False, 193: True, 194: False, 195: False, 196: True, 197: True, 198: True, 199: False, 200: False, 201: True, 202: True, 203: False, 204: True, 205: True, 206: True, 207: True, 208: True, 209: True, 210: False, 211: False, 212: True, 213: True, 214: False, 215: False, 216: True, 217: True, 218: True, 219: True, 220: False, 221: True, 222: True, 223: True, 224: True, 225: False, 226: False, 227: False, 228: True, 229: True, 230: False, 231: True, 232: False, 233: False, 234: True, 235: True, 236: False, 237: True, 238: False, 239: False, 240: False, 241: True, 242: True, 243: True, 244: True, 245: True, 246: False, 247: True, 248: True, 249: False, 250: True, 251: True, 252: True, 253: True, 254: True, 255: True, 256: True, 257: False, 258: False, 259: False, 260: False, 261: False, 262: True, 263: True, 264: False, 265: False, 266: False, 267: True, 268: False, 269: False, 270: False, 271: False, 272: True, 273: True, 274: False, 275: True, 276: True, 277: False, 278: True, 279: False, 280: False, 281: True, 282: True, 283: False, 284: False, 285: False, 286: True, 287: True, 288: False, 289: False, 290: True, 291: True, 292: True, 293: True, 294: True, 295: False, 296: False, 297: True, 298: False, 299: False, 300: True, 301: False, 302: True, 303: False, 304: True, 305: True, 306: False, 307: True, 308: True, 309: True, 310: True, 311: False, 312: True, 313: True, 314: True, 315: True, 316: True, 317: True, 318: False, 319: True, 320: True, 321: False, 322: True, 323: True, 324: True, 325: False, 326: True, 327: False, 328: False, 329: False, 330: False, 331: False, 332: False, 333: True, 334: False, 335: True, 336: True, 337: False, 338: False, 339: True, 340: True, 341: True, 342: True, 343: True, 344: True, 345: True, 346: False, 347: True, 348: True, 349: False, 350: True, 351: False, 352: True, 353: False, 354: False, 355: False, 356: True, 357: True, 358: False, 359: True, 360: True, 361: True, 362: False, 363: True, 364: True, 365: True, 366: False, 367: True, 368: False, 369: False, 370: False, 371: True, 372: True, 373: True, 374: True, 375: True, 376: True, 377: True, 378: False, 379: False, 380: True, 381: True, 382: False, 383: False, 384: False, 385: True, 386: True, 387: True, 388: False, 389: True, 390: True, 391: True, 392: True, 393: True, 394: True, 395: False, 396: True, 397: False, 398: True, 399: True, 400: False, 401: False, 402: True, 403: True, 404: False, 405: True, 406: True, 407: True, 408: True, 409: False, 410: True, 411: False, 412: True, 413: True, 414: False, 415: True, 416: True, 417: True, 418: False, 419: True, 420: False, 421: False, 422: True, 423: True, 424: True, 425: False, 426: False, 427: False, 428: True, 429: True, 430: False, 431: True, 432: True, 433: False, 434: False, 435: False, 436: False, 437: False, 438: False, 439: True, 440: False, 441: True, 442: True, 443: True, 444: False, 445: False, 446: False, 447: True, 448: True, 449: False, 450: True, 451: True, 452: True, 453: True, 454: False, 455: False, 456: True, 457: True, 458: False, 459: False, 460: True, 461: True, 462: True, 463: False, 464: True, 465: True, 466: True, 467: True, 468: False, 469: True, 470: False, 471: True, 472: False, 473: False, 474: True, 475: True, 476: True, 477: True, 478: True, 479: False, 480: True, 481: False, 482: True, 483: True, 484: True, 485: True, 486: True, 487: True, 488: True, 489: False, 490: False, 491: True, 492: True, 493: True, 494: True, 495: True, 496: False, 497: False, 498: False, 500: True}
IS_A_LESS = {12: False, 13: False, 14: True, 15: True, 16: True, 17: False, 18: False, 19: True, 20: False, 21: True, 22: True, 23: False, 24: False, 25: True, 26: False, 27: False, 28: False, 29: True, 30: True, 31: False, 32: False, 33: True, 34: False, 35: False, 36: True, 37: True, 38: False, 39: False, 40: True, 41: False, 42: False, 43: True, 44: True, 45: True, 46: False, 47: False, 48: True, 49: True, 50: True, 51: True, 52: False, 53: True, 54: True, 55: False, 56: False, 57: True, 58: False, 59: False, 60: True, 61: True, 62: False, 63: True, 64: True, 65: True, 66: False, 67: False, 68: True, 69: True, 70: False, 71: False, 72: True, 73: False, 74: True, 75: True, 76: False, 77: False, 78: True, 79: False, 80: False, 81: True, 82: False, 83: True, 84: True, 85: False, 86: False, 87: False, 88: False, 89: False, 90: True, 91: True, 92: True, 93: True, 94: False, 95: False, 96: False, 97: False, 98: True, 99: True, 100: False, 101: True, 102: True, 103: False, 104: True, 105: True, 106: False, 107: False, 108: False, 109: False, 110: True, 111: False, 112: True, 113: True, 114: True, 115: False, 116: True, 117: True, 118: False, 119: False, 120: False, 121: True, 122: True, 123: False, 124: True, 125: True, 126: True, 127: True, 128: True, 129: False, 130: True, 131: True, 132: True, 133: True, 134: True, 135: False, 136: True, 137: False, 138: True, 139: False, 140: True, 141: False, 142: False, 143: True, 144: False, 145: True, 146: True, 147: False, 148: False, 149: False, 150: True, 151: True, 152: True, 153: False, 154: True, 155: False, 156: True, 157: True, 158: False, 159: False, 160: False, 161: True, 162: False, 163: True, 164: True, 165: True, 166: False, 167: True, 168: True, 169: True, 170: True, 171: False, 172: True, 173: True, 174: True, 175: True, 176: False, 177: True, 178: True, 179: True, 180: True, 181: False, 182: False, 183: True, 184: False, 185: True, 186: True, 187: True, 188: False, 189: False, 190: False, 191: False, 192: False, 193: True, 194: False, 195: False, 196: True, 197: True, 198: True, 199: False, 200: False, 201: True, 202: True, 203: False, 204: True, 205: True, 206: True, 207: True, 208: True, 209: True, 210: False, 211: False, 212: True, 213: True, 214: False, 215: False, 216: True, 217: True, 218: True, 219: True, 220: False, 221: True, 222: True, 223: True, 224: False, 225: False, 226: False, 227: False, 228: True, 229: True, 230: False, 231: True, 232: False, 233: False, 234: True, 235: True, 236: False, 237: True, 238: False, 239: False, 240: False, 241: False, 242: True, 243: True, 244: True, 245: True, 246: False, 247: False, 248: True, 249: False, 250: True, 251: True, 252: True, 253: True, 254: True, 255: True, 256: True, 257: False, 258: False, 259: False, 260: False, 261: False, 262: True, 263: False, 264: False, 265: False, 266: False, 267: True, 268: False, 269: False, 270: False, 271: False, 272: True, 273: True, 274: False, 275: True, 276: True, 277: False, 278: True, 279: False, 280: False, 281: True, 282: True, 283: False, 284: False, 285: False, 286: True, 287: True, 288: False, 289: False, 290: True, 291: True, 292: True, 293: True, 294: True, 295: False, 296: False, 297: True, 298: False, 299: False, 300: True, 301: False, 302: True, 303: False, 304: True, 305: True, 306: False, 307: True, 308: True, 309: True, 310: False, 311: False, 312: True, 313: True, 314: True, 315: True, 316: True, 317: True, 318: False, 319: True, 320: True, 321: False, 322: True, 323: True, 324: True, 325: False, 326: True, 327: False, 328: False, 329: False, 330: False, 331: False, 332: True, 333: True, 334: False, 335: True, 336: True, 337: False, 338: False, 339: True, 340: True, 341: True, 342: True, 343: True, 344: False, 345: True, 346: False, 347: True, 348: True, 349: False, 350: True, 351: False, 352: True, 353: False, 354: False, 355: False, 356: True, 357: True, 358: False, 359: True, 360: True, 361: True, 362: False, 363: False, 364: True, 365: True, 366: False, 367: True, 368: False, 369: False, 370: False, 371: True, 372: True, 373: True, 374: True, 375: True, 376: True, 377: True, 378: False, 379: False, 380: True, 381: True, 382: False, 383: False, 384: False, 385: True, 386: True, 387: True, 388: False, 389: False, 390: True, 391: True, 392: True, 393: True, 394: True, 395: False, 396: True, 397: False, 398: True, 399: True, 400: False, 401: False, 402: True, 403: True, 404: False, 405: True, 406: True, 407: True, 408: True, 409: False, 410: True, 411: False, 412: True, 413: True, 414: False, 415: True, 416: True, 417: True, 418: False, 419: True, 420: False, 421: False, 422: True, 423: True, 424: True, 425: False, 426: False, 427: False, 428: True, 429: True, 430: False, 431: True, 432: True, 433: False, 434: False, 435: False, 436: False, 437: False, 438: False, 439: True, 440: False, 441: True, 442: True, 443: True, 444: False, 445: False, 446: False, 447: True, 448: True, 449: False, 450: True, 451: True, 452: True, 453: True, 454: False, 455: False, 456: True, 457: True, 458: False, 459: False, 460: True, 461: True, 462: True, 463: False, 464: True, 465: True, 466: True, 467: True, 468: False, 469: True, 470: False, 471: True, 472: False, 473: False, 474: True, 475: True, 476: True, 477: True, 478: False, 479: False, 480: True, 481: False, 482: True, 483: True, 484: True, 485: True, 486: True, 487: True, 488: True, 489: False, 490: False, 491: False, 492: True, 493: True, 494: True, 495: True, 496: False, 497: False, 498: False, 500: True}
