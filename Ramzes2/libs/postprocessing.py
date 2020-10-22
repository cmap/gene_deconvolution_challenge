import numpy as np
import pandas as pd

def getResultTable(data, predict):
    data['res_geneLowPos'] = predict[:, 0] * (data['max_value'] - data['min_value']) + data['min_value']
    data['res_geneHighPos'] = predict[:, 1] * (data['max_value'] - data['min_value']) + data['min_value']

    resLowPart = data.pivot(index='gene_low', columns='col', values='res_geneLowPos')
    resLowPart.insert(0, 'id', resLowPart.index.values)

    resHighPart = data.pivot(index='gene_high', columns='col', values='res_geneHighPos')
    resHighPart.insert(0, 'id', resHighPart.index.values)

    res = pd.concat([resLowPart, resHighPart])

    return res


def saveToFile(dir, testName, data):
    if dir[-1] != '/':
        dir += '/'

    outFileName = dir + testName + '.gct'

    out = open(outFileName, 'w')
    out.write('#1.3\n')
    out.write(str(data.shape[0]) + '\t' + str(data.shape[1] - 1) + '\t0\t0\n')

    data.to_csv(out, sep='\t', index=False)
    out.close()
