import pandas as pd
import numpy as np
import os
import sys
import multiprocessing
from functools import partial
from sklearn.cluster import KMeans
import time
# python main.py 1 input\\DPK.CP001_A549_24H_X1_B42  1 output 1 0 1 DPK.CP001_A549_24H_X1_B42
# python main.py input\\LITMUS.KD017_A549_96H_X1_B42 output 0 LITMUS.KD017_A549_96H_X1_B42


def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))

def compute_counter(x):
    d = {}
    el = x[0]
    d[el] = 0
    for i in x:
        if el == i:
            d[i] += 1
        else:
            el = i
            d[i] = 1
    return d

def compute_answer(input_files, path_input, gr_barcode):
    res_dict = {}
    for input_file in input_files:
        path = os.path.join(path_input, input_file)
        experiment_df  = pd.read_csv(path, sep = '\t')
        gp_experiment_df = experiment_df.groupby('barcode_id').FI.apply(lambda x: x.values).to_dict()
        tmp_dict = {}
        for bar_id, x in gr_barcode.items():
            vals = np.sort(gp_experiment_df[bar_id])
            vals = vals[vals.tolist().count(0):]

            out_k = 9
            shape_check = round(vals.shape[0] / out_k)
            check_div = vals[:-1] / vals[1:]
            check = np.where(check_div[:shape_check] < 0.8)[0]
            if check.shape[0] != 0:
                vals = vals[1 + check[-1]:]
            shape_check = round(vals.shape[0] / out_k)

            check_div = vals[:-1] / vals[1:]
            check = np.where(check_div[(out_k - 1) * shape_check:] < 0.8)[0]
            if check.shape[0] != 0:
                vals = vals[:(out_k - 1) * shape_check + check[0]]

            max_v = max(vals)
            min_v = min(vals)
            if vals.shape[0] % 2 == 0:
                div = (vals[int(vals.shape[0] / 2)] + vals[int(vals.shape[0] / 2) - 1]) / 2 / np.mean(vals)
            else:
                div = vals[int(vals.shape[0] / 2) + 1] / np.mean(vals)

            k_check = 1.2
            if div / k_check > 2 or div < 1 / k_check:
                num_bins = int(vals.shape[0] / 3)
                split = (max_v - min_v) / num_bins
                vals_tmp = vals - min_v
                vals_tmp = vals_tmp // split

                perc_25 = int(vals.shape[0] * 0.25)
                tmp = vals_tmp[:-1] - vals_tmp[1:]

                tmp = np.where(tmp[perc_25:] < -1)[0]
                if tmp.shape[0] != 0:
                    middle = tmp[0] + perc_25

                    if middle > vals.shape[0] / 2:
                        tmp_dict[x[0]] = vals[int(middle / 2)]
                        tmp_dict[x[1]] = vals[middle + int((vals.shape[0] - middle) / 2)]
                    else:
                        tmp_dict[x[1]] = vals[int(middle / 2)]
                        tmp_dict[x[0]] = vals[middle + int((vals.shape[0] - middle) / 2)]
                else:

                    km = KMeans(2, random_state=228, max_iter=1, n_init=1)
                    res = km.fit_predict(vals.reshape(-1, 1))

                    tmp = res[:-1] - res[1:]
                    middle = np.where(tmp != 0)[0][0] + 1
                    if middle > vals.shape[0] / 2:
                        tmp_dict[x[0]] = vals[int(middle / 2)]
                        tmp_dict[x[1]] = vals[middle + int((vals.shape[0] - middle) / 2)]
                    else:
                        tmp_dict[x[1]] = vals[int(middle / 2)]
                        tmp_dict[x[0]] = vals[middle + int((vals.shape[0] - middle) / 2)]

            else:
                num_bins = int(vals.shape[0] / 3)
                split = (max_v - min_v) / num_bins
                vals_tmp = vals - min_v
                vals_tmp = vals_tmp // split

                perc_15 = int(vals.shape[0] * 0.15)
                tmp = vals_tmp[:-1] - vals_tmp[1:]
                min_tmp = np.argmin(tmp[perc_15: -perc_15])

                if tmp[min_tmp + perc_15] < -2:
                    middle = min_tmp + perc_15

                    if middle > vals.shape[0] / 2:
                        tmp_dict[x[0]] = vals[int(middle / 2)]
                        tmp_dict[x[1]] = vals[middle + int((vals.shape[0] - middle) / 2)]
                    else:
                        tmp_dict[x[1]] = vals[int(middle / 2)]
                        tmp_dict[x[0]] = vals[middle + int((vals.shape[0] - middle) / 2)]

                else:

                    km = KMeans(2, random_state=228, max_iter=1, n_init=1)
                    res = km.fit_predict(vals.reshape(-1, 1))

                    tmp = res[:-1] - res[1:]
                    middle = np.where(tmp != 0)[0][0] + 1
                    if middle > vals.shape[0] / 2:
                        tmp_dict[x[0]] = vals[int(middle / 2)]
                        tmp_dict[x[1]] = vals[middle + int((vals.shape[0] - middle) / 2)]
                    else:
                        tmp_dict[x[1]] = vals[int(middle / 2)]
                        tmp_dict[x[0]] = vals[middle + int((vals.shape[0] - middle) / 2)]

        res_dict[input_file[-7:-4]] = tmp_dict
    return res_dict

def main():
    path_input = sys.argv[2]
    path_output = sys.argv[4]
    output_name = sys.argv[8] + '.gct'
    st_t = time.time()
    input_files = os.listdir(path_input)

    barcode_ids = pd.read_csv('barcode_to_gene_map.txt', sep='\t')
    barcode_ids = barcode_ids[~barcode_ids.barcode_id.isin([11, 499])]
    gr_barcode = barcode_ids.groupby('barcode_id').gene_id.apply(list).to_dict()

    new_df = pd.DataFrame(index = barcode_ids.gene_id, columns = [x[-7:-4] for x in input_files])
    num_proc = multiprocessing.cpu_count()

    subList = split(input_files, num_proc)
    with multiprocessing.Pool(processes=(num_proc)) as pool:
        results = pool.map(partial(compute_answer, path_input = path_input,
                                   gr_barcode = gr_barcode, ), subList)

    for res in results:
        for name_col in res:
            new_df[name_col] = [res[name_col][x] for x in new_df.index.values]

    med_dict = new_df.median(axis=1).to_dict()
    dict_id = dict(zip(new_df.index.values, range(new_df.shape[0])))
    new_df1 = new_df.copy()
    new_df1_vals = new_df1.values

    for k, i in gr_barcode.items():
        x0 = med_dict[i[0]]
        x1 = med_dict[i[1]]
        i0 = dict_id[i[0]]
        i1 = dict_id[i[1]]
        if x0 > x1:
            swap_check = new_df1_vals[i0] > new_df1_vals[i1]
        else:
            swap_check = new_df1_vals[i1] > new_df1_vals[i0]
        l1 = []
        l2 = []
        for y0, y1, c in zip(new_df1_vals[i0], new_df1_vals[i1], swap_check):
            if c:
                l1 += [y0]
                l2 += [y1]
            else:
                l1 += [y1]
                l2 += [y0]
        new_df1_vals[i0] = l1
        new_df1_vals[i1] = l2
    new_df1 = pd.DataFrame(data=new_df1_vals, index=barcode_ids.gene_id, columns=[x[-7:-4] for x in input_files])

    new_df1['avg'] = new_df1.mean(axis=1)
    new_df1_vals = new_df1.values
    for x in range(new_df1.shape[0]):
        aver = new_df1_vals[x][-1]
        tmp = new_df1_vals[x][:-1]
        new_df1_vals[x] = np.where(abs(tmp - aver) > aver, aver, tmp).tolist() + [aver]
    new_df1 = pd.DataFrame(data=new_df1_vals, index=barcode_ids.gene_id, columns=[x[-7:-4] for x in input_files] + ['avg'])
    new_df1 = new_df1.drop(['avg'], axis=1)

    new_df1.index.name = 'id'
    new_df1.to_csv('out.csv', sep='\t')
    with open(os.path.join(path_output, output_name), 'w') as f:
        f.write('#1.3\n')
        f.write(str(new_df.shape[0]) + '\t' + str(new_df.shape[1]) + '\t' + str(0) + '\t' + str(0) + '\n')
        with open('out.csv') as infile:
            for line in infile:
                f.write(line)

if __name__ == "__main__":
    main()