import numpy as np
from tqdm import tqdm


def generate_fake_normal_mixture(count, hist_count=32):
    resX = []
    resY = []
    resV = []
    for i in tqdm(range(count), desc="fake data"):
        elements_total = int(np.random.uniform(50, 100))
        elements_low = int(np.random.uniform(0.3, 0.36)*elements_total)
        elements_high = elements_total - elements_low

        mu1 = np.random.uniform(0, 200)
        sigma1 = np.exp(np.random.uniform(np.log(0.5), np.log(100)))
        mu2 = np.random.uniform(0, 200)
        sigma2 = np.exp(np.random.uniform(np.log(0.5), np.log(100)))

        v = np.concatenate([
            np.random.normal(mu1, sigma1, elements_low),
            np.random.normal(mu2, sigma2, elements_high),
        ])

        min_v = np.min(v)
        max_v = np.max(v)

        resY.append([(mu1-min_v)/(max_v-min_v), (mu2-min_v)/(max_v-min_v)])

        h, _ = np.histogram(v, bins=hist_count)
        h = np.array(h) / np.max(h)

        resX.append(h)
        resV.append((v-min_v)/(max_v-min_v))

    resX = np.vstack(resX)
    resY = np.vstack(resY)

    return resX, resY, resV
