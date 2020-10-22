import os
import sys
import numpy as np

from math import log
#from sklearn.mixture import GaussianMixture

from gene import em, conf, speedup


def median(values):
    return values[int(len(values) / 2)]


def estimate1(values):
    # Overall accuracy score 1e6 * COR * AUC = 583132.125461057
    # Overall accuracy score 1e6 * COR * AUC = 667946.727840407
    # OVERALL SCORE = 393097.953488719
    #
    # With hack:
    # Overall accuracy score 1e6 * COR * AUC = 666632.626106542
    # Overall accuracy score 1e6 * COR * AUC = 788683.655631056
    # OVERALL SCORE = 533210.483424636

    values = sorted(values)

    i1 = round(len(values) * 1 / 3)
    cluster1A = values[:i1-5]
    cluster1B = values[i1+5:]

    value1A = median(cluster1A)
    value1B = median(cluster1B)
    var1 = np.var(cluster1A) + np.var(cluster1B)
    
    
    i2 = round(len(values) * 2 / 3)
    cluster2A = values[:i2-5]
    cluster2B = values[i2+5:]

    value2A = median(cluster2A)
    value2B = median(cluster2B)
    var2 = np.var(cluster2A) + np.var(cluster2B)

    if var1 < var2:
        return value1A, value1B, None
    else:
        return value2B, value2A, None


def estimate4(X, ratio1to2):
    n = len(X)
    pivot = int(0.33*n if ratio1to2 else 0.67*n)

    estA = median(X[:pivot])
    estB = median(X[pivot:])

    if not ratio1to2:
        estA, estB = estB, estA

    return estA, estB, None


#classif = GaussianMixture(n_components=2, covariance_type='diag')

def estimate2(values):
    # Overall accuracy score 1e6 * COR * AUC = 570716.740514427
    # Overall accuracy score 1e6 * COR * AUC = 653222.204752789
    # OVERALL SCORE = 376208.423342753

    #classif = GaussianMixture(n_components=2, covariance_type='diag')
    classif.fit(np.array(values).reshape(-1, 1))
    #classif.fit(values.reshape(-1, 1))
    #classif.fit(values)

    weights = classif.weights_
    means = classif.means_
    
    if weights[0] < weights[1]:
        return means[0][0], means[1][0]
    else:
        return means[1][0], means[0][0]


def medianw(X, Q):
    cs = np.cumsum(Q)
    cs /= cs[-1]

    for i, x in enumerate(cs):
        if x > 0.5:
            break

    if i == 0:
        return X[0]

    d = cs[i] - cs[i-1]
    c1 = 1 - ((0.5 - cs[i-1]) / d)
    c2 = 1 - ((cs[i] - 0.5) / d)

    return c1 * X[i-1] + c2 * X[i]


def refine(X, w, mu, sigma):
    n = len(X)
    #w = [0.33, 0.67]
    p = np.zeros((n, 2))
    Q = np.zeros((n, 2))

    speedup.pdf(p, X, 2, mu, sigma)
    speedup.posterior(Q, p, w)

    #valuesA = X[Q[:,0] > 0.5]
    #valuesB = X[Q[:,0] < 0.5]

    #print("size", len(valuesA), len(valuesB))
    estA = medianw(X, Q[:,0])
    estB = medianw(X, Q[:,1])

    if not (0 < (estA + estB) < 100000):
        estA = median(X)
        estB = estA

    return estA, estB
    #return np.median(valuesA), np.median(valuesB)


def bic(n, k, L):
    return log(n) * k - 2*L


def estimate3(X, noise=.02, single_thold=2):
    X = np.array(X)
    #X = X[X > 10]
    X = X[5:-5]

    n = len(X) - 1

    #print("n vals =", len(values))

    c1 = X[int(.05*n)]
    c2 = X[int(.95*n)]
    c3 = c2 - c1

    def norm(X, shift=c1, scale=c3):
        return (X - shift) / scale

    def denorm(X, shift=c1, scale=c3):
        return (X * scale) + shift

    X = norm(X)

    #scale = max(X)
    #X = X / scale

    #X = np.log2(X)

    #valMax = 1.0
    #valMin = min(X)

    #noise = .1

    w1, mu1, sigma1, L1 = em.em(X, 2, [X[int(0.1*n)], X[int(0.9*n)]], [.2, .2], [0.33, 0.67], 50, 0.001, noise=noise)
    w2, mu2, sigma2, L2 = em.em(X, 2, [X[int(0.9*n)], X[int(0.1*n)]], [.2, .2], [0.33, 0.67], 50, 0.001, noise=noise)
    #w3, mu3, sigma3, L3 = em.em_with_p(X, 2, [X[int(0.1*n)], X[int(0.9*n)]], [.2, .2], [0.5, 0.5], 100, 1e-4, noise=noise)
    #w4, mu4, sigma4, L4 = em.em(X, 1, [1.0], [.3], [1.0], 100, 1e0, noise=noise)

    L1 = bic(n, 6, L1)
    L2 = bic(n, 6, L2)
    #L4 = bic(n, 3, L4)

    #mu1 = 2**mu1
    #mu2 = 2**mu2

    mu1 = denorm(mu1)
    mu2 = denorm(mu2)
    #mu3 = denorm(mu3)
    #mu4 = denorm(mu4)

    sigma1 *= c3
    sigma2 *= c3
    #sigma3 *= c3
    #sigma4 *= c3


    X = denorm(X)

    if conf.DEBUG and 1:
        print()
        print("EST1: [{:.0f} {:.0f}] [{:.0f} {:.0f}] [{:.2f} {:.2f}] {:.2f}".format(mu1[0], mu1[1], sigma1[0], sigma1[1], w1[0], w1[1], L1))
        print("EST2: [{:.0f} {:.0f}] [{:.0f} {:.0f}] [{:.2f} {:.2f}] {:.2f}".format(mu2[0], mu2[1], sigma2[0], sigma2[1], w2[0], w2[1], L2))
        #print("EST3: [{:.0f} {:.0f}] [{:.0f} {:.0f}] [{:.2f} {:.2f}] {:.2f}".format(mu3[0], mu3[1], sigma3[0], sigma3[1], w3[0], w3[1], L3))
        #print("EST4: [{:.0f}] [{:.0f}] {:.2f}".format(mu4[0], sigma4[0], L4))

    if 0 and L4 < L1 and L4 < L2:
        w = [w4[0], w4[0]]
        mu = [mu4[0], mu4[0]]
        sigma = [sigma4[0], sigma4[0]]
    elif L1 < L2:
        w = w1
        mu = mu1
        sigma = sigma1
    else:
        w = w2
        mu = mu2
        sigma = sigma2


    if (abs(mu[0] - mu[1])) < single_thold * (sum(sigma)):
        mu = X[n // 2]
        mu = [mu, mu]
        sigma = [1, 1]


    #values = np.array(list(sorted(values)))
    #Q = em2.posterior(values, 2, [.33, .67], mu, sigma)

    #if abs((L1 - L2) / L1) < 1e-2:
    #    p = np.sum(Q, axis=0) / len(Q)
    #    if p[0] > p[1]:
    #        Q = em2.posterior(values, 2, [.67, .33], mu, sigma)

    #mu[0] = np.median(values[Q[:,0] > .5])
    #mu[1] = np.median(values[Q[:,1] > .5])

    #print("MED: [{:.0f} {:.0f}]".format(mu[0], mu[1]))

    #return refine(values, w, mu, sigma)
    return mu[0], mu[1], dict(w=w, mu=mu, sigma=sigma)

    #if L1 < L2 and L1 < L3:
    #    return mu1[0], mu1[1]
    #elif L2 < L1 and L2 < L3:
    #    return mu2[1], mu2[0]
    #else:
    #    return mu3[1], mu3[0]

    #return mu1[0], mu1[1]

