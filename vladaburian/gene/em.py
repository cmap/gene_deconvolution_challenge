
import numpy as np
#from scipy.stats import norm
#from math import log, sqrt

from gene import speedup


def posterior(X, k, w, mu, sigma, noise=.0):
    Q = np.zeros((len(X), k))

    for i in range(k):
        Q[:, i] = w[i] * norm.pdf(X, mu[i], sigma[i])

    Q /= (np.sum(Q, axis=1) + noise)[:, np.newaxis]

    return Q


def posterior2(p, w, noise):
    Q = p * w
    Q /= (Q.sum(axis=1) + noise)[:, np.newaxis]
    return Q


def likelihood(X, k, p, mu, sigma, Q):
    #L = np.zeros((len(X), k))

    L = speedup.logpdf(X, k, mu, sigma)
    #for i in range(k):
    #    speedup.logpdf(L[:,i], X, mu[i], sigma[i])
    #    #[:, i] = norm.logpdf(X, mu[i], sigma[i])

    L *= Q
    L = np.sum(L)

    L += np.sum(np.sum(Q, axis=0) * np.log(p))

    return L


def em(X, k, mu, sigma, w, steps, thold=1e-3, noise=.0):
    w = np.array(w)
    mu = np.array(mu)
    sigma = np.array(sigma)

    w, mu, sigma, L = speedup.em(X, k, w, mu, sigma, steps, thold, noise)
    return w, mu, sigma, L


    Lprev = 0

    Lresult = None
    result = None

    n = len(X)

    p = np.zeros((n, k))
    Q = np.zeros((n, k))

    for i in range(steps):
        # E step
        #speedup.pdf(p, X, k, mu, sigma)
        #p[:,0] = norm.pdf(X, mu[0], sigma[0])
        #p[:,1] = norm.pdf(X, mu[1], sigma[1])
        #Q = posterior2(p, w, noise)
        #speedup.posterior(Q, p, w, noise)

        Q = speedup.posterior(X, k, w, mu, sigma, noise)

        #Q[p[:, 0] < (0.2 / sigma[0]), 0] = 0.0
        #Q[p[:, 1] < (0.2 / sigma[1]), 1] = 0.0

        # M step
        mu = np.sum(Q*X[:,np.newaxis], axis=0) / np.sum(Q,axis=0)
        sigma = np.sum(Q * np.power(X[:, np.newaxis] - mu, 2), axis=0) / np.sum(Q, axis=0)
        sigma = np.power(sigma, 0.5)

        #speedup.logpdf(p, X, k, mu, sigma)
        #L = speedup.likelihood(p, w, Q)

        L = speedup.likelihood(X, k, w, mu, sigma, Q)
        #print(i, L)

        if (Lresult is None) or (Lresult < L):
            Lresult = L
            result = [mu, sigma, L]

        #result = [mu, sigma, L]

        if abs((L - Lprev) / L) < thold:
            break

        Lprev = L

    w = np.sum(Q, axis=0) / len(Q)
    #result[2] = log(len(X)) * (2*k) - 2 * result[2]

    return (w, *result)


def em_with_p(X, k, mu, sigma, w, steps, thold=1e-3, noise=.05):
    Lprev = 0

    for i in range(steps):
        # E step
        Q = posterior(X, k, w, mu, sigma, noise)

        # M step
        w = np.sum(Q, axis=0) / len(Q)
        mu = np.sum(Q*X[:,np.newaxis], axis=0) / np.sum(Q,axis=0)
        sigma = np.sum(Q * np.power(X[:, np.newaxis] - mu, 2), axis=0) / np.sum(Q, axis=0)
        sigma = np.power(sigma, 0.5)

        L = likelihood(X, k, w, mu, sigma, Q)
        #print(i, L)

        if abs((L - Lprev) / L) < thold:
            break

        Lprev = L

    return w, mu, sigma, L



if __name__ == '__main__':
    import matplotlib.pyplot as plt

    pi = np.array([0.33, 0.67])
    NSamples = 1000.0

    mus = [0.0, 3]
    sigmas = [1.0, 1.0]

    x1 = np.random.normal(mus[0], sigmas[0], int(pi[0] * NSamples))
    x2 = np.random.normal(mus[1], sigmas[1], int(pi[1] * NSamples))

    X = np.hstack([x1, x2])

    p, estMu, estSigma, L = em(X, 2, [-1, 10], [1, 1], [0.33, 0.67], 100)
    print(p, estMu, estSigma, L)

    p, estMu, estSigma, L = em(X, 2, [10, -1], [1, 1], [0.33, 0.67], 100)
    print(p, estMu, estSigma, L)

    plt.hist(X, 50, density=True)
    plt.show()
