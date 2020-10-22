
#include <cmath>
#include <iostream>

#include <Eigen/Core>


#define CHECK(_cond) \
	do { \
		if (!(_cond)) { \
			throw std::runtime_error("Failed condition: '" #_cond "' in " __FILE__ " @ " + std::to_string(__LINE__)); \
		} \
	} while (0)



using ArrayXd = Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using ArrayXdRef = Eigen::Ref<ArrayXd>;

using Vector = Eigen::Array<double, 1, Eigen::Dynamic>;
using VectorRef = Eigen::Ref<Vector>;


void test(const Vector &v)
{
}

ArrayXd pdf(const Vector &X, int k, const Vector &mu, const Vector &sigma)
{
	// 1 / (sqrt(2*pi) * sigma) * exp(-1/2 * ((x - mu) / sigma)^2)
	ArrayXd res(X.size(), k);

	for (int j = 0; j < k; ++j)
	{
		double c1 = 1 / (sqrt(2*M_PI) * sigma(j));
		double c2 = -0.5 / (sigma(j) * sigma(j));

		double mu_j = mu(j);

		for (int i = 0, n = X.size(); i < n; ++i)
		{
			double arg = (X(i) - mu_j);
			arg = arg * arg * c2;

			res(i,j) = c1 * exp(arg);
		}
	}

	return res;
}


ArrayXd logpdf(const VectorRef X, int k, const VectorRef mu, const VectorRef sigma)
{
#if 0
	auto c1 = 1 / (-2 * sigma * sigma);
	auto c2 = Eigen::log((2 * M_PI) * sigma * sigma) / 2;

	auto err = (X - mu);
	auto err2 = err * err;

	res = err2 * c1 - c2;
#else

	ArrayXd res(X.size(), k);

	for (int j = 0; j < k; ++j) {
		double sigma_j_sqr = sigma(j) * sigma(j);
		double mu_j = mu[j];

		double c1 = 1 / (-2 * sigma_j_sqr);
		double c2 = (log(2 * M_PI * sigma_j_sqr) / 2);

		for (int i = 0, n = X.size(); i < n; ++i) {
			double err = X(i) - mu_j;
			res(i,j) = err * err * c1 - c2;
		}
	}

	return res;

#endif
}


double likelihood(const VectorRef X, int k, const VectorRef w, const VectorRef mu, const VectorRef sigma, const ArrayXdRef Q)
{
	ArrayXd p = logpdf(X, k, mu, sigma);
	return (p * Q).sum() + (Q.colwise().sum() * w.log()).sum();
}


ArrayXd posterior(const VectorRef X, int k, const VectorRef w, const VectorRef mu, const VectorRef sigma, double noise)
{
	ArrayXd Q = pdf(X, k, mu, sigma);

	Q.rowwise() *= w;
	Q.colwise() /= (Q.rowwise().sum() + noise);

	return Q;
}

struct GMM
{
	GMM() {}
	GMM(int k): w(k), mu(k), sigma(k) {}
	GMM(Vector w, Vector mu, Vector sigma): w(w), mu(mu), sigma(sigma) {}

	Vector w;
	Vector mu;
	Vector sigma;
};


GMM estimate(VectorRef X, ArrayXdRef Q)
{
	//CHECK(k == gmm.w.size());
	//CHECK(k == gmm.mu.size());
	//CHECK(k == gmm.sigma.size());
	//CHECK(Q.rows() == k);
	//CHECK(Q.cols() = X.size());

	int k = Q.cols();

	GMM gmm(k);

	for (int j = 0; j < k; ++j)
	{
		double Q_j_sum = Q.col(j).sum();
		double mu_j = (Q.col(j).transpose() * X).sum() / Q_j_sum;

		gmm.w(j) = Q_j_sum;
		gmm.mu(j) = mu_j;

		auto err = X - mu_j;
		double sigma_j = (Q.col(j).transpose() * (err * err)).sum() / Q_j_sum;

		gmm.sigma(j) = sqrt(sigma_j);
	}

	gmm.w /= gmm.w.sum();

	return gmm;
}


auto em(const VectorRef &X, int k, const VectorRef &w, const VectorRef &mu, const VectorRef &sigma, int steps, double thold, double noise)
{
	int n = X.size();

	GMM model(w, mu, sigma);
	GMM result(k);

	double Lresult = -INFINITY;
	double Lprev = -INFINITY;
	double L = 0.0;

	ArrayXd Q;

	for (int i = 0; i < steps; ++i)
	{
		// E step
		Q = posterior(X, k, w, model.mu, model.sigma, noise);

		// M step
		model = estimate(X, Q);
		L = likelihood(X, k, w, model.mu, model.sigma, Q);


		if (Lresult < L || i == 0) {
			result = model;
			Lresult = L;
		}

		if (abs((L - Lprev) / L) < thold) {
			break;
		}

		Lprev = L;
	}

	return std::make_tuple(result, Lresult);
}

