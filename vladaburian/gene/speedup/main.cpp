#include "io.cpp"
#include "em.cpp"
#include "code2gene.cpp"


std::pair<double, double> detectPeaks(std::vector<double> &values, double noise=0.02)
{
	//std::cout << values.size() << "\n";
	//Vector X = Eigen::Map<Vector>(values.data(), values.size());

	Vector X(values.size());

	for (int i = 0, n = values.size(); i < n; ++i) {
		X[i] = values[i];
	}

	//return {0, 0};

	//X = X[5:-5]

	int n = X.size() - 1;

	double c1 = X[int(0.05*n)];
	double c2 = X[int(0.95*n)];
	double c3 = c2 - c1;

	//std::cout << c1 << "  " << c2 << "\n";


	auto norm = [=](auto x) {
		return (x - c1) / c3;
	};

	auto denorm = [=](auto x) {
		return (x * c3) + c1;
	};

	X = norm(X);


	GMM gmm1, gmm2;
	double L1, L2;

	{
		Vector mu(2); mu << X[int(0.1*n)], X[int(0.9*n)];
		Vector sigma(2); sigma << 0.2, 0.2;
		Vector w(2); w << 0.33, 0.67;

		std::tie(gmm1, L1) = em(X, 2, mu, sigma, w, 50, 0.001, noise);
	}

	{
		Vector mu(2); mu << X[int(0.9*n)], X[int(0.1*n)];
		Vector sigma(2); sigma << 0.2, 0.2;
		Vector w(2); w << 0.33, 0.67;

		std::tie(gmm2, L2) = em(X, 2, mu, sigma, w, 50, 0.001, noise);
	}

	GMM gmm;

	//mu1 = denorm(mu1)
	//mu2 = denorm(mu2)

	//sigma1 *= c3
	//sigma2 *= c3

	//X = denorm(X);

	//if conf.DEBUG and 1:
	//    print()
	//    print("EST1: [{:.0f} {:.0f}] [{:.0f} {:.0f}] [{:.2f} {:.2f}] {:.2f}".format(mu1[0], mu1[1], sigma1[0], sigma1[1], w1[0], w1[1], L1))
	//    print("EST2: [{:.0f} {:.0f}] [{:.0f} {:.0f}] [{:.2f} {:.2f}] {:.2f}".format(mu2[0], mu2[1], sigma2[0], sigma2[1], w2[0], w2[1], L2))

	//std::cout << L1 << "..." << L2 << "\n";

	if (L1 > L2) {
		gmm = gmm1;
	} else {
		gmm = gmm2;
	}

	gmm.mu = denorm(gmm.mu);
	//    w = w1
	//    mu = mu1
	//    sigma = sigma1
	//else:
	//    w = w2
	//    mu = mu2
	//    sigma = sigma2

	//return mu[0], mu[1], dict(w=w, mu=mu, sigma=sigma)

	//return {0, 0};
	return {gmm.mu[0], gmm.mu[1]};
}


int main()
{
	std::string root = "/home/vladimir/workspaces/workspace-topcoder/cmap/competitor_pack_v2/input/";
	std::string plate = "DPK.CP001_A549_24H_X1_B42";

	auto dataset = load_lbx(root + plate, plate);

	std::map<int, std::map<std::string, double>> result;
	std::mutex mtx;

	tbb::task_group tg;

	for (const auto &kvp: dataset) {
		int barcode = kvp.first;

		tg.run([&,barcode]() {
			int geneA, geneB;

			{
				auto it = CODE2GENE.find(barcode);

				if (it == CODE2GENE.end()) {
					return;
					//continue;
				}

				std::tie(geneA, geneB) = it->second;
			}

			//std::cout << "geneA = " << geneA << "\ngeneB = " << geneB << "\n";

			std::map<std::string, double> resultA;
			std::map<std::string, double> resultB;

			for (auto &kvp: dataset.at(barcode)) {
				const std::string &well = kvp.first;

				double A, B;

				std::tie(A, B) = detectPeaks(kvp.second);

				resultA[well] = A;
				resultB[well] = B;
			}

			std::unique_lock lock(mtx);
			result[geneA] = std::move(resultA);
			result[geneB] = std::move(resultB);
		});
	}

	tg.wait();

	save_gct(root + plate + ".gct", result);

	return 0;
}
