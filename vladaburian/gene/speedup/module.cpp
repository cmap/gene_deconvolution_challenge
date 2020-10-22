
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include "em.cpp"
#include "io.cpp"


PYBIND11_MODULE(speedup, m)
{
	m.def("load_lbx", &load_lbx);
	m.def("save_gct", &save_gct);

    m.def("pdf", &pdf);
    m.def("logpdf", &logpdf);
    m.def("likelihood", &likelihood);
    m.def("posterior", &posterior);

    m.def("em", [](const VectorRef &X, int k, const VectorRef &w, const VectorRef &mu, const VectorRef &sigma, int steps, double thold, double noise)
    {
    	auto [gmm, L] = em(X, k, w, mu, sigma, steps, thold, noise);
    	return std::make_tuple(gmm.w, gmm.mu, gmm.sigma, L);
    });

}
