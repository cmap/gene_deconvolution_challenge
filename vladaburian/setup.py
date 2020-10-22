
from setuptools import setup, Extension
from distutils.core import setup


ext_modules = [
    Extension(
        'gene.speedup',
        #['gene/speedup/em.cpp', 'gene/speedup/io.cpp'],
        ['gene/speedup/module.cpp'],
        include_dirs=[
            "/usr/include/eigen3/",
            "/tmp/pybind11-2.2.4/include/",
            # Path to pybind11 headers
            #get_pybind_include(),
            #get_pybind_include(user=True)
        ],
        extra_compile_args=['--std=c++17', '-flto'],
        extra_link_args=['-lstdc++fs', '-ltbb'],
        language='c++'
    ),
]


setup(
    name='gene',
    packages=['gene'],
    install_requires=['numpy', 'scipy'],
    ext_modules=ext_modules,
)
