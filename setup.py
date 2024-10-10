from setuptools import setup, find_packages
from pybind11.setup_helpers import Pybind11Extension, build_ext

ext_modules = [
    Pybind11Extension(
        "galaxy_sed.bindings.pybind_wrapper",
        ["galaxy_sed/bindings/pybind_wrapper.cpp"],
        include_dirs=["galaxy_sed/cpp"],
        extra_compile_args=["-std=c++11"]
    ),
]

setup(
    name="galaxy_sed",
    version="0.1",
    author="Ryusei Kano",
    packages=find_packages(),
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    install_requires=["numpy", "scipy", "pybind11"]
)
