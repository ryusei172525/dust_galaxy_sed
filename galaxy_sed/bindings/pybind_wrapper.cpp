#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "sed_calculator.hpp"

namespace py = pybind11;

PYBIND11_MODULE(pybind_wrapper, m) {
    py::class_<SEDCalculator>(m, "SEDCalculator")
        .def(py::init<>())
        .def("calculate", &SEDCalculator::calculate);
}