/**
 * @file
 * @brief emodlib Python bindings.
*/

#include "pybind11/pybind11.h"

#include "emodlib/version.hpp"

#include "malaria.cpp"


namespace py = pybind11;
namespace em = emodlib;


PYBIND11_MODULE(_emodlib_py, m)
{
    m.doc() = "A collection of algorithms for disease-transmission modeling.";

    m.attr("__version__") = em::version::version_str;

    py::module malaria_m = m.def_submodule("malaria", "The malaria intra-host module of emodlib");
    add_malaria_bindings(malaria_m);
}
