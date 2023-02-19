/**
 * @file
 * @brief emodlib Python bindings.
*/

#include "pybind11/pybind11.h"

#include "malaria.cpp"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

PYBIND11_MODULE(_emodlib_py, m)
{
    m.doc() = "A collection of algorithms for disease-transmission modeling.";

    py::module malaria_m = m.def_submodule("malaria", "The malaria intra-host module of emodlib");
    add_malaria_bindings(malaria_m);

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
