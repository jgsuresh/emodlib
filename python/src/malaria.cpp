/**
 * @file
 * @brief emodlib malaria Python bindings.
*/

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "emodlib/malaria/IntrahostComponent.h"

namespace py = pybind11;
namespace emm = emodlib::malaria;


void add_malaria_bindings(py::module& m) {

    using namespace emm;

    // ==== Binding of the intrahost component ==== //
    py::class_<IntrahostComponent> (m, "IntrahostComponent")
    
//        .def(py::init<>())
    
        .def("create", &IntrahostComponent::Create)
        .def("update", &IntrahostComponent::Update)
        .def("challenge", &IntrahostComponent::Challenge)
        .def("treat", &IntrahostComponent::Treat)
    
        .def_property_readonly("parasite_density", &IntrahostComponent::GetParasiteDensity)
        .def_property_readonly("gametocyte_density", &IntrahostComponent::GetGametocyteDensity)
        .def_property_readonly("fever_temperature", &IntrahostComponent::GetFeverTemperature);

}
