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
    using namespace py::literals;


    // ==== Binding of the intrahost component ==== //
    py::class_<IntrahostComponent> (m, "IntrahostComponent")

        .def_static("create", &IntrahostComponent::Create)

        .def_static("configure",
                    &IntrahostComponent::params::Configure,
                    "Configure the component from a ParamSet dictionary",
                    "pset"_a)

        .def("update",
             &IntrahostComponent::Update,
             "Update the intrahost model state by dt",
             "dt"_a)

        .def("challenge",
             &IntrahostComponent::Challenge,
             "Challenge with a new infection")

        .def("treat",
             &IntrahostComponent::Treat,
             "Treat and clear all infections")

        .def_property_readonly("n_infections", &IntrahostComponent::GetNumInfections)

        .def_property_readonly("parasite_density", &IntrahostComponent::GetParasiteDensity)
        .def_property_readonly("gametocyte_density", &IntrahostComponent::GetGametocyteDensity)
        .def_property_readonly("fever_temperature", &IntrahostComponent::GetFeverTemperature);

    // TODO: emodlib#9 (readwrite for init) + emodlib#11 (readonly for testing)
    // py::class_<Infection>
    // py::class_<Susceptibility>
    // py::class_<IMalariaAntibody>

}
