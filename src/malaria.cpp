/**
 * @file
 * @brief emodlib malaria Python bindings.
*/

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "emodlib/malaria/IntrahostComponent.h"
#include "emodlib/malaria/MalariaAntibody.h"

namespace py = pybind11;
namespace emm = emodlib::malaria;


// Following template trampoline pattern for combining virtual functions + inheritance:
// https://pybind11.readthedocs.io/en/stable/advanced/classes.html#combining-virtual-functions-and-inheritance

template <class IAntibodyBase = emm::IMalariaAntibody>
class PyIMalariaAntibody : public IAntibodyBase {
public:
     using IAntibodyBase::IAntibodyBase; // Inherit constructors
     void IncreaseAntigenCount(int64_t antigenCount) override {
          PYBIND11_OVERRIDE_PURE(void, IAntibodyBase, IncreaseAntigenCount, antigenCount); }
     int64_t GetAntigenCount() const override {
          PYBIND11_OVERRIDE_PURE(int64_t, IAntibodyBase, GetAntigenCount, ); }
     float GetAntibodyCapacity() const override {
          PYBIND11_OVERRIDE_PURE(float, IAntibodyBase, GetAntibodyCapacity, ); }
     float GetAntibodyConcentration() const override {
          PYBIND11_OVERRIDE_PURE(float, IAntibodyBase, GetAntibodyConcentration, ); }
};

template <class MalariaAntibodyBase = emm::MalariaAntibody>
class PyMalariaAntibody : public PyIMalariaAntibody<MalariaAntibodyBase> {
public:
     using PyIMalariaAntibody<MalariaAntibodyBase>::PyIMalariaAntibody; // Inherit constructors
     void IncreaseAntigenCount(int64_t antigenCount) override {
          PYBIND11_OVERRIDE(void, MalariaAntibodyBase, IncreaseAntigenCount, antigenCount); }
     int64_t GetAntigenCount() const override {
          PYBIND11_OVERRIDE(int64_t, MalariaAntibodyBase, GetAntigenCount, ); }
     float GetAntibodyCapacity() const override {
          PYBIND11_OVERRIDE(float, MalariaAntibodyBase, GetAntibodyCapacity, ); }
     float GetAntibodyConcentration() const override {
          PYBIND11_OVERRIDE(float, MalariaAntibodyBase, GetAntibodyConcentration, ); }
};


void add_malaria_bindings(py::module& m) {

    using namespace emm;
    using namespace py::literals;


    // ==== Binding of the intrahost component ==== //
    py::class_<IntrahostComponent> (m, "IntrahostComponent")

        .def_static("create", &IntrahostComponent::Create)

        .def_static("_configure_from_params",
                    &IntrahostComponent::params::Configure,
                    "Configure the IntrahostComponent params from a ParamSet dictionary",
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
        .def_property_readonly("fever_temperature", &IntrahostComponent::GetFeverTemperature)

        .def_property_readonly("infectiousness", &IntrahostComponent::GetInfectiousness)

        .def_property_readonly("susceptibility", &IntrahostComponent::GetSusceptibility)
        .def_property_readonly("infections", &IntrahostComponent::GetInfections);


    // TODO: emodlib#9 (readwrite for init) + emodlib#11 (readonly for testing)
    // py::class_<Infection>
    // py::class_<Susceptibility>
    // py::class_<IMalariaAntibody>


    py::class_<Susceptibility> (m, "Susceptibility")

          .def_static("create", &Susceptibility::Create)

          .def_static("configure",
                      &Susceptibility::params::Configure,
                      "Configure the Susceptibility params from a ParamSet dictionary",
                      "pset"_a)

          .def("update",
               &Susceptibility::Update,
               "Update the susceptibility state by dt",
               "dt"_a)

          .def_property("age", &Susceptibility::get_age, &Susceptibility::set_age)

          .def_property("maternal_antibody_strength",
                        &Susceptibility::get_maternal_antibody_strength,
                        &Susceptibility::set_maternal_antibody_strength)

          .def_property("pyrogenic_threshold",
                        &Susceptibility::get_pyrogenic_threshold,
                        &Susceptibility::set_pyrogenic_threshold)

          .def_property("fever_kill_rate",
                        &Susceptibility::get_fever_kill_rate,
                        &Susceptibility::set_fever_kill_rate);


     py::class_<Infection> (m, "Infection")

          .def_static("create", &Infection::Create,
               "Create an Infection object with pointer to Susceptibility",
               "susceptibility"_a, "hepatocytes"_a=1)

          .def_static("configure",
                      &Infection::params::Configure,
                      "Configure the Infection params from a ParamSet dictionary",
                      "pset"_a)

          .def("update",
               &Infection::Update,
               "Update the infection state by dt",
               "dt"_a)

          .def_property_readonly("msp_type", &Infection::get_msp_type)

          .def_property_readonly("pfemp1_major_types", &Infection::get_pfemp1_major_types)

          .def_property_readonly("msp_antibody", &Infection::get_msp_antibody);


     py::class_<IMalariaAntibody, PyIMalariaAntibody<>> (m, "IMalariaAntibody");

     py::class_<MalariaAntibody, IMalariaAntibody, PyMalariaAntibody<>> (m, "MalariaAntibody")
          .def_property_readonly("antigen_count", &MalariaAntibody::GetAntigenCount)
          .def_property_readonly("antibody_capacity", &MalariaAntibody::GetAntibodyCapacity)
          .def_property_readonly("antibody_concentration", &MalariaAntibody::GetAntibodyConcentration);


     py::class_<MalariaAntibodyMSP, MalariaAntibody, PyMalariaAntibody<MalariaAntibodyMSP>> (m, "MalariaAntibodyMSP");

}
