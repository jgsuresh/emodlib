/**
 * @file
 * @brief emodlib malaria Python bindings.
*/

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "emodlib/malaria/MalariaConfig.h"
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
     // EMOD 2.22: Updated signature with currentTime and dt
     void IncreaseAntigenCount(int64_t antigenCount, float currentTime, float dt) override {
          PYBIND11_OVERRIDE_PURE(void, IAntibodyBase, IncreaseAntigenCount, antigenCount, currentTime, dt); }
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
     // EMOD 2.22: Updated signature with currentTime and dt
     void IncreaseAntigenCount(int64_t antigenCount, float currentTime, float dt) override {
          PYBIND11_OVERRIDE(void, MalariaAntibodyBase, IncreaseAntigenCount, antigenCount, currentTime, dt); }
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


    // ==== Binding of the MalariaConfig ==== //
    py::class_<MalariaConfig, std::shared_ptr<MalariaConfig>> (m, "MalariaConfig", py::dynamic_attr())

        .def(py::init<>())

        .def("configure", &MalariaConfig::Configure,
             "Configure from a ParamSet dictionary",
             "pset"_a)

        .def_static("from_params", &MalariaConfig::FromParamSet,
                    "Create a configured MalariaConfig from a ParamSet dictionary",
                    "pset"_a)

        // IntrahostComponent-level params
        .def_readwrite("random_seed", &MalariaConfig::randomSeed)
        .def_readwrite("max_ind_inf", &MalariaConfig::max_ind_inf)
        .def_readwrite("falciparum_MSP_variants", &MalariaConfig::falciparumMSPVars)
        .def_readwrite("falciparum_nonspecific_types", &MalariaConfig::falciparumNonSpecTypes)
        .def_readwrite("falciparum_PfEMP1_variants", &MalariaConfig::falciparumPfEMP1Vars)
        .def_readwrite("base_gametocyte_mosquito_survival", &MalariaConfig::base_gametocyte_mosquito_survival)
        .def_readwrite("cytokine_gametocyte_inactivation", &MalariaConfig::cytokine_gametocyte_inactivation)

        // Susceptibility params
        .def_readwrite("memory_level", &MalariaConfig::memory_level)
        .def_readwrite("hyperimmune_decay_rate", &MalariaConfig::hyperimmune_decay_rate)
        .def_readwrite("MSP1_antibody_growthrate", &MalariaConfig::MSP1_antibody_growthrate)
        .def_readwrite("antibody_stimulation_c50", &MalariaConfig::antibody_stimulation_c50)
        .def_readwrite("antibody_capacity_growthrate", &MalariaConfig::antibody_capacity_growthrate)
        .def_readwrite("minimum_adapted_response", &MalariaConfig::minimum_adapted_response)
        .def_readwrite("non_specific_growth", &MalariaConfig::non_specific_growth)
        .def_readwrite("antibody_csp_decay_days", &MalariaConfig::antibody_csp_decay_days)
        .def_readwrite("antibody_days_to_long_term_decay", &MalariaConfig::antibody_days_to_long_term_decay)
        .def_readwrite("antibody_long_term_decay_days", &MalariaConfig::antibody_long_term_decay_days)
        .def_readwrite("maternal_antibody_decay_rate", &MalariaConfig::maternal_antibody_decay_rate)
        .def_readwrite("PfHRP2_boost_rate", &MalariaConfig::PfHRP2_boost_rate)
        .def_readwrite("PfHRP2_decay_rate", &MalariaConfig::PfHRP2_decay_rate)
        .def_readwrite("pyrogenic_threshold", &MalariaConfig::pyrogenic_threshold)
        .def_readwrite("fever_IRBC_killrate", &MalariaConfig::fever_IRBC_killrate)
        .def_readwrite("erythropoiesis_anemia_effect", &MalariaConfig::erythropoiesis_anemia_effect)

        // Infection params
        .def_readwrite("incubation_period", &MalariaConfig::incubation_period)
        .def_readwrite("antibody_IRBC_killrate", &MalariaConfig::antibody_IRBC_killrate)
        .def_readwrite("non_specific_antigenicity", &MalariaConfig::non_specific_antigenicity)
        .def_readwrite("MSP1_merozoite_kill", &MalariaConfig::MSP1_merozoite_kill)
        .def_readwrite("gametocyte_stage_survival", &MalariaConfig::gametocyte_stage_survival)
        .def_readwrite("base_gametocyte_sexratio", &MalariaConfig::base_gametocyte_sexratio)
        .def_readwrite("base_gametocyte_production", &MalariaConfig::base_gametocyte_production)
        .def_readwrite("antigen_switch_rate", &MalariaConfig::antigen_switch_rate)
        .def_readwrite("merozoites_per_hepatocyte", &MalariaConfig::merozoites_per_hepatocyte)
        .def_readwrite("merozoites_per_schizont", &MalariaConfig::merozoites_per_schizont)
        .def_readwrite("RBC_destruction_multiplier", &MalariaConfig::RBC_destruction_multiplier)
        .def_readwrite("n_asexual_cycles_wo_gametocytes", &MalariaConfig::n_asexual_cycles_wo_gametocytes);


    // ==== Binding of the DebugAntigens struct (for forcing explicit antigens) ==== //
    py::class_<DebugAntigens> (m, "DebugAntigens")
        .def(py::init<>())
        .def_readwrite("msp_type", &DebugAntigens::msp_type)
        .def_readwrite("nonspec_type", &DebugAntigens::nonspec_type)
        .def_readwrite("irbc_types", &DebugAntigens::irbc_types)
        .def_readwrite("minor_epitope_types", &DebugAntigens::minor_epitope_types);


    // ==== Binding of the intrahost component ==== //
    py::class_<IntrahostComponent> (m, "IntrahostComponent")

        .def_static("create", &IntrahostComponent::Create,
                    "Create an IntrahostComponent with the given configuration",
                    "config"_a)

        .def("update",
             &IntrahostComponent::Update,
             "Update the intrahost model state by dt",
             "dt"_a)

        .def("challenge",
             &IntrahostComponent::Challenge,
             "Challenge with a new infection")

        .def("challenge_with_antigens",
             &IntrahostComponent::ChallengeWithAntigens,
             "Challenge with a new infection using explicit antigens (for debugging)",
             "antigens"_a)

        .def("treat",
             &IntrahostComponent::Treat,
             "Treat and clear all infections")

        .def_property_readonly("n_infections", &IntrahostComponent::GetNumInfections)

        .def_property_readonly("parasite_density", &IntrahostComponent::GetParasiteDensity)
        .def_property_readonly("gametocyte_density", &IntrahostComponent::GetGametocyteDensity)
        .def_property_readonly("fever_temperature", &IntrahostComponent::GetFeverTemperature)

        .def_property_readonly("infectiousness", &IntrahostComponent::GetInfectiousness)

        .def_property_readonly("susceptibility", &IntrahostComponent::GetSusceptibility)
        .def_property_readonly("infections", &IntrahostComponent::GetInfections)

        // Aggregate infection state (for validation without raw pointer access)
        .def_property_readonly("total_irbc_counts", &IntrahostComponent::get_total_irbc_counts)
        .def_property_readonly("total_hepatocytes", &IntrahostComponent::get_total_hepatocytes)
        .def_property_readonly("total_male_gametocytes", &IntrahostComponent::get_total_male_gametocytes)
        .def_property_readonly("total_female_gametocytes", &IntrahostComponent::get_total_female_gametocytes)
        .def_property_readonly("total_asexual_cycles", &IntrahostComponent::get_total_asexual_cycles);


    // TODO: emodlib#9 (readwrite for init) + emodlib#11 (readonly for testing)
    // py::class_<Infection>
    // py::class_<Susceptibility>
    // py::class_<IMalariaAntibody>


    py::class_<Susceptibility> (m, "Susceptibility")

          .def_static("create", &Susceptibility::Create,
                      "Create a Susceptibility with the given configuration",
                      "config"_a)

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
                        &Susceptibility::set_fever_kill_rate)

          // Existing readonly properties
          .def_property_readonly("rbc_count", &Susceptibility::get_RBC_count)
          .def_property_readonly("inv_microliters_blood", &Susceptibility::get_inv_microliters_blood)
          .def_property_readonly("rbc_availability", &Susceptibility::get_RBC_availability)
          .def_property_readonly("fever", &Susceptibility::get_fever)
          .def_property_readonly("fever_celsius", &Susceptibility::get_fever_celsius)
          .def_property_readonly("cytokines", &Susceptibility::get_cytokines)
          .def_property_readonly("fever_killing_rate", &Susceptibility::get_fever_killing_rate)
          .def_property_readonly("parasite_density", &Susceptibility::get_parasite_density)
          .def_property_readonly("maternal_antibodies", &Susceptibility::get_maternal_antibodies)

          // Additional getters for detailed state exposure (for validation)
          .def_property_readonly("rbc_capacity", &Susceptibility::get_rbc_capacity)
          .def_property_readonly("csp_antibody", &Susceptibility::get_csp_antibody)
          .def_property_readonly("active_msp_antibodies", &Susceptibility::get_active_msp_antibodies)
          .def_property_readonly("active_pfemp1_minor_antibodies", &Susceptibility::get_active_pfemp1_minor_antibodies)
          .def_property_readonly("active_pfemp1_major_antibodies", &Susceptibility::get_active_pfemp1_major_antibodies);


     py::class_<Infection> (m, "Infection")

          .def_static("create", &Infection::Create,
               "Create an Infection object with pointer to Susceptibility and configuration",
               "susceptibility"_a, "config"_a, "hepatocytes"_a=1)

          .def("update",
               &Infection::Update,
               "Update the infection state by dt",
               "dt"_a)

          .def_property_readonly("msp_type", &Infection::get_msp_type)
          .def_property_readonly("pfemp1_major_types", &Infection::get_pfemp1_major_types)
          .def_property_readonly("msp_antibody", &Infection::get_msp_antibody)

          // Existing getters
          .def("get_male_gametocytes", &Infection::get_MaleGametocytes, "stage"_a)
          .def("get_female_gametocytes", &Infection::get_FemaleGametocytes, "stage"_a)
          .def_property_readonly("asexual_density", &Infection::get_asexual_density)
          .def_property_readonly("mature_gametocyte_density", &Infection::get_mature_gametocyte_density)
          .def_property_readonly("is_cleared", &Infection::IsCleared)

          // Additional getters for detailed state exposure (for validation)
          .def_property_readonly("irbc_counts", &Infection::get_irbc_counts)
          .def_property_readonly("hepatocyte_count", &Infection::get_hepatocyte_count)
          .def_property_readonly("liver_stage_timer", &Infection::get_liver_stage_timer)
          .def_property_readonly("asexual_phase", &Infection::get_asexual_phase)
          .def_property_readonly("asexual_cycle_timer", &Infection::get_asexual_cycle_timer)
          .def_property_readonly("asexual_cycle_count", &Infection::get_asexual_cycle_count)
          .def_property_readonly("all_male_gametocytes", &Infection::get_all_male_gametocytes)
          .def_property_readonly("all_female_gametocytes", &Infection::get_all_female_gametocytes)
          .def_property_readonly("minor_epitope_types", &Infection::get_minor_epitope_types);


     py::class_<IMalariaAntibody, PyIMalariaAntibody<>> (m, "IMalariaAntibody")
          .def_property_readonly("antigen_count", &IMalariaAntibody::GetAntigenCount)
          .def_property_readonly("antibody_capacity", &IMalariaAntibody::GetAntibodyCapacity)
          .def_property_readonly("antibody_concentration", &IMalariaAntibody::GetAntibodyConcentration)
          // EMOD 2.22: Time tracking for deferred decay
          .def_property_readonly("time_last_active", &IMalariaAntibody::GetTimeLastActive)
          .def_property_readonly("active_index", &IMalariaAntibody::GetActiveIndex);

     py::class_<MalariaAntibody, IMalariaAntibody, PyMalariaAntibody<>> (m, "MalariaAntibody")
          .def_property_readonly("antigen_count", &MalariaAntibody::GetAntigenCount)
          .def_property_readonly("antibody_capacity", &MalariaAntibody::GetAntibodyCapacity)
          .def_property_readonly("antibody_concentration", &MalariaAntibody::GetAntibodyConcentration)
          .def_property_readonly("time_last_active", &MalariaAntibody::GetTimeLastActive)
          .def_property_readonly("active_index", &MalariaAntibody::GetActiveIndex);


     py::class_<MalariaAntibodyMSP, MalariaAntibody, PyMalariaAntibody<MalariaAntibodyMSP>> (m, "MalariaAntibodyMSP");

}
