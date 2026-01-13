/**
 * @file
 * @brief emodlib vector Python bindings.
 */

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "emodlib/vector/VectorConfig.h"
#include "emodlib/vector/VectorCohort.h"
#include "emodlib/vector/VectorPopulation.h"
#include "emodlib/vector/TransmissionSimulation.h"

namespace py = pybind11;
namespace emv = emodlib::vector;


void add_vector_bindings(py::module& m) {

    using namespace emv;
    using namespace py::literals;


    // ==== VectorState enum ==== //
    py::enum_<VectorState>(m, "VectorState")
        .value("EGG", VectorState::EGG)
        .value("LARVA", VectorState::LARVA)
        .value("IMMATURE", VectorState::IMMATURE)
        .value("ADULT", VectorState::ADULT)
        .value("INFECTED", VectorState::INFECTED)
        .value("INFECTIOUS", VectorState::INFECTIOUS)
        .value("MALE", VectorState::MALE)
        .export_values();


    // ==== VectorConfig ==== //
    py::class_<VectorConfig>(m, "VectorConfig")

        .def(py::init<>())

        // Arrhenius parameters
        .def_readwrite("aquatic_arrhenius1", &VectorConfig::aquatic_arrhenius1)
        .def_readwrite("aquatic_arrhenius2", &VectorConfig::aquatic_arrhenius2)
        .def_readwrite("infected_arrhenius1", &VectorConfig::infected_arrhenius1)
        .def_readwrite("infected_arrhenius2", &VectorConfig::infected_arrhenius2)
        .def_readwrite("cycle_arrhenius1", &VectorConfig::cycle_arrhenius1)
        .def_readwrite("cycle_arrhenius2", &VectorConfig::cycle_arrhenius2)

        // Life expectancy and mortality
        .def_readwrite("adult_life_expectancy", &VectorConfig::adult_life_expectancy)
        .def_readwrite("male_life_expectancy", &VectorConfig::male_life_expectancy)
        .def_readwrite("aquatic_mortality_rate", &VectorConfig::aquatic_mortality_rate)
        .def_readwrite("immature_duration", &VectorConfig::immature_duration)
        .def_readwrite("enable_age_dependent_mortality", &VectorConfig::enable_age_dependent_mortality)

        // Feeding and host-seeking
        .def_readwrite("anthropophily", &VectorConfig::anthropophily)
        .def_readwrite("indoor_feeding_fraction", &VectorConfig::indoor_feeding_fraction)
        .def_readwrite("days_between_feeds", &VectorConfig::days_between_feeds)

        // Reproduction
        .def_readwrite("egg_batch_size", &VectorConfig::egg_batch_size)
        .def_readwrite("infected_egg_batch_factor", &VectorConfig::infected_egg_batch_factor)
        .def_readwrite("egg_survival_rate", &VectorConfig::egg_survival_rate)
        .def_readwrite("egg_hatch_duration", &VectorConfig::egg_hatch_duration)

        // Transmission
        .def_readwrite("transmission_rate", &VectorConfig::transmission_rate)
        .def_readwrite("acquire_modifier", &VectorConfig::acquire_modifier)
        .def_readwrite("sporozoites_per_bite", &VectorConfig::sporozoites_per_bite)

        // Population
        .def_readwrite("carrying_capacity", &VectorConfig::carrying_capacity)
        .def_readwrite("initial_infected_fraction", &VectorConfig::initial_infected_fraction)
        .def_readwrite("enable_lifecycle", &VectorConfig::enable_lifecycle)

        // Computed values (read-only)
        .def_property_readonly("adult_mortality_rate", &VectorConfig::GetAdultMortalityRate)
        .def_property_readonly("male_mortality_rate", &VectorConfig::GetMaleMortalityRate)
        .def_property_readonly("biting_rate", &VectorConfig::GetBitingRate)
        .def_property_readonly("daily_emergence", &VectorConfig::GetDailyEmergence)

        // Helper methods
        .def("get_larval_progress", &VectorConfig::GetLarvalProgress,
             "Calculate larval development progress for one timestep",
             "temp_celsius"_a, "dt"_a)
        .def("get_infected_progress", &VectorConfig::GetInfectedProgress,
             "Calculate infected (EIP) progress for one timestep",
             "temp_celsius"_a, "dt"_a)
        .def("get_survival_probability", &VectorConfig::GetSurvivalProbability,
             "Get combined survival probability for one timestep",
             "age_days"_a, "base_mortality_rate"_a)
        .def_static("get_age_dependent_mortality", &VectorConfig::GetAgeDependentMortality,
                    "Calculate age-dependent mortality rate (Styer et al. 2007)",
                    "age_days"_a);


    // ==== VectorCohort ==== //
    py::class_<VectorCohort>(m, "VectorCohort")

        .def(py::init<>())
        .def(py::init<float, VectorState, float>(),
             "pop"_a, "state"_a, "age_days"_a = 0.0f)

        // Core state
        .def_readwrite("population", &VectorCohort::population)
        .def_readwrite("age", &VectorCohort::age)
        .def_readwrite("progress", &VectorCohort::progress)
        .def_readwrite("state", &VectorCohort::state)
        .def_readwrite("days_since_infection", &VectorCohort::days_since_infection)

        // State queries
        .def("is_susceptible", &VectorCohort::IsSusceptible)
        .def("is_exposed", &VectorCohort::IsExposed)
        .def("is_infectious", &VectorCohort::IsInfectious)
        .def("is_aquatic", &VectorCohort::IsAquatic)
        .def("is_adult", &VectorCohort::IsAdult)
        .def("is_empty", &VectorCohort::IsEmpty)
        .def("is_progress_complete", &VectorCohort::IsProgressComplete)

        // State transitions
        .def("apply_mortality", &VectorCohort::ApplyMortality,
             "Apply mortality - reduce population by survival probability",
             "survival_prob"_a)
        .def("increment_age", &VectorCohort::IncrementAge,
             "Increment age by timestep",
             "dt"_a)
        .def("add_progress", &VectorCohort::AddProgress,
             "Add development progress",
             "delta_progress"_a)
        .def("reset_progress", &VectorCohort::ResetProgress)
        .def("transition_to", &VectorCohort::TransitionTo,
             "Transition to next lifecycle state",
             "new_state"_a)
        .def("increment_infection_days", &VectorCohort::IncrementInfectionDays,
             "Increment infection days (for EIP tracking)",
             "dt"_a = 1)

        // Splitting and merging
        .def("split", &VectorCohort::Split,
             "Split off a fraction of this cohort",
             "fraction"_a)
        .def("split_count", &VectorCohort::SplitCount,
             "Split off a specific number from this cohort",
             "count"_a)
        .def("merge", &VectorCohort::Merge,
             "Merge another cohort into this one",
             "other"_a);


    // ==== VectorPopulation ==== //
    py::class_<VectorPopulation>(m, "VectorPopulation")

        .def(py::init<>())

        .def("initialize",
             py::overload_cast<const VectorConfig&>(&VectorPopulation::Initialize),
             "Initialize population with configuration",
             "config"_a)

        .def("initialize",
             py::overload_cast<const VectorConfig&, uint32_t>(&VectorPopulation::Initialize),
             "Initialize with specific random seed",
             "config"_a, "seed"_a)

        .def("update", &VectorPopulation::Update,
             "Update population for one timestep",
             "dt"_a, "temperature"_a, "human_infectiousness"_a)

        .def("get_sporozoite_challenges", &VectorPopulation::GetSporozoiteChallenges,
             "Get sporozoite challenges for each human",
             "n_humans"_a)

        // Population queries
        .def_property_readonly("total_population", &VectorPopulation::GetTotalPopulation)
        .def_property_readonly("adult_population", &VectorPopulation::GetAdultPopulation)
        .def_property_readonly("susceptible_count", &VectorPopulation::GetSusceptibleCount)
        .def_property_readonly("exposed_count", &VectorPopulation::GetExposedCount)
        .def_property_readonly("infectious_count", &VectorPopulation::GetInfectiousCount)
        .def_property_readonly("infectious_fraction", &VectorPopulation::GetInfectiousFraction)
        .def_property_readonly("infected_fraction", &VectorPopulation::GetInfectedFraction)
        .def_property_readonly("egg_count", &VectorPopulation::GetEggCount)
        .def_property_readonly("larva_count", &VectorPopulation::GetLarvaCount)
        .def_property_readonly("immature_count", &VectorPopulation::GetImmatureCount)
        .def_property_readonly("current_day", &VectorPopulation::GetCurrentDay);


    // ==== TransmissionSimulation ==== //
    py::class_<TransmissionSimulation>(m, "TransmissionSimulation")

        .def(py::init<>())

        .def("initialize", &TransmissionSimulation::Initialize,
             "Initialize simulation with configurations",
             "malaria_config"_a, "vector_config"_a, "n_humans"_a, "seed"_a = 42)

        .def("update", &TransmissionSimulation::Update,
             "Update simulation for one timestep",
             "dt"_a, "temperature"_a)

        .def("seed_infections", &TransmissionSimulation::SeedInfections,
             "Seed initial infections in humans",
             "n_infections"_a)

        .def("challenge_human", &TransmissionSimulation::ChallengeHuman,
             "Challenge a specific human with sporozoites",
             "human_index"_a)

        // Human population metrics
        .def_property_readonly("prevalence", &TransmissionSimulation::GetPrevalence)
        .def_property_readonly("mean_parasite_density", &TransmissionSimulation::GetMeanParasiteDensity)
        .def_property_readonly("mean_gametocyte_density", &TransmissionSimulation::GetMeanGametocyteDensity)
        .def_property_readonly("mean_infectiousness", &TransmissionSimulation::GetMeanInfectiousness)
        .def_property_readonly("infected_count", &TransmissionSimulation::GetInfectedCount)

        // Vector population metrics
        .def_property_readonly("vector_population", &TransmissionSimulation::GetVectorPopulation)
        .def_property_readonly("vector_infectious_fraction", &TransmissionSimulation::GetVectorInfectiousFraction)
        .def_property_readonly("vector_infected_fraction", &TransmissionSimulation::GetVectorInfectedFraction)
        .def_property_readonly("vector_susceptible_count", &TransmissionSimulation::GetVectorSusceptibleCount)
        .def_property_readonly("vector_exposed_count", &TransmissionSimulation::GetVectorExposedCount)
        .def_property_readonly("vector_infectious_count", &TransmissionSimulation::GetVectorInfectiousCount)

        // State
        .def_property_readonly("current_day", &TransmissionSimulation::GetCurrentDay)
        .def_property_readonly("n_humans", &TransmissionSimulation::GetNumHumans);

}
