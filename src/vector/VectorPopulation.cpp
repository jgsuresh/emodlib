/**
 * @file VectorPopulation.cpp
 *
 * @brief Implementation of vector population dynamics.
 */

#include "emodlib/vector/VectorPopulation.h"

#include <cmath>
#include <numeric>

namespace emodlib
{
    namespace vector
    {
        // =====================================================
        // Construction and initialization
        // =====================================================

        VectorPopulation::VectorPopulation()
            : m_current_day(0)
        {
        }

        void VectorPopulation::Initialize(const VectorConfig& config)
        {
            Initialize(config, 42);  // Default seed
        }

        void VectorPopulation::Initialize(const VectorConfig& config, uint32_t seed)
        {
            m_config = config;
            m_rng.seed(seed);
            m_current_day = 0;

            // Clear any existing cohorts
            m_eggs.clear();
            m_larvae.clear();
            m_immature.clear();
            m_adults.clear();
            m_infected.clear();
            m_infectious.clear();

            // Initialize with equilibrium population
            InitializeEquilibriumPopulation();
        }

        void VectorPopulation::InitializeEquilibriumPopulation()
        {
            // At equilibrium, number of mosquitoes of age a is:
            // N(a) = N(0) * survival^a
            // where N(0) = daily_emergence

            float daily_emergence = m_config.GetDailyEmergence();
            float base_mortality = m_config.GetAdultMortalityRate();

            // Create cohorts up to ~3 lifespans (captures >95% of population)
            int max_age = static_cast<int>(3 * m_config.adult_life_expectancy);

            for (int age = 0; age < max_age; ++age)
            {
                // Calculate expected population at this age
                float survival_to_age = std::pow(
                    m_config.GetSurvivalProbability(static_cast<float>(age), base_mortality),
                    static_cast<float>(age)
                );
                float count = daily_emergence * survival_to_age;

                if (count < 0.1f)
                {
                    break;  // Stop when negligible
                }

                // Handle initial infected fraction if specified
                if (m_config.initial_infected_fraction > 0.0f && age >= 10)
                {
                    // Some old mosquitoes start infectious
                    float infected_count = count * m_config.initial_infected_fraction;
                    float susceptible_count = count - infected_count;

                    if (susceptible_count > 0.1f)
                    {
                        VectorCohort cohort(susceptible_count, VectorState::ADULT,
                                           static_cast<float>(age));
                        m_adults.push_back(cohort);
                    }

                    if (infected_count > 0.1f)
                    {
                        VectorCohort cohort(infected_count, VectorState::INFECTIOUS,
                                           static_cast<float>(age));
                        cohort.days_since_infection = 10;  // Already past EIP
                        m_infectious.push_back(cohort);
                    }
                }
                else
                {
                    VectorCohort cohort(count, VectorState::ADULT, static_cast<float>(age));
                    m_adults.push_back(cohort);
                }
            }
        }

        // =====================================================
        // Main update
        // =====================================================

        void VectorPopulation::Update(float dt, float temperature,
                                      const std::vector<float>& human_infectiousness)
        {
            // Update order follows EMOD convention:
            // infectious -> infected -> adult -> immature -> larval -> egg -> laying

            Update_Infectious_Queue(dt, temperature);
            Update_Infected_Queue(dt, temperature);
            Update_Adult_Queue(dt, temperature, human_infectiousness);
            Update_Immature_Queue(dt, temperature);
            Update_Larval_Queue(dt, temperature);
            Update_Egg_Queue(dt, temperature);

            // Add new eggs from adult oviposition
            AddEmergence();

            // Remove empty cohorts
            CleanupEmptyCohorts();

            m_current_day++;
        }

        // =====================================================
        // Population queries
        // =====================================================

        float VectorPopulation::GetTotalInQueue(const std::vector<VectorCohort>& queue) const
        {
            float total = 0.0f;
            for (const auto& cohort : queue)
            {
                total += cohort.population;
            }
            return total;
        }

        float VectorPopulation::GetTotalPopulation() const
        {
            return GetEggCount() + GetLarvaCount() + GetImmatureCount() +
                   GetAdultPopulation();
        }

        float VectorPopulation::GetAdultPopulation() const
        {
            return GetSusceptibleCount() + GetExposedCount() + GetInfectiousCount();
        }

        float VectorPopulation::GetSusceptibleCount() const
        {
            return GetTotalInQueue(m_adults);
        }

        float VectorPopulation::GetExposedCount() const
        {
            return GetTotalInQueue(m_infected);
        }

        float VectorPopulation::GetInfectiousCount() const
        {
            return GetTotalInQueue(m_infectious);
        }

        float VectorPopulation::GetEggCount() const
        {
            return GetTotalInQueue(m_eggs);
        }

        float VectorPopulation::GetLarvaCount() const
        {
            return GetTotalInQueue(m_larvae);
        }

        float VectorPopulation::GetImmatureCount() const
        {
            return GetTotalInQueue(m_immature);
        }

        float VectorPopulation::GetInfectiousFraction() const
        {
            float adult_pop = GetAdultPopulation();
            if (adult_pop <= 0.0f)
            {
                return 0.0f;
            }
            return GetInfectiousCount() / adult_pop;
        }

        float VectorPopulation::GetInfectedFraction() const
        {
            float adult_pop = GetAdultPopulation();
            if (adult_pop <= 0.0f)
            {
                return 0.0f;
            }
            return (GetExposedCount() + GetInfectiousCount()) / adult_pop;
        }

        // =====================================================
        // Transmission interface
        // =====================================================

        std::vector<float> VectorPopulation::GetSporozoiteChallenges(int n_humans)
        {
            std::vector<float> challenges(n_humans, 0.0f);

            if (n_humans <= 0)
            {
                return challenges;
            }

            float infectious_count = GetInfectiousCount();
            float total_adults = GetAdultPopulation();

            if (infectious_count <= 0.0f || total_adults <= 0.0f)
            {
                return challenges;
            }

            // Bites per human per day
            float biting_rate = m_config.GetBitingRate();
            float bites_per_human = total_adults * biting_rate *
                                    m_config.anthropophily / n_humans;

            float infectious_fraction = infectious_count / total_adults;
            float infectious_bites_per_human = bites_per_human * infectious_fraction;

            // Poisson draw for infectious bites per human
            std::poisson_distribution<int> bite_dist(infectious_bites_per_human);
            std::bernoulli_distribution transmit_dist(m_config.transmission_rate);

            for (int i = 0; i < n_humans; ++i)
            {
                int n_infectious_bites = bite_dist(m_rng);

                for (int b = 0; b < n_infectious_bites; ++b)
                {
                    if (transmit_dist(m_rng))
                    {
                        challenges[i] += m_config.sporozoites_per_bite;
                    }
                }
            }

            return challenges;
        }

        // =====================================================
        // Cleanup
        // =====================================================

        void VectorPopulation::CleanupEmptyCohorts()
        {
            auto remove_empty = [](std::vector<VectorCohort>& queue)
            {
                queue.erase(
                    std::remove_if(queue.begin(), queue.end(),
                                   [](const VectorCohort& c) { return c.IsEmpty(); }),
                    queue.end()
                );
            };

            remove_empty(m_eggs);
            remove_empty(m_larvae);
            remove_empty(m_immature);
            remove_empty(m_adults);
            remove_empty(m_infected);
            remove_empty(m_infectious);
        }

        // =====================================================
        // Stub implementations (filled in chunks 4-7)
        // =====================================================

        void VectorPopulation::ApplyMortality(float temperature)
        {
            float base_adult_mortality = m_config.GetAdultMortalityRate();

            // Aquatic mortality (eggs, larvae, immature)
            float aquatic_survival = 1.0f - m_config.aquatic_mortality_rate;
            float egg_survival = m_config.egg_survival_rate;

            for (auto& cohort : m_eggs)
            {
                cohort.ApplyMortality(egg_survival);
            }

            for (auto& cohort : m_larvae)
            {
                cohort.ApplyMortality(aquatic_survival);
            }

            for (auto& cohort : m_immature)
            {
                cohort.ApplyMortality(aquatic_survival);
            }

            // Adult mortality (age-dependent via Styer et al.)
            for (auto& cohort : m_adults)
            {
                float survival = m_config.GetSurvivalProbability(
                    cohort.age, base_adult_mortality);
                cohort.ApplyMortality(survival);
            }

            for (auto& cohort : m_infected)
            {
                float survival = m_config.GetSurvivalProbability(
                    cohort.age, base_adult_mortality);
                cohort.ApplyMortality(survival);
            }

            for (auto& cohort : m_infectious)
            {
                float survival = m_config.GetSurvivalProbability(
                    cohort.age, base_adult_mortality);
                cohort.ApplyMortality(survival);
            }
        }

        void VectorPopulation::UpdateDevelopment(float dt, float temperature)
        {
            // Larval development progress (Arrhenius equation)
            float larval_progress = m_config.GetLarvalProgress(temperature, dt);

            for (auto& cohort : m_larvae)
            {
                cohort.AddProgress(larval_progress);
            }

            // Immature development (fixed rate, not temperature-dependent)
            float immature_progress = dt / m_config.immature_duration;

            for (auto& cohort : m_immature)
            {
                cohort.AddProgress(immature_progress);
            }

            // Infected (EIP) development progress (Arrhenius equation)
            float infected_progress = m_config.GetInfectedProgress(temperature, dt);

            for (auto& cohort : m_infected)
            {
                cohort.AddProgress(infected_progress);
                cohort.IncrementInfectionDays(static_cast<int>(dt));
            }
        }

        void VectorPopulation::ProcessTransitions()
        {
            // TODO: Chunk 6
        }

        void VectorPopulation::ProcessTransmission(const std::vector<float>& human_infectiousness)
        {
            // TODO: Chunk 7
        }

        void VectorPopulation::AddEmergence()
        {
            // Add new adult cohort from daily emergence
            float daily_emergence = m_config.GetDailyEmergence();
            if (daily_emergence > 0.1f)
            {
                VectorCohort new_adults(daily_emergence, VectorState::ADULT, 0.0f);
                m_adults.push_back(new_adults);
            }
        }

        void VectorPopulation::Update_Egg_Queue(float dt, float temperature)
        {
            // TODO: Chunk 6 - egg hatching
        }

        void VectorPopulation::Update_Larval_Queue(float dt, float temperature)
        {
            // Temperature-dependent larval development (Arrhenius equation)
            // progress = A1 * exp(-A2 / (T + 273.15)) * dt
            float larval_progress = m_config.GetLarvalProgress(temperature, dt);

            // Aquatic mortality
            float aquatic_survival = 1.0f - m_config.aquatic_mortality_rate;

            // Track cohorts that complete development
            std::vector<VectorCohort> new_immature;

            for (auto& cohort : m_larvae)
            {
                // Apply mortality first
                cohort.ApplyMortality(aquatic_survival);

                // Add development progress
                cohort.AddProgress(larval_progress);
                cohort.IncrementAge(dt);

                // Check for transition to immature stage
                if (cohort.IsProgressComplete())
                {
                    VectorCohort transitioning = cohort;
                    transitioning.TransitionTo(VectorState::IMMATURE);
                    new_immature.push_back(transitioning);
                    cohort.population = 0.0f;  // Mark for removal
                }
            }

            // Add newly transitioned cohorts to immature queue
            for (auto& cohort : new_immature)
            {
                m_immature.push_back(cohort);
            }
        }

        void VectorPopulation::Update_Immature_Queue(float dt, float temperature)
        {
            // Immature development (fixed duration, not temperature-dependent)
            // progress = dt / immature_duration
            float immature_progress = dt / m_config.immature_duration;

            // Aquatic mortality
            float aquatic_survival = 1.0f - m_config.aquatic_mortality_rate;

            // Track cohorts that complete development (become adults)
            std::vector<VectorCohort> new_adults;

            for (auto& cohort : m_immature)
            {
                // Apply mortality first
                cohort.ApplyMortality(aquatic_survival);

                // Add development progress
                cohort.AddProgress(immature_progress);
                cohort.IncrementAge(dt);

                // Check for transition to adult stage
                if (cohort.IsProgressComplete())
                {
                    VectorCohort transitioning = cohort;
                    transitioning.TransitionTo(VectorState::ADULT);
                    new_adults.push_back(transitioning);
                    cohort.population = 0.0f;  // Mark for removal
                }
            }

            // Add newly emerged adults to adult queue
            for (auto& cohort : new_adults)
            {
                m_adults.push_back(cohort);
            }
        }

        void VectorPopulation::Update_Adult_Queue(float dt, float temperature,
                                                  const std::vector<float>& human_infectiousness)
        {
            float base_mortality = m_config.GetAdultMortalityRate();

            for (auto& cohort : m_adults)
            {
                // Age-dependent mortality (Styer et al.)
                float survival = m_config.GetSurvivalProbability(cohort.age, base_mortality);
                cohort.ApplyMortality(survival);

                // Increment age
                cohort.IncrementAge(dt);
            }

            // Transmission from humans to vectors handled in Chunk 7
        }

        void VectorPopulation::Update_Infected_Queue(float dt, float temperature)
        {
            // Temperature-dependent EIP progression (Arrhenius equation)
            // progress = A1 * exp(-A2 / (T + 273.15)) * dt
            // When progress >= 1.0, mosquito becomes infectious
            float infected_progress = m_config.GetInfectedProgress(temperature, dt);
            float base_mortality = m_config.GetAdultMortalityRate();

            // Track cohorts that complete EIP (become infectious)
            std::vector<VectorCohort> new_infectious;

            for (auto& cohort : m_infected)
            {
                // Age-dependent mortality (Styer et al.)
                float survival = m_config.GetSurvivalProbability(cohort.age, base_mortality);
                cohort.ApplyMortality(survival);

                // Add EIP progress
                cohort.AddProgress(infected_progress);
                cohort.IncrementAge(dt);
                cohort.IncrementInfectionDays(static_cast<int>(dt));

                // Check for transition to infectious stage
                if (cohort.IsProgressComplete())
                {
                    VectorCohort transitioning = cohort;
                    transitioning.TransitionTo(VectorState::INFECTIOUS);
                    new_infectious.push_back(transitioning);
                    cohort.population = 0.0f;  // Mark for removal
                }
            }

            // Add newly infectious cohorts to infectious queue
            for (auto& cohort : new_infectious)
            {
                m_infectious.push_back(cohort);
            }
        }

        void VectorPopulation::Update_Infectious_Queue(float dt, float temperature)
        {
            float base_mortality = m_config.GetAdultMortalityRate();

            for (auto& cohort : m_infectious)
            {
                // Age-dependent mortality (Styer et al.)
                float survival = m_config.GetSurvivalProbability(cohort.age, base_mortality);
                cohort.ApplyMortality(survival);

                // Increment age and infection days
                cohort.IncrementAge(dt);
                cohort.IncrementInfectionDays(static_cast<int>(dt));
            }
        }

    }  // namespace vector

}  // namespace emodlib
