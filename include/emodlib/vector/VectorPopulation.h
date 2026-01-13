/**
 * @file VectorPopulation.h
 *
 * @brief Vector population dynamics with cohort-based lifecycle tracking.
 *
 * Manages mosquito population through lifecycle stages:
 * Egg -> Larva -> Immature -> Adult -> Infected -> Infectious
 *
 * Temperature-dependent development via Arrhenius equations.
 * Age-dependent mortality via Styer et al. formula.
 */

#pragma once

#include <vector>
#include <random>
#include <algorithm>

#include "VectorConfig.h"
#include "VectorCohort.h"

namespace emodlib
{
    namespace vector
    {
        /**
         * @brief Manages a population of mosquito cohorts with transmission dynamics.
         *
         * Update cycle each timestep:
         * 1. Apply mortality to all cohorts
         * 2. Progress development (temperature-dependent)
         * 3. Handle stage transitions
         * 4. Process transmission (human <-> vector)
         * 5. Add new emergence (eggs from adults)
         */
        class VectorPopulation
        {
        public:
            // =====================================================
            // Construction and initialization
            // =====================================================

            VectorPopulation();

            /**
             * @brief Initialize population with configuration.
             *
             * Creates equilibrium age distribution of adults.
             *
             * @param config Vector configuration parameters
             */
            void Initialize(const VectorConfig& config);

            /**
             * @brief Initialize with specific random seed.
             */
            void Initialize(const VectorConfig& config, uint32_t seed);

            // =====================================================
            // Main update
            // =====================================================

            /**
             * @brief Update population for one timestep.
             *
             * @param dt Timestep in days (typically 1.0)
             * @param temperature Temperature in Celsius
             * @param human_infectiousness Infectiousness of each human [0,1]
             */
            void Update(float dt, float temperature,
                       const std::vector<float>& human_infectiousness);

            // =====================================================
            // Transmission interface
            // =====================================================

            /**
             * @brief Get sporozoite challenges for each human.
             *
             * @param n_humans Number of humans
             * @return Vector of sporozoite counts per human
             */
            std::vector<float> GetSporozoiteChallenges(int n_humans);

            // =====================================================
            // Population queries
            // =====================================================

            float GetTotalPopulation() const;
            float GetAdultPopulation() const;
            float GetSusceptibleCount() const;
            float GetExposedCount() const;
            float GetInfectiousCount() const;

            float GetInfectiousFraction() const;
            float GetInfectedFraction() const;

            // Lifecycle stage counts
            float GetEggCount() const;
            float GetLarvaCount() const;
            float GetImmatureCount() const;

            int GetCurrentDay() const { return m_current_day; }

        private:
            // =====================================================
            // Cohort queues by lifecycle stage
            // =====================================================

            std::vector<VectorCohort> m_eggs;
            std::vector<VectorCohort> m_larvae;
            std::vector<VectorCohort> m_immature;
            std::vector<VectorCohort> m_adults;      // Susceptible adults
            std::vector<VectorCohort> m_infected;    // In EIP
            std::vector<VectorCohort> m_infectious;  // Can transmit

            // =====================================================
            // Configuration and state
            // =====================================================

            VectorConfig m_config;
            std::mt19937 m_rng;
            int m_current_day;

            // =====================================================
            // Update methods (implemented in chunks 4-7)
            // =====================================================

            void ApplyMortality(float temperature);
            void UpdateDevelopment(float dt, float temperature);
            void ProcessTransitions();
            void ProcessTransmission(const std::vector<float>& human_infectiousness);
            void AddEmergence();
            void AddEggLaying();
            void CleanupEmptyCohorts();

            // =====================================================
            // Lifecycle stage updates
            // =====================================================

            void Update_Egg_Queue(float dt, float temperature);
            void Update_Larval_Queue(float dt, float temperature);
            void Update_Immature_Queue(float dt, float temperature);
            void Update_Adult_Queue(float dt, float temperature,
                                   const std::vector<float>& human_infectiousness);
            void Update_Infected_Queue(float dt, float temperature);
            void Update_Infectious_Queue(float dt, float temperature);

            // =====================================================
            // Helpers
            // =====================================================

            void InitializeEquilibriumPopulation();
            float GetTotalInQueue(const std::vector<VectorCohort>& queue) const;
        };

    }  // namespace vector

}  // namespace emodlib
