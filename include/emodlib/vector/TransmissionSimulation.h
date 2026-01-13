/**
 * @file TransmissionSimulation.h
 *
 * @brief Coordinates human intrahost models with vector population for
 *        complete malaria transmission cycle simulation.
 *
 * This class enables validation of pymalaria's transmission implementation
 * by running the same dynamics in C++ emodlib.
 */

#pragma once

#include <vector>
#include <memory>
#include <random>

#include "emodlib/malaria/MalariaConfig.h"
#include "emodlib/malaria/IntrahostComponent.h"
#include "emodlib/vector/VectorConfig.h"
#include "emodlib/vector/VectorPopulation.h"

namespace emodlib
{
    namespace vector
    {
        /**
         * @brief Coordinates human-vector transmission dynamics.
         *
         * Update cycle each timestep:
         * 1. Update all human intrahost models
         * 2. Collect human infectiousness values
         * 3. Update vector population with human infectiousness
         * 4. Get sporozoite challenges from vectors
         * 5. Apply challenges to humans
         */
        class TransmissionSimulation
        {
        public:
            // =====================================================
            // Construction and initialization
            // =====================================================

            TransmissionSimulation();

            /**
             * @brief Initialize simulation with configurations.
             *
             * @param malaria_config Configuration for human intrahost models
             * @param vector_config Configuration for vector population
             * @param n_humans Number of humans in simulation
             * @param seed Random seed
             */
            void Initialize(
                std::shared_ptr<malaria::MalariaConfig> malaria_config,
                const VectorConfig& vector_config,
                int n_humans,
                uint32_t seed = 42);

            // =====================================================
            // Main update
            // =====================================================

            /**
             * @brief Update simulation for one timestep.
             *
             * @param dt Timestep in days (typically 1.0)
             * @param temperature Temperature in Celsius
             */
            void Update(float dt, float temperature);

            // =====================================================
            // Intervention interface
            // =====================================================

            /**
             * @brief Seed initial infections in humans.
             *
             * @param n_infections Number of humans to infect
             */
            void SeedInfections(int n_infections);

            /**
             * @brief Challenge a specific human with sporozoites.
             *
             * @param human_index Index of human to challenge
             */
            void ChallengeHuman(int human_index);

            // =====================================================
            // Output queries
            // =====================================================

            // Human population metrics
            float GetPrevalence() const;
            float GetMeanParasiteDensity() const;
            float GetMeanGametocyteDensity() const;
            float GetMeanInfectiousness() const;
            int GetInfectedCount() const;

            // Vector population metrics
            float GetVectorPopulation() const;
            float GetVectorInfectiousFraction() const;
            float GetVectorInfectedFraction() const;
            float GetVectorSusceptibleCount() const;
            float GetVectorExposedCount() const;
            float GetVectorInfectiousCount() const;

            // Simulation state
            int GetCurrentDay() const { return m_current_day; }
            int GetNumHumans() const { return static_cast<int>(m_humans.size()); }

            // Access to individual components (for detailed inspection)
            const std::vector<malaria::IntrahostComponent>& GetHumans() const { return m_humans; }
            const VectorPopulation& GetVectors() const { return m_vectors; }

        private:
            // =====================================================
            // Internal state
            // =====================================================

            std::vector<malaria::IntrahostComponent> m_humans;
            VectorPopulation m_vectors;

            std::shared_ptr<malaria::MalariaConfig> m_malaria_config;
            VectorConfig m_vector_config;

            std::mt19937 m_rng;
            int m_current_day;

            // =====================================================
            // Internal helpers
            // =====================================================

            std::vector<float> CollectHumanInfectiousness() const;
            void ApplySporozoiteChallenges(const std::vector<float>& challenges);
        };

    }  // namespace vector

}  // namespace emodlib
