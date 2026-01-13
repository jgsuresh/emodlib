/**
 * @file TransmissionSimulation.cpp
 *
 * @brief Implementation of coordinated human-vector transmission simulation.
 */

#include "emodlib/vector/TransmissionSimulation.h"

#include <algorithm>
#include <cmath>

namespace emodlib
{
    namespace vector
    {
        // =====================================================
        // Construction and initialization
        // =====================================================

        TransmissionSimulation::TransmissionSimulation()
            : m_current_day(0)
        {
        }

        void TransmissionSimulation::Initialize(
            std::shared_ptr<malaria::MalariaConfig> malaria_config,
            const VectorConfig& vector_config,
            int n_humans,
            uint32_t seed)
        {
            m_malaria_config = malaria_config;
            m_vector_config = vector_config;
            m_rng.seed(seed);
            m_current_day = 0;

            // Create human intrahost models
            m_humans.clear();
            m_humans.reserve(n_humans);
            for (int i = 0; i < n_humans; ++i)
            {
                m_humans.push_back(malaria::IntrahostComponent::Create(m_malaria_config));
            }

            // Initialize vector population
            m_vectors.Initialize(m_vector_config, seed + 1000);
        }

        // =====================================================
        // Main update
        // =====================================================

        void TransmissionSimulation::Update(float dt, float temperature)
        {
            // 1. Update all human intrahost models
            for (auto& human : m_humans)
            {
                human.Update(dt);
            }

            // 2. Collect human infectiousness values
            std::vector<float> human_infectiousness = CollectHumanInfectiousness();

            // 3. Update vector population with human infectiousness
            m_vectors.Update(dt, temperature, human_infectiousness);

            // 4. Get sporozoite challenges from vectors
            std::vector<float> challenges = m_vectors.GetSporozoiteChallenges(
                static_cast<int>(m_humans.size()));

            // 5. Apply challenges to humans
            ApplySporozoiteChallenges(challenges);

            m_current_day++;
        }

        // =====================================================
        // Intervention interface
        // =====================================================

        void TransmissionSimulation::SeedInfections(int n_infections)
        {
            // Randomly select humans to infect
            std::vector<int> indices(m_humans.size());
            for (size_t i = 0; i < indices.size(); ++i)
            {
                indices[i] = static_cast<int>(i);
            }

            std::shuffle(indices.begin(), indices.end(), m_rng);

            int to_infect = std::min(n_infections, static_cast<int>(m_humans.size()));
            for (int i = 0; i < to_infect; ++i)
            {
                m_humans[indices[i]].Challenge();
            }
        }

        void TransmissionSimulation::ChallengeHuman(int human_index)
        {
            if (human_index >= 0 && human_index < static_cast<int>(m_humans.size()))
            {
                m_humans[human_index].Challenge();
            }
        }

        // =====================================================
        // Output queries - Human population
        // =====================================================

        float TransmissionSimulation::GetPrevalence() const
        {
            if (m_humans.empty())
            {
                return 0.0f;
            }

            int infected = GetInfectedCount();
            return static_cast<float>(infected) / static_cast<float>(m_humans.size());
        }

        float TransmissionSimulation::GetMeanParasiteDensity() const
        {
            if (m_humans.empty())
            {
                return 0.0f;
            }

            float total = 0.0f;
            for (const auto& human : m_humans)
            {
                total += human.GetParasiteDensity();
            }
            return total / static_cast<float>(m_humans.size());
        }

        float TransmissionSimulation::GetMeanGametocyteDensity() const
        {
            if (m_humans.empty())
            {
                return 0.0f;
            }

            float total = 0.0f;
            for (const auto& human : m_humans)
            {
                total += human.GetGametocyteDensity();
            }
            return total / static_cast<float>(m_humans.size());
        }

        float TransmissionSimulation::GetMeanInfectiousness() const
        {
            if (m_humans.empty())
            {
                return 0.0f;
            }

            float total = 0.0f;
            for (const auto& human : m_humans)
            {
                total += human.GetInfectiousness();
            }
            return total / static_cast<float>(m_humans.size());
        }

        int TransmissionSimulation::GetInfectedCount() const
        {
            int count = 0;
            for (const auto& human : m_humans)
            {
                if (human.GetNumInfections() > 0)
                {
                    count++;
                }
            }
            return count;
        }

        // =====================================================
        // Output queries - Vector population
        // =====================================================

        float TransmissionSimulation::GetVectorPopulation() const
        {
            return m_vectors.GetAdultPopulation();
        }

        float TransmissionSimulation::GetVectorInfectiousFraction() const
        {
            return m_vectors.GetInfectiousFraction();
        }

        float TransmissionSimulation::GetVectorInfectedFraction() const
        {
            return m_vectors.GetInfectedFraction();
        }

        float TransmissionSimulation::GetVectorSusceptibleCount() const
        {
            return m_vectors.GetSusceptibleCount();
        }

        float TransmissionSimulation::GetVectorExposedCount() const
        {
            return m_vectors.GetExposedCount();
        }

        float TransmissionSimulation::GetVectorInfectiousCount() const
        {
            return m_vectors.GetInfectiousCount();
        }

        // =====================================================
        // Internal helpers
        // =====================================================

        std::vector<float> TransmissionSimulation::CollectHumanInfectiousness() const
        {
            std::vector<float> infectiousness(m_humans.size());
            for (size_t i = 0; i < m_humans.size(); ++i)
            {
                infectiousness[i] = m_humans[i].GetInfectiousness();
            }
            return infectiousness;
        }

        void TransmissionSimulation::ApplySporozoiteChallenges(
            const std::vector<float>& challenges)
        {
            for (size_t i = 0; i < m_humans.size() && i < challenges.size(); ++i)
            {
                // Challenge if sporozoites were delivered
                // Each sporozoite has a chance to establish infection
                if (challenges[i] > 0.0f)
                {
                    // Poisson number of successful hepatocyte invasions
                    // For simplicity, challenge once if any sporozoites arrived
                    // (EMOD uses more complex sporozoite survival calculation)
                    m_humans[i].Challenge();
                }
            }
        }

    }  // namespace vector

}  // namespace emodlib
