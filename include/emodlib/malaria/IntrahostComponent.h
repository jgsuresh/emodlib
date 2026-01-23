/**
 * @file IntrahostComponent.h
 *
 * @brief Malaria intrahost component interface
 */

#pragma once

#include <list>
#include <memory>

#include "MalariaConfig.h"
#include "InfectionMalaria.h"
#include "SusceptibilityMalaria.h"


namespace emodlib
{

    namespace malaria
    {

        class IntrahostComponent
        {

        public:

            static IntrahostComponent* Create(std::shared_ptr<MalariaConfig> config);

            void Update(float dt);

            void Challenge();
            void ChallengeWithAntigens(const DebugAntigens& antigens);  // Debug: explicit antigens
            bool ChallengeWithSporozoites(int n_sporozoites);
            bool ChallengeWithBites(int n_bites = 1);
            void Treat();

            int GetNumInfections() const;

            float GetParasiteDensity() const;
            float GetGametocyteDensity() const;
            float GetFeverTemperature() const;

            float GetInfectiousness() const;

            Susceptibility* GetSusceptibility() const;
            std::list<Infection*> GetInfections() const;

            std::shared_ptr<MalariaConfig> GetConfig() const;

            // Aggregate infection state getters (for validation without raw pointer access)
            std::vector<int64_t> get_total_irbc_counts() const;
            int64_t get_total_hepatocytes() const;
            std::vector<int64_t> get_total_male_gametocytes() const;
            std::vector<int64_t> get_total_female_gametocytes() const;
            int32_t get_total_asexual_cycles() const;

        private:

            std::shared_ptr<MalariaConfig> m_config;
            Susceptibility* susceptibility;
            std::list<Infection*> infections;


            IntrahostComponent(std::shared_ptr<MalariaConfig> config);

        };

    }

}
