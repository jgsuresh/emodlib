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
            void Treat();

            int GetNumInfections() const;

            float GetParasiteDensity() const;
            float GetGametocyteDensity() const;
            float GetFeverTemperature() const;

            float GetInfectiousness() const;

            Susceptibility* GetSusceptibility() const;
            std::list<Infection*> GetInfections() const;

            std::shared_ptr<MalariaConfig> GetConfig() const;

        private:

            std::shared_ptr<MalariaConfig> m_config;
            Susceptibility* susceptibility;
            std::list<Infection*> infections;


            IntrahostComponent(std::shared_ptr<MalariaConfig> config);

        };

    }

}
