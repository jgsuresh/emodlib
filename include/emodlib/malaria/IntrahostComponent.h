/**
 * @file IntrahostComponent.h
 *
 * @brief Malaria intrahost component interface
 */

#pragma once

#include <list>

#include "emodlib/ParamSet.h"

#include "InfectionMalaria.h"
#include "SusceptibilityMalaria.h"


namespace emodlib
{

    namespace malaria
    {
    
        class IntrahostComponent
        {

        public:
            struct params
            {
                static void Configure(const ParamSet& pset);
            };
            
            static IntrahostComponent* Create();

            void Update();

            void Challenge();
            void Treat();

            float GetParasiteDensity() const;
            float GetGametocyteDensity() const;
            float GetFeverTemperature() const;

        private:
            IntrahostComponent();

            Susceptibility* susceptibility;
            std::list<Infection*> infections;
        };

    }

}
