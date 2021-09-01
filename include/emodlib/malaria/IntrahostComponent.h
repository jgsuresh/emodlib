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
            static float increment_parasite;
            static float increment_gametocyte;
            static float increment_fever;

            static void Configure(const ParamSet& pset);
            
            static IntrahostComponent* Create();
            virtual ~IntrahostComponent() {}

            void Update();

            void Challenge();
            void Treat();

            float GetParasiteDensity();
            float GetGametocyteDensity();
            float GetFeverTemperature();

        private:
            IntrahostComponent();

            SusceptibilityMalaria* susceptibility;
            std::list<InfectionMalaria*> infections;
        };

    }

}
