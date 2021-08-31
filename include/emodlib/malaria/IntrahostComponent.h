/**
 * @file IntrahostComponent.h
 *
 * @brief Malaria intrahost component interface
 */

#pragma once

#include <list>

#include "InfectionMalaria.h"
#include "SusceptibilityMalaria.h"

namespace emodlib
{

    namespace malaria
    {

        class IntrahostComponent
        {

        public:
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
