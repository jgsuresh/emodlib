/**
 * @file SusceptibilityMalaria.h
 *
 * @brief Malaria susceptibility interface
 */

#pragma once

#include "emodlib/ParamSet.h"

namespace emodlib
{

    namespace malaria
    {

        class Susceptibility
        {

        public:
            struct params
            {
                static float increment_fever;
                static void Configure(const ParamSet& pset);
            };
            
            Susceptibility();

            void Update();

            float GetFeverTemperature() const;

        private:

            float fever_temperature;
        };

    }

}
