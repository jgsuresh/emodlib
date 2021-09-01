/**
 * @file InfectionMalaria.h
 *
 * @brief Malaria infection interface
 */

#pragma once

#include "emodlib/ParamSet.h"

namespace emodlib
{

    namespace malaria
    {

        class Infection
        {

        public:
            struct params
            {
                static float increment_parasite;
                static float increment_gametocyte;
                static void Configure(const ParamSet& pset);
            };
            
            Infection();

            void Update();

            float GetParasiteDensity() const;
            float GetGametocyteDensity() const;

        private:

            float parasite_density;
            float gametocyte_density;
        };

    }

}
