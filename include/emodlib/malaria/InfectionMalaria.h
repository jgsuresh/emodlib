/**
 * @file InfectionMalaria.h
 *
 * @brief Malaria infection interface
 */

#pragma once

namespace emodlib
{

    namespace malaria
    {

        class InfectionMalaria
        {

        public:
            InfectionMalaria();
            virtual ~InfectionMalaria() {}

            void Update();

            float GetParasiteDensity();
            float GetGametocyteDensity();

        private:

            float parasite_density;
            float gametocyte_density;
        };

    }

}
