/**
 * @file InfectionMalaria.cpp
 *
 * @brief Malaria infection implementation
 */

#include "InfectionMalaria.h"

namespace emodlib
{

    namespace malaria
    {

        InfectionMalaria::InfectionMalaria()
            : parasite_density(0)
            , gametocyte_density(0)
        {

        }

        void InfectionMalaria::Update()
        {
            // DUMMY LOGIC
            parasite_density += 1.0f;
            gametocyte_density += 1.0f;
        }

        float InfectionMalaria::GetParasiteDensity()
        {
            return parasite_density;
        }

        float InfectionMalaria::GetGametocyteDensity()
        {
            return gametocyte_density;
        }

    }

}
