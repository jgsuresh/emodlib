/**
 * @file InfectionMalaria.cpp
 *
 * @brief Malaria infection implementation
 */

#include "InfectionMalaria.h"

#include "IntrahostComponent.h"  // for static params


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
            parasite_density += IntrahostComponent::increment_parasite;
            gametocyte_density += IntrahostComponent::increment_gametocyte;
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
