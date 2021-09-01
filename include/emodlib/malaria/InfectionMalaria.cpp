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

        float Infection::params::increment_parasite = 1.0f;
        float Infection::params::increment_gametocyte = 1.0f;
    
        void Infection::params::Configure(const ParamSet& pset)
        {
            increment_parasite = pset["increment_parasite"].cast<float>();
            increment_gametocyte = pset["increment_gametocyte"].cast<float>();
        }
    
        Infection::Infection()
            : parasite_density(0)
            , gametocyte_density(0)
        {

        }

        void Infection::Update()
        {
            // DUMMY LOGIC
            parasite_density += params::increment_parasite;
            gametocyte_density += params::increment_gametocyte;
        }

        float Infection::GetParasiteDensity() const
        {
            return parasite_density;
        }

        float Infection::GetGametocyteDensity() const
        {
            return gametocyte_density;
        }

    }

}
