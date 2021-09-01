/**
 * @file SusceptibilityMalaria.cpp
 *
 * @brief Malaria susceptibility implementation
 */

#include "SusceptibilityMalaria.h"


namespace emodlib
{

    namespace malaria
    {
    
        float Susceptibility::params::increment_fever = 1.0f;

        void Susceptibility::params::Configure(const ParamSet& pset)
        {
            increment_fever = pset["increment_fever"].cast<float>();
        }
    
        Susceptibility::Susceptibility()
            : fever_temperature(37.0f)
        {

        }

        void Susceptibility::Update()
        {
            // DUMMY LOGIC
            fever_temperature += params::increment_fever;
        }

        float Susceptibility::GetFeverTemperature() const
        {
            return fever_temperature;
        }

    }

}
