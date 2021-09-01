/**
 * @file SusceptibilityMalaria.cpp
 *
 * @brief Malaria susceptibility implementation
 */

#include "SusceptibilityMalaria.h"

#include "IntrahostComponent.h"  // for static params


namespace emodlib
{

    namespace malaria
    {

        SusceptibilityMalaria::SusceptibilityMalaria()
            : fever_temperature(37.0f)
        {

        }

        void SusceptibilityMalaria::Update()
        {
            // DUMMY LOGIC
            fever_temperature += IntrahostComponent::increment_fever;
        }

        float SusceptibilityMalaria::GetFeverTemperature()
        {
            return fever_temperature;
        }

    }

}
