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

        SusceptibilityMalaria::SusceptibilityMalaria()
            : fever_temperature(37.0f)
        {

        }

        void SusceptibilityMalaria::Update()
        {
            // DUMMY LOGIC
            fever_temperature += 1.0f;
        }

        float SusceptibilityMalaria::GetFeverTemperature()
        {
            return fever_temperature;
        }

    }

}
