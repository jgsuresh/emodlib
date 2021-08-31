/**
 * @file SusceptibilityMalaria.h
 *
 * @brief Malaria susceptibility interface
 */

#pragma once

namespace emodlib
{

    namespace malaria
    {

        class SusceptibilityMalaria
        {

        public:
            SusceptibilityMalaria();
            virtual ~SusceptibilityMalaria() {}

            void Update();

            float GetFeverTemperature();

        private:

            float fever_temperature;
        };

    }

}
