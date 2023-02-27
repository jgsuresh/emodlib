/**
 * @file IntrahostComponent.h
 *
 * @brief Malaria intrahost component interface
 */

#pragma once

#include <list>
#include <memory>

#include "emodlib/ParamSet.h"
#include "emodlib/utils/RANDOM.h"

#include "InfectionMalaria.h"
#include "SusceptibilityMalaria.h"


namespace emodlib
{

    namespace malaria
    {

        class IntrahostComponent
        {

        public:

            struct params
            {
                static int randomSeed;

                static int max_ind_inf;

                static int falciparumMSPVars;
                static int falciparumNonSpecTypes;
                static int falciparumPfEMP1Vars;

                static void Configure(const ParamSet& pset);
            };


            static std::shared_ptr<RANDOMBASE> p_rng;

            static IntrahostComponent* Create();

            void Update(float dt);

            void Challenge();
            void Treat();

            int GetNumInfections() const;

            float GetParasiteDensity() const;
            float GetGametocyteDensity() const;
            float GetFeverTemperature() const;

            Susceptibility* GetSusceptibility() const;
            std::list<Infection*> GetInfections() const;

        private:

            Susceptibility* susceptibility;
            std::list<Infection*> infections;


            IntrahostComponent();

        };

    }

}
