/**
 * @file InfectionMalaria.h
 *
 * @brief Malaria infection interface
 */

#pragma once

#include "emodlib/ParamSet.h"

#include "MalariaEnums.h"
//#include "MalariaContexts.h"
//#include "IMalariaAntibody.h"

namespace emodlib
{

    namespace malaria
    {

        class Infection
        {

        public:
            struct params
            {
                static ParasiteSwitchType::Enum parasite_switch_type;
                static MalariaStrains::Enum     malaria_strains;

                static float antibody_IRBC_killrate;
                static float MSP1_merozoite_kill;
                static float gametocyte_stage_survival;
                static float base_gametocyte_sexratio;
                static float base_gametocyte_production;
                static float antigen_switch_rate;
                static float merozoites_per_hepatocyte;
                static float merozoites_per_schizont;
                static float non_specific_antigenicity;
                static float RBC_destruction_multiplier;
                static int   n_asexual_cycles_wo_gametocytes;
                
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
