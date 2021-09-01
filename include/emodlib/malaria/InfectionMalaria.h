/**
 * @file InfectionMalaria.h
 *
 * @brief Malaria infection interface
 */

#pragma once

#include "emodlib/ParamSet.h"
#include "emodlib/utils/suids.hpp"

#include "Malaria.h"
#include "MalariaEnums.h"
#include "IMalariaAntibody.h"


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
            
            
            static suids::distributed_generator infectionSuidGenerator;
            
            
            static Infection *Create(int initial_hepatocytes=1);
            
            void Update();

            suids::suid GetSuid() const;
            float GetParasiteDensity() const;
            float GetGametocyteDensity() const;

        private:

            suids::suid suid; // unique id of this infection within the system
            
            float duration;   // local timer
            float total_duration;
            float incubation_timer;
            float infectious_timer;
            
            double m_IRBCtimer;
            int32_t m_hepatocytes;
            AsexualCycleStatus::Enum m_asexual_phase;
            int32_t m_asexual_cycle_count;
            
            int32_t m_MSPtype;        // allow variation in MSP from clone to clone
            int32_t m_nonspectype;    // what is the set of minor_epitope_types
            int32_t m_minor_epitope_type[CLONAL_PfEMP1_VARIANTS];
            int32_t m_IRBCtype[CLONAL_PfEMP1_VARIANTS];
            
            IMalariaAntibody* m_MSP_antibody;
            std::vector< pfemp1_antibody_t > m_PfEMP1_antibodies;

            std::vector<int64_t> m_IRBC_count;
            int64_t m_malegametocytes[GametocyteStages::Count];
            int64_t m_femalegametocytes[GametocyteStages::Count];

            // govern distribution of next merozoites -- TODO: why are these not in Params?
            double m_gametorate;
            double m_gametosexratio;
            
            double m_inv_microliters_blood;   // tracks blood volume based on age
            
            float parasite_density;
            float gametocyte_density;
            
            
            Infection();
            void Initialize(int initial_hepatocytes);
        };

    }

}
