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

        class Susceptibility;
    
        class Infection
        {

        public:
            
            struct params
            {
                // TODO: emodlib#8 (boost + enums)
                // static ParasiteSwitchType::Enum parasite_switch_type;
                // static MalariaStrains::Enum     malaria_strains;

                static float incubation_period;
                static float antibody_IRBC_killrate;
                static float non_specific_antigenicity;
                static float MSP1_merozoite_kill;
                static float gametocyte_stage_survival;
                static float base_gametocyte_sexratio;
                static float base_gametocyte_production;
                static float antigen_switch_rate;
                static float merozoites_per_hepatocyte;
                static float merozoites_per_schizont;
                static float RBC_destruction_multiplier;
                static int   n_asexual_cycles_wo_gametocytes;
                
                static void Configure(const ParamSet& pset);
            };
            
            
            static suids::distributed_generator infectionSuidGenerator;
            
            
            static Infection *Create(Susceptibility* _susceptibility, int initial_hepatocytes=1);
            
            void Update(float dt);

            suids::suid GetSuid() const;
            int64_t get_MaleGametocytes(int stage) const;
            int64_t get_FemaleGametocytes(int stage) const;
            float get_asexual_density() const;
            float get_mature_gametocyte_density() const;
            bool IsCleared() const;

            
        private:

            suids::suid suid; // unique id of this infection within the system
            
            float m_liver_stage_timer;
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

            // placeholders for infection-level variation in merozoite-to-gametocyte dynamics
            double m_gametorate;
            double m_gametosexratio;
                        
            Susceptibility* immunity;
            
            
            Infection();
            void Initialize(Susceptibility* _susceptibility, int initial_hepatocytes);
            
            void malariaProcessHepatocytes(float dt);
            void processEndOfAsexualCycle();
            void malariaIRBCAntigenSwitch(double merozoitesurvival = 1.0);
            void malariaCycleGametocytes(double merozoitesurvival = 1.0);
            void malariaImmuneStimulation(float dt);
            void malariaImmunityIRBCKill(float dt);
            void malariaImmunityGametocyteKill(float dt);
            void malariaCheckInfectionStatus(float dt);  // TODO: emodlib#3 (InfectionStateChange::Cleared)
            void apply_MatureGametocyteKillProbability(float pkill);

        };

    }

}
