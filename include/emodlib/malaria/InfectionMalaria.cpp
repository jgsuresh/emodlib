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

        ParasiteSwitchType::Enum Infection::params::parasite_switch_type = ParasiteSwitchType::RATE_PER_PARASITE_7VARS;
        MalariaStrains::Enum     Infection::params::malaria_strains = MalariaStrains::FALCIPARUM_RANDOM_STRAIN;

        float Infection::params::antibody_IRBC_killrate = 0.0f;
        float Infection::params::MSP1_merozoite_kill = 0.0f;
        float Infection::params::gametocyte_stage_survival = 0.0f;
        float Infection::params::base_gametocyte_sexratio = 0.0f;
        float Infection::params::base_gametocyte_production = 0.0f;
        float Infection::params::antigen_switch_rate = 0.0f;
        float Infection::params::merozoites_per_hepatocyte = 0.0f;
        float Infection::params::merozoites_per_schizont = 0.0f;
        float Infection::params::non_specific_antigenicity = 0.0f;
        float Infection::params::RBC_destruction_multiplier = 0.0f;
        int   Infection::params::n_asexual_cycles_wo_gametocytes = 0;
    
        void Infection::params::Configure(const ParamSet& pset)
        {

        }
    
    
        suids::distributed_generator Infection::infectionSuidGenerator(0, 0);

    
        Infection::Infection()
            : suid(suids::nil_suid())
    
            , duration(0.0f)
            , total_duration(0.0f)
            , incubation_timer(0.0f)
            , infectious_timer(0.0f)
    
            , m_IRBCtimer(0.0)
            , m_hepatocytes(0)
            , m_asexual_phase(AsexualCycleStatus::NoAsexualCycle)
            , m_asexual_cycle_count(0)
    
            , m_MSPtype(0)
            , m_nonspectype(0)
            , m_minor_epitope_type()
            , m_IRBCtype()
    
            , m_MSP_antibody(nullptr)
            , m_PfEMP1_antibodies(CLONAL_PfEMP1_VARIANTS)
    
            , m_IRBC_count(CLONAL_PfEMP1_VARIANTS)
            , m_malegametocytes()
            , m_femalegametocytes()
    
            , m_gametorate(0.0)
            , m_gametosexratio(0.0)
    
            , m_inv_microliters_blood(INV_MICROLITERS_BLOOD_ADULT)
    
            , parasite_density(0)
            , gametocyte_density(0)
        {

        }

        Infection* Infection::Create(int initial_hepatocytes)
        {
            Infection *newinfection = new Infection();
            newinfection->Initialize(initial_hepatocytes);
            
            return newinfection;
        }
    
        void Infection::Initialize(int initial_hepatocytes)
        {
            suid = infectionSuidGenerator();  // next suid from generator
            m_hepatocytes = initial_hepatocytes;
            
            // TODO: SetParameters() does all the antigenic variation random draws
        }
    
        void Infection::Update()
        {

        }

        suids::suid Infection::GetSuid() const
        {
            return suid;
            
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
