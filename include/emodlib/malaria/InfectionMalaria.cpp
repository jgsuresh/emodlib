/**
 * @file InfectionMalaria.cpp
 *
 * @brief Malaria infection implementation
 */

#include "InfectionMalaria.h"

#include "IntrahostComponent.h"
#include "SusceptibilityMalaria.h"


namespace emodlib
{

    namespace malaria
    {

        ParasiteSwitchType::Enum Infection::params::parasite_switch_type = ParasiteSwitchType::RATE_PER_PARASITE_7VARS;
//        MalariaStrains::Enum     Infection::params::malaria_strains = MalariaStrains::FALCIPARUM_RANDOM_STRAIN;

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
    
    
        suids::distributed_generator Infection::infectionSuidGenerator(0, 0);

    
        void Infection::params::Configure(const ParamSet& pset)
        {

        }
    
    
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
    
            , susceptibility(nullptr)
    
            , parasite_density(0)
            , gametocyte_density(0)
        {

        }

        Infection* Infection::Create(Susceptibility* _susceptibility, int initial_hepatocytes)
        {
            Infection *newinfection = new Infection();
            newinfection->Initialize(_susceptibility, initial_hepatocytes);
            
            return newinfection;
        }
    
        void Infection::Initialize(Susceptibility* _susceptibility, int initial_hepatocytes)
        {
            suid = infectionSuidGenerator();  // next suid from generator
            m_hepatocytes = initial_hepatocytes;
            
            // Here we set the antigenic repertoire of the infection
            // Can be completely distinct strains, or partially overlapping repertoires of antigens
            // Bull, P. C., B. S. Lowe, et al. (1998). "Parasite antigens on the infected red cell surface are targets for naturally acquired immunity to malaria." Nat Med 4(3): 358-360.
            // Recker, M., S. Nee, et al. (2004). "Transient cross-reactive immune responses can orchestrate antigenic variation in malaria." Nature 429(6991): 555-558.
            // In our model, not all antigens are expressed at the same time, but switching occurs.  This just sets the total repertoire
                        
            auto rng = IntrahostComponent::p_rng;
            
            m_MSPtype = rng->uniformZeroToN16(IntrahostComponent::params::falciparumMSPVars);
            m_nonspectype = rng->uniformZeroToN16(IntrahostComponent::params::falciparumNonSpecTypes);

            for (int i = 0; i < CLONAL_PfEMP1_VARIANTS; i++)
            {
                m_IRBCtype[i] = rng->uniformZeroToN16(IntrahostComponent::params::falciparumPfEMP1Vars);
                m_minor_epitope_type[i] = rng->uniformZeroToN16(MINOR_EPITOPE_VARS_PER_SET) + MINOR_EPITOPE_VARS_PER_SET * m_nonspectype;
            }
            
            susceptibility = _susceptibility;
            
            m_MSP_antibody = susceptibility->RegisterAntibody(MalariaAntibodyType::MSP1, m_MSPtype);
            
            for( int ivariant = 0; ivariant < m_PfEMP1_antibodies.size(); ivariant++ )
            {
                m_PfEMP1_antibodies[ivariant].major = nullptr;
                m_PfEMP1_antibodies[ivariant].minor = nullptr;

                if ( m_IRBC_count[ivariant] > 0 )
                {
                    m_PfEMP1_antibodies[ivariant].minor  = susceptibility->RegisterAntibody(MalariaAntibodyType::PfEMP1_minor, m_minor_epitope_type[ivariant]);
                    m_PfEMP1_antibodies[ivariant].major  = susceptibility->RegisterAntibody(MalariaAntibodyType::PfEMP1_major, m_IRBCtype[ivariant]);
                }
            }
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
