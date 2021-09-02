/**
 * @file SusceptibilityMalaria.cpp
 *
 * @brief Malaria susceptibility implementation
 */

#include "SusceptibilityMalaria.h"

#include <iostream>

#include "emodlib/utils/Common.h"
#include "emodlib/utils/Sigmoid.h"
#include "Malaria.h"

#include "MalariaAntibody.h"


namespace emodlib
{

    namespace malaria
    {
    
        float  Susceptibility::params::memory_level                      = 0.0f;
        float  Susceptibility::params::hyperimmune_decay_rate            = 0.0f;
        float  Susceptibility::params::MSP1_antibody_growthrate          = 0.0f;
        float  Susceptibility::params::antibody_stimulation_c50          = 0.0f;
        float  Susceptibility::params::antibody_capacity_growthrate      = 0.0f;
        float  Susceptibility::params::minimum_adapted_response          = 0.0f;
        float  Susceptibility::params::non_specific_growth               = 0.0f;
        float  Susceptibility::params::antibody_csp_decay_days           = 0.0f;
    
        bool   Susceptibility::params::enable_maternal_antibodies_transmission  = false;
        MaternalAntibodiesType::Enum Susceptibility::params::maternal_antibodies_type = MaternalAntibodiesType::OFF;
        float  Susceptibility::params::maternal_antibody_protection      = 0.0f;
        float  Susceptibility::params::maternal_antibody_decay_rate      = 0.0f;
        InnateImmuneVariationType::Enum Susceptibility::params::innate_immune_variation_type = InnateImmuneVariationType::NONE;
        float  Susceptibility::params::base_gametocyte_mosquito_survival = 1.0f;
        float  Susceptibility::params::cytokine_gametocyte_inactivation  = 1.0f;
        float  Susceptibility::params::erythropoiesis_anemia_effect      = 0.0f;
        float  Susceptibility::params::pyrogenic_threshold               = 0.0f;
        float  Susceptibility::params::fever_IRBC_killrate               = 0.0f;
    
        void Susceptibility::params::Configure(const ParamSet& pset)
        {

        }
    
    
        Susceptibility::Susceptibility()
            : age(0)
    
            , m_antigenic_flag(0)
            , m_maternal_antibody_strength(0)
            , m_CSP_antibody(nullptr)  // IMalariaAntibody* and vector<IMalariaAntibody*> assigned in Initialize()
            , m_active_MSP_antibodies()
            , m_active_PfEMP1_minor_antibodies()
            , m_active_PfEMP1_major_antibodies()
    
            , m_RBC(0)
            , m_RBCcapacity(0)
            , m_RBCproduction(0)
            , m_inv_microliters_blood(0.0f)  // assigned in Initialize() as function of age
    
            , m_cytokines(0.0f)
            , m_ind_pyrogenic_threshold(0.0f)  // TODO: m_ind_* include individual heterogeneity
            , m_ind_fever_kill_rate(0.0f)
            , m_cytokine_stimulation(0.0f)
            , m_parasite_density(0.0f)
        {

        }

        Susceptibility* Susceptibility::Create()
        {
            Susceptibility *newsusceptibility = new Susceptibility();
            newsusceptibility->Initialize();
            
            return newsusceptibility;
        }
    
        void Susceptibility::Initialize()
        {
            age = 20 * DAYSPERYEAR;  // TODO: update + access from elsewhere; similarly age-dependent biting risk
            
            recalculateBloodCapacity(age);
            m_RBC = m_RBCcapacity;

            // Track individual pyrogenic thresholds + fever killing rates as instance variables
            m_ind_pyrogenic_threshold = params::pyrogenic_threshold;
            m_ind_fever_kill_rate = params::fever_IRBC_killrate;
            
            m_CSP_antibody = MalariaAntibodyCSP::CreateAntibody(0);

            // MSP + PfEMP1 antibodies are added upon infection
        }
    
        IMalariaAntibody* Susceptibility::RegisterAntibody(MalariaAntibodyType::Enum type, int variant, float capacity)
        {
            std::vector<IMalariaAntibody*> *variant_vector;
            IMalariaAntibody* (*typed_create_antibody)(int,float);

            switch( type )
            {
            case MalariaAntibodyType::CSP:
                return m_CSP_antibody; // only one CSP variant, so ignore second argument for now.

            case MalariaAntibodyType::MSP1:
                variant_vector = &m_active_MSP_antibodies;
                typed_create_antibody = MalariaAntibodyMSP::CreateAntibody;
                break;

            case MalariaAntibodyType::PfEMP1_minor:
                variant_vector = &m_active_PfEMP1_minor_antibodies;
                typed_create_antibody = MalariaAntibodyPfEMP1Minor::CreateAntibody;
                break;

            case MalariaAntibodyType::PfEMP1_major:
                variant_vector = &m_active_PfEMP1_major_antibodies;
                typed_create_antibody = MalariaAntibodyPfEMP1Major::CreateAntibody;
                break;

            default:
                std::cout << "Unknown MalariaAntibodyType enum used";
                throw;
            }
            
            IMalariaAntibody* antibody = nullptr;
            for (auto tmp_antibody : *variant_vector)
            {
                if ( tmp_antibody->GetAntibodyVariant() == variant )
                {
                    antibody = tmp_antibody;
                    break;
                }
            }

            if (antibody == nullptr) // make a new antibody if it hasn't been created yet
            {
                antibody = typed_create_antibody(variant, capacity);
                variant_vector->push_back(antibody);
            }

            return antibody;
        }
    
        void Susceptibility::Update(float dt)
        {
            age += dt;
            
            recalculateBloodCapacity(age);
            
            // Red blood cell dynamics
            if (Susceptibility::params::erythropoiesis_anemia_effect > 0)
            {
                // This is the amount of "erythropoietin", assume absolute amounts of erythropoietin correlate linearly with absolute increases in hemoglobin
                float anemia_erythropoiesis_multiplier = exp( Susceptibility::params::erythropoiesis_anemia_effect * (1 - get_RBC_availability()) );
                m_RBC = int64_t(m_RBC - (m_RBC * .00833 - m_RBCproduction * anemia_erythropoiesis_multiplier) * dt); // *.00833 ==/120 (AVERAGE_RBC_LIFESPAN)
            }
            else
            {
                m_RBC = int64_t(m_RBC - (m_RBC * .00833 - m_RBCproduction) * dt); // *.00833 ==/120 (AVERAGE_RBC_LIFESPAN)
            }
            
            // Cytokines decay with time constant of 12 hours
            m_cytokines -= (m_cytokines * 2 * dt);
            if (m_cytokines < 0) { m_cytokines = 0; }
            
            // Reset parasite density
            m_parasite_density = 0; // this is accumulated in updateImmunityPfEMP1Minor

            // decay maternal antibodies
            m_maternal_antibody_strength -= dt * m_maternal_antibody_strength * Susceptibility::params::maternal_antibody_decay_rate;
            if ( m_maternal_antibody_strength < 0 ) { m_maternal_antibody_strength = 0; }

            // antibody capacities increase and antibodies released if antigen present, only process if antigens are present at all
            // concept of antibody stimulation threshold seen in other models--(Molineaux, Diebner et al. 2001; Paget-McNicol, Gatton et al. 2002; Dietz, Raddatz et al. 2006)
            // first CSP, then rest (but only process rest if there is an active infection), have to process CSP every time step
            updateImmunityCSP(dt);
            
            // now all other antigens
            if ( !m_antigenic_flag )
            {
                // NO ANTIGENS.  All antibodies decay to zero and all antibody_capacities decay towards 0.3
                decayAllAntibodies(dt);
            }
            else
            {
                // Update antigen-antibody reactions for MSP and PfEMP1 minor/major epitopes, including cytokine stimulation
                float temp_cytokine_stimulation = 0; // used to track total stimulation of cytokines due to rupturing schizonts
                updateImmunityMSP(dt, temp_cytokine_stimulation);
                updateImmunityPfEMP1Minor(dt);
                updateImmunityPfEMP1Major(dt);
                
                // inflammatory immune response--Stevenson, M. M. and E. M. Riley (2004). "Innate immunity to malaria." Nat Rev Immunol 4(3): 169-180.
                // now let cytokine be increased in response to IRBCs and ruptured schizonts, if any
                // pyrogenic threshold similar to previous models--(Molineaux, Diebner et al. 2001; Paget-McNicol, Gatton et al. 2002; Maire, Smith et al. 2006)
                m_cytokines = float(m_cytokines + CYTOKINE_STIMULATION_SCALE * Sigmoid::basic_sigmoid(m_ind_pyrogenic_threshold, m_cytokine_stimulation) * dt * 2);//12-hour time constant
                m_cytokines = float(m_cytokines + CYTOKINE_STIMULATION_SCALE * Sigmoid::basic_sigmoid(m_ind_pyrogenic_threshold, temp_cytokine_stimulation));//one time spike for rupturing schizonts
                m_cytokine_stimulation = 0; // and reset for next time step

                // reset antigenic presence and IRBC counters
                m_antigenic_flag = 0;
                for (auto antibody : m_active_MSP_antibodies)
                {
                    antibody->ResetCounters();
                }

                for (auto antibody : m_active_PfEMP1_minor_antibodies)
                {
                    antibody->ResetCounters();
                }

                for (auto antibody : m_active_PfEMP1_major_antibodies)
                {
                    antibody->ResetCounters();
                }
            }
        }

        void Susceptibility::recalculateBloodCapacity( float _age )
        {
            // How many RBCs a person should have determined by age.
            // This sets the daily production of red blood cells for adults to maintain
            // standard equilibrium RBC concentrations given RBC lifetime
            if ( _age > (20 * DAYSPERYEAR) )
            {
                // 2.0*10^11 (RBCs/day)*(120 days)=2.4x10^13 RBCs ~= 5 liters * 5x10^6 RBCs/microliter
                m_RBCproduction         = ADULT_RBC_PRODUCTION;
                m_inv_microliters_blood = float(1 / ( (0.225 * (7300/DAYSPERYEAR) + 0.5) * 1e6 ));
            }
            else
            {
                // Sets daily production of red blood cells for children to set their equilibrium RBC concentrations and blood volume given an RBC lifetime
                // Only approximate due to linear increase in blood volume from 0.5 to 5 liters with age, a better growth model would be nonlinear
                m_RBCproduction         = int64_t(INFANT_RBC_PRODUCTION + (_age * .000137) * (ADULT_RBC_PRODUCTION - INFANT_RBC_PRODUCTION)); //*.000137==/(20*DAYSPERYEAR)
                m_inv_microliters_blood = float(1 / ( (0.225 * (_age/DAYSPERYEAR) + 0.5 ) * 1e6 ));
            }
            
            m_RBCcapacity = m_RBCproduction * AVERAGE_RBC_LIFESPAN;  // Health equilibrium of RBC is production*lifetime.  This is the total number of RBC per human
        }
    
        void Susceptibility::updateImmunityCSP( float dt )
        {
            
        }
        
        void Susceptibility::updateImmunityMSP( float dt, float& temp_cytokine_stimulation )
        {
            
        }
        
        void Susceptibility::updateImmunityPfEMP1Minor( float dt )
        {
            
        }
        
        void Susceptibility::updateImmunityPfEMP1Major( float dt )
        {
            
        }
        
        void Susceptibility::decayAllAntibodies( float dt )
        {
            
        }
    
        void Susceptibility::SetAntigenPresent()
        {
            m_antigenic_flag = 1;
        }
    
        long long Susceptibility::get_RBC_count() const
        {
            return m_RBC;
        }

        double Susceptibility::get_RBC_availability() const
        {
            return (m_RBCcapacity > 0) ? (double(m_RBC) / m_RBCcapacity) : 0.0;
        }
    
        // Fever tracks the level of cytokines
        // This changes a limited cytokine range to more closely match the range of fevers experienced by patients
        float Susceptibility::get_fever() const
        {
            return FEVER_DEGREES_CELSIUS_PER_UNIT_CYTOKINES * m_cytokines;
        }
    
        float Susceptibility::get_fever_celsius() const
        {
            return 37.0f + get_fever();
        }

    }

}
