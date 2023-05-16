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

        float  Susceptibility::params::memory_level                      = 0.2f;
        float  Susceptibility::params::hyperimmune_decay_rate            = 0.0f;
        float  Susceptibility::params::MSP1_antibody_growthrate          = 0.02f;
        float  Susceptibility::params::antibody_stimulation_c50          = 10.0f;
        float  Susceptibility::params::antibody_capacity_growthrate      = 0.1f;
        float  Susceptibility::params::minimum_adapted_response          = 0.02f;
        float  Susceptibility::params::non_specific_growth               = 0.5f;
        float  Susceptibility::params::antibody_csp_decay_days           = DEFAULT_ANTIBODY_CSP_DECAY_DAYS;

        // TODO: emodlib#9 (maternal antibody init) + emodlib#8 (boost + enums)
        // bool   Susceptibility::params::enable_maternal_antibodies_transmission  = false;
        // MaternalAntibodiesType::Enum Susceptibility::params::maternal_antibodies_type = MaternalAntibodiesType::OFF;
        // float  Susceptibility::params::maternal_antibody_protection      = 0.1f;
        float  Susceptibility::params::maternal_antibody_decay_rate      = 0.01f;

        // TODO: emodlib#9 (innate heterogeneity init) + emodlib#8 (boost + enums)
        // InnateImmuneVariationType::Enum Susceptibility::params::innate_immune_variation_type = InnateImmuneVariationType::NONE;
        float  Susceptibility::params::pyrogenic_threshold               = 1000.0f;
        float  Susceptibility::params::fever_IRBC_killrate               = DEFAULT_FEVER_IRBC_KILL_RATE;

        float  Susceptibility::params::erythropoiesis_anemia_effect      = 3.5f;


        void Susceptibility::params::Configure(const ParamSet& pset)
        {
            memory_level = pset["Antibody_Memory_Level"].cast<float>();
            hyperimmune_decay_rate = -log((0.4f - memory_level) / (1.0f - memory_level)) / 120.0f;  // This sets the decay rate towards memory level so that the decay from antibody levels of 1 to levels of 0.4 is consistent
            MSP1_antibody_growthrate = pset["Max_MSP1_Antibody_Growthrate"].cast<float>();
            antibody_stimulation_c50 = pset["Antibody_Stimulation_C50"].cast<float>();
            antibody_capacity_growthrate = pset["Antibody_Capacity_Growth_Rate"].cast<float>();
            minimum_adapted_response = pset["Min_Adapted_Response"].cast<float>();
            non_specific_growth = pset["Nonspecific_Antibody_Growth_Rate_Factor"].cast<float>();
            antibody_csp_decay_days = pset["Antibody_CSP_Decay_Days"].cast<float>();

            maternal_antibody_decay_rate = pset["Maternal_Antibody_Decay_Rate"].cast<float>();

            pyrogenic_threshold = pset["Pyrogenic_Threshold"].cast<float>();
            fever_IRBC_killrate = pset["Fever_IRBC_Kill_Rate"].cast<float>();

            erythropoiesis_anemia_effect = pset["Erythropoiesis_Anemia_Effect"].cast<float>();
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
            , m_ind_pyrogenic_threshold(0.0f)
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
            age = 20 * DAYSPERYEAR;  // TODO: emodlib#10 (demographic components)

            // TODO: emodlib#10 (transmission components)
            // m_age_dependent_biting_risk = BitingRiskAgeFactor(age);

            recalculateBloodCapacity(age);
            m_RBC = m_RBCcapacity;

            // Track individual pyrogenic thresholds + fever killing rates as instance variables
            // TODO: emodlib#9 (innate heterogeneity init)
            m_ind_pyrogenic_threshold = params::pyrogenic_threshold;
            m_ind_fever_kill_rate = params::fever_IRBC_killrate;

            // TODO: emodlib#9 (maternal antibody init)

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

        void Susceptibility::UpdateActiveAntibody( pfemp1_antibody_t &pfemp1_variant, int minor_variant, int major_variant )
        {
            if(pfemp1_variant.minor == nullptr)
            {
                pfemp1_variant.minor = RegisterAntibody(MalariaAntibodyType::PfEMP1_minor, minor_variant);
            }

            if(pfemp1_variant.major == nullptr)
            {
                pfemp1_variant.major = RegisterAntibody(MalariaAntibodyType::PfEMP1_major, major_variant);
            }
        }

        void Susceptibility::remove_RBCs(int64_t infectedAsexual, int64_t infectedGametocytes, double RBC_destruction_multiplier)
        {
            m_RBC -= ( int64_t(infectedAsexual*RBC_destruction_multiplier) + infectedGametocytes );
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
            if ( !m_CSP_antibody->GetAntigenicPresence() )
            {
                m_CSP_antibody->Decay( dt );
                return;
            }

            // Hyper-immune response (could potentially keep this as part of the update in ExposeToInfectivity)
            if (m_CSP_antibody->GetAntibodyCapacity() > 0.4)
            {
                m_CSP_antibody->UpdateAntibodyCapacityByRate( dt, 0.33f );
            }

            m_CSP_antibody->UpdateAntibodyConcentration( dt );
        }

        void Susceptibility::updateImmunityMSP( float dt, float& temp_cytokine_stimulation )
        {
            // Merozoite-specific immunity
            // Blackman, M. J., H. G. Heidrich, et al. (1990).
            // "A single fragment of a malaria merozoite surface protein remains on the parasite during
            //  red cell invasion and is the target of invasion-inhibiting antibodies."
            // J Exp Med 172(1): 379-382.

            for (auto antibody : m_active_MSP_antibodies)
            {
                if ( !antibody->GetAntigenicPresence() )
                {
                    antibody->Decay( dt );
                    continue;
                }

                // Temporary cytokines stimulated by spikes in MSP antigenic presence after schizont bursts
                temp_cytokine_stimulation += antibody->StimulateCytokines( dt, m_inv_microliters_blood );

                antibody->UpdateAntibodyCapacity( dt, m_inv_microliters_blood );
                antibody->UpdateAntibodyConcentration( dt );
            }
        }

        void Susceptibility::updateImmunityPfEMP1Minor( float dt )
        {
            // Minor epitope IRBC antigens
            // Recker, M., S. Nee, et al. (2004).
            // "Transient cross-reactive immune responses can orchestrate antigenic variation in malaria."
            // Nature 429(6991): 555-558.

            for (auto antibody : m_active_PfEMP1_minor_antibodies)
            {
                if ( !antibody->GetAntigenicPresence() )
                {
                    antibody->Decay( dt );
                    continue;
                }

                antibody->UpdateAntibodyCapacity( dt, m_inv_microliters_blood );
                antibody->UpdateAntibodyConcentration( dt );

                // Accumulate parasite density
                m_parasite_density += float(antibody->GetAntigenCount()) * m_inv_microliters_blood;
            }
        }

        void Susceptibility::updateImmunityPfEMP1Major( float dt )
        {
            for (auto antibody : m_active_PfEMP1_major_antibodies)
            {
                if ( !antibody->GetAntigenicPresence() )
                {
                    antibody->Decay( dt );
                    continue;
                }

                // Cytokines released at low antibody concentration (if capacity hasn't switched into high proliferation rate yet)
                if ( antibody->GetAntibodyCapacity() <= 0.4 )
                {
                    m_cytokine_stimulation += antibody->StimulateCytokines( dt, m_inv_microliters_blood );
                }

                antibody->UpdateAntibodyCapacity( dt, m_inv_microliters_blood );
                antibody->UpdateAntibodyConcentration( dt );
            }
        }

        void Susceptibility::decayAllAntibodies( float dt )
        {
            // CSP handled outside check for any active infection

            for (auto antibody : m_active_MSP_antibodies)
            {
                antibody->Decay( dt );
            }

            for (auto antibody : m_active_PfEMP1_minor_antibodies)
            {
                antibody->Decay( dt );
            }

            for (auto antibody : m_active_PfEMP1_major_antibodies)
            {
                antibody->Decay( dt );
            }
        }

        void Susceptibility::SetAntigenPresent()
        {
            m_antigenic_flag = 1;
        }

        long long Susceptibility::get_RBC_count() const
        {
            return m_RBC;
        }

        float Susceptibility::get_inv_microliters_blood() const
        {
            return m_inv_microliters_blood;
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

        float Susceptibility::get_cytokines() const
        {
            return m_cytokines;
        }

        float Susceptibility::get_fever_killing_rate() const
        {
            return m_ind_fever_kill_rate;
        }

        float Susceptibility::get_parasite_density() const
        {
            return m_parasite_density;
        }

        float Susceptibility::get_maternal_antibodies() const
        {
            return m_maternal_antibody_strength;
        }

        float Susceptibility::get_age() const
        {
            return age;
        }

        void Susceptibility::set_age(float _age)
        {
            age = _age;
        }

        float Susceptibility::get_maternal_antibody_strength() const
        {
            return m_maternal_antibody_strength;
        }

        void Susceptibility::set_maternal_antibody_strength(float _matAb)
        {
            m_maternal_antibody_strength = _matAb;
        }

        float Susceptibility::get_pyrogenic_threshold() const
        {
            return m_ind_pyrogenic_threshold;
        }

        void Susceptibility::set_pyrogenic_threshold(float _threshold)
        {
            m_ind_pyrogenic_threshold = _threshold;
        }

        float Susceptibility::get_fever_kill_rate() const
        {
            return m_ind_fever_kill_rate;
        }

        void Susceptibility::set_fever_kill_rate(float _rate)
        {
            m_ind_fever_kill_rate = _rate;
        }
    }

}
