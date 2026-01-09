#include "MalariaAntibody.h"

#include <cmath>
#include <cfloat>

#include "emodlib/utils/Sigmoid.h"
#include "MalariaConfig.h"


#define NON_TRIVIAL_ANTIBODY_THRESHOLD  (0.0000001)
#define TWENTY_DAY_DECAY_CONSTANT       (0.05f)
#define B_CELL_PROLIFERATION_THRESHOLD  (0.4)
#define B_CELL_PROLIFERATION_CONSTANT   (0.33f)
#define ANTIBODY_RELEASE_THRESHOLD      (0.3)
#define ANTIBODY_RELEASE_FACTOR         (4)


namespace emodlib
{

    namespace malaria
    {

        MalariaAntibody::MalariaAntibody(MalariaConfig* config)
            : m_config(config)
            , m_antibody_capacity(0.0f)
            , m_antibody_concentration(0.0f)
            , m_antigen_count(0)
            , m_antibody_type(MalariaAntibodyType::CSP)
            , m_antibody_variant(0)
            , m_active_index(-1)
            , m_time_last_active(-1.0f)
        {
        }

        void MalariaAntibody::Initialize( MalariaAntibodyType::Enum type, int variant, float capacity, float concentration )
        {
            m_antibody_type          = type;
            m_antibody_variant       = variant;
            m_antibody_capacity      = capacity;
            m_antibody_concentration = concentration;
        }

        // EMOD 2.22: Decay now takes total inactive time, not per-timestep dt
        // Called from IncreaseAntigenCount when antibody becomes active
        void MalariaAntibody::Decay( float decay_time )
        {
            // decay_time - time (in days) the antibody has spent being inactive
            // this is only called when antibodies are activated

            // allow the decay of anti-CSP concentrations greater than unity (e.g. after boosting by vaccine)
            if( (m_antibody_type == MalariaAntibodyType::CSP) && (m_antibody_concentration > m_antibody_capacity) )
            {
                m_antibody_concentration -= m_antibody_concentration * decay_time / m_config->antibody_csp_decay_days;
            }
            else
            {
                // otherwise do the normal behavior of decaying antibody concentration based on capacity

                // don't do multiplication and subtraction unless antibody levels non-trivial
                if(m_antibody_concentration > NON_TRIVIAL_ANTIBODY_THRESHOLD)
                {
                    // EMOD 2.22: exponential decay instead of linear
                    m_antibody_concentration = m_antibody_concentration * std::exp( -decay_time * TWENTY_DAY_DECAY_CONSTANT );
                }

                if(m_antibody_capacity > m_config->memory_level)
                {
                    // antibody capacity decays to a medium value (.3) dropping below .4 in ~120 days from 1.0
                    // EMOD 2.22: exponential decay formula
                    m_antibody_capacity = ( m_antibody_capacity - m_config->memory_level )
                        * std::exp( -decay_time * m_config->hyperimmune_decay_rate )
                        + m_config->memory_level;
                } // stays around memory level until antibody_days_to_long_term_decay kicks in

                // --------------------------------------------------------------------------------
                // --- EMOD 2.22: Long-term decay
                // --- If the antibody has been dormant for a long time, start a gradual decay.
                // --- This is to help reduce the issue where older people have too much immunity.
                // --------------------------------------------------------------------------------
                if(( decay_time >= m_config->antibody_days_to_long_term_decay ) &&
                    ( m_antibody_capacity > FLT_EPSILON ))
                {
                    float delta_time = decay_time - m_config->antibody_days_to_long_term_decay;
                    m_antibody_capacity = m_antibody_capacity * std::exp( -delta_time / m_config->antibody_long_term_decay_days );
                }
            }
        }

        float MalariaAntibody::StimulateCytokines( float dt, float inv_uL_blood )
        {
            // Cytokines released at low antibody concentration (if capacity hasn't switched into high proliferation rate yet)
            return ( 1 - m_antibody_concentration ) * float(m_antigen_count) * inv_uL_blood;
        }

        // Let's use the MSP version of antibody growth in the base class ...
        void MalariaAntibody::UpdateAntibodyCapacity( float dt, float inv_uL_blood )
        {
            float growth_rate = m_config->MSP1_antibody_growthrate;
            float threshold   = m_config->antibody_stimulation_c50;

            m_antibody_capacity += growth_rate  * (1.0f - m_antibody_capacity) * float(Sigmoid::basic_sigmoid( threshold, float(m_antigen_count) * inv_uL_blood));

            // rapid B cell proliferation above a threshold given stimulation
            if (m_antibody_capacity > B_CELL_PROLIFERATION_THRESHOLD)
            {
                m_antibody_capacity += ( 1.0f - m_antibody_capacity ) * B_CELL_PROLIFERATION_CONSTANT * dt;
            }

            if (m_antibody_capacity > 1.0)
            {
                m_antibody_capacity = 1.0;
            }
        }

        // Different arguments used by CSP update called directly from IndividualHumanMalaria::ExposeToInfectivity
        // and also in SusceptibilityMalaria::updateImmunityCSP
        void MalariaAntibody::UpdateAntibodyCapacityByRate( float dt, float growth_rate )
        {
            m_antibody_capacity += growth_rate * dt * (1 - m_antibody_capacity);

            if (m_antibody_capacity > 1.0)
            {
                m_antibody_capacity = 1.0;
            }
        }

        // The minor PfEMP1 version is similar but not exactly the same...
        void MalariaAntibodyPfEMP1Minor::UpdateAntibodyCapacity( float dt, float inv_uL_blood )
        {
            float min_stimulation = m_config->antibody_stimulation_c50 * m_config->minimum_adapted_response;
            float growth_rate     = m_config->antibody_capacity_growthrate * m_config->non_specific_growth;
            float threshold       = m_config->antibody_stimulation_c50;

            if (m_antibody_capacity <= B_CELL_PROLIFERATION_THRESHOLD)
            {
                m_antibody_capacity += growth_rate * dt * (1.0f - m_antibody_capacity) * float(Sigmoid::basic_sigmoid(threshold, float(m_antigen_count) * inv_uL_blood + min_stimulation));
            }
            else
            {
                //rapid B cell proliferation above a threshold given stimulation
                m_antibody_capacity += (1.0f - m_antibody_capacity) * B_CELL_PROLIFERATION_CONSTANT * dt;
            }

            if (m_antibody_capacity > 1.0)
            {
                m_antibody_capacity = 1.0;
            }
        }

        // The major PfEMP1 version is slightly different again...
        void MalariaAntibodyPfEMP1Major::UpdateAntibodyCapacity( float dt, float inv_uL_blood )
        {
            float min_stimulation = m_config->antibody_stimulation_c50 * m_config->minimum_adapted_response;
            float growth_rate     = m_config->antibody_capacity_growthrate;
            float threshold       = m_config->antibody_stimulation_c50;

            if (m_antibody_capacity <= B_CELL_PROLIFERATION_THRESHOLD)
            {
                //ability and number of B-cells to produce antibodies, with saturation
                m_antibody_capacity += growth_rate * dt * (1.0f - m_antibody_capacity) * float(Sigmoid::basic_sigmoid(threshold, float(m_antigen_count) * inv_uL_blood + min_stimulation));

                // check for antibody capacity out of range
                if (m_antibody_capacity > 1.0)
                {
                    m_antibody_capacity = 1.0;
                }
            }
            else
            {
                //rapid B cell proliferation above a threshold given stimulation
                m_antibody_capacity += (1.0f - m_antibody_capacity) * B_CELL_PROLIFERATION_CONSTANT * dt;
            }
        }

        void MalariaAntibody::UpdateAntibodyConcentration( float dt )
        {
            // release of antibodies and effect of B cell proliferation on capacity
            // antibodies released after capacity passes 0.3
            // detection and proliferation in lymph nodes, etc...
            // and circulating memory cells
            if ( m_antibody_capacity > ANTIBODY_RELEASE_THRESHOLD )
            {
                m_antibody_concentration += ( m_antibody_capacity - m_antibody_concentration ) * ANTIBODY_RELEASE_FACTOR * dt;
            }

            if ( m_antibody_concentration > m_antibody_capacity )
            {
                m_antibody_concentration = m_antibody_capacity;
            }
        }

        void MalariaAntibodyCSP::UpdateAntibodyConcentration( float dt )
        {
            // allow the decay of anti-CSP concentrations greater than unity (e.g. after boosting by vaccine)
            if ( m_antibody_concentration > m_antibody_capacity )
            {
                m_antibody_concentration -= m_antibody_concentration * dt / m_config->antibody_csp_decay_days;
            }
            else
            {
                // otherwise do the normal behavior of incrementing antibody concentration based on capacity
                MalariaAntibody::UpdateAntibodyConcentration( dt );
            }
        }

        void MalariaAntibody::ResetCounters()
        {
            m_antigen_count = 0;
            m_active_index = -1;
        }

        // EMOD 2.22: IncreaseAntigenCount now handles time-based decay
        void MalariaAntibody::IncreaseAntigenCount( int64_t antigenCount, float currentTime, float dt )
        {
            m_antigen_count += antigenCount;

            // --------------------------------------------------------------------------------------------
            // --- currentTime here is not sim currentTime but currentTime + X*infectious_timestep(i.e. dt)
            // --- so if this was updated the previous infectious_timestep then it was just active
            // --- and should not decay.   The difference between the currentTime and the last time active
            // --- needs to be greater than the infectious_timestep so we know that it spent sometime
            // --- unactive and needs to be decayed.
            // --- We subtract the dt because we don't decay for the current time step.
            // --------------------------------------------------------------------------------------------
            float decay_time = currentTime - m_time_last_active - dt;
            if( (m_time_last_active >= 0.0f) && (decay_time > 0.0f) )
            {
                Decay( decay_time );
            }
            m_time_last_active = currentTime;
        }

        void MalariaAntibody::SetAntigenicPresence( bool antigenPresent )
        {
            m_antigen_count = (antigenPresent ? 1 : 0);
        }

        int64_t MalariaAntibody::GetAntigenCount() const
        {
            return m_antigen_count;
        }

        bool MalariaAntibody::GetAntigenicPresence() const
        {
            return (m_antigen_count > 0);
        }

        float MalariaAntibody::GetAntibodyCapacity() const
        {
            return m_antibody_capacity;
        }

        float MalariaAntibody::GetAntibodyConcentration() const
        {
            return m_antibody_concentration;
        }

        void MalariaAntibody::SetAntibodyCapacity( float antibody_capacity )
        {
            m_antibody_capacity = antibody_capacity;
        }

        void MalariaAntibody::SetAntibodyConcentration( float antibody_concentration )
        {
            m_antibody_concentration = antibody_concentration;
        }

        MalariaAntibodyType::Enum MalariaAntibody::GetAntibodyType() const
        {
            return m_antibody_type;
        }

        int MalariaAntibody::GetAntibodyVariant() const
        {
            return m_antibody_variant;
        }

        // EMOD 2.22: Time tracking methods
        void MalariaAntibody::SetTimeLastActive( float time )
        {
            m_time_last_active = time;
        }

        float MalariaAntibody::GetTimeLastActive() const
        {
            return m_time_last_active;
        }

        void MalariaAntibody::SetActiveIndex( int32_t index )
        {
            m_active_index = index;
        }

        int32_t MalariaAntibody::GetActiveIndex() const
        {
            return m_active_index;
        }

        //------------------------------------------------------------------

        IMalariaAntibody* MalariaAntibodyCSP::CreateAntibody( MalariaConfig* config, int variant, float capacity )
        {
            MalariaAntibodyCSP * antibody = new MalariaAntibodyCSP(config);
            antibody->Initialize( MalariaAntibodyType::CSP, variant, capacity );

            return antibody;
        }

        IMalariaAntibody* MalariaAntibodyMSP::CreateAntibody( MalariaConfig* config, int variant, float capacity )
        {
            MalariaAntibodyMSP * antibody = new MalariaAntibodyMSP(config);
            antibody->Initialize( MalariaAntibodyType::MSP1, variant, capacity );

            return antibody;
        }

        IMalariaAntibody* MalariaAntibodyPfEMP1Minor::CreateAntibody( MalariaConfig* config, int variant, float capacity )
        {
            MalariaAntibodyPfEMP1Minor * antibody = new MalariaAntibodyPfEMP1Minor(config);
            antibody->Initialize( MalariaAntibodyType::PfEMP1_minor, variant, capacity );

            return antibody;
        }

        IMalariaAntibody* MalariaAntibodyPfEMP1Major::CreateAntibody( MalariaConfig* config, int variant, float capacity )
        {
            MalariaAntibodyPfEMP1Major * antibody = new MalariaAntibodyPfEMP1Major(config);
            antibody->Initialize( MalariaAntibodyType::PfEMP1_major, variant, capacity );

            return antibody;
        }

    }

}
