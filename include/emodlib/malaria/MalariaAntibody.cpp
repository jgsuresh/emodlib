#include "MalariaAntibody.h"

#include "emodlib/utils/Sigmoid.h"
#include "SusceptibilityMalaria.h"


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

        MalariaAntibody::MalariaAntibody()
            : m_antigen_count(0)
            , m_antigen_present(false)
        {
        }

        void MalariaAntibody::Initialize( MalariaAntibodyType::Enum type, int variant, float capacity, float concentration )
        {
            m_antibody_type          = type;
            m_antibody_variant       = variant;
            m_antibody_capacity      = capacity;
            m_antibody_concentration = concentration;
        }

        void MalariaAntibody::Decay( float dt )
        {
            // don't do multiplication and subtraction unless antibody levels non-trivial
            if ( m_antibody_concentration > NON_TRIVIAL_ANTIBODY_THRESHOLD )
            {
                m_antibody_concentration -= m_antibody_concentration * TWENTY_DAY_DECAY_CONSTANT * dt;  //twenty day decay constant
            }

            // antibody capacity decays to a medium value (.3) dropping below .4 in ~120 days from 1.0
            if ( m_antibody_capacity > Susceptibility::params::memory_level )
            {
                m_antibody_capacity -= ( m_antibody_capacity - Susceptibility::params::memory_level) * Susceptibility::params::hyperimmune_decay_rate * dt;
            }
        }

        void MalariaAntibodyCSP::Decay( float dt )
        {
            // allow the decay of anti-CSP concentrations greater than unity (e.g. after boosting by vaccine)
            if ( m_antibody_concentration > m_antibody_capacity )
            {
                m_antibody_concentration -= m_antibody_concentration * dt / Susceptibility::params::antibody_csp_decay_days;
            }
            else
            {
                // otherwise do the normal behavior of decaying antibody concentration based on capacity
                MalariaAntibody::Decay( dt );
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
            float growth_rate = Susceptibility::params::MSP1_antibody_growthrate;
            float threshold   = Susceptibility::params::antibody_stimulation_c50;

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
            float min_stimulation = Susceptibility::params::antibody_stimulation_c50 * Susceptibility::params::minimum_adapted_response;
            float growth_rate     = Susceptibility::params::antibody_capacity_growthrate * Susceptibility::params::non_specific_growth;
            float threshold       = Susceptibility::params::antibody_stimulation_c50;

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
            float min_stimulation = Susceptibility::params::antibody_stimulation_c50 * Susceptibility::params::minimum_adapted_response;
            float growth_rate     = Susceptibility::params::antibody_capacity_growthrate;
            float threshold       = Susceptibility::params::antibody_stimulation_c50;

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
                m_antibody_concentration -= m_antibody_concentration * dt / Susceptibility::params::antibody_csp_decay_days;
            }
            else
            {
                // otherwise do the normal behavior of incrementing antibody concentration based on capacity
                MalariaAntibody::UpdateAntibodyConcentration( dt );
            }
        }

        void MalariaAntibody::ResetCounters()
        {
            m_antigen_present = false;
            m_antigen_count   = 0;
        }

        void MalariaAntibody::IncreaseAntigenCount( int64_t antigenCount )
        {
            if( antigenCount > 0 )
            {
                m_antigen_count += antigenCount;
                m_antigen_present = true;
            }
        }

        void MalariaAntibody::SetAntigenicPresence( bool antigenPresent )
        {
            m_antigen_present = antigenPresent;
        }

        int64_t MalariaAntibody::GetAntigenCount() const
        {
            return m_antigen_count;
        }

        bool MalariaAntibody::GetAntigenicPresence() const
        {
            return m_antigen_present;
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

        //------------------------------------------------------------------

        IMalariaAntibody* MalariaAntibodyCSP::CreateAntibody( int variant, float capacity )
        {
            MalariaAntibodyCSP * antibody = new MalariaAntibodyCSP();
            antibody->Initialize( MalariaAntibodyType::CSP, variant, capacity );

            return antibody;
        }

        IMalariaAntibody* MalariaAntibodyMSP::CreateAntibody( int variant, float capacity )
        {
            MalariaAntibodyMSP * antibody = new MalariaAntibodyMSP();
            antibody->Initialize( MalariaAntibodyType::MSP1, variant, capacity );

            return antibody;
        }

        IMalariaAntibody* MalariaAntibodyPfEMP1Minor::CreateAntibody( int variant, float capacity )
        {
            MalariaAntibodyPfEMP1Minor * antibody = new MalariaAntibodyPfEMP1Minor();
            antibody->Initialize( MalariaAntibodyType::PfEMP1_minor, variant, capacity );

            return antibody;
        }

        IMalariaAntibody* MalariaAntibodyPfEMP1Major::CreateAntibody( int variant, float capacity )
        {
            MalariaAntibodyPfEMP1Major * antibody = new MalariaAntibodyPfEMP1Major();
            antibody->Initialize( MalariaAntibodyType::PfEMP1_major, variant, capacity );

            return antibody;
        }

    }

}
