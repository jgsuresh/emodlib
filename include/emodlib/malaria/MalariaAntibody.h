#pragma once

#include "IMalariaAntibody.h"

namespace emodlib
{

    namespace malaria
    {
    
        class MalariaAntibody : public IMalariaAntibody
        {
        public:

            virtual ~MalariaAntibody() {}

        public:
            // IMalariaAntibody methods
            virtual void  Decay( float dt ) override;
            virtual float StimulateCytokines( float dt, float inv_uL_blood ) override;
            virtual void  UpdateAntibodyCapacity( float dt, float inv_uL_blood ) override;
            virtual void  UpdateAntibodyCapacityByRate( float dt, float growth_rate ) override;
            virtual void  UpdateAntibodyConcentration( float dt ) override;
            virtual void  ResetCounters() override;

            virtual void  IncreaseAntigenCount( int64_t antigenCount ) override;
            virtual void  SetAntigenicPresence( bool antigenPresent ) override;

            virtual int64_t GetAntigenCount() const override;
            virtual bool    GetAntigenicPresence() const override;
            virtual float   GetAntibodyCapacity() const override;
            virtual float   GetAntibodyConcentration() const override;

            virtual void    SetAntibodyCapacity( float antibody_capacity ) override;
            virtual void    SetAntibodyConcentration(float antibody_concentration) override;

            virtual MalariaAntibodyType::Enum GetAntibodyType() const override;
            virtual int GetAntibodyVariant() const override;

        protected:
            float   m_antibody_capacity;
            float   m_antibody_concentration;
            int64_t m_antigen_count;
            bool    m_antigen_present;

            MalariaAntibodyType::Enum m_antibody_type;
            int m_antibody_variant;

            MalariaAntibody();
            void Initialize( MalariaAntibodyType::Enum type, int variant, float capacity = 0, float concentration = 0 );
        };

        // -----------------------------------------------------------

        class MalariaAntibodyCSP : public MalariaAntibody
        {
        public:
            static IMalariaAntibody* CreateAntibody( int variant, float capacity=0.0f );
            virtual void UpdateAntibodyConcentration( float dt ) override;
            virtual void Decay( float dt ) override;
        };

        class MalariaAntibodyMSP : public MalariaAntibody
        {
        public:
            static IMalariaAntibody* CreateAntibody( int variant, float capacity=0.0f );
        };

        class MalariaAntibodyPfEMP1Minor : public MalariaAntibody
        {
        public:
            static IMalariaAntibody* CreateAntibody( int variant, float capacity=0.0f );
            virtual void UpdateAntibodyCapacity( float dt, float inv_uL_blood ) override;
        };

        class MalariaAntibodyPfEMP1Major : public MalariaAntibody
        {
        public:
            static IMalariaAntibody* CreateAntibody( int variant, float capacity=0.0f );
            virtual void UpdateAntibodyCapacity( float dt, float inv_uL_blood ) override;
        };
    }

}
