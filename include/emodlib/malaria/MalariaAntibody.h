#pragma once

#include "IMalariaAntibody.h"

namespace emodlib
{

    namespace malaria
    {

        struct MalariaConfig;  // forward declaration

        class MalariaAntibody : public IMalariaAntibody
        {
        public:

            virtual ~MalariaAntibody() {}

        public:
            // IMalariaAntibody methods
            virtual void  Decay( float decay_time ) override;  // decay_time = total inactive time
            virtual float StimulateCytokines( float dt, float inv_uL_blood ) override;
            virtual void  UpdateAntibodyCapacity( float dt, float inv_uL_blood ) override;
            virtual void  UpdateAntibodyCapacityByRate( float dt, float growth_rate ) override;
            virtual void  UpdateAntibodyConcentration( float dt ) override;
            virtual void  ResetCounters() override;

            // EMOD 2.22: IncreaseAntigenCount now takes currentTime and dt for time-based decay
            virtual void  IncreaseAntigenCount( int64_t antigenCount, float currentTime, float dt ) override;
            virtual void  SetAntigenicPresence( bool antigenPresent ) override;

            virtual int64_t GetAntigenCount() const override;
            virtual bool    GetAntigenicPresence() const override;
            virtual float   GetAntibodyCapacity() const override;
            virtual float   GetAntibodyConcentration() const override;

            virtual void    SetAntibodyCapacity( float antibody_capacity ) override;
            virtual void    SetAntibodyConcentration(float antibody_concentration) override;

            virtual MalariaAntibodyType::Enum GetAntibodyType() const override;
            virtual int GetAntibodyVariant() const override;

            // EMOD 2.22: Time tracking for deferred decay
            virtual void    SetTimeLastActive( float time ) override;
            virtual float   GetTimeLastActive() const override;
            virtual void    SetActiveIndex( int32_t index ) override;
            virtual int32_t GetActiveIndex() const override;

        protected:
            MalariaConfig* m_config;  // non-owning pointer to configuration

            float   m_antibody_capacity;
            float   m_antibody_concentration;
            int64_t m_antigen_count;

            MalariaAntibodyType::Enum m_antibody_type;
            int m_antibody_variant;

            // EMOD 2.22: Time tracking for deferred decay
            int32_t m_active_index;
            float   m_time_last_active;

            MalariaAntibody(MalariaConfig* config);
            void Initialize( MalariaAntibodyType::Enum type, int variant, float capacity = 0, float concentration = 0 );
        };

        // -----------------------------------------------------------

        class MalariaAntibodyCSP : public MalariaAntibody
        {
        public:
            MalariaAntibodyCSP(MalariaConfig* config) : MalariaAntibody(config) {}
            static IMalariaAntibody* CreateAntibody( MalariaConfig* config, int variant, float capacity=0.0f );
            virtual void UpdateAntibodyConcentration( float dt ) override;
            // Note: Decay() now handles CSP-specific logic in base class via type check
        };

        class MalariaAntibodyMSP : public MalariaAntibody
        {
        public:
            MalariaAntibodyMSP(MalariaConfig* config) : MalariaAntibody(config) {}
            static IMalariaAntibody* CreateAntibody( MalariaConfig* config, int variant, float capacity=0.0f );
        };

        class MalariaAntibodyPfEMP1Minor : public MalariaAntibody
        {
        public:
            MalariaAntibodyPfEMP1Minor(MalariaConfig* config) : MalariaAntibody(config) {}
            static IMalariaAntibody* CreateAntibody( MalariaConfig* config, int variant, float capacity=0.0f );
            virtual void UpdateAntibodyCapacity( float dt, float inv_uL_blood ) override;
        };

        class MalariaAntibodyPfEMP1Major : public MalariaAntibody
        {
        public:
            MalariaAntibodyPfEMP1Major(MalariaConfig* config) : MalariaAntibody(config) {}
            static IMalariaAntibody* CreateAntibody( MalariaConfig* config, int variant, float capacity=0.0f );
            virtual void UpdateAntibodyCapacity( float dt, float inv_uL_blood ) override;
        };
    }

}
