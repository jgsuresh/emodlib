/**
 * @file InfectionMalaria.cpp
 *
 * @brief Malaria infection implementation
 */

#include "InfectionMalaria.h"

#include <iostream>
#include <numeric>

#include "emodlib/utils/Common.h"
#include "emodlib/utils/Sigmoid.h"

#include "IntrahostComponent.h"
#include "SusceptibilityMalaria.h"


namespace emodlib
{

    namespace malaria
    {

        // TODO: emodlib#8 (boost + enums)
        // ParasiteSwitchType::Enum Infection::params::parasite_switch_type = ParasiteSwitchType::RATE_PER_PARASITE_7VARS;
        // MalariaStrains::Enum     Infection::params::malaria_strains = MalariaStrains::FALCIPARUM_RANDOM_STRAIN;

        float Infection::params::incubation_period = 7.0f; // liver stage duration

        float Infection::params::antibody_IRBC_killrate = DEFAULT_ANTIBODY_IRBC_KILLRATE;
        float Infection::params::non_specific_antigenicity = DEFAULT_NON_SPECIFIC_ANTIGENICITY;
        float Infection::params::MSP1_merozoite_kill = DEFAULT_MSP1_MEROZOITE_KILL;
        float Infection::params::gametocyte_stage_survival = DEFAULT_GAMETOCYTE_STAGE_SURVIVAL;
        float Infection::params::base_gametocyte_sexratio = DEFAULT_BASE_GAMETOCYTE_SEX_RATIO;
        float Infection::params::base_gametocyte_production = DEFAULT_BASE_GAMETOCYTE_PRODUCTION;
        float Infection::params::antigen_switch_rate = DEFAULT_ANTIGEN_SWITCH_RATE;
        float Infection::params::merozoites_per_hepatocyte = DEFAULT_MEROZOITES_PER_HEPATOCYTE;
        float Infection::params::merozoites_per_schizont = DEFAULT_MEROZOITES_PER_SCHIZONT;
        float Infection::params::RBC_destruction_multiplier = DEFAULT_RBC_DESTRUCTION_MULTIPLIER;
        int   Infection::params::n_asexual_cycles_wo_gametocytes = DEFAULT_ASEXUAL_CYCLES_WITHOUT_GAMETOCYTES;


        suids::distributed_generator Infection::infectionSuidGenerator(0, 0);


        void Infection::params::Configure(const ParamSet& pset)
        {
            incubation_period = pset["Base_Incubation_Period"].cast<float>();  // TODO: emodlib#6 (gaussian distribution)

            antibody_IRBC_killrate = pset["Antibody_IRBC_Kill_Rate"].cast<float>();
            non_specific_antigenicity = pset["Nonspecific_Antigenicity_Factor"].cast<float>();
            MSP1_merozoite_kill = pset["MSP1_Merozoite_Kill_Fraction"].cast<float>();
            gametocyte_stage_survival = pset["Gametocyte_Stage_Survival_Rate"].cast<float>();
            base_gametocyte_sexratio = pset["Base_Gametocyte_Fraction_Male"].cast<float>();
            base_gametocyte_production = pset["Base_Gametocyte_Production_Rate"].cast<float>();
            antigen_switch_rate = pset["Antigen_Switch_Rate"].cast<float>();
            merozoites_per_hepatocyte = pset["Merozoites_Per_Hepatocyte"].cast<float>();
            merozoites_per_schizont = pset["Merozoites_Per_Schizont"].cast<float>();
            RBC_destruction_multiplier = pset["RBC_Destruction_Multiplier"].cast<float>();
            n_asexual_cycles_wo_gametocytes = pset["Number_Of_Asexual_Cycles_Without_Gametocytes"].cast<int>();
        }


        Infection::Infection()
            : suid(suids::nil_suid())

            , m_liver_stage_timer(0.0f)
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

            , immunity(nullptr)
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

            immunity = _susceptibility;

            m_MSP_antibody = immunity->RegisterAntibody(MalariaAntibodyType::MSP1, m_MSPtype);

            for( int ivariant = 0; ivariant < m_PfEMP1_antibodies.size(); ivariant++ )
            {
                m_PfEMP1_antibodies[ivariant].major = nullptr;
                m_PfEMP1_antibodies[ivariant].minor = nullptr;

                if ( m_IRBC_count[ivariant] > 0 )
                {
                    m_PfEMP1_antibodies[ivariant].minor  = immunity->RegisterAntibody(MalariaAntibodyType::PfEMP1_minor, m_minor_epitope_type[ivariant]);
                    m_PfEMP1_antibodies[ivariant].major  = immunity->RegisterAntibody(MalariaAntibodyType::PfEMP1_major, m_IRBCtype[ivariant]);
                }
            }
        }

        void Infection::Update(float dt)
        {
            m_liver_stage_timer += dt;  // increment latent period

            if (m_hepatocytes > 0)
            {
                malariaProcessHepatocytes(dt);
            }

            if (m_asexual_phase > AsexualCycleStatus::NoAsexualCycle)
            {
                // do not decrement timer if it was just set by the hepatocytes this time step (asexual_phase==2), or else the timer gets decreased one timestep too many
                if (m_asexual_phase == AsexualCycleStatus::HepatocyteRelease)
                {
                    m_asexual_phase = AsexualCycleStatus::AsexualCycle;
                }
                else
                {
                    m_IRBCtimer -= dt;
                }

                // process end of asexual cycle events if appropriate
                if (m_IRBCtimer <= 0)
                {
                    processEndOfAsexualCycle();
                }

                // check for death due to death of all RBCs
                if (immunity->get_RBC_count() < 1)
                {
                    std::cout << "Individual has no more red-blood cells";
                    throw;  // TODO: emodlib#3 (InfectionStateChange::Killed)
                }

                // Immune Interaction
                // Infection Effect on Immune System--
                // in Susceptibility object update, antibody capacities increase and antibodies produced in response to antigenic-specific parasite load, tolerance--lack of inflamatory response-- develops
                malariaImmuneStimulation(dt);

                // Immune and Drug Effects on Infection
                malariaImmunityIRBCKill(dt);

                // Immune and Drug Effects on Gametocytes
                malariaImmunityGametocyteKill(dt);

                //make sure MSP type generates antibodies during an ongoing infection, not just during the short time of IRBC rupturing, since the stimulation may persist
                m_MSP_antibody->IncreaseAntigenCount(1);
                immunity->SetAntigenPresent(); // NOTE: this has an interesting behavior in that it continues to update MSP capacity AFTER there are no IRBC (only gametocytes)
            }

            // check for death, clearance, and take care of end-of-timestep bookkeeping
            malariaCheckInfectionStatus(dt);
        }

        void Infection::malariaProcessHepatocytes(float dt)
        {
            // check for valid inputs
            if (dt > 0 && immunity && m_hepatocytes > 0)
            {
                // TODO: emodlib#4 (hepatocyte drug killing)

                // ----------------------------------------------------------------------------------------------------------------------
                // --- latency in hepatocyte phase Collins, W. E. and G. M. Jeffery (1999).
                // --- "A retrospective examination of sporozoite- and trophozoite-induced infections with Plasmodium falciparum:
                // --- development of parasitologic and clinical immunity during primary infection." Am J Trop Med Hyg 61(1 Suppl): 4-19.
                // --- process start of asexual phase if the incubation period is over and there are still hepatocytes
                // ----------------------------------------------------------------------------------------------------------------------
                if (m_asexual_phase == AsexualCycleStatus::NoAsexualCycle &&
                     m_liver_stage_timer >= Infection::params::incubation_period)
                {
                    m_IRBC_count.assign(CLONAL_PfEMP1_VARIANTS, 0);

                    // testing starting with multiple antigens, which reduces the probability of a single first variant being cleared by a pre-existing antibody response
                    // picked starting with 5 variants after exploring different options in work developing Intrahost model
                    const int INITIAL_PFEMP1_VARIANTS = 5;
                    #pragma loop(hint_parallel(8))
                    for ( int i=0; i<INITIAL_PFEMP1_VARIANTS; i++ )
                    {
                        m_IRBC_count[i] = int64_t(m_hepatocytes * Infection::params::merozoites_per_hepatocyte / INITIAL_PFEMP1_VARIANTS);
                        immunity->UpdateActiveAntibody( m_PfEMP1_antibodies[i], m_minor_epitope_type[i], m_IRBCtype[i] ); // insert into set of antigens the immune system has ever "seen"
                    }

                    // now back to normal
                    m_hepatocytes   = 0;
                    m_IRBCtimer     = 2.0;  // P. falciparum has a 2-day asexual cycle
                    m_asexual_phase = AsexualCycleStatus::HepatocyteRelease;    // HepatocyteRelease means the asexual cycle is just beginning, so that the IRBCtimer is not decremented in the same time step
                }
            }
        }

        void Infection::processEndOfAsexualCycle()
        {
            // Merozoite-specific antibodies can limit merozoite success--Blackman, M. J., H. G. Heidrich, et al. (1990).
            // "A single fragment of a malaria merozoite surface protein remains on the parasite during red cell invasion
            // and is the target of invasion-inhibiting antibodies." J Exp Med 172(1): 379-382.
            double RBCavailability = immunity->get_RBC_availability();

            // Merozoite survival limited at very low density according to density-dependent probability-of-success formula
            double merozoitesurvival = std::max(0.0, (1.0 - Infection::params::MSP1_merozoite_kill * m_MSP_antibody->GetAntibodyConcentration() ) * EXPCDF(-RBCavailability / MEROZOITE_LIMITING_RBC_THRESHOLD));

            // How many rupture for this infection handed to suscept object for total stimulation calculations
            int64_t totalIRBC = 0;
            totalIRBC = std::accumulate( m_IRBC_count.begin(), m_IRBC_count.end(), totalIRBC );
            m_MSP_antibody->IncreaseAntigenCount( totalIRBC );

            // Move immature gametocytes forward a stage and create initial stage gametocytes from previous merozoites
            // This is the last function to use m_IRBC_count from the previous cycle
            malariaCycleGametocytes(merozoitesurvival);

            // Calculate antigenic switching and create asexual IRBCss for next asexual cycle
            // After this function, m_IRBC_count will have been updated
            malariaIRBCAntigenSwitch(merozoitesurvival);

            totalIRBC = 0;
            //std::accumulate( m_IRBC_count.begin(), m_IRBC_count.end(), totalIRBC );
            #pragma loop(hint_parallel(8))
            for ( int j = 0; j < CLONAL_PfEMP1_VARIANTS; j++ )
            {
                if ( m_IRBC_count[j] > 0 )
                {
                    totalIRBC += m_IRBC_count[j];
                    immunity->UpdateActiveAntibody( m_PfEMP1_antibodies[j], m_minor_epitope_type[j], m_IRBCtype[j] ); // insert into set of antigens the immune system has ever "seen"
                }
            }

            // Uninfected RBC killing diminishing in proportion to RBC availability
            double destruction_factor_ = std::max(1.0, Infection::params::RBC_destruction_multiplier * EXPCDF(-RBCavailability / MEROZOITE_LIMITING_RBC_THRESHOLD) );
            immunity->remove_RBCs( totalIRBC, m_malegametocytes[0] + m_femalegametocytes[0], destruction_factor_ );

            // reset timer for next asexual cycle
            m_IRBCtimer = 2.0;

            // increment counter of completed asexual cycles
            m_asexual_cycle_count++;
        }

        // Calculates the antigenic switching when an asexual cycle completes and creates next generation of IRBC's
        void Infection::malariaIRBCAntigenSwitch(double merozoitesurvival)
        {
            int64_t switchingIRBC[SWITCHING_IRBC_VARIANT_COUNT];
            std::vector<int64_t> tmpIRBCcount(CLONAL_PfEMP1_VARIANTS);

            // check for valid range of input, and only create next cycle if valid
            if (merozoitesurvival < 0)
            {
                std::cout << "merozoitesurvival should not be negative";
                throw;
            }

            // Several antigen switching mechanisms are supported
            #pragma loop(hint_parallel(8))
            for (int j = 0; j < CLONAL_PfEMP1_VARIANTS; j++)
            {
                // parasite switching studied in Paget-McNicol, S., M. Gatton, et al. (2002). "The Plasmodium falciparum var gene switching rate, switching mechanism and patterns of parasite recrudescence described by mathematical modelling." Parasitology 124(Pt 3): 225-235.
                // experimental studies in Horrocks, P., R. Pinches, et al. (2004). "Variable var transition rates underlie antigenic variation in malaria." Proceedings of the National Academy of Sciences of the United States of America 101(30): 11129-11134.
                // review in Horrocks, P., S. A. Kyes, et al. (2005). Molecular Aspects of Antigenic Variation in Plasmodium falciparum. Molecular Approaches to Malaria. I. W. Sherman. Washington DC, ASM Press: 399-415.

                if ( m_IRBC_count[j] <= 0 ) continue; // no IRBC means no contribution to next time step

                int64_t temp_sum_IRBC = 0;
                if (Infection::params::antigen_switch_rate > 0)
                {
                    #pragma loop(hint_parallel(8))
                    for ( int iswitch = 0; iswitch < SWITCHING_IRBC_VARIANT_COUNT; iswitch++ )
                    {
                        switchingIRBC[iswitch] = (iswitch < 7) ? IntrahostComponent::p_rng->Poisson(Infection::params::antigen_switch_rate * m_IRBC_count[j]) : 0;
                    }

                    // now test to see if these add up to more than 100 percent
                    temp_sum_IRBC = std::accumulate(switchingIRBC, switchingIRBC + SWITCHING_IRBC_VARIANT_COUNT, temp_sum_IRBC);

                    // if more than 100 percent minus those switching to gametocyte production, scale down in multiplicative way
                    if (temp_sum_IRBC > ((1.0 - m_gametorate)*m_IRBC_count[j]))
                    {
                        #pragma loop(hint_parallel(8))
                        for (int iswitch = 0; iswitch < SWITCHING_IRBC_VARIANT_COUNT; iswitch++)
                            switchingIRBC[iswitch] = int64_t(switchingIRBC[iswitch] * ((1.0f - m_gametorate) * m_IRBC_count[j] / temp_sum_IRBC));

                        temp_sum_IRBC = int64_t((1.0 - m_gametorate) * m_IRBC_count[j]);
                    }
                }

                // Now switch to next stages based on predetermined number of switching IRBC's
                tmpIRBCcount[j] = int64_t(tmpIRBCcount[j] + ((1.0 - m_gametorate) * m_IRBC_count[j] - temp_sum_IRBC) * Infection::params::merozoites_per_schizont * merozoitesurvival);
                if (Infection::params::antigen_switch_rate > 0)
                {
                    #pragma loop(hint_parallel(8))
                    for ( int iswitch = 0; iswitch < SWITCHING_IRBC_VARIANT_COUNT; iswitch++)
                    {
                        tmpIRBCcount[(j + iswitch + 1) % CLONAL_PfEMP1_VARIANTS]  = int64_t(tmpIRBCcount[(j + iswitch + 1) % CLONAL_PfEMP1_VARIANTS] + switchingIRBC[iswitch] * Infection::params::merozoites_per_schizont * merozoitesurvival);
                    }
                }
            }

            m_IRBC_count.swap(tmpIRBCcount); // swap temporarily accumulated vector of next time step into data member
        }

        // Moves all falciparum gametocytes forward a development stage when an asexual cycle completes, and creates the stage 0 immature gametocytes
        void Infection::malariaCycleGametocytes(double merozoitesurvival)
        {
            // set gametocyte production rate for next cycle
            if ( m_asexual_cycle_count >= Infection::params::n_asexual_cycles_wo_gametocytes )
            {
                m_gametorate     = double(Infection::params::base_gametocyte_production); // gametocyte production used by all switching calculations, here is where factors modifying production would go
                m_gametosexratio = double(Infection::params::base_gametocyte_sexratio);
            }

            // check for valid range of input, and only create next cycle if valid
            if (merozoitesurvival >= 0)
            {
                #pragma loop(hint_parallel(8))
                //process gametocytes--5 stages--Sinden, R. E., G. A. Butcher, et al. (1996). "Regulation of Infectivity of Plasmodium to the Mosquito Vector." Advances in Parasitology 38: 53-117.
                for (int j = GametocyteStages::Mature; j > 0; j--) // move developing gametocytes forward a class, moving backwards through stages to not override next stage's values
                {
                    m_malegametocytes[j] = int64_t(m_malegametocytes[j] + m_malegametocytes[j - 1] * Infection::params::gametocyte_stage_survival);
                    m_malegametocytes[j - 1] = 0;

                    if (m_malegametocytes[j] < 1)
                        m_malegametocytes[j] = 0;

                    m_femalegametocytes[j] = int64_t(m_femalegametocytes[j] + (m_femalegametocytes[j - 1] * Infection::params::gametocyte_stage_survival));
                    m_femalegametocytes[j - 1] = 0;

                    if (m_femalegametocytes[j] < 1)
                         m_femalegametocytes[j] = 0;
                }

                #pragma loop(hint_parallel(8))
                // Now create the new stage 1 gametocytes based on production ratios and the prevIRBC counts
                for (int j = 0; j < CLONAL_PfEMP1_VARIANTS; j++)
                {
                    // review of production rates and sex ratios in Sinden, R. E., G. A. Butcher, et al. (1996). "Regulation of Infectivity of Plasmodium to the Mosquito Vector." Advances in Parasitology 38: 53-117.
                    // each factor may be variable, but here we leave it constant at the moment, conservatively not including the possible senescence of transmission in late infection
                    m_malegametocytes[GametocyteStages::Stage0]   = int64_t(m_malegametocytes[GametocyteStages::Stage0]   + m_IRBC_count[j] * m_gametorate * m_gametosexratio * merozoitesurvival * Infection::params::merozoites_per_schizont);
                    m_femalegametocytes[GametocyteStages::Stage0] = int64_t(m_femalegametocytes[GametocyteStages::Stage0] + m_IRBC_count[j] * m_gametorate * (1.0 - m_gametosexratio) * merozoitesurvival * Infection::params::merozoites_per_schizont);
                }
            }
        }

        // Calculates stimulation of immune system by malaria infection
        void Infection::malariaImmuneStimulation(float dt)
        {
            // check for valid inputs
            if ( dt <= 0 || immunity == nullptr )
            {
                std::cout << "Invalid input to malariaImmuneStimulation" << std::endl;
                return;
            }

            // antibody capacity for MSP 1 and MSP-2 are above
            // antibody capacity for the different RBC surface variants
            // transfer total IRBC to array owned by Susceptibility_Malaria, which then calculates total immune stimulation by all concurrent infections
            #pragma loop(hint_parallel(8))
            for (int i = 0; i < CLONAL_PfEMP1_VARIANTS; i++)
            {
                if (m_IRBC_count[i] < 0)
                {
                    m_IRBC_count[i] = 0;
                    std::cout << "malariaImmuneStimulation() IRBC count at index" << i << "should not be negative" << std::endl;
                }

                // only update if there are actually IRBCs
                if (m_IRBC_count[i] > 0)
                {
                    // PfEMP-1 major epitopes
                    m_PfEMP1_antibodies[i].major->IncreaseAntigenCount(m_IRBC_count[i]);

                    // PfEMP-1 minor epitopes
                    m_PfEMP1_antibodies[i].minor->IncreaseAntigenCount(m_IRBC_count[i]);

                    // Notify susceptibility that there is antigen present
                    immunity->SetAntigenPresent();
                }
            }
        }

        // Calculates the IRBC killing from drugs and immune action
        void Infection::malariaImmunityIRBCKill(float dt)
        {
            // check for valid inputs
            if (dt > 0 && immunity)
            {
                // inflammatory response--Stevenson, M. M. and E. M. Riley (2004).
                // "Innate immunity to malaria." Nat Rev Immunol 4(3): 169-180.

                // Offset basic sigmoid: effect rises as basic sigmoid beginning from a fever of MIN_FEVER_DEGREES_KILLING
                double fever_cytokine_killrate = (immunity->get_fever() > MIN_FEVER_DEGREES_KILLING) ? immunity->get_fever_killing_rate() * Sigmoid::basic_sigmoid(1.0, immunity->get_fever() - MIN_FEVER_DEGREES_KILLING) : 0.0;

                // TODO: emodlib#4 (asexual-stage drug killing)
                double drug_killrate = 0;

                #pragma loop(hint_parallel(8))
                for (int i = 0; i < CLONAL_PfEMP1_VARIANTS; i++)
                {
                    if ( m_IRBC_count[i] == 0 ) continue; // don't need to estimate killing if there are no IRBC of this variant to kill!

                    // total = antibodies (major, minor, maternal) + fever + drug
                    double pkill = EXPCDF(-dt * ( (m_PfEMP1_antibodies[i].major->GetAntibodyConcentration() + Infection::params::non_specific_antigenicity * m_PfEMP1_antibodies[i].minor->GetAntibodyConcentration() + immunity->get_maternal_antibodies() ) * Infection::params::antibody_IRBC_killrate + fever_cytokine_killrate + drug_killrate));

                    // Now here there is an interesting issue: to save massive amounts of computational time, can use a Gaussian approximation for the true binomial, but this returns a float
                    // This is fine for large numbers of killed IRBC's, but an issue arises for small numbers
                    // big question, is 1.5 killed IRBC's 1 or 2 killed?

                    double tempval1 = m_IRBC_count[i] * pkill;
                    if ( tempval1 > 0 ) // don't need to smear the killing by a random number if it is going to be zero
                        tempval1 = IntrahostComponent::p_rng->eGauss() * sqrt(tempval1 * (1.0 - pkill)) + tempval1;

                    if (tempval1 < 0.5)
                        tempval1 = 0;


                    // so add a continuity correction 0.5, and then convert to integer
                    m_IRBC_count[i] -= int64_t(tempval1 + 0.5);

                    if (m_IRBC_count[i] < 1)
                        m_IRBC_count[i] = 0;   // check for too large a time step

                }
            }

        }

        // Calculates immature gametocyte killing from drugs and immune action
        void Infection::malariaImmunityGametocyteKill(float dt)
        {
            // check for valid inputs
            if (dt > 0 && immunity)
            {
                #pragma loop(hint_parallel(8))
                for (int i = 0; i < GametocyteStages::Mature; i++)
                {
                    // Currently have fever and inflammatory cytokines limiting infectivity,
                    // rather than killing gametocytes.  See IndividualHumanMalaria::DepositInfectiousnessFromGametocytes()
                    // We leave this variable here at zero incase we change DepositInfectiousnessFromGametocytes()
                    double fever_cytokine_killrate = 0; // 0 = don't kill due to fever

                    // TODO: emodlib#4 (early- and late-stage gametocyte drug killing)
                    double drug_killrate = 0;

                    // no randomness in gametocyte killing, but a continuity correction
                    double gametocyte_kill_fraction = EXPCDF( -dt * (fever_cytokine_killrate + drug_killrate) );

                    m_malegametocytes[i] -= int64_t( 0.5 + m_malegametocytes[i] * gametocyte_kill_fraction );
                    if (m_malegametocytes[i] < 1)
                        m_malegametocytes[i] = 0;

                    m_femalegametocytes[i] -= int64_t( 0.5 + m_femalegametocytes[i] * gametocyte_kill_fraction );
                    if (m_femalegametocytes[i] < 1)
                        m_femalegametocytes[i] = 0;
                }

                // TODO: emodlib#5 (mature gametocyte decay)

                float drugGametocyteKill = 0;  // TODO: emodlib#4 (mature gametocyte drug killing)
                double pkill = EXPCDF( -dt * (0.277 + drugGametocyteKill) ); // half-life of 2.5 days corresponds to a decay time constant of 3.6 days, 0.277 = 1/3.6
                apply_MatureGametocyteKillProbability( pkill );
            }
        }

        int64_t ApplyKillProbability( float pkill, int64_t initial_gc, double eGauss )
        {
            double numkilled = (eGauss * sqrt( pkill * initial_gc * (1.0 - pkill) ) + pkill * initial_gc);
            numkilled = std::max( 0.0, numkilled ); //can't add by killing
            int64_t new_gc = int64_t( initial_gc - numkilled );
            return std::max( (int64_t)0, new_gc );
        }

        void Infection::apply_MatureGametocyteKillProbability(float pkill)
        {
            // Gaussian approximation of binomial errors for male and female mature gametocytes
            m_femalegametocytes[ GametocyteStages::Mature ] = ApplyKillProbability( pkill, m_femalegametocytes[ GametocyteStages::Mature ], IntrahostComponent::p_rng->eGauss() );
            m_malegametocytes[ GametocyteStages::Mature ] = ApplyKillProbability(   pkill, m_malegametocytes[   GametocyteStages::Mature ], IntrahostComponent::p_rng->eGauss() );
        }

        void Infection::malariaCheckInfectionStatus(float dt)
        {
            // TODO: emodlib#3 (InfectionStateChange::Cleared)
            // if (hepatocytes + IRBC + gametocytes) = 0
        }

        suids::suid Infection::GetSuid() const
        {
            return suid;
        }

        int64_t Infection::get_MaleGametocytes(int stage) const
        {
            return m_malegametocytes[stage];
        }

        int64_t Infection::get_FemaleGametocytes(int stage) const
        {
            return m_femalegametocytes[stage];
        }

        float Infection::get_asexual_density() const
        {
            int64_t totalIRBC = 0;
            totalIRBC = std::accumulate(m_IRBC_count.begin(), m_IRBC_count.end(), totalIRBC);
            return totalIRBC * immunity->get_inv_microliters_blood();
        }

        float Infection::get_mature_gametocyte_density() const
        {
            int64_t mature_female_gametocytes = get_FemaleGametocytes(GametocyteStages::Mature);
            return mature_female_gametocytes * immunity->get_inv_microliters_blood();
        }

        bool Infection::IsCleared() const {

            int64_t totalIRBC = 0;
            totalIRBC = std::accumulate( m_IRBC_count.begin(), m_IRBC_count.end(), totalIRBC );

            int64_t totalgametocytes = 0;
            for (int i = 0; i <= GametocyteStages::Mature; i++)
            {
                totalgametocytes += m_malegametocytes[i] + m_femalegametocytes[i];
            }

            return (totalIRBC + m_hepatocytes + totalgametocytes) < 1;
        }

        int32_t Infection::get_msp_type() const
        {
            return m_MSPtype;
        }

        std::vector<int32_t> Infection::get_pfemp1_major_types() const
        {
            std::vector<int32_t> vi;
            vi.assign(m_IRBCtype, m_IRBCtype + CLONAL_PfEMP1_VARIANTS);
            return vi;
        }

        IMalariaAntibody* Infection::get_msp_antibody() const
        {
            return m_MSP_antibody;
        }


    }

}
