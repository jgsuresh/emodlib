/**
 * @file MalariaConfig.cpp
 *
 * @brief Malaria configuration implementation
 */

#include "MalariaConfig.h"

#include <cmath>


namespace emodlib
{

    namespace malaria
    {

        MalariaConfig::MalariaConfig()
            : rng(std::make_shared<PSEUDO_DES>(0, 256))  // Initialize with default seed
            , suidGenerator(0, 0)
        {
        }


        void MalariaConfig::Configure(const ParamSet& pset)
        {
            // IntrahostComponent-level params
            randomSeed = pset["Run_Number"].cast<int>();
            rng = std::make_shared<PSEUDO_DES>(randomSeed, 256);

            max_ind_inf = pset["Max_Individual_Infections"].cast<int>();

            falciparumMSPVars = pset["Falciparum_MSP_Variants"].cast<int>();
            falciparumNonSpecTypes = pset["Falciparum_Nonspecific_Types"].cast<int>();
            falciparumPfEMP1Vars = pset["Falciparum_PfEMP1_Variants"].cast<int>();

            base_gametocyte_mosquito_survival = pset["Base_Gametocyte_Mosquito_Survival_Rate"].cast<float>();
            cytokine_gametocyte_inactivation = pset["Cytokine_Gametocyte_Inactivation"].cast<float>();

            // Susceptibility params (from nested susceptibility_params)
            const ParamSet& susc_pset = pset["susceptibility_params"];

            memory_level = susc_pset["Antibody_Memory_Level"].cast<float>();
            hyperimmune_decay_rate = -log((0.4f - memory_level) / (1.0f - memory_level)) / 120.0f;
            MSP1_antibody_growthrate = susc_pset["Max_MSP1_Antibody_Growthrate"].cast<float>();
            antibody_stimulation_c50 = susc_pset["Antibody_Stimulation_C50"].cast<float>();
            antibody_capacity_growthrate = susc_pset["Antibody_Capacity_Growth_Rate"].cast<float>();
            minimum_adapted_response = susc_pset["Min_Adapted_Response"].cast<float>();
            non_specific_growth = susc_pset["Nonspecific_Antibody_Growth_Rate_Factor"].cast<float>();
            antibody_csp_decay_days = susc_pset["Antibody_CSP_Decay_Days"].cast<float>();

            maternal_antibody_decay_rate = susc_pset["Maternal_Antibody_Decay_Rate"].cast<float>();
            if (susc_pset.contains("Maternal_Antibody_Protection")) {
                maternal_antibody_protection = susc_pset["Maternal_Antibody_Protection"].cast<float>();
            }

            pyrogenic_threshold = susc_pset["Pyrogenic_Threshold"].cast<float>();
            fever_IRBC_killrate = susc_pset["Fever_IRBC_Kill_Rate"].cast<float>();

            erythropoiesis_anemia_effect = susc_pset["Erythropoiesis_Anemia_Effect"].cast<float>();

            // Infection params (from nested infection_params)
            const ParamSet& inf_pset = pset["infection_params"];

            incubation_period = inf_pset["Base_Incubation_Period"].cast<float>();
            antibody_IRBC_killrate = inf_pset["Antibody_IRBC_Kill_Rate"].cast<float>();
            non_specific_antigenicity = inf_pset["Nonspecific_Antigenicity_Factor"].cast<float>();
            MSP1_merozoite_kill = inf_pset["MSP1_Merozoite_Kill_Fraction"].cast<float>();
            gametocyte_stage_survival = inf_pset["Gametocyte_Stage_Survival_Rate"].cast<float>();
            base_gametocyte_sexratio = inf_pset["Base_Gametocyte_Fraction_Male"].cast<float>();
            base_gametocyte_production = inf_pset["Base_Gametocyte_Production_Rate"].cast<float>();
            antigen_switch_rate = inf_pset["Antigen_Switch_Rate"].cast<float>();
            merozoites_per_hepatocyte = inf_pset["Merozoites_Per_Hepatocyte"].cast<float>();
            merozoites_per_schizont = inf_pset["Merozoites_Per_Schizont"].cast<float>();
            RBC_destruction_multiplier = inf_pset["RBC_Destruction_Multiplier"].cast<float>();
            n_asexual_cycles_wo_gametocytes = inf_pset["Number_Of_Asexual_Cycles_Without_Gametocytes"].cast<int>();
        }


        std::shared_ptr<MalariaConfig> MalariaConfig::FromParamSet(const ParamSet& pset)
        {
            auto config = std::make_shared<MalariaConfig>();
            config->Configure(pset);
            return config;
        }

    }

}
