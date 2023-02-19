/**
 * @file MalariaEnums.h
 *
 * @brief Malaria enumerations
 */

#pragma once

namespace emodlib
{

    namespace malaria
    {

        // ENUM defs for MALARIA_STRAINS
        // FALCIPARUM_RANDOM_STRAIN is default
        namespace MalariaStrains {
            enum Enum {
                FALCIPARUM_NONRANDOM_STRAIN = 11,
                FALCIPARUM_RANDOM50_STRAIN = 12,
                FALCIPARUM_RANDOM_STRAIN = 13,
                FALCIPARUM_STRAIN_GENERATOR = 20,
            };
        }
        
        // ENUM defs for PARASITE_SWITCH_TYPE
        // RATE_PER_PARASITE_7VARS is default
        namespace ParasiteSwitchType {
            enum Enum {
                CONSTANT_SWITCH_RATE_2VARS = 0,
                RATE_PER_PARASITE_7VARS = 1,
                RATE_PER_PARASITE_5VARS_DECAYING = 2,
            };
        }
    
        // ENUM defs for MATERNAL_ANTIBODIES_TYPE
        // SIMPLE_WANING draws a PfEMP1 antibody strength is initialized
        //               at a fraction of the mother's level and wanes exponentially
        // CONSTANT_INITIAL_IMMUNITY ignores the mother's level (or mean adult population level) of immunity
        //               and gives some constant initial immunity level to all children, which wanes exponentially.
        namespace MaternalAntibodiesType {
            enum Enum {
                OFF = 0,
                SIMPLE_WANING = 1,
                CONSTANT_INITIAL_IMMUNITY = 2,
            };
        }
    
        namespace InnateImmuneVariationType {
            enum Enum {
                NONE = 0,
                PYROGENIC_THRESHOLD = 1,
                CYTOKINE_KILLING = 2,
                PYROGENIC_THRESHOLD_VS_AGE = 3,
            };
        }
    
        namespace AsexualCycleStatus {
            enum Enum {
                NoAsexualCycle = 0,
                AsexualCycle = 1,
                HepatocyteRelease = 2,
            };
        }
    
        // CSP:    Circumsporozoite protein
        // MSP1:   Merozoite surface protein
        // PfEMP1: Plasmodium falciparum erythrocyte membrane protein (minor non-specific epitopes)
        //                                                            (major epitopes)

        namespace MalariaAntibodyType {
            enum Enum {
                CSP = 0,
                MSP1 = 1,
                PfEMP1_minor = 2,
                PfEMP1_major = 3,
                N_MALARIA_ANTIBODY_TYPES = 4,
            };
        }
    
        // 5 stages of development and mature gametocytes
        namespace GametocyteStages {
            enum Enum {
                Stage0 = 0,
                Stage1 = 1,
                Stage2 = 2,
                Stage3 = 3,
                Stage4 = 4,
                Mature = 5,
                Count = 6,
            };
        }
    
    }

}
