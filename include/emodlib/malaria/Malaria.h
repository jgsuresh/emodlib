#pragma once


#define CLONAL_PfEMP1_VARIANTS (50)
#define MINOR_EPITOPE_VARS_PER_SET 5

#define INV_MICROLITERS_BLOOD_ADULT (1.0/5e6) // 5 liters of blood/adult (http://hypertextbook.com/facts/1998/LanNaLee.shtml)

// Typical red blood cell (RBC) concentration: 5x10^6 RBC/microliter
// Use 5 liters for average adult blood volume
// 5x10^6 x 5x10^6 = 2.5x10^13 RBC/adult
// With average RBC lifetime of 120 days, daily adult RBC production is
// 2.5x10^13 / 120 ~= 2x10^11 RBC/day

#define ADULT_RBC_PRODUCTION    (2e11)      // RBC production per day: Williams Hematology p.231

// Use 1.2 liters for average infant blood volume
// 5x10^6 x 1.2x10^6 = 6x10^12 RBC/infant
// With average RBC lifetime of 120 days, daily infant RBC production is
// 6x10^12 / 120 = 5x10^10

#define INFANT_RBC_PRODUCTION   ( 1.5e10)     // Redefine down from 5e10 so that infants start with physiologic anemia = 11-12g/dl, then increase to adult with increasing blood volume

#define AVERAGE_RBC_LIFESPAN        (120)     // days
#define AVERAGE_RBC_CONCENTRATION   (5e6)     // RBC/microliter

#define FEVER_DEGREES_CELSIUS_PER_UNIT_CYTOKINES (4)
#define CYTOKINE_STIMULATION_SCALE (1.0)
#define MIN_FEVER_DEGREES_KILLING (1.5)

#define MEROZOITE_LIMITING_RBC_THRESHOLD (0.2)
#define SWITCHING_IRBC_VARIANT_COUNT    (10)

#define DEFAULT_MSP_VARIANTS 100
#define DEFAULT_NONSPECIFIC_TYPES 20
#define DEFAULT_PFEMP1_VARIANTS 1000

#define DEFAULT_ANTIBODY_IRBC_KILLRATE      (2.0)
#define DEFAULT_NON_SPECIFIC_ANTIGENICITY   (0.2)
#define DEFAULT_MSP1_MEROZOITE_KILL         (0.5)
#define DEFAULT_GAMETOCYTE_STAGE_SURVIVAL   (1.0)
#define DEFAULT_BASE_GAMETOCYTE_SEX_RATIO   (0.2f)
#define DEFAULT_BASE_GAMETOCYTE_PRODUCTION  (0.02f)
#define DEFAULT_ANTIGEN_SWITCH_RATE         (2.0e-9f)
#define DEFAULT_MEROZOITES_PER_HEPATOCYTE   (15000)
#define DEFAULT_MEROZOITES_PER_SCHIZONT     (16)
#define DEFAULT_RBC_DESTRUCTION_MULTIPLIER  (9.5)
#define DEFAULT_ASEXUAL_CYCLES_WITHOUT_GAMETOCYTES (1)

#define DEFAULT_ANTIBODY_CSP_DECAY_DAYS (90)
#define DEFAULT_FEVER_IRBC_KILL_RATE    (0.15f)

// TODO: emodlib#7 (infectiousness calculations)
#define DEFAULT_BASE_GAMETOCYTE_MOSQUITO_SURVIVAL   (0.01f)
#define DEFAULT_CYTOKINE_GAMETOCYTE_INACTIVATION    (0.02f)

#define MICROLITERS_PER_BLOODMEAL (2.0)  // from Vector.h
