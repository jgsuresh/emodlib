#pragma once


#define CLONAL_PfEMP1_VARIANTS (50)

#define DEFAULT_MSP_VARIANTS 100
#define DEFAULT_NONSPECIFIC_TYPES 20
#define DEFAULT_PFEMP1_VARIANTS 1000
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
