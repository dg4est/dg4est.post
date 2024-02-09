/**
 * File:   analysis_reynolds_stresses.h
 * Author: akirby
 *
 * Created on November 30, 2020, 4:53 PM
 */

#ifndef ANALYSIS_REYNOLDS_STRESSES_H
#define ANALYSIS_REYNOLDS_STRESSES_H

/* system header files */
#include <string.h>

/* header files */
#include "analysis_var_inputs.hxx"
#include "analysis_var_extracts.hxx"
#include "analysis_read_inputs.hxx"
#include "analysis_temporal_mean.hxx"

#ifdef __cplusplus
extern "C" {
#endif

void reynolds_stresses(inputs_t *data,char standalone);
void compute_extract_stresses(extracts_t *stresses,extracts_t *means,inputs_t *data,int ex);
void compute_fluctuations_mean(extracts_t *stresses,extracts_t *means,char *file_name,Real frac);
void compute_resolved_turbulent_kinetic_energy(extracts_t *stresses);
void output_stresses(extracts_t *ext,inputs_t *inputs,int sindex,int eindex,int gl_mode,int mode,char standalone);

#ifdef __cplusplus
}
#endif
#endif /* ANALYSIS_REYNOLDS_STRESSES_H */