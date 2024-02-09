/**
 * File:   analysis_les_subgridscales.h
 * Author: akirby
 *
 * Created on December 1, 2020, 8:13 PM
 */

#ifndef ANALYSIS_LES_SUBGRIDSCALES_H
#define ANALYSIS_LES_SUBGRIDSCALES_H

/* system header files */
#include <math.h>
#include <string.h>
#include <assert.h>

/* header files */
#include "analysis_var_inputs.hxx"
#include "analysis_var_extracts.hxx"
#include "analysis_read_inputs.hxx"
#include "analysis_temporal_mean.hxx"
#include "analysis_reynolds_stresses.hxx"

#ifdef __cplusplus
extern "C" {
#endif

void les_subgridscales(inputs_t *data,char standalone);
void compute_tke(extracts_t *tke,extracts_t *stresses,extracts_t *mean,inputs_t *inputs,char *file_name);
void output_sgs_quantities(extracts_t *ext,inputs_t *inputs,
                           int sindex,int timestep,int mode,char standalone);

#ifdef __cplusplus
}
#endif
#endif /* ANALYSIS_LES_SUBGRIDSCALES_H */