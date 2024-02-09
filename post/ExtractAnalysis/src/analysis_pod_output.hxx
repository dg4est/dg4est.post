/**
 * File:   analysis_pod_output.hxx
 * Author: akirby
 *
 * Created on February 8, 2024, 2:01 PM
 */

#ifndef ANALYSIS_POD_OUTPUT_HXX
#define ANALYSIS_POD_OUTPUT_HXX

/* header files */
#include "analysis_read_inputs.hxx"
#include "analysis_var_inputs.hxx"
#include "analysis_var_extracts.hxx"

#ifdef __cplusplus
extern "C" {
#endif

void output_POD_energies(extracts_t *ext,inputs_t *inputs,int sindex,int neig);
void output_ROM_quantities(extracts_t *ext,inputs_t *inputs,int sindex,int timestep,int mode,char standalone);
void output_TVC_quantities(extracts_t *ext,inputs_t *inputs,int sindex);

#ifdef __cplusplus
}
#endif
#endif /* ANALYSIS_POD_OUTPUT_HXX */