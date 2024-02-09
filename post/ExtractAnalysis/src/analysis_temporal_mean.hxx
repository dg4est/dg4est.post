/**
 * File:   analysis_temporal_mean.h
 * Author: akirby
 *
 * Created on November 21, 2020, 12:23 PM
 */

#ifndef ANALYSIS_TEMPORAL_MEAN_H
#define ANALYSIS_TEMPORAL_MEAN_H

/* system header files */
#include <math.h>
#include <string.h>

/* header files */
#include "analysis_read_inputs.hxx"
#include "analysis_var_inputs.hxx"
#include "analysis_var_extracts.hxx"

#ifdef __cplusplus
extern "C" {
#endif

void temporal_mean(inputs_t *inputs,char standalone);
void compute_extract_mean(extracts_t *mean,inputs_t *data,int ex);
Real compute_extract_stddev(extracts_t *stats,inputs_t *data,int ex);
void compute_mean(extracts_t *mean,inputs_t *inputs,char *file_name,Real frac);
void compute_stddev(extracts_t *stats,char *file_name,Real frac);
void compute_stddev_final(extracts_t *stats);
void output_mean(extracts_t *ext,inputs_t *inputs,int sindex,int eindex,int gl_mode,int mode,char standalone);
void output_statistics(extracts_t *ext,inputs_t *inputs,int sindex,int eindex,int gl_mode,int mode,char standalone,Real frac);

#ifdef __cplusplus
}
#endif
#endif /* ANALYSIS_TEMPORAL_MEAN_H */