/**
 * File:   analysis_binary2ascii.hxx
 * Author: akirby
 *
 * Created on February 9, 2024, 12:50 PM
 */

#ifndef ANALYSIS_BINARY2ASCII_HXX
#define ANALYSIS_BINARY2ASCII_HXX

/* system header files */
#include <math.h>
#include <string.h>
#include <assert.h>

/* header files */
#include "analysis_var_inputs.hxx"
#include "analysis_var_extracts.hxx"
#include "analysis_read_inputs.hxx"

#ifdef __cplusplus
extern "C" {
#endif

void binary2ascii(inputs_t *data);
void convert_extract_files(extracts_t *stats,inputs_t *data,int ex);
void read_binary_data(extracts_t *stats,char *file_name);
void write_ascii_data(extracts_t *ext,inputs_t *inputs,int sindex,int istep,int mode);

#ifdef __cplusplus
}
#endif
#endif /* ANALYSIS_BINARY2ASCII_HXX */