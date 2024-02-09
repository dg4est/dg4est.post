/**
 * File:   analysis_read_inputs.h
 * Author: akirby
 *
 * Created on November 21, 2020, 10:55 PM
 */

#ifndef READ_INPUTS_H
#define READ_INPUTS_H

/* header files */
#include "analysis_var_inputs.hxx"

#ifndef DOXYGEN_IGNORE
#  include <math.h>
#  include <ctype.h>
#  include <unistd.h>
#  include <string.h>
#  include <stdlib.h>
#  include <stdio.h>
#  include <errno.h>
#  include <sys/stat.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

inputs_t* read_inputs(int argc,char **argv);

/* read utility functions */
int form_file_name(inputs_t *data,int ex,int i,char **filename,char disp,const char *msg);
void create_directory(const char dir[]);
char find_keyword_integer(char *filename,const char *keyword,int *integer,int required);
char find_keyword_two_integers(char *filename,const char *keyword,int *int1,int *int2,int required);
char find_keyword_three_integers(char *filename,const char *keyword,int *int1,int *int2,int *int3,int required);
char find_keyword_real(char *filename,const char *keyword,Real *dbl,int required);
char find_keyword_two_reals(char *filename,const char *keyword,Real *dbl1,Real *dbl2,int required);
char find_keyword_three_reals(char *filename,const char *keyword,Real *dbl1,Real *dbl2,Real *dbl3,int required);
char find_keyword_int_real(char *filename,const char *keyword,int *int1,Real *dbl1,int required);
char find_keyword_string(char *filename,const char *keyword,char *string,int required);
char find_keyword_string_caps(char *filename,char upper_flag,const char *keyword,char *string,int required);

#ifdef __cplusplus
}
#endif
#endif /* READ_INPUTS_H */