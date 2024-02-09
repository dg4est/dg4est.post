/**
 * File:   analysis_var_inputs.h
 * Author: akirby
 *
 * Created on November 30, 2020, 12:24 PM
 */

#ifndef ANALYSIS_VAR_INPUTS_H
#define ANALYSIS_VAR_INPUTS_H

/* system header files */

#ifdef _OPENMP
#  include <omp.h>
#endif

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

/* header files */
#include "precision_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define eC "\x1B[0m"
#define GC "\x1B[1;32m"
#define gC "\x1B[0;32m"
#define rC "\x1B[0;31m"
#define BUFF_SIZE 1024      /**< @brief string buffer size */

#define OPTIONAL 0
#define REQUIRED 1

#define TO_UPPER 1

#define DISPLAY_OFF 0
#define DISPLAY_ON  1

/* POD variable indicator bit masks */
#define POD_MASK_DENSITY  0b00000001 /**< fluid density POD flag */
#define POD_MASK_U_VEL    0b00000010 /**< fluid x-velocity POD flag */
#define POD_MASK_V_VEL    0b00000100 /**< fluid y-velocity POD flag */
#define POD_MASK_W_VEL    0b00001000 /**< fluid z-velocity POD flag */
#define POD_MASK_PRESSURE 0b00010000 /**< fluid pressure POD flag */
#define MAX_NPODVAR 5                /**< max number of POD Variables */

static
char PODVAR_MASKS[MAX_NPODVAR] = {
    POD_MASK_DENSITY,
    POD_MASK_U_VEL,
    POD_MASK_V_VEL,
    POD_MASK_W_VEL,
    POD_MASK_PRESSURE,
};

static
char PODVAR_NAMES[MAX_NPODVAR][BUFF_SIZE] = {
    "RHO",
    "U",
    "V",
    "W",
    "P"
};

typedef enum {
  ANALYSIS_MEAN     = 1, /**< temporal mean analysis */
  ANALYSIS_REYNOLDS = 2, /**< Reynolds stress analsis */
  ANALYSIS_LES_SGS  = 3, /**< LES SGS mesh resolution analysis */
  ANALYSIS_POD_TRAD = 4, /**< Proper Orthogonal Decomposition analysis: Traditional */
  ANALYSIS_POD_SNAP = 5  /**< Proper Orthogonal Decomposition analysis: Snapshot */
}
analysis_types;

#define TIMER(time,...)      \
  do {                       \
      double t_start,t_end;  \
      t_start = MPI_Wtime(); \
      __VA_ARGS__            \
      t_end = MPI_Wtime();   \
      time += t_end-t_start;  \
  } while(0)

typedef struct {
    char args_used;
    char input_file[BUFF_SIZE];         /**< Input file name */
    char input_data_path[BUFF_SIZE];    /**< Input data file path */
    char output_data_path[BUFF_SIZE];   /**< Output file path */

    char input_base_name[BUFF_SIZE];    /**< Residual statistics output file */
    char input_extension_str[5];        /**< Residual statistics output file */
    int input_extension;

    int plot_tec_sind;
    int plot_tec_eind;
    int plot_tec_ind_skip;

    int ext_ind_skip;
    int ext_sind;
    int ext_eind;
    int file_ind_skip;
    int file_sind;
    int file_eind;

    int multiple_extract_mean;
    int Cartesian2Polar;
    int output_format;
    int output_ascii;
    int output_binary;
    int analysis_type;

    char pod_vars_mask;         /**< POD variables bit-wise masking */
    int npod;                   /**< Number of POD variables */
    int pod_nmode;              /**< Number of POD modes */
    int pod_plot_nmodes;        /**< Number of POD plot modes */
    int pod_plot_ROM;           /**< POD plot reduced order model flag */
    int pod_plot_TVC;           /**< POD plot time-varying modes flag */
    int nstep;                  /**< Number of steps */
}
inputs_t;

static
inputs_t *inputs_new(){
    inputs_t *inputs = (inputs_t *) malloc(sizeof(inputs_t));
    inputs->args_used = 0;
    inputs->input_extension = 1;
    inputs->ext_ind_skip  = 1;
    inputs->ext_sind = 0;
    inputs->ext_eind = 0;
    inputs->file_ind_skip = 1;
    inputs->file_sind = 0;
    inputs->file_eind = 0;
    inputs->multiple_extract_mean = 0;
    inputs->Cartesian2Polar  = 0;
    inputs->output_format = 0;
    inputs->output_ascii = 0;
    inputs->output_binary = 0;
    inputs->analysis_type = ANALYSIS_MEAN;
    inputs->npod = 0;
    inputs->pod_nmode = 0;
    inputs->pod_plot_nmodes = 0;
    inputs->pod_plot_ROM = 0;
    inputs->pod_plot_TVC = 0;
    return inputs;
}

static
void inputs_destroy(inputs_t *inputs){
    if(inputs) free(inputs);
}

static
const char* analysis2string(int type){
    switch(type){
        case(ANALYSIS_MEAN):     return "Temporal Mean";
        case(ANALYSIS_REYNOLDS): return "Reynolds Stress";
        case(ANALYSIS_LES_SGS):  return "LES-SGS Mesh Resolution";
        case(ANALYSIS_POD_TRAD): return "Proper Orthogonal Decomposition - Traditional";
        case(ANALYSIS_POD_SNAP): return "Proper Orthogonal Decomposition - Snapshot";
        default: return "ANALYSIS TYPE NOT RECONGNIZED!";
    }
}

static
void inputs_dump(inputs_t *inputs){
    printf("\n+============== Input Data =============+\n");
    printf(" input_file: %s\n",inputs->input_file);
    printf(" output_data_path: %s\n",inputs->output_data_path);
    printf(" input_data_path: %s\n",inputs->input_data_path);
    printf("\n");
    printf(" extract_index_skip: %d\n",inputs->ext_ind_skip);
    printf(" extract_start_index: %d\n",inputs->ext_sind);
    printf(" extract_end_index: %d\n",inputs->ext_eind);
    printf(" file_index_skip: %d\n",inputs->file_ind_skip);
    printf(" file_start_index: %d\n",inputs->file_sind);
    printf(" file_end_index: %d\n",inputs->file_eind);
    printf(" > Number of Time Steps: %d\n",inputs->nstep);
    printf("\n");
    printf(" analysis_type: %s\n",analysis2string(inputs->analysis_type));
    printf(" average_multiple_extracts: %d\n",inputs->multiple_extract_mean);
    printf(" transform_Cartesian2Polar: %d\n",inputs->Cartesian2Polar);
    printf("+=======================================+\n\n");
}

static
void inputs_example(){
    printf("===============================================================================\n");
    printf("#                       DG4EST EXTRACT ANALYSIS INPUTS                        #\n");
    printf("===============================================================================\n");
    printf("---------------\n");
    printf("# Output Info #\n");
    printf("---------------\n");
    printf("output_format:  1  # Analysis Output File Format:\n");
    printf("                   #   [1] ascii only (tecplot)\n");
    printf("                   #   [2] binary only (analysis)\n");
    printf("                   #   [3] both formats\n");
    printf("plot_tec_start_index: 0\n");
    printf("plot_tec_end_index: 100\n");
    printf("plot_tec_index_skip:  1\n");
    printf("\n");
    printf("--------------\n");
    printf("# Input Info #\n");
    printf("--------------\n");
    printf("input_data_path: WRK/solution/extracts/lines/bin/\n");
    printf("\n");
    printf("extract_start_index: 0  # extract starting id, .e.g. dg4est_extract000_<?>.bin --> 0\n");
    printf("extract_end_index:   5  # extract ending id,   .e.g. dg4est_extract005_<?>.bin --> 5\n");
    printf("extract_index_skip:  1  # interval between extract indices\n");
    printf("\n");
    printf("file_start_index:    0  # file starting id, e.g. dg4est_extract<?>_0000000.bin --> 0\n");
    printf("file_end_index:    100  # file ending id,   e.g. dg4est_extract<?>_0000100.bin --> 100\n");
    printf("file_index_skip:     1  # interval between file indices\n");
    printf("\n");
    printf("------------\n");
    printf("# Analysis #\n");
    printf("------------\n");
    printf("analysis_type: 0  # Analysis Options:\n");
    printf("                  #   [1] Temporal Mean\n");
    printf("                  #   [2] Reynolds Stresses\n");
    printf("                  #   [3] LES SGS Kinetic Energy\n");
    printf("                  #   [4] POD/PCA Traditional\n");
    printf("                  #   [5] POD/PCA Snapshot\n");
    printf("average_multiple_extracts: 0  # Combine and average multiples extracts: e.g. radial averaging\n");
    printf("transform_Cartesian2Polar: 0  # Perform coordinate transformation flag\n");
    printf("\n");
    printf("POD_vars: \"rho,u,v,p\"  # POD Field Options:\n");
    printf("                         #  [rho] fluid density\n");
    printf("                         #  [u,v,w] x/y/z-velocity components\n");
    printf("                         #  [p] fluid pressure\n");
    printf("POD_nmode: 50            # POD Number of reconstruction modes (<1 set to all)\n");
    printf("POD_plot_nmode: 10       # POD Number of plot modes\n");
    printf("POD_plot_ROM: 0          # POD plot reduce order model flag: [0] OFF, [1] ON\n");
    printf("POD_plot_TVC: 0          # POD Number of time-varying coefficients plot modes\n");
}
#ifdef __cplusplus
}
#endif
#endif /* ANALYSIS_VAR_INPUTS_H */