/**
 * \file    main.c
 * \ingroup utils_group
 * \author  akirby
 * \brief   dg4est data extracts analysis.
 */

/* header files */
#include "analysis_var_inputs.hxx"
#include "analysis_read_inputs.hxx"
#include "analysis_temporal_mean.hxx"
#include "analysis_reynolds_stresses.hxx"
#include "analysis_les_subgridscales.hxx"
#include "analysis_pod_traditional.hxx"
#include "analysis_pod_snapshot.hxx"

int main(int argc,char **argv) {
    MPI_Init(&argc, &argv);
    inputs_t *data = read_inputs(argc,argv);

    switch(data->analysis_type){
        case(ANALYSIS_MEAN):     temporal_mean(data,1); break;
        case(ANALYSIS_REYNOLDS): reynolds_stresses(data,1); break;
        case(ANALYSIS_LES_SGS):  les_subgridscales(data,1); break;
        case(ANALYSIS_POD_TRAD): pod_traditional(data,1); break;
        case(ANALYSIS_POD_SNAP): pod_snapshot(data,1); break;
    }
    inputs_destroy(data);
    return 0;
}