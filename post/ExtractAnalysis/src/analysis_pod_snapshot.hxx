/**
 * File:   analysis_pod_snapshot.h
 * Author: akirby
 *
 * Created on January 21, 2024, 5:13 PM
 */

#ifndef ANALYSIS_POD_SNAPSHOT_H
#define ANALYSIS_POD_SNAPSHOT_H

/* system header files */
#include <math.h>
#include <string.h>
#include <assert.h>

/* header files */
#include "analysis_read_inputs.hxx"
#include "analysis_var_inputs.hxx"
#include "analysis_var_extracts.hxx"
#include "analysis_temporal_mean.hxx"
#include "analysis_reynolds_stresses.hxx"
#include "analysis_pod_output.hxx"

#ifdef __cplusplus
extern "C" {
#endif

/* ================== */
/* External Functions */
/* ================== */
/** External covariance calculation function.
 *
 * @param [in]    nfield            number of variables at each point
 * @param [in]    npt               number of points
 * @param [in]    oneOnStepMone     statistics scaling factor: 1/(nStep-1)
 * @param [in]    dU                [nfield][npt]: fluctuations
 * @param [inout] R_tot             [nfield][npt][npt]: covariance matrix
 */
void POD_snap_covariance_(const int *pod_var,const int *nfield,const int *npt,
                          const Real *dU,Real *Y_data);

/** External eigenvalue and eigenvector solve function.
 *
 * @param npt
 * @param nstep
 * @param nStepMone
 * @param oneOnStepMone
 * @param KT
 * @param Y
 * @param R
 * @param eigval
 * @param energy_total
 */
void POD_snap_solve_(const int *npt,const int *nstep,
                     const Real *nStepMone,const Real *oneOnStepMone,
                     Real *KT,Real *Y,Real *R,Real *eigval,Real *energy_total);

void POD_snap_reconstruct_solution_(const int *iextract,const int *index,
                                    const int *pod_var,const int *nfield,
                                    const int *plot_dim,const int *npt,
                                    const int *nx,const int *ny,const int *nz,
                                    const int *nstep,const int *nmode,
                                    const int *plot_nmode,const Real *energy_total,
                                    const Real *eigenvalues,const Real *coordinates,
                                    const Real *U,const Real *Ubar,const Real *dU,
                                    const Real *R_tot,Real *R_a_p,Real *TVC,
                                    Real *modes, Real *UROM,
                                    const int *plot_TVC,const int *plot_ROM);
/* ================ */
/* Public Functions */
/* ================ */
void pod_snapshot(inputs_t *data,char standalone);
void compute_extract_pod_snap(extracts_t *fluctuations,extracts_t *means,inputs_t *data,int ex);
void compute_pod_snap_reconstruction(extracts_t *analysis,extracts_t *means,inputs_t *data,int ex,int index,char *filename);
void compute_snap_fluctuations(extracts_t *fluctuations,extracts_t *means,char *file_name);
void compute_snap_reconstruction(extracts_t *analysis,inputs_t *data,int iextract,int istep);
void compute_snap_covariance(extracts_t *analysis,int istep);
void compute_snap_eigenproblem(extracts_t *analysis);

#ifdef __cplusplus
}
#endif
#endif /* ANALYSIS_POD_SNAPSHOT_H */