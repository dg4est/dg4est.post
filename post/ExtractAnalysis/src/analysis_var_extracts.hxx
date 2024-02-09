/**
 * File:   analysis_var_extracts.h
 * Author: akirby
 *
 * Created on November 21, 2020, 12:25 PM
 */

#ifndef ANALYSIS_VAR_EXTRACTS_H
#define ANALYSIS_VAR_EXTRACTS_H

/* system header files */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

/* header files */
#include "precision_types.h"
#include "memory.hxx"
#include "analysis_var_inputs.hxx"

#ifdef __cplusplus
extern "C" {
#endif

/* FORTRAN Base Index */
#define FBASE 1

/* Output Mode */
#define LOCAL_MODE 0
#define GLOBAL_MODE 1

/* Extract Mode */
#define PT_MODE 0
#define LINE_MODE 1
#define PLANE_MODE 2
#define VOLUME_MODE 3

/* Extracts indexes */
#define EX_X     0
#define EX_Y     1
#define EX_Z     2
#define EX_RHO   3
#define EX_U     4
#define EX_V     5
#define EX_W     6
#define EX_PRES  7
#define EX_MACH  8
#define EX_VORT  9
#define EX_QCRT  10
#define EX_IBLK  11
#define EX_MUSGS 12
#define EX_LSGS  13
#define EX_KESGS 14
#define EX_NUM   15

/* Extracts indexes */
#define MEAN_RHO   (EX_RHO  -3)
#define MEAN_U     (EX_U    -3)
#define MEAN_V     (EX_V    -3)
#define MEAN_W     (EX_W    -3)
#define MEAN_PRES  (EX_PRES -3)
#define MEAN_MACH  (EX_MACH -3)
#define MEAN_VORT  (EX_VORT -3)
#define MEAN_QCRT  (EX_QCRT -3)
#define MEAN_IBLK  (EX_IBLK -3)
#define MEAN_MUSGS (EX_MUSGS-3)
#define MEAN_LSGS  (EX_LSGS -3)
#define MEAN_KESGS (EX_KESGS-3)
#define MEAN_NUM   (EX_NUM  -3)

/* Reynolds Stress indexes */
#define RE_UU 0 /* u'u' */
#define RE_VV 1 /* u'v' */
#define RE_WW 2 /* w'w' */
#define RE_UV 3 /* u'v' */
#define RE_UW 4 /* u'w' */
#define RE_VW 5 /* v'w' */
#define RE_KE 6 /* kinetic energy */
#define RE_IBLK 7 /* iblank */
#define RE_NUM 8

/* Fluctuation indexes */
#define FL_RHO  0 /* density' */
#define FL_U    1 /* u' */
#define FL_V    2 /* v' */
#define FL_W    3 /* w' */
#define FL_PRES 4 /* pressure' */
#define FL_NUM  5

/* TKE indexes */
#define TKE_RES  0 /* resolved tke */
#define TKE_SGS  1 /* subgrid scale tke */
#define TKE_PERC 2 /* percentage: (tke_res) / (tke_res + tke_sgs) */
#define TKE_PRSF 3 /* pressure fluctuation */
#define TKE_PRMS 4 /* root-mean-square pressure fluctation */
#define TKE_IBLK 5 /* iblank */
#define TKE_NUM  6

typedef struct {
    char pod_vars_mask;         /* POD variables bit-wise masking */
    int npod;                           /* number of POD variables read */
    int pod_vars[MAX_NPODVAR];          /* POD variable list */
    char *pod_vars_names[MAX_NPODVAR];  /* POD variable list */
    Real pod_total_energy[MAX_NPODVAR]; /* POD total energies */

    dg::memory<Real> R_tot;
    dg::memory<Real> TVC;      /* [POD_nModes,nStep] */
    dg::memory<Real> UROM;     /* [nElem] */
    dg::memory<Real> K_T;      /* [nStep,nStep] */
    dg::memory<Real> R_a_p;    /* Traditional: [nElem], Snap: [nStep] */
    dg::memory<Real> Y_data;   /* [nElem,nStep] */
    dg::memory<Real> lambda;   /* [nStep] */
    dg::memory<Real> mode;     /* [nElem] */
    dg::memory<Real> eig;      /* [nElem] */
    dg::memory<Real> energies; /* [npod] */
}
pod_t;

typedef struct {
    int extract_dim;/**< extraction dimension */
    int npts;       /**< number of points */
    int ncf;        /**< number of coordinates and fields per point */
    int nf;         /**< number of fields per point */
    int nx;         /**< number x-direction points */
    int ny;         /**< number y-direction points */
    int nz;         /**< number z-direction points */
    int nStep;      /**< number of time steps */
    Real oneOnStepMone; /**< 1.0/(nStep-1) */

    dg::memory<Real> U;      /**< extract solution data */
    dg::memory<Real> UM;     /**< extract mean solution data */
    dg::memory<Real> dU;     /**< extract fluctuation data */
    dg::memory<Real> data;     /**< extract data */
    dg::memory<Real> tke;      /**< extract tke data */
    dg::memory<Real> mean;     /**< extract mean data */
    dg::memory<Real> stddev;   /**< extract standard deviation data */
    dg::memory<Real> stresses; /**< extract reynolds stresses data */
    dg::memory<Real> coordinates; /**< extract coordinates */
    dg::memory<Real> coordinates_polar; /**< extract polar coordinates */
    pod_t pod;  /**< POD data */
}
extracts_t;

static
extracts_t extracts_new(){
    extracts_t ext;
    ext.ncf = 0;
    ext.nf = 0;
    ext.nx = 0;
    ext.ny = 0;
    ext.nz = 0;
    ext.npts = 0;
    ext.extract_dim = 0;
    return ext;
}

static
void extracts_destroy(extracts_t *ext){
    if(ext){
        ext->U.free();
        ext->UM.free();
        ext->dU.free();

        ext->data.free();
        ext->tke.free();
        ext->mean.free();
        ext->stddev.free();
        ext->stresses.free();
        ext->coordinates.free();
        ext->coordinates_polar.free();

        // pod variables
        ext->pod.R_tot.free();
        ext->pod.TVC.free();
        ext->pod.UROM.free();
        ext->pod.K_T.free();
        ext->pod.R_a_p.free();
        ext->pod.Y_data.free();
        ext->pod.lambda.free();
        ext->pod.mode.free();
        ext->pod.eig.free();
        ext->pod.energies.free();
    }
}

#ifdef __cplusplus
}
#endif
#endif /* ANALYSIS_VAR_EXTRACTS_H */