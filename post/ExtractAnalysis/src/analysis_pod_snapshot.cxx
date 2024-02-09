/**
 * \file    analysis_pod_snapshot.c
 * \ingroup utils_group
 * \author  akirby
 *
 * \brief   Computes the Proper Orthogonal Decomposition using the snapshot
 *          method. This function should be used when npts < nsteps.
 */

/* header files */
#include "analysis_pod_snapshot.hxx"

void pod_snapshot(inputs_t *data,char standalone){
    char *filename;
    int ex,i;
    Real aveTime = 0.0;
    Real podTime = 0.0;

    /* setup data structures */
    extracts_t mean_ = extracts_new(); extracts_t *mean = &mean_;
    extracts_t analysis_ = extracts_new(); extracts_t *analysis = &analysis_;

    /* set fields to calculate pod */
    analysis->pod.npod = data->npod;
    analysis->pod.pod_vars_mask = data->pod_vars_mask;

    int pod_count=0;
    for (i = 0; i < MAX_NPODVAR; i++) {
        if(data->pod_vars_mask & PODVAR_MASKS[i]){
            analysis->pod.pod_vars[pod_count] = i;
            analysis->pod.pod_vars_names[pod_count] = PODVAR_NAMES[i];
            pod_count++;
        }
    }

    TIMER(aveTime,
    if (data->multiple_extract_mean) {
        /* compute mean flow over all extracts */
        for (ex = data->ext_sind; ex <= data->ext_eind; ex+=data->ext_ind_skip) {
            compute_extract_mean(mean,data,ex);
        }
    }
    );

    TIMER(podTime,
    for (ex = data->ext_sind; ex <= data->ext_eind; ex+=data->ext_ind_skip) {
        if(!data->multiple_extract_mean) TIMER(aveTime,compute_extract_mean(mean,data,ex););

        /* ====================== */
        /* eigenvalue calculation */
        /* ====================== */
        compute_extract_pod_snap(analysis,mean,data,ex);

        /* write tec POD energy files */
        output_POD_energies(analysis,data,ex,analysis->nStep);

        /* construct Reduce Order Model and/or Time Varying Coefficients */
        if(data->pod_plot_ROM || data->pod_plot_TVC){
            for (i = data->plot_tec_sind; i <= data->plot_tec_eind; i+=data->plot_tec_ind_skip) {
                if(!form_file_name(data,ex,i,&filename,DISPLAY_OFF,"POD RECON")) continue;

                /* reconstruct POD/flow-field */
                compute_pod_snap_reconstruction(analysis,mean,data,ex,i,filename);

                /* write ROM tec file */
                if(data->pod_plot_ROM) output_ROM_quantities(analysis,data,ex,i,
                                                             analysis->extract_dim,
                                                             standalone);
                /* clean up filename variable */
                if(filename) free(filename); filename = NULL;
            }

            /* write TVC file */
            if(data->pod_plot_TVC) output_TVC_quantities(analysis,data,ex);

            /* reset extract arrays */
            if(!data->multiple_extract_mean) mean->mean.fill(0.0);
            analysis->pod.R_tot.fill(0.0);
        }
    }
    );

    /* 3. display wall-clock times */
    printf("+============================================+\n");
    printf(" [POD TRAD] Performance Statistics\n"
           "  Matrix Dimensions:\n"
           "              fields = %d\n"
           "     spatial samples = %d\n"
           "    temporal samples = %d\n"
           "  Wall-Clock Timing: \n"
           "    Soln. Mean [%2.d fields]: %f seconds\n"
           "    Soln. POD  [%2.d fields]: %f seconds\n",
            mean->nf,mean->npts,analysis->nStep,
            mean->nf,aveTime,
            data->npod,podTime);
    printf("+============================================+\n");

    /* 4. clean up memory */
    extracts_destroy(mean);
    extracts_destroy(analysis);

}

void compute_extract_pod_snap(extracts_t *analysis,extracts_t *means,inputs_t *data,int ex){
    Real flucTime=0.0;
    Real covaTime=0.0;
    Real eignTime=0.0;
    char *filename;
    int istep;
    int i;

    /* compute unbiased factor */
    analysis->nStep = (int) (data->file_eind-data->file_sind+data->file_ind_skip)/(Real)(data->file_ind_skip);
    analysis->oneOnStepMone = 1.0 / ((Real)analysis->nStep  - 1.0);

    /* ============================ */
    /* auto-correlation calculation */
    /* ============================ */
#ifdef _OPENMP
    char *omp_nthread = getenv("OMP_NUM_THREADS");
    if(omp_nthread){
        printf("Number of OpenMP Threads: %d\n",atoi(getenv("OMP_NUM_THREADS")));
    }
#endif
    for (istep=0,i = data->file_sind; i <= data->file_eind; i+=data->file_ind_skip,istep++) {
        if(!form_file_name(data,ex,i,&filename,DISPLAY_OFF,"FLUCTUATION")) continue;

        /* display status */
        Real progress = (istep/(Real) analysis->nStep)*100.0;
        printf("\r[COVARIANCE] Calculation in progress (%d vars): %.1f%%",analysis->pod.npod,progress);

        TIMER(flucTime,compute_snap_fluctuations(analysis,means,filename););
        TIMER(covaTime,compute_snap_covariance(analysis,istep););

        /* clean up filename variable */
        if(filename) free(filename); filename = NULL;
    }
    printf("\n");

    /* ========================== */
    /* compute eigenvalue problem */
    /* ========================== */
    TIMER(eignTime,compute_snap_eigenproblem(analysis););
    printf("[POD] w-time (seconds):\n"
           "   <fluctuation> %f sec,\n"
           "    <covariance> %f sec,\n"
           "  <eigenproblem> %f sec \n",
            flucTime,covaTime,eignTime);
}

void compute_pod_snap_reconstruction(extracts_t *analysis,extracts_t *means,inputs_t *data,
                                     int ex,int index,char *filename){
    /* re-compute solution fluctuation */
    compute_snap_fluctuations(analysis,means,filename);

    /* reconstruct flow solution */
    compute_snap_reconstruction(analysis,data,ex,index);
}

void compute_snap_fluctuations(extracts_t *analysis,extracts_t *means,char *file_name){
    int fill_coordinates = 0;
    const int nhead = 6;
    int sizes[nhead];
    int pt;

    /* Header Info:
     *   sizes[0] = plot_dim; // plotting dimension
     *   sizes[1] = PDSIZE;   // number of fields per point
     *   sizes[2] = npts;     // number of points
     *   sizes[3] = nx;
     *   sizes[4] = ny;
     *   sizes[5] = nz;
     */

    /* ================= */
    /* read extract data */
    /* ================= */
    FILE *fpb = fopen(file_name, "rb");
    size_t ret;

    /* read header info */
    ret = fread(sizes,sizeof(int),nhead,fpb);

    analysis->extract_dim = sizes[0];
    analysis->ncf = sizes[1];
    analysis->npts = sizes[2];
    analysis->nx = sizes[3];
    analysis->ny = sizes[4];
    analysis->nz = sizes[5];
    analysis->nf = FL_NUM;

    /* =============== */
    /* allocate memory */
    /* =============== */
    /* field data */
    analysis->data.newmalloc(analysis->ncf*analysis->npts);
    fill_coordinates = analysis->coordinates.newmalloc(3*analysis->npts);

    /* fluctuation data */
    analysis->U.newmalloc(FL_NUM*analysis->npts);
    analysis->UM.newmalloc(FL_NUM*analysis->npts);
    analysis->dU.newmalloc(FL_NUM*analysis->npts,0.0);

    /* ===================== */
    /* read point field data */
    /* ===================== */
    ret = fread(analysis->data.ptr(),sizeof(Real),analysis->ncf*analysis->npts,fpb);

    /* fill coordinates */
    if(fill_coordinates) {
        for (pt = 0; pt < analysis->npts; ++pt) {
            Real *co = &analysis->coordinates[3*pt];
            Real *u  = &analysis->data[analysis->ncf*pt];

            co[0] = u[0]; // x
            co[1] = u[1]; // y
            co[2] = u[2]; // z
        }
    }

    /* close file */
    fclose(fpb);

    /* =========================== */
    /* calculate fluctuation means */
    /* =========================== */
    if(analysis->ncf != EX_NUM){
        printf("ERROR: ncf(%d)!=EX_NUM(%d) [%s, line %d]\n",
                analysis->ncf,EX_NUM,
                __FILE__, __LINE__);
    }
    for (pt = 0; pt < analysis->npts; ++pt) {
        Real *U  = &analysis->data[EX_NUM*pt];
        Real *UM = &means->mean[MEAN_NUM*pt];

        Real *dU = &analysis->dU[FL_NUM*pt];
        Real *ULoc = &analysis->U[FL_NUM*pt];
        Real *UAve = &analysis->UM[FL_NUM*pt];

        /* stash U values for fluctuating data */
        ULoc[FL_RHO]  = U[EX_RHO];
        ULoc[FL_U]    = U[EX_U];
        ULoc[FL_V]    = U[EX_V];
        ULoc[FL_W]    = U[EX_W];
        ULoc[FL_PRES] = U[EX_PRES];

        UAve[FL_RHO]  = UM[MEAN_RHO];
        UAve[FL_U]    = UM[MEAN_U];
        UAve[FL_V]    = UM[MEAN_V];
        UAve[FL_W]    = UM[MEAN_W];
        UAve[FL_PRES] = UM[MEAN_PRES];

        /* NOTE: U/UM = [rho,u,v,w,p] */
        dU[FL_RHO]  = U[EX_RHO]  - UM[MEAN_RHO]; // rho' = rho - <rho>
        dU[FL_U]    = U[EX_U]    - UM[MEAN_U];   // u' = u - <u>
        dU[FL_V]    = U[EX_V]    - UM[MEAN_V];   // v' = v - <v>
        dU[FL_W]    = U[EX_W]    - UM[MEAN_W];   // w' = w - <w>
        dU[FL_PRES] = U[EX_PRES] - UM[MEAN_PRES];// p' = p - <p>
    }
}

void compute_snap_reconstruction(extracts_t *analysis,inputs_t *data,int iextract,int istep){
    pod_t *pod = &analysis->pod;
    Real * const U  = analysis->U.ptr();
    Real * const UM = analysis->UM.ptr();
    Real * const dU = analysis->dU.ptr();
    Real * const coordinates = analysis->coordinates.ptr();

    /* check pod_nmode input variable */
    if(data->pod_nmode < 1 || data->pod_nmode > analysis->nStep){
        data->pod_nmode = analysis->nStep;
    }

    const int NPOD = pod->npod;
    const int NPTS = analysis->npts;
    const int NSTEP = analysis->nStep;
    const int NMODES = data->pod_nmode;
    const int PLOT_NMODES = data->pod_plot_nmodes;
    const int nfield = FL_NUM;
    int i;

    /* Reconstruct all energies at this time state
     * time varying modal coefficients: a_p(nModes)
     * NOTE: These are the time varying coefficients because
     *       we are projecting onto the i-th (in time) instantaneous
     *       flow solution. Plot the first m-coefficients in time by
     *       appending the time varying coefficients per mode to a
     *       file for each instantaneous flow solution.
     */

    /* allocate memory (1st pass) */
    pod->mode.newmalloc(NPTS*NPOD);
    pod->UROM.newmalloc(NPTS*NPOD);
    pod->R_a_p.newmalloc(NSTEP*NPOD);
    pod->TVC.newmalloc(NMODES*NSTEP*NPOD,0.0);

    pod->mode.fill(0.0);
    pod->UROM.fill(0.0);
    pod->R_a_p.fill(0.0);
    for (i = 0; i < NPOD; i++) {
        const Real * const R      = &pod->R_tot[NPTS*NSTEP*i];
        const Real * const eigvec = &pod->eig[NSTEP*i];
        Real * const Rap          = &pod->R_a_p[NSTEP*i];
        Real * const modes        = &pod->mode[NPTS*i];
        Real * const UROM         = &pod->UROM[NPTS*i];
        Real * const TVC_STEPS    = &pod->TVC[NMODES*NSTEP*i]; //@(:,:,i)
        Real * const TVC          = &TVC_STEPS[NMODES*istep];  //@(:,istep)

        const Real *energy = &pod->pod_total_energy[pod->pod_vars[i]];
        const int pod_var = pod->pod_vars[i]+FBASE; // FORTRAN indexing (1)

        // Rap, TVC, UROM, mode
        if(data->pod_plot_ROM) printf("[ROM-%3s]: Step %d",
                                      pod->pod_vars_names[i],istep);
        POD_snap_reconstruct_solution_(&iextract,&istep,&pod_var,&nfield,
                                       &analysis->extract_dim,&NPTS,
                                       &analysis->nx,&analysis->ny,&analysis->nz,
                                       &NSTEP,&NMODES,&PLOT_NMODES,
                                       energy,eigvec,coordinates,
                                       U,UM,dU,
                                       R,Rap,TVC,
                                       modes,UROM,
                                       &data->pod_plot_TVC,
                                       &data->pod_plot_ROM);
    }
}

void compute_snap_covariance(extracts_t *analysis,int istep){
    pod_t *pod = &analysis->pod;
    Real * const dU = analysis->dU.ptr();

    const int NPOD = pod->npod;
    const int NPTS = analysis->npts;
    const int NSTEP = analysis->nStep;
    const int nfield = FL_NUM;
    int i;

    /* allocate memory (1st pass) */
    pod->Y_data.newmalloc(NPTS*NSTEP*NPOD);

    /* =========================== */
    /* copy fluctuations to Y_data */
    /* =========================== */
    for (i = 0; i < NPOD; i++) {
        Real * const Y_STEPS = &pod->Y_data[NPTS*NSTEP*i];
        Real * const Y = &Y_STEPS[NPTS*istep];
        const int pod_var = pod->pod_vars[i]+FBASE; // add FORTRAN Base-1 index

        // copy operator
        POD_snap_covariance_(&pod_var,&nfield,&NPTS,dU,Y);
    }
}

void compute_snap_eigenproblem(extracts_t *analysis){
    pod_t *pod = &analysis->pod;
    const int NPOD = pod->npod;
    const int NPTS = analysis->npts;
    const int NSTEP = analysis->nStep;
    const Real nStepMone = (Real)NSTEP-1.0;
    const Real oneOnStepMone = 1.0/nStepMone;
    int i;

    /* allocate memory (1st pass) */
    pod->R_tot.newmalloc(NPTS*NSTEP*NPOD,0.0);
    pod->eig.newmalloc(NSTEP*NPOD,0.0);
    pod->K_T.newmalloc(NSTEP*NSTEP,0.0); // workspace -- not NPOD

    /* ==================================== */
    /* calculation eigenvalues and vectors  */
    /* ==================================== */
    for (i = 0; i < NPOD; i++) {
        Real *eigval = &pod->eig[NSTEP*i];   // get address to POD variable
        Real *R = &pod->R_tot[NPTS*NSTEP*i]; // get address to POD variable
        Real *Y = &pod->Y_data[NPTS*NSTEP*i];// get address to POD variable
        Real *KT = pod->K_T.ptr();            // K_T workspace address

        int pod_var = pod->pod_vars[i];
        Real *energy_total = &pod->pod_total_energy[pod_var];

        POD_snap_solve_(&NPTS,&NSTEP,&nStepMone,&oneOnStepMone,KT,Y,R,eigval,energy_total);
        printf("TOTAL ENERGY[%s]=%f\n",pod->pod_vars_names[i],*energy_total);
    }
}