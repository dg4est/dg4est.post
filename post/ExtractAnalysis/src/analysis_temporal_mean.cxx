/**
 * \file    analysis_temporal_mean.c
 * \ingroup utils_group
 * \author  akirby
 *
 * \brief   Computes the temporal mean field from dg4est extract data.
 */

/* header files */
#include "analysis_temporal_mean.hxx"

void temporal_mean(inputs_t *data,char standalone){
    Real oneOts;
    int ex;

    /* setup data structures */
    extracts_t stats_ = extracts_new(); extracts_t *stats = &stats_;

    if (data->multiple_extract_mean) {
        /* compute mean flow over all extracts */
        for (ex = data->ext_sind; ex <= data->ext_eind; ex+=data->ext_ind_skip) {
            compute_extract_mean(stats,data,ex);
        }
    }

    for (ex = data->ext_sind; ex <= data->ext_eind; ex+=data->ext_ind_skip) {
        if(!data->multiple_extract_mean) compute_extract_mean(stats,data,ex);
        oneOts = compute_extract_stddev(stats,data,ex);

        /* write out extract mean file */
        if (!data->multiple_extract_mean) {
            output_mean(stats,data,ex,ex,LOCAL_MODE,stats->extract_dim,standalone);
            output_statistics(stats,data,ex,ex,LOCAL_MODE,stats->extract_dim,standalone,oneOts);

            /* reset extract arrays */
            stats->mean.fill(0.0);
            stats->stddev.fill(0.0);
        }
    }

    /* write out extract mean file */
    if (data->multiple_extract_mean) {
        output_mean(stats,data,data->ext_sind,data->ext_eind,GLOBAL_MODE,stats->extract_dim,standalone);
        output_statistics(stats,data,data->ext_sind,data->ext_eind,GLOBAL_MODE,stats->extract_dim,standalone,oneOts);
    }

    /* 3. clean up memory */
    extracts_destroy(stats);
}

void compute_extract_mean(extracts_t *stats,inputs_t *data,int ex){
    char *filename;
    int i;

    int nfile = (data->file_eind - data->file_sind + data->file_ind_skip) / data->file_ind_skip;
    int nextr = (data->ext_eind - data->ext_sind + data->ext_ind_skip) / data->ext_ind_skip;

    Real oneOts = (data->multiple_extract_mean) ? (1.0 / ((Real) (nfile*nextr))):
                                                  (1.0 / ((Real) (nfile)));
    /* =================== */
    /* average calculation */
    /* =================== */
    for (i = data->file_sind; i <= data->file_eind; i+=data->file_ind_skip) {
        if(!form_file_name(data,ex,i,&filename,DISPLAY_OFF,"AVERAGE")) continue;

        compute_mean(stats,data,filename,oneOts);

        /* clean up filename variable */
        if(filename) free(filename); filename = NULL;
    }
}

Real compute_extract_stddev(extracts_t *stats,inputs_t *data,int ex){
    char *filename;
    int i;

    int nfile = (data->file_eind - data->file_sind + data->file_ind_skip) / data->file_ind_skip;
    int nextr = (data->ext_eind - data->ext_sind + data->ext_ind_skip) / data->ext_ind_skip;

    Real oneOts = (data->multiple_extract_mean) ? (1.0 / ((Real) (nfile*nextr))):
                                                  (1.0 / ((Real) (nfile)));

    /* ======================= */
    /* fluctuation calculation */
    /* ======================= */
    for (i = data->file_sind; i <= data->file_eind; i+=data->file_ind_skip) {
        if(!form_file_name(data,ex,i,&filename,DISPLAY_OFF,"STDDEV")) continue;

        compute_stddev(stats,filename,oneOts);

        /* clean up filename variable */
        if(filename) free(filename); filename = NULL;
    }
    compute_stddev_final(stats);
    return oneOts;
}

void compute_mean(extracts_t *stats,inputs_t *inputs,char *file_name,Real frac){
    int fill_coordinates = 0;
    const int nhead = 6;
    int sizes[nhead];
    int f,pt;

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
    ret = fread(sizes, sizeof(int), nhead, fpb);

    stats->extract_dim = sizes[0];
    stats->ncf = sizes[1];
    stats->npts = sizes[2];
    stats->nx = sizes[3];
    stats->ny = sizes[4];
    stats->nz = sizes[5];
    stats->nf = MEAN_NUM; // subtract number of coordinates
#if 0
    printf("Header Info:\n"
           "  Dim : %d\n"
           "  ncf : %d\n"
           "  npts: %d\n"
           "  nx  : %d\n"
           "  ny  : %d\n"
           "  nz  : %d\n",
           sizes[0],sizes[1],sizes[2],sizes[3],sizes[4],sizes[5]);
#endif

    /* allocate array data */
    stats->data.newmalloc(stats->ncf*stats->npts);
    stats->mean.newmalloc(stats->nf*stats->npts,0.0);
    fill_coordinates = stats->coordinates.newmalloc(3*stats->npts);

    /* read point data */
    ret = fread(stats->data.ptr(),sizeof(Real),stats->ncf*stats->npts,fpb);

    /* fill coordinates */
    if(fill_coordinates) {
        for (pt = 0; pt < stats->npts; ++pt) {
            Real *co = &stats->coordinates[3*pt];
            Real *u  = &stats->data[stats->ncf*pt];

            co[0] = u[0]; // x
            co[1] = u[1]; // y
            co[2] = u[2]; // z
        }
    }

    /* close bin file */
    fclose(fpb);

    if (fill_coordinates && inputs->Cartesian2Polar) {
        stats->coordinates_polar.newmalloc(3*stats->npts);

        for (pt = 0; pt < stats->npts; ++pt) {
            Real *co = &stats->coordinates[3*pt];
            Real *pc = &stats->coordinates_polar[3*pt];

            // FIXME: r, theta, z
            pc[0] = co[0]; // x
            pc[1] = co[1]; // y
            pc[2] = co[2]; // z
        }
    }

    /* calculate mean */
    for (pt = 0; pt < stats->npts; ++pt) {
        Real *u  = &stats->data[stats->ncf*pt+3]; // +3 to skip coordinates;
        Real *m  = &stats->mean[stats->nf*pt];
        for (f = 0; f < stats->nf; ++f) {
            m[f] += frac*u[f];
        }
    }
}

void compute_stddev(extracts_t *stats,char *file_name,Real frac){
    int fill_coordinates = 0;
    const int nhead = 6;
    int sizes[nhead];
    int pt,f;

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

    stats->extract_dim = sizes[0];
    stats->ncf = sizes[1];
    stats->npts = sizes[2];
    stats->nx = sizes[3];
    stats->ny = sizes[4];
    stats->nz = sizes[5];
    stats->nf = MEAN_NUM;  // std. dev. for each field

    /* =============== */
    /* allocate memory */
    /* =============== */
    stats->data.newmalloc(stats->ncf*stats->npts);
    stats->stddev.newmalloc(MEAN_NUM*stats->npts,0.0);
    fill_coordinates = stats->coordinates.newmalloc(3*stats->npts);

    /* ===================== */
    /* read point field data */
    /* ===================== */
    ret = fread(stats->data.ptr(),sizeof(Real),stats->ncf*stats->npts,fpb);

    /* fill coordinates */
    if(fill_coordinates) {
        for (pt = 0; pt < stats->npts; ++pt) {
            Real *co = &stats->coordinates[3*pt];
            Real *u  = &stats->data[stats->ncf*pt];

            co[0] = u[0]; // x
            co[1] = u[1]; // y
            co[2] = u[2]; // z
        }
    }

    /* close file */
    fclose(fpb);

    /* ============================ */
    /* calculate standard deviation */
    /* ============================ */
    for (pt = 0; pt < stats->npts; ++pt) {
        Real *U  = &stats->data[EX_NUM*pt+3]; // +3 to skip coordinates (x,y,z)
        Real *UM = &stats->mean[MEAN_NUM*pt];
        Real *stddev = &stats->stddev[MEAN_NUM*pt];

        for (f = 0; f < MEAN_NUM; ++f) {
            Real fluc = U[f] - UM[f];
            stddev[f] += frac*(fluc*fluc); // sum((u'u')/N)
        }
    }
    // sqrt is applied after all extracts
}

void compute_stddev_final(extracts_t *stats){
    int pt,f;

    /* =========================== */
    /* finalize standard deviation */
    /* =========================== */
    for (pt = 0; pt < stats->npts; ++pt) {
        Real *stddev = &stats->stddev[MEAN_NUM*pt];

        for (f = 0; f < MEAN_NUM; ++f) {
            stddev[f] = sqrt(stddev[f]);
        }
    }
}

void output_mean(extracts_t *ext,inputs_t *inputs,int sindex,int eindex,int gl_mode,int mode,char standalone){
    int nx,ny,nz;
    int plot_dim;
    int pt,i;

    /* ================== */
    /* write data to file */
    /* ================== */
    char file_ascii[48] = { '\0' };
    char file_binary[48] = { '\0' };

    if (gl_mode == LOCAL_MODE) {
        snprintf(file_ascii,43, "/tec/extract%03d.Mean.tec",sindex);
        snprintf(file_binary,43,"/bin/extract%03d.Mean.bin",sindex);
    } else {
        snprintf(file_ascii,48, "/tec/extracts%03d-%03d.Mean.tec",sindex,eindex);
        snprintf(file_binary,48,"/bin/extracts%03d-%03d.Mean.bin",sindex,eindex);
    }

    char *file_name_ascii = (char*) malloc(strlen(inputs->output_data_path) + strlen(file_ascii) + 1);
    char *file_name_binary = (char*) malloc(strlen(inputs->output_data_path) + strlen(file_binary) + 1);

    strcpy(file_name_ascii,inputs->output_data_path); strcpy(file_name_binary,inputs->output_data_path);
    strcat(file_name_ascii,file_ascii);               strcat(file_name_binary,file_binary);

    /* fill array */
    for (pt = 0; pt < ext->npts; pt++) {
        Real *arr = &ext->data[ext->ncf*pt];
        Real *co = &ext->coordinates[3*pt];
        Real *p = &ext->mean[ext->nf*pt];

        arr[0] = co[0]; arr[1] = co[1]; arr[2] = co[2];
        for(i = 0; i < ext->nf; i++) arr[i+3] = p[i];
    }

    /* ---------------- */
    /* write ascii data */
    /* ---------------- */
    if (inputs->output_ascii || !standalone) {
        FILE *fpa = fopen(file_name_ascii, "w");
        fprintf(fpa,"TITLE = \"DG4EST Extractions\"\n");
        fprintf(fpa,"VARIABLES = \"x\", \"y\", \"z\", "
                    "\"rho\", \"u\", \"v\", \"w\", "
                    "\"pressure\", \"mach\", "
                    "\"vorticity\", \"qcriterion\", \"iblank\","
                    "\"mu_sgs\", \"l_sgs\", \"tke_sgs\"\n");

        nx = ny = nz = 1;
        if (mode == PT_MODE || mode == LINE_MODE) {
            plot_dim = ext->extract_dim;
            nx = ext->nx;
            fprintf(fpa,"ZONE T=\"Zone 1\", I=%d, F=POINT\n",ext->npts);
        } else if (mode == PLANE_MODE) {
            plot_dim = 2;
            nx = ext->nx;
            ny = ext->ny;
            fprintf(fpa,"ZONE T=\"Zone 1\", I=%d, J=%d, F=POINT\n",nx,ny);
        } else if (mode == VOLUME_MODE) {
            plot_dim = 3;
            nx = ext->nx;
            ny = ext->ny;
            nz = ext->nz;
            fprintf(fpa,"ZONE T=\"Zone 1\", I=%d, J=%d, K=%d, F=POINT\n",nx,ny,nz);
        }

        /* write point data */
        for (pt = 0; pt < ext->npts; pt++) {
            Real *p = &ext->data[ext->ncf*pt];
            for (i = 0; i < ext->ncf; i++) {
                (p[i]<0.0) ? fprintf(fpa,"%.15f ",p[i]):fprintf(fpa," %.15f ",p[i]);
            }
            fprintf(fpa,"\n");
        }
        /* close tec file */
        fclose(fpa);
        DGOUT(stdout,"[analysis] Output tecplot ascii mean file: %s\n",file_name_ascii);
    }

    /* ----------------- */
    /* write binary data */
    /* ----------------- */
    if (inputs->output_binary || !standalone) {
        FILE *fpb = fopen(file_name_binary, "wb");

        const int nhead = 6;
        int sizes[nhead];

        sizes[0] = plot_dim;
        sizes[1] = ext->ncf;
        sizes[2] = ext->npts;
        sizes[3] = ext->nx;
        sizes[4] = ext->ny;
        sizes[5] = ext->nz;

        /* write size data */
        fwrite(sizes,sizeof(int),nhead,fpb);

        /* write point data */
        fwrite(ext->data.ptr(),sizeof(Real),ext->ncf*ext->npts,fpb);

        /* close bin file */
        fclose(fpb);
        DGOUT(stdout,"[analysis] Output binary mean file: %s\n",file_name_binary);
    }

    /* clean up memory */
    if(file_name_ascii) {free(file_name_ascii);} file_name_ascii = NULL;
    if(file_name_binary) {free(file_name_binary);} file_name_binary = NULL;
}

void output_statistics(extracts_t *ext,inputs_t *inputs,int sindex,int eindex,int gl_mode,int mode,char standalone,Real frac){
    int nx,ny,nz;
    int plot_dim;
    int pt,i;

    /* ================== */
    /* write data to file */
    /* ================== */
    char file_ascii[49] = { '\0' };

    if (gl_mode == LOCAL_MODE) {
        snprintf(file_ascii,44, "/tec/extract%03d.Stats.tec",sindex);
    } else {
        snprintf(file_ascii,49, "/tec/extracts%03d-%03d.Stats.tec",sindex,eindex);
    }
    char *file_name_ascii = (char*) malloc(strlen(inputs->output_data_path) + strlen(file_ascii) + 1);
    strcpy(file_name_ascii,inputs->output_data_path);
    strcat(file_name_ascii,file_ascii);

    /* ---------------- */
    /* write ascii data */
    /* ---------------- */
    if (inputs->output_ascii || !standalone) {
        FILE *fpa = fopen(file_name_ascii, "w");
        fprintf(fpa,"TITLE = \"DG4EST Extraction Statistics: mean, standard deviation, std. error of mean\"\n");
        fprintf(fpa,"VARIABLES = \"x\", \"y\", \"z\", "
                "\"rho_mean\", \"u_mean\", \"v_mean\", \"w_mean\", \"pressure_mean\", \"mach_mean\", \"vorticity_mean\", \"qcriterion_mean\", \"iblank_mean\", \"mu_sgs_mean\", \"l_sgs_mean\", \"tke_sgs_mean\", "
                "\"rho_sd\", \"u_sd\", \"v_sd\", \"w_sd\", \"pressure_sd\", \"mach_sd\", \"vorticity_sd\", \"qcriterion_sd\", \"iblank_sd\", \"mu_sgs_sd\", \"l_sgs_sd\", \"tke_sgs_sd\", "
                "\"rho_err\", \"u_err\", \"v_err\", \"w_err\", \"pressure_err\", \"mach_err\", \"vorticity_err\", \"qcriterion_err\", \"iblank_err\", \"mu_sgs_err\", \"l_sgs_err\", \"tke_sgs_err\""
                "\n");

        nx = ny = nz = 1;
        if (mode == PT_MODE || mode == LINE_MODE) {
            plot_dim = ext->extract_dim;
            nx = ext->nx;
            fprintf(fpa,"ZONE T=\"Zone 1\", I=%d, F=POINT\n",ext->npts);
        } else if (mode == PLANE_MODE) {
            plot_dim = 2;
            nx = ext->nx;
            ny = ext->ny;
            fprintf(fpa,"ZONE T=\"Zone 1\", I=%d, J=%d, F=POINT\n",nx,ny);
        } else if (mode == VOLUME_MODE) {
            plot_dim = 3;
            nx = ext->nx;
            ny = ext->ny;
            nz = ext->nz;
            fprintf(fpa,"ZONE T=\"Zone 1\", I=%d, J=%d, K=%d, F=POINT\n",nx,ny,nz);
        }

        Real z = 1.960; //z-score for 95% confidence
        Real oneOsqrtN = sqrt(frac); // frac = 1/N

        printf("Number of Samples: %f\n",1/frac);

        /* write point data */
        for (pt = 0; pt < ext->npts; pt++) {
            Real *co = &ext->coordinates[3*pt];
            Real *UM = &ext->mean[MEAN_NUM*pt];
            Real *SD = &ext->stddev[MEAN_NUM*pt];

            // coordinates
            for (i = 0; i < 3; i++) {
                (co[i]<0.0) ? fprintf(fpa,"%.15f ",co[i]):fprintf(fpa," %.15f ",co[i]);
            }
            // mean
            for (i = 0; i < MEAN_NUM; i++) {
                (UM[i]<0.0) ? fprintf(fpa,"%.15f ",UM[i]):fprintf(fpa," %.15f ",UM[i]);
            }
            // std. deviation
            for (i = 0; i < MEAN_NUM; i++) {
                (SD[i]<0.0) ? fprintf(fpa,"%.15f ",SD[i]):fprintf(fpa," %.15f ",SD[i]);
            }
            // 95% confidence interval
            for (i = 0; i < MEAN_NUM; i++) {
                (z*oneOsqrtN*SD[i]<0.0) ? fprintf(fpa,"%.15f ",z*oneOsqrtN*SD[i]):fprintf(fpa," %.15f ",z*oneOsqrtN*SD[i]);
            }
            fprintf(fpa,"\n");
        }
        /* close tec file */
        fclose(fpa);
        DGOUT(stdout,"[analysis] Output tecplot ascii statistics file: %s\n",file_name_ascii);
    }

    /* clean up memory */
    if(file_name_ascii) {free(file_name_ascii);} file_name_ascii = NULL;
}