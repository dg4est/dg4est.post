/**
 * \file    analysis_reynolds_stresses.c
 * \ingroup utils_group
 * \author  akirby
 *
 * \brief   Computes the Reynolds Stresses of flow fields from dg4est extract data.
 */

/* header files */
#include "analysis_reynolds_stresses.hxx"

void reynolds_stresses(inputs_t *data,char standalone){
    int ex;

    /* setup data structures */
    extracts_t stresses_ = extracts_new(); extracts_t *stresses = &stresses_;
    extracts_t mean_ = extracts_new(); extracts_t *mean = &mean_;

    if (data->multiple_extract_mean) {
        /* compute mean flow over all extracts */
        for (ex = data->ext_sind; ex <= data->ext_eind; ex+=data->ext_ind_skip) {
            compute_extract_mean(mean,data,ex);
        }
    }

    /* =========================== */
    /* Reynolds Stress calculation */
    /* =========================== */
    for (ex = data->ext_sind; ex <= data->ext_eind; ex+=data->ext_ind_skip) {
        if(!data->multiple_extract_mean) compute_extract_mean(mean,data,ex);
        compute_extract_stresses(stresses,mean,data,ex);

        /* write out extract mean file */
        if (!data->multiple_extract_mean) {
            output_stresses(stresses,data,ex,ex,LOCAL_MODE,stresses->extract_dim,standalone);

            /* reset extract arrays */
            mean->mean.fill(0.0);
            stresses->stresses.fill(0.0);
        }
    }

    /* write out extract mean file */
    if (data->multiple_extract_mean) {
        output_stresses(stresses,data,data->ext_sind,data->ext_eind,GLOBAL_MODE,stresses->extract_dim,standalone);
    }

    /* 3. clean up memory */
    extracts_destroy(mean);
    extracts_destroy(stresses);
}

void compute_extract_stresses(extracts_t *stresses,extracts_t *means,inputs_t *data,int ex){
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
        if(!form_file_name(data,ex,i,&filename,DISPLAY_OFF,NULL)) continue;

        compute_fluctuations_mean(stresses,means,filename,oneOts);

        /* clean up filename variable */
        if(filename) free(filename); filename = NULL;
    }
    compute_resolved_turbulent_kinetic_energy(stresses);
}

void compute_fluctuations_mean(extracts_t *stresses,extracts_t *means,char *file_name,Real frac){
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

    stresses->extract_dim = sizes[0];
    stresses->ncf = sizes[1];
    stresses->npts = sizes[2];
    stresses->nx = sizes[3];
    stresses->ny = sizes[4];
    stresses->nz = sizes[5];
    stresses->nf = RE_NUM; // Reynolds Stresses + pressure_fluc_square + resolved TKE

    /* =============== */
    /* allocate memory */
    /* =============== */
    stresses->data.newmalloc(stresses->ncf*stresses->npts);
    stresses->stresses.newmalloc(RE_NUM*stresses->npts,0.0);
    fill_coordinates = stresses->coordinates.newmalloc(3*stresses->npts);

    /* ===================== */
    /* read point field data */
    /* ===================== */
    ret = fread(stresses->data.ptr(),sizeof(Real),stresses->ncf*stresses->npts,fpb);

    /* fill coordinates */
    if(fill_coordinates) {
        for (pt = 0; pt < stresses->npts; ++pt) {
            Real *co = &stresses->coordinates[3*pt];
            Real *u  = &stresses->data[stresses->ncf*pt];

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
    for (pt = 0; pt < stresses->npts; ++pt) {
        Real *U  = &stresses->data[EX_NUM*pt+3]; // +3 to skip coordinates
        Real *UM = &means->mean[MEAN_NUM*pt];
        Real *Re = &stresses->stresses[RE_NUM*pt];

        /* NOTE: U/UM = [rho,u,v,w,p] */
        Real uf = U[EX_U] - UM[MEAN_U];       // u' = u - <u>
        Real vf = U[EX_V] - UM[MEAN_V];       // v' = v - <v>
        Real wf = U[EX_W] - UM[MEAN_W];       // w' = w - <w>

        Re[RE_UU] += frac*(uf*uf); // <u'u'>
        Re[RE_VV] += frac*(vf*vf); // <v'v'>
        Re[RE_WW] += frac*(wf*wf); // <w'w'>
        Re[RE_UV] += frac*(uf*vf); // <u'v'>
        Re[RE_UW] += frac*(uf*wf); // <u'w'>
        Re[RE_VW] += frac*(vf*wf); // <v'w'>
        Re[RE_IBLK] = UM[MEAN_IBLK];
    }
}

void compute_resolved_turbulent_kinetic_energy(extracts_t *stresses){
    int pt;
    /* =========================== */
    /* calculate fluctuation means */
    /* =========================== */
    for (pt = 0; pt < stresses->npts; ++pt) {
        Real *Re = &stresses->stresses[RE_NUM*pt];
        Real  uu =  Re[RE_UU];
        Real  vv =  Re[RE_VV];
        Real  ww =  Re[RE_WW];
        Real *ke = &Re[RE_KE];

        ke[0] = 0.5*(uu + vv + ww);
    }
}

void output_stresses(extracts_t *ext,inputs_t *inputs,int sindex,int eindex,int gl_mode,int mode,char standalone){
    int nx,ny,nz;
    int plot_dim;
    int pt,i;

    /* ================== */
    /* write data to file */
    /* ================== */
    char file_ascii[48] = { '\0' };
    char file_binary[48] = { '\0' };

    if (gl_mode == LOCAL_MODE) {
        snprintf(file_ascii,37, "/tec/extract%03d.ReynoldsStresses.tec",sindex);
        snprintf(file_binary,37,"/bin/extract%03d.ReynoldsStresses.bin",sindex);
    } else {
        snprintf(file_ascii,48, "/tec/extracts%03d-%03d.ReynoldsStresses.tec",sindex,eindex);
        snprintf(file_binary,48,"/bin/extracts%03d-%03d.ReynoldsStresses.bin",sindex,eindex);
    }

    char *file_name_ascii = (char *) malloc(strlen(inputs->output_data_path) + strlen(file_ascii) + 1);
    char *file_name_binary = (char *) malloc(strlen(inputs->output_data_path) + strlen(file_binary) + 1);

    strcpy(file_name_ascii,inputs->output_data_path); strcpy(file_name_binary,inputs->output_data_path);
    strcat(file_name_ascii,file_ascii);               strcat(file_name_binary,file_binary);

    /* fill array */
    for (pt = 0; pt < ext->npts; pt++) {
        Real *arr = &ext->data[(ext->nf+3)*pt]; // reuse contiguous buffer
        Real *co = &ext->coordinates[3*pt];
        Real *p = &ext->stresses[ext->nf*pt];

        arr[0] = co[0]; arr[1] = co[1]; arr[2] = co[2];
        for(i = 0; i < ext->nf; i++) arr[i+3] = p[i];
    }

    /* ---------------- */
    /* write ascii data */
    /* ---------------- */
    if (inputs->output_ascii && standalone) {
        FILE *fpa = fopen(file_name_ascii, "w");
        fprintf(fpa,"TITLE = \"DG4EST Reynolds Stresses\"\n");
        fprintf(fpa,"VARIABLES = \"x\", \"y\", \"z\", "
                    "\"u'u'\", \"v'v'\", \"w'w'\", \"u'v'\", \"u'w'\", \"v'w'\", \"ke_resolved\", \"iblank\"\n");

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
            Real *p = &ext->data[(ext->nf+3)*pt];
            for (i = 0; i < (ext->nf+3); i++) {
                (p[i]<0.0) ? fprintf(fpa,"%.15f ",p[i]):fprintf(fpa," %.15f ",p[i]);
            }
            fprintf(fpa,"\n");
        }
        /* close tec file */
        fclose(fpa);
        DGOUT(stdout,"[analysis] Output tecplot ascii Reynolds Stresses file: %s\n",file_name_ascii);
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
        fwrite(ext->data.ptr(),sizeof(Real),(ext->nf+3)*ext->npts,fpb);

        /* close bin file */
        fclose(fpb);
        DGOUT(stdout,"[analysis] Output binary Reynolds Stresses file: %s\n",file_name_binary);
    }

    /* clean up memory */
    if(file_name_ascii) {free(file_name_ascii);} file_name_ascii = NULL;
    if(file_name_binary) {free(file_name_binary);} file_name_binary = NULL;
}