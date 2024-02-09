/**
 * \file    analysis_les_subgridscales.c
 * \ingroup utils_group
 * \author  akirby
 *
 * \brief   Computes the ratio of resolved turbulent kinetic energy to sub-grid scale
 *          turbulent kinetic energy from dg4est extract data.
 */

/* header files */
#include "analysis_les_subgridscales.hxx"

void les_subgridscales(inputs_t *data,char standalone){
    char *filename;
    int ex,i;

    /* setup data structures */
    extracts_t mean_ = extracts_new(); extracts_t *mean = &mean_;
    extracts_t stresses_ = extracts_new(); extracts_t *stresses = &stresses_;
    extracts_t tke_ = extracts_new(); extracts_t *tke = &tke_;

    for (ex = data->ext_sind; ex <= data->ext_eind; ex+=data->ext_ind_skip) {
        compute_extract_mean(mean,data,ex);
        compute_extract_stresses(stresses,mean,data,ex);

        /* ========================== */
        /* kinetic energy calculation */
        /* ========================== */
        for (i = data->plot_tec_sind; i <= data->plot_tec_eind; i+=data->plot_tec_ind_skip) {
            if(!form_file_name(data,ex,i,&filename,DISPLAY_OFF,"SGS TKE")) continue;

            compute_tke(tke,stresses,mean,data,filename);

            /* write tec file */
            output_sgs_quantities(tke,data,ex,i,tke->extract_dim,standalone);

            /* clean up filename variable */
            if(filename) free(filename); filename = NULL;
        }

        /* reset extract arrays */
        mean->mean.fill(0.0);
        stresses->stresses.fill(0.0);
    }
    /* 3. clean up memory */
    extracts_destroy(mean);
    extracts_destroy(stresses);
    extracts_destroy(tke);
}

void compute_tke(extracts_t *tke,extracts_t *stresses,extracts_t *means,inputs_t *inputs,char *file_name){
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
    ret = fread(sizes, sizeof(int), nhead, fpb);

    tke->extract_dim = sizes[0];
    tke->ncf = sizes[1];
    tke->npts = sizes[2];
    tke->nx = sizes[3];
    tke->ny = sizes[4];
    tke->nz = sizes[5];
    tke->nf = tke->ncf-3; // subtract number of coordinates
#if 1
    printf("Header Info: %s\n"
           "  Dim : %d\n"
           "  ncf : %d\n"
           "  npts: %d\n"
           "  nx  : %d\n"
           "  ny  : %d\n"
           "  nz  : %d\n",
           file_name,sizes[0],sizes[1],sizes[2],sizes[3],sizes[4],sizes[5]);
#endif
    assert(tke->ncf == EX_NUM);

    /* allocate array data */
    tke->data.newmalloc(EX_NUM*tke->npts);
    tke->tke.newmalloc(TKE_NUM*tke->npts,0.0);
    fill_coordinates = tke->coordinates.newmalloc(3*tke->npts);

    /* read point data */
    ret = fread(tke->data.ptr(),sizeof(Real),EX_NUM*tke->npts,fpb);

    /* fill coordinates */
    if (fill_coordinates) {
        for (pt = 0; pt < tke->npts; ++pt) {
            Real *co = &tke->coordinates[3*pt];
            Real *u  = &tke->data[EX_NUM*pt];

            co[0] = u[EX_X]; // x
            co[1] = u[EX_Y]; // y
            co[2] = u[EX_Z]; // z
        }
    }

    /* close bin file */
    fclose(fpb);

    if (fill_coordinates && inputs->Cartesian2Polar) {
        tke->coordinates_polar.newmalloc(3*tke->npts);

        for (pt = 0; pt < tke->npts; ++pt) {
            Real *co = &tke->coordinates[3*pt];
            Real *pc = &tke->coordinates_polar[3*pt];

            // FIXME: r, theta, z
            pc[0] = co[0]; // x
            pc[1] = co[1]; // y
            pc[2] = co[2]; // z
        }
    }

    /* calculate sgs quantities */
    for (pt = 0; pt < tke->npts; ++pt) {
        /* tke->array contains instantaneous flow from extracts */
        Real *U = &tke->data[EX_NUM*pt];
        Real *UM = &means->mean[MEAN_NUM*pt];
        Real *REY = &stresses->stresses[RE_NUM*pt];
        Real *TKE = &tke->tke[TKE_NUM*pt];

        Real ke_res = REY[RE_KE];
        Real ke_sgs = U[EX_KESGS];
        Real pf = U[EX_PRES] - UM[MEAN_PRES]; // p' = p - <p>

        TKE[TKE_RES] = ke_res;
        TKE[TKE_SGS] = ke_sgs;
        TKE[TKE_PERC] = (ke_res)/(ke_res + ke_sgs);
        TKE[TKE_PRSF] = pf;
        TKE[TKE_PRMS] = sqrt(pf*pf);
        TKE[TKE_IBLK] = U[EX_IBLK];
    }
}

void output_sgs_quantities(extracts_t *ext,inputs_t *inputs,
                           int sindex,int timestep,int mode,char standalone){
    int nx,ny,nz;
    int plot_dim;
    int pt,i;

    /* ================== */
    /* write data to file */
    /* ================== */
    char file_ascii[49] = { '\0' };

    snprintf(file_ascii,49, "/tec/sgs/extract%03d_%07d.SGS.tec",sindex,timestep);

    char *file_name_ascii = (char *) malloc(strlen(inputs->output_data_path) + strlen(file_ascii) + 1);

    strcpy(file_name_ascii,inputs->output_data_path);
    strcat(file_name_ascii,file_ascii);
    printf("%s\n",file_name_ascii);

    /* ---------------- */
    /* write ascii data */
    /* ---------------- */
    if (inputs->output_ascii || !standalone) {
        FILE *fpa = fopen(file_name_ascii, "w");
        fprintf(fpa,"TITLE = \"DG4EST Sub-Grid Scale Kinetic Energy\"\n");
        fprintf(fpa,"VARIABLES = \"x\", \"y\", \"z\", "
                    "\"ke_resolved\", \"ke_sgs\", \"ke_percentage\","
                    "\"p'\", \"p'_rms\", \"iblank\"\n");

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
            Real *co = &ext->coordinates[3*pt];
            Real *p = &ext->tke[TKE_NUM*pt];

            for (i = 0; i < 3; i++) {
                (co[i]<0.0) ? fprintf(fpa,"%.15f ",co[i]):fprintf(fpa," %.15f ",co[i]);
            }
            for (i = 0; i < TKE_NUM; i++) {
                (p[i]<0.0) ? fprintf(fpa,"%.15f ",p[i]):fprintf(fpa," %.15f ",p[i]);
            }
            fprintf(fpa,"\n");
        }
        /* close tec file */
        fclose(fpa);
        DGOUT(stdout,"[analysis] Output tecplot ascii tke msgs file: %s\n",file_name_ascii);
    }

    /* clean up memory */
    if(file_name_ascii) {free(file_name_ascii);} file_name_ascii = NULL;
}