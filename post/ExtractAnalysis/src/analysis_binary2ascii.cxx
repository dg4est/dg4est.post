/**
 * \file    analysis_binary2ascii.cxx
 * \ingroup utils_group
 * \author  akirby
 *
 * \brief   Convert binary-formatted data into ascii format for visualization.
 */

/* header files */
#include "analysis_binary2ascii.hxx"

void binary2ascii(inputs_t *data){
    int ex;

    /* setup data structures */
    extracts_t solution_ = extracts_new(); extracts_t *solution = &solution_;

    for (ex = data->ext_sind; ex <= data->ext_eind; ex+=data->ext_ind_skip) {
        convert_extract_files(solution,data,ex);
    }

    /* clean up memory */
    extracts_destroy(solution);
}

void convert_extract_files(extracts_t *stats,inputs_t *data,int ex){
    char *filename;
    int i;

    /* ============================ */
    /* read binary and output ascii */
    /* ============================ */
    for (i = data->file_sind; i <= data->file_eind; i+=data->file_ind_skip) {
        if(!form_file_name(data,ex,i,&filename,DISPLAY_OFF,"AVERAGE")) continue;

        read_binary_data(stats,filename);
        write_ascii_data(stats,data,ex,i,stats->extract_dim);

        /* clean up filename variable */
        if(filename) free(filename); filename = NULL;
    }
}

void read_binary_data(extracts_t *stats,char *file_name){
    const int nhead = 6;
    int sizes[nhead];

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

    /* read point data */
    ret = fread(stats->data.ptr(),sizeof(Real),stats->ncf*stats->npts,fpb);

    /* close bin file */
    fclose(fpb);
}

void write_ascii_data(extracts_t *ext,inputs_t *inputs,int sindex,int istep,int mode){
    int nx,ny,nz;
    int plot_dim;
    int pt,i;

    /* ================== */
    /* write data to file */
    /* ================== */
    char file_ascii[49] = { '\0' };
    snprintf(file_ascii,49, "/tec/extract%03d_%07d.tec",sindex,istep);

    char *file_name_ascii = (char*) malloc(strlen(inputs->output_data_path) +
                                           strlen(file_ascii) + 1);

    strcpy(file_name_ascii,inputs->output_data_path);
    strcat(file_name_ascii,file_ascii);

    /* ---------------- */
    /* write ascii data */
    /* ---------------- */
    if (inputs->output_ascii) {
        FILE *fpa = fopen(file_name_ascii, "w");
        fprintf(fpa,"TITLE = \"DG4EST Extract: %s\"\n",file_ascii);
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
            Real *p = &ext->data[EX_NUM*pt];
            for (i = 0; i < ext->ncf; i++) {
                (p[i]<0.0) ? fprintf(fpa,"%.15f ",p[i]):fprintf(fpa," %.15f ",p[i]);
            }
            fprintf(fpa,"\n");
        }
        /* close tec file */
        fclose(fpa);
        DGOUT(stdout,"[analysis] Output tecplot ascii file: %s\n",file_name_ascii);
    }

    /* clean up memory */
    if(file_name_ascii) {free(file_name_ascii);} file_name_ascii = NULL;
}