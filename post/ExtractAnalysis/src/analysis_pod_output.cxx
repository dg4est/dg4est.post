/**
 * \file    analysis_pod_output.cxx
 * \ingroup utils_group
 * \author  akirby
 */

/* header files */
#include "analysis_pod_output.hxx"

void output_POD_energies(extracts_t *ext,inputs_t *inputs,int sindex,int neig){
    pod_t *pod = &ext->pod;
    const int NPOD = pod->npod;
    int var,i;

    /* ================== */
    /* write data to file */
    /* ================== */
    for (var = 0; var < NPOD; var++) {
        char file_ascii[120] = { '\0' };
        char file_extract_name[80] = { '\0' };

        snprintf(file_ascii,120, "/tec/pod/extract%03d.POD.energies.%s.tec",sindex,pod->pod_vars_names[var]);
        snprintf(file_extract_name,80, "extract%03d.POD.energies.%s.tec",sindex,pod->pod_vars_names[var]);
        char *file_name_ascii = (char *) malloc(strlen(inputs->output_data_path) + strlen(file_ascii) + 1);

        strcpy(file_name_ascii,inputs->output_data_path);
        strcat(file_name_ascii,file_ascii);

        Real *eigval = &pod->eig[neig*var]; // get address to POD variable
        Real energy_total = pod->pod_total_energy[var];

        /* ---------------- */
        /* write ascii data */
        /* ---------------- */
        if (inputs->output_ascii) {
            FILE *fpa = fopen(file_name_ascii, "w");
            fprintf(fpa,"TITLE = \"Wake Analysis plot\"\n");
            fprintf(fpa,"VARIABLES = \"Mode Number\", \"Mode Energy %%\", \"Energy Sum %%\", \"Total Energy\"\n");
            fprintf(fpa,"%25d  %15.12f  %15.12f %18.15f\n",0,0.0,0.0,energy_total);

            Real energy_sum = 0.0;
            for (i = neig-1; i >= 0; i--) {
                energy_sum += pod->eig[i]*100.0;
                fprintf(fpa,"%25d  %15.12f  %15.12f %18.15f\n",
                             neig-i,eigval[i]*100.0,energy_sum,energy_total);
            }
            /* close tec file */
            fclose(fpa);
            DGOUT(stdout,"[analysis] Output tecplot ascii POD energy file: %s\n",file_name_ascii);

            /* make symbolic link to energy file for automated tecplot layout */
	    int ret;
            ret = symlink(file_extract_name,"analysis/tec/pod/energy.tec");
        }
        /* clean up memory */
        if(file_name_ascii) {free(file_name_ascii);} file_name_ascii = NULL;
    }
}

void output_ROM_quantities(extracts_t *ext,inputs_t *inputs,
                           int sindex,int timestep,int mode,
                           char standalone){
    pod_t *pod = &ext->pod;
    const int NPTS = ext->npts;
    int nx,ny,nz;
    int plot_dim;
    int pt,i;

    /* ================== */
    /* write data to file */
    /* ================== */
    char file_ascii[53] = { '\0' };
    snprintf(file_ascii,53, "/tec/pod/extract%03d_%07d.POD.ROM.tec",sindex,timestep);
    char *file_name_ascii = (char *) malloc(strlen(inputs->output_data_path) + strlen(file_ascii) + 1);

    strcpy(file_name_ascii,inputs->output_data_path);
    strcat(file_name_ascii,file_ascii);

    /* ---------------- */
    /* write ascii data */
    /* ---------------- */
    if (inputs->output_ascii || !standalone) {
        FILE *fpa = fopen(file_name_ascii, "w");
        fprintf(fpa,"TITLE = \"DG4EST POD ROM\"\n");
        fprintf(fpa,"VARIABLES = \"x\", \"y\", \"z\"");
        for(i = 0; i < pod->npod; i++) fprintf(fpa,", \"%s\"",pod->pod_vars_names[i]);
        fprintf(fpa,"\n");

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

            for (i = 0; i < 3; i++) {
                (co[i]<0.0) ? fprintf(fpa,"%.15f ",co[i]):fprintf(fpa," %.15f ",co[i]);
            }
            for (i = 0; i < pod->npod; i++) {
                Real podvar = pod->UROM[NPTS*i+pt];
                (podvar<0.0) ? fprintf(fpa,"%.15f ",podvar):fprintf(fpa," %.15f ",podvar);
            }
            fprintf(fpa,"\n");
        }
        /* close tec file */
        fclose(fpa);
      //DGOUT(stdout,"[analysis] Output tecplot ascii POD file: %s\n",file_name_ascii);
    }

    /* clean up memory */
    if(file_name_ascii) {free(file_name_ascii);} file_name_ascii = NULL;
}

void output_TVC_quantities(extracts_t *ext,inputs_t *inputs,int sindex){
    pod_t *pod = &ext->pod;
    const int NPOD = pod->npod;
    const int NSTEP = ext->nStep;
    const int NMODES = inputs->pod_nmode;
    int var,istep,mode;

    /* ================== */
    /* write data to file */
    /* ================== */
    for (var = 0; var < NPOD; var++) {
        char file_ascii[120] = { '\0' };

        snprintf(file_ascii,120, "/tec/pod/extract%03d.POD.TVC.%s.tec",sindex,pod->pod_vars_names[var]);
        char *file_name_ascii = (char *) malloc(strlen(inputs->output_data_path) + strlen(file_ascii) + 1);

        strcpy(file_name_ascii,inputs->output_data_path);
        strcat(file_name_ascii,file_ascii);

        /* ---------------- */
        /* write ascii data */
        /* ---------------- */
        if (inputs->output_ascii) {
            FILE *fpa = fopen(file_name_ascii, "w");
            fprintf(fpa,"TITLE = \"DG4EST POD Time-Varying Coefficients\"\n");
            fprintf(fpa,"VARIABLES = \"Time\"");
            for (mode = 0; mode < inputs->pod_plot_TVC; mode++) fprintf(fpa,", \"mode%d\"",mode+1);
            fprintf(fpa,"\n");

            Real * const TVC_STEPS = &pod->TVC[NMODES*NSTEP*var]; //@(:,:,i)
            for (istep = 0; istep < NSTEP; istep++) {
                Real * const TVC = &TVC_STEPS[NMODES*istep]; //@(:,istep)
                fprintf(fpa,"%6d  ",istep);
                for (mode = 0; mode < inputs->pod_plot_TVC; mode++) {
                    fprintf(fpa,"%15.12f  ",TVC[mode]);
                }
                fprintf(fpa,"\n");
            }
            /* close tec file */
            fclose(fpa);
            DGOUT(stdout,"[analysis] Output tecplot ascii POD TVC file: %s\n",file_name_ascii);
        }
        /* clean up memory */
        if(file_name_ascii) {free(file_name_ascii);} file_name_ascii = NULL;
    }
}
