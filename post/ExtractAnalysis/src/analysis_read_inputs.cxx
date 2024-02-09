/**
 * \file    analysis_read_inputs.c
 * \author  akirby
 *
 * \brief   Input file reader and utilities.
 */

/* header files */
#include "analysis_read_inputs.hxx"

inputs_t* read_inputs(int argc,char **argv){
    inputs_t *data = inputs_new();
    char analysis_read = 0;
    int i;

    /* ====================================================================== */
    /* read command line */
    if (argc < 2) {
        DGOUT(stdout,"[USAGE]: ./Extract.Analysis < input.file\n");
        DGOUT(stdout,rC "\tNo Input File Found!\n" eC);

        printf("Enter 'y' to see a sample input file: ");
        char c = getchar();
        if(c == 'y' || c == 'Y') inputs_example();
        exit(EXIT_SUCCESS);
    } else {
        DGOUT(stdout," >> Input file: %s\n",argv[1]);
        strcpy(data->input_file,argv[1]);

        if (argc == 3) {
            data->analysis_type = atoi(argv[2]);
            analysis_read = 1;
        }
        if (argc == 4) {
            data->ext_sind = atoi(argv[2]);
            data->ext_eind = atoi(argv[3]);
            data->args_used = 1;
            printf("Reading Extract indices from command line: %d %d\n",
                    data->ext_sind,data->ext_eind);
        }
    }
    /* ====================================================================== */

    /* ====================================================================== */
    find_keyword_integer(data->input_file,"output_format:",&data->output_format,OPTIONAL);
    find_keyword_integer(data->input_file,"plot_tec_start_index:",&data->plot_tec_sind,OPTIONAL);
    find_keyword_integer(data->input_file,"plot_tec_end_index:",&data->plot_tec_eind,OPTIONAL);
    find_keyword_integer(data->input_file,"plot_tec_index_skip:",&data->plot_tec_ind_skip,OPTIONAL);

    find_keyword_string(data->input_file,"input_data_path:",data->input_data_path,REQUIRED);

    /* extracts info */
    find_keyword_integer(data->input_file,"extract_index_skip:",&data->ext_ind_skip,OPTIONAL);
    if (!data->args_used) {
        find_keyword_integer(data->input_file,"extract_start_index:",&data->ext_sind,REQUIRED);
        find_keyword_integer(data->input_file,"extract_end_index:",&data->ext_eind,REQUIRED);
    }

    find_keyword_integer(data->input_file,"file_index_skip:",&data->file_ind_skip,OPTIONAL);
    find_keyword_integer(data->input_file,"file_start_index:",&data->file_sind,REQUIRED);
    find_keyword_integer(data->input_file,"file_end_index:",&data->file_eind,REQUIRED);

    /* analysis info */
    if(!analysis_read) find_keyword_integer(data->input_file,"analysis_type:",&data->analysis_type,REQUIRED);
    find_keyword_integer(data->input_file,"average_multiple_extracts:",&data->multiple_extract_mean,OPTIONAL);
    find_keyword_integer(data->input_file,"transform_Cart2Polar:",&data->Cartesian2Polar,OPTIONAL);

    /* POD variables */
    data->npod = 0;
    data->pod_vars_mask = 0; // initialize to 0
    if(data->analysis_type == ANALYSIS_POD_TRAD || data->analysis_type == ANALYSIS_POD_SNAP){
        char podvar_string[BUFF_SIZE] = { '\0' };

        char notfound = find_keyword_string_caps(data->input_file,TO_UPPER,
                                                 "POD_vars:",podvar_string,
                                                 REQUIRED);
        for (i = 0; i < MAX_NPODVAR; i++) {
            if(!notfound && strstr(podvar_string,PODVAR_NAMES[i])) {
              //printf("FOUND POD VARIABLE: %s\n",PODVAR_NAMES[i]);
                data->pod_vars_mask |= PODVAR_MASKS[i];
                data->npod++;
            }
        }
        find_keyword_integer(data->input_file,"POD_nmode:",&data->pod_nmode,OPTIONAL);
        find_keyword_integer(data->input_file,"POD_plot_nmodes:",&data->pod_plot_nmodes,OPTIONAL);
        find_keyword_integer(data->input_file,"POD_plot_ROM:",&data->pod_plot_ROM,OPTIONAL);
        find_keyword_integer(data->input_file,"POD_plot_TVC:",&data->pod_plot_TVC,OPTIONAL);
    }
    /* ====================================================================== */

    /* set input base name */
    strcpy(data->output_data_path,"analysis");
    strcpy(data->input_base_name,"dg4est_ext");

    /* check input values */
    if(data->plot_tec_ind_skip < 1) data->plot_tec_ind_skip = 1;
    if(data->ext_ind_skip < 1) data->ext_ind_skip = 1;
    if(data->file_ind_skip < 1) data->file_ind_skip = 1;

    if(data->plot_tec_sind < data->file_sind) data->plot_tec_sind = data->file_sind;
    if(data->plot_tec_eind > data->file_eind) data->plot_tec_eind = data->file_eind;

    /* set number of steps */
    data->nstep = (int) (data->file_eind-data->file_sind+data->file_ind_skip)/(Real)(data->file_ind_skip);

    switch(data->output_format){
        case(1): data->output_ascii = 1; break;
        case(2): data->output_binary = 1; break;
        case(3): data->output_ascii = data->output_binary = 1; break;
    }

    switch(data->input_extension){
        case(1): strcpy(data->input_extension_str,".bin"); break;
        case(2): strcpy(data->input_extension_str,".tec"); break;
        default: strcpy(data->input_extension_str,".bin"); break;
    }

    char PODFLAG = (data->analysis_type==ANALYSIS_POD_TRAD) ||
                   (data->analysis_type==ANALYSIS_POD_SNAP);
    /* build output directory and display read inputs */
    create_directory("analysis");
    create_directory("analysis/bin");
    create_directory("analysis/tec");

    if(data->analysis_type==ANALYSIS_LES_SGS)  create_directory("analysis/tec/sgs");
    if(PODFLAG) create_directory("analysis/tec/pod");
    if(PODFLAG) create_directory("analysis/tec/pod/modes");

    /* make symbolic links to tecplot templates */
    int ret;
    if(PODFLAG) ret=symlink("../../../tec.templates/energy.lay","analysis/tec/pod/energy.lay");
    inputs_dump(data);
    return data;
}

int form_file_name(inputs_t *data,int ex,int i,char **filename,char disp,const char *msg){
    char lastc = data->input_data_path[strlen(data->input_data_path)-1];
    int slash_missing = (lastc == '/') ? 0:1;

    char file_timestep[16] = {'\0'};
    snprintf(file_timestep,16,"%03d_%07d.bin",ex,i);

    if (slash_missing) {
        *filename = (char *) malloc(strlen(data->input_data_path) +
                                    1 + // slash
                                    strlen(data->input_base_name) +
                                    strlen(file_timestep) + 1);
        strcpy(*filename,data->input_data_path);
        strcat(*filename,"/");
        strcat(*filename,data->input_base_name);
        strcat(*filename,file_timestep);
    } else {
        *filename = (char *) malloc(strlen(data->input_data_path) +
                                    strlen(data->input_base_name) +
                                    strlen(file_timestep) + 1);
        strcpy(*filename,data->input_data_path);
        strcat(*filename,data->input_base_name);
        strcat(*filename,file_timestep);
    }

    /* ==================== */
    /* check file existence */
    /* ==================== */
    const char *out = (msg==NULL) ? "Analysis":msg;
    int file_exists = (access(*filename, F_OK) != -1) ? 1:0;
    if (!file_exists) {
        printf(rC "[%s] File does not exist (skipping): %s\n" eC,out,*filename);
        fflush(stdout);
    } else {
        if(disp) printf(gC "[%s] File found: %s\n" eC,out,*filename);
        if(disp) fflush(stdout);
    }
    return file_exists;
}

void create_directory(const char dir[]){
    DGOUT(stdout,"Creating %s directory\n",dir);
    mkdir(dir, S_IRWXU | S_IRWXG | S_IRWXO);
}

char find_keyword_integer(char *filename,const char *keyword,
                          int *integer,int required){
    FILE *fp;
    char buff[BUFF_SIZE] = {'\0'};
    fp = fopen(filename,"r");

    unsigned long length = strlen(keyword);

    while (fp != NULL && fgets(buff,sizeof(buff),fp) != NULL) {
        if (strstr(buff,keyword)) {

            *integer = atoi(&buff[length]);

            fclose(fp);
            return EXIT_SUCCESS;
        }
    }
    /* otherwise */
    if (required) {
        DGOUT(stdout,"\033[1;31m>>>>> Keyword [%s] not found! Using default.\033[0m\n",
                keyword);
    }

    fclose(fp);
    return EXIT_FAILURE;
}

char find_keyword_two_integers(char *filename,const char *keyword,
                               int *integer1,int *integer2,int required){
    FILE *fp;
    char buff[BUFF_SIZE];
    fp = fopen(filename,"r");

    unsigned long length = strlen(keyword);

    while (fp != NULL && fgets(buff,sizeof(buff),fp) != NULL) {
        if (strstr(buff,keyword)) {

          sscanf(&buff[length],"%d %d", integer1,integer2);
          sscanf(&buff[length],"%d, %d",integer1,integer2);

          fclose(fp);
          return EXIT_SUCCESS;
        }
    }
    /* otherwise */
    if (required) {
        DGOUT(stdout,"\033[1;31m>>>>> Keyword [%s] not found! Using default.\033[0m\n",
                keyword);
    }

    fclose(fp);
    return EXIT_FAILURE;
}

char find_keyword_three_integers(char *filename,const char *keyword,
                                 int *integer1,int *integer2,int *integer3,int required){
    FILE *fp;
    char buff[BUFF_SIZE] = {'\0'};
    fp = fopen(filename,"r");

    unsigned long length = strlen(keyword);

    while (fp != NULL && fgets(buff,sizeof(buff),fp) != NULL) {
        if (strstr(buff,keyword)) {

            sscanf(&buff[length],"%d %d %d",  integer1,integer2,integer3);
            sscanf(&buff[length],"%d, %d, %d",integer1,integer2,integer3);
            sscanf(&buff[length],"%d, %d %d", integer1,integer2,integer3);
            sscanf(&buff[length],"%d %d, %d", integer1,integer2,integer3);

            fclose(fp);
            return EXIT_SUCCESS;
        }
    }
    /* otherwise */
    if (required) {
        DGOUT(stdout,"\033[1;31m>>>>> Keyword [%s] not found! Using default.\033[0m\n",
                keyword);
    }

    fclose(fp);
    return EXIT_FAILURE;
}

char find_keyword_real(char *filename,const char *keyword,
                       Real *dbl,int required){
    FILE *fp;
    char buff[BUFF_SIZE];
    fp = fopen(filename,"r");

    unsigned long length = strlen(keyword);

    while (fp != NULL && fgets(buff,sizeof(buff),fp) != NULL) {
        if (strstr(buff,keyword)) {

            sscanf(&buff[length],RealFormat,dbl);

            fclose(fp);
            return EXIT_SUCCESS;
        }
    }
    /* otherwise */
    if (required) {
        DGOUT(stdout,"\033[1;31m>>>>> Keyword [%s] not found! Using default.\033[0m\n",
                keyword);
    }

    fclose(fp);
    return EXIT_FAILURE;
}

char find_keyword_two_reals(char *filename,const char *keyword,
                            Real *dbl1,Real *dbl2,int required){
    FILE *fp;
    char buff[BUFF_SIZE] = {'\0'};
    char pref[BUFF_SIZE] = {'\0'};
    fp = fopen(filename,"r");

    unsigned long length = strlen(keyword);

    while (fp !=NULL && fgets(buff,sizeof(buff),fp) != NULL) {
        if (strstr(buff,keyword)) {
            sprintf(pref,"%s %s",RealFormat,RealFormat);
            sscanf(&buff[length],pref,dbl1,dbl2);

            sprintf(pref,"%s, %s",RealFormat,RealFormat);
            sscanf(&buff[length],pref,dbl1,dbl2);

            sprintf(pref,"%s,%s",RealFormat,RealFormat);
            sscanf(&buff[length],pref,dbl1,dbl2);

            fclose(fp);
            return EXIT_SUCCESS;
        }
    }
    /* otherwise */
    if (required) {
        DGOUT(stdout,"\033[1;31m>>>>> Keyword [%s] not found! Using default.\033[0m\n",
                keyword);
    }

    fclose(fp);
    return EXIT_FAILURE;
}

char find_keyword_three_reals(char *filename,const char *keyword,
                              Real *dbl1,Real *dbl2,Real *dbl3,int required){
    FILE *fp;
    char buff[BUFF_SIZE] = {'\0'};
    char pref[BUFF_SIZE] = {'\0'};
    fp = fopen(filename,"r");

    unsigned long length = strlen(keyword);

    while (fp !=NULL && fgets(buff,sizeof(buff),fp) != NULL) {
        if (strstr(buff,keyword)) {
            sprintf(pref,"%s %s %s",RealFormat,RealFormat,RealFormat);
            sscanf(&buff[length],pref,  dbl1,dbl2,dbl3);

            sprintf(pref,"%s, %s, %s",RealFormat,RealFormat,RealFormat);
            sscanf(&buff[length],pref,dbl1,dbl2,dbl3);

            sprintf(pref,"%s %s, %s",RealFormat,RealFormat,RealFormat);
            sscanf(&buff[length],pref, dbl1,dbl2,dbl3);

            sprintf(pref,"%s, %s %s",RealFormat,RealFormat,RealFormat);
            sscanf(&buff[length],pref, dbl1,dbl2,dbl3);

            fclose(fp);
            return EXIT_SUCCESS;
        }
    }
    /* otherwise */
    if (required) {
        DGOUT(stdout,"\033[1;31m>>>>> Keyword [%s] not found! Using default.\033[0m\n",
                keyword);
    }

    fclose(fp);
    return EXIT_FAILURE;
}

char find_keyword_int_real(char *filename,const char *keyword,
                           int *int1,Real *dbl1,int required){
    FILE *fp;
    char buff[BUFF_SIZE] = {'\0'};
    char pref[BUFF_SIZE] = {'\0'};
    fp = fopen(filename,"r");

    unsigned long length = strlen(keyword);

    while (fp !=NULL && fgets(buff,sizeof(buff),fp) != NULL) {
        if (strstr(buff,keyword)) {
            sprintf(pref,"%s %s",IntFormat,RealFormat);
            sscanf(&buff[length],pref,int1,dbl1);

            sprintf(pref,"%s, %s",IntFormat,RealFormat);
            sscanf(&buff[length],pref,int1,dbl1);

            sprintf(pref,"%s ,%s",IntFormat,RealFormat);
            sscanf(&buff[length],pref,int1,dbl1);

            fclose(fp);
            return EXIT_SUCCESS;
        }
    }
    /* otherwise */
    if (required) {
        DGOUT(stdout,"\033[1;31m>>>>> Keyword [%s] not found! Using default.\033[0m\n",
                keyword);
    }

    fclose(fp);
    return EXIT_FAILURE;
}

char find_keyword_string(char *filename,const char *keyword,
                         char *string,int required){
    FILE *fp;
    char buff[BUFF_SIZE] = {'\0'};
    unsigned long i,length;

    fp = fopen(filename,"r");
    length = strlen(keyword);

    while(fp != NULL && fgets(buff,sizeof(buff),fp) != NULL) {
        if (strstr(buff,keyword)) {

            int ind=0;
            for (i = length; i < BUFF_SIZE; ++i) {
                if (buff[i] != ' ') {
                    string[ind] = buff[i];
                    ind = ind+1;
            }
            }

            string[strcspn(string, "\n")] = 0;

            if (strlen(string) >= BUFF_SIZE) {
                DGOUT(stdout,"[CartDG] string length is greater than allowable");
                exit(EXIT_FAILURE);
            }

    fclose(fp);
            return EXIT_SUCCESS;
        }
    }
    /* otherwise */
    if (required) {
        DGOUT(stdout,"\033[1;31m>>>>> Keyword [%s] not found! Using default.\033[0m\n",
                keyword);
    }

    fclose(fp);
    return EXIT_FAILURE;
}

char find_keyword_string_caps(char *filename,char upper_flag,
                              const char *keyword,char *string,int required){
    FILE *fp;
    char buff[BUFF_SIZE] = {'\0'};
    unsigned long i,length;

    fp = fopen(filename,"r");
    length = strlen(keyword);

    while(fp != NULL && fgets(buff,sizeof(buff),fp) != NULL) {
        if (strstr(buff,keyword)) {

            int ind=0;
            for (i = length; i < BUFF_SIZE; ++i) {
                string[ind++] = (upper_flag) ? toupper(buff[i]):buff[i];
            }

            /* cutoff comments after # sign */
            char *ptr = strchr(string, '#');
            if(ptr != NULL) *ptr = '\0';

            string[strcspn(string, "\n")] = 0;

            if (strlen(string) >= BUFF_SIZE) {
                DGOUT(stdout,"[CartDG] string length is greater than allowable");
                exit(EXIT_FAILURE);
            }
            fclose(fp);
            return EXIT_SUCCESS;
        }
    }
    /* otherwise */
    if (required) {
        DGOUT(stdout,"\033[1;31m>>>>> Keyword [%s] not found! Using default.\033[0m\n",
                keyword);
    }

    fclose(fp);
    return EXIT_FAILURE;
}
