/*!
\file  
\brief Parsing of command-line arguments

This file parses the command line arguments

\date   Started 6/16/14
\author George, David
 */


#include "orion.h"

/*-------------------------------------------------------------------
 * Command-line options 
 *-------------------------------------------------------------------*/
static struct gk_option long_options[] = {
        {"mode",              1,      0,      CMD_MODE},
        {"simtype",           1,      0,      CMD_SIMTYPE},

        {"rowmodel",          1,      0,      CMD_ROWMODEL},
        {"minspan",           1,      0,      CMD_MINSPAN},
        {"niters",            1,      0,      CMD_NITERS},
        {"mintp",             1,      0,      CMD_MINTP},
        {"n2frac",            1,      0,      CMD_N2FRAC},
        {"scost",             1,      0,      CMD_SCOST},
        {"diff",              1,      0,      CMD_DIFFMODE},
        {"simt",              1,      0,      CMD_DIFFSIMT},
        {"ldim",              1,      0,      CMD_DIFFLDIM},
        {"dg",                0,      0,      CMD_DIFFGLOBAL},
        {"svdmaxit",          1,      0,      CMD_SVDMAXIT},
        {"seed",              1,      0,      CMD_SEED},
        {"of",                1,      0,      CMD_OFILE},

        {"w",                 0,      0,      CMD_WRITEPATHS},
        {"wftrs",             0,      0,      CMD_WRITEFTRS},

        {"nthreads",          1,      0,      CMD_NUMTHREADS},
        {"dbglvl",            1,      0,      CMD_DBGLVL},
        {"help",              0,      0,      CMD_HELP},
        {0,                   0,      0,      0}
};

/*-------------------------------------------------------------------
 * Mappings for the various parameter values
 *-------------------------------------------------------------------*/
static gk_StringMap_t simtype_options[] = {
        {"cos",                SIMTYPE_COS},
        {"ejc",                SIMTYPE_EJC},
        {"sqe",                SIMTYPE_SQE},
        {NULL,                 0}
};

static gk_StringMap_t mode_options[] = {
        {"orion",              MODE_ORION},
        {"diff",               MODE_DIFF},
        {"otm",                MODE_OTM},
        {NULL,                 0}
};

static gk_StringMap_t diffmode_options[] = {
        {"pw",                 DIFFMODE_PW},
        {"tpw",                DIFFMODE_TPW},
        {"dr",                 DIFFMODE_DR},
        {"tdr",                DIFFMODE_TDR},
        {"prw",                DIFFMODE_PRW},
        {NULL,                 0}
};


static gk_StringMap_t rowmodel_options[] = {
        {"none",               ROWMODEL_NONE},
        {"log",                ROWMODEL_LOG},
        {"sqrt",               ROWMODEL_SQRT},
        {NULL,                 0}
};


/*-------------------------------------------------------------------
 * Mini help
 *-------------------------------------------------------------------*/
static char helpstr[][100] =
{
        " ",
        "Usage: orion+ [options] fstem nprotos",
        " ",
        " ",
        " Options",
        "  -mode=text",
        "     Specifies the algorithm that should be applied.",
        "     Possible values are:",
        "        orion   Original orion algorithm",
        "        diff    Orion with diffused usage vectors [Default]",
        "        otm     Orion with micro-proto topic modeling",
        " ",
        "  -diff=text",
        "     Diffusion model (for 'diff' algorithm mode only). .",
        "     Possible values are:",
        "        pw      pairwise similarities, S = X*X^T",
        "        ptw     thresholded pairwise similarities, S = L2AP(S, t) [Default]",
        "        dr      Dimensionality reduction",
        "        tdr     Dimensionality reduction + thresholded pairwise similarities",
        "        prw     Personalized random walk",
        " ",
        "  -simtype=text",
        "     Specifies the similarity/distance function to use for clustering.",
        "     Possible values are:",
        "        cos     Cosine similarity",
        "        ejc     Extended Jacquard coefficient",
        "        sqe     Squered error [Default]",
        " ",
        "  -rowmodel=text",
        "     Specifies how the values will be scaled.",
        "     Possible values are:",
        "        none    No scaling.",
        "        log     Take the log of the raw values [Default]",
        "        sqrt    Take the square-root of the raw values",
        " ",
        "  -minspan=int",
        "     Specifies the minimum number of consecutive weeks that a proto must",
        "     cover a user's time-series.",
        "     Default value is 5.",
        " ",
        "  -niters=int",
        "     Specifies the maximum number of clustering iterations.",
        "     Default value is 20.",
        " ",
        "  -scost=float",
        "     A per segment cost as a fraction of the error.",
        "     Default value is .01.",
        " ",
        "  -simt=float",
        "     Minimum similarity for thresholded pairwise similarity computation in DIFF mode.",
        "     Default value is .9.",
        " ",
        "  -ldim=int",
        "     Number of latent dimensions for dimensionality reduction in DIFF mode.",
        "     Default value is 20.",
        " ",
        "  -dg",
        "     Use a global diffusion matrix in DIFF mode.",
        " ",
        "  -mintp=float",
        "     Specifies the minimum P2P transition probability to be analyzed.",
        "     Default value is .20.",
        " ",
        "  -svdmaxit=int",
        "     Maximum number of iterations for SVD.",
        "     Default value is 500.",
        " ",
        "  -n2frac=float",
        "     Specifies the maximum fraction of a vector to be analyzed for features",
        "     Default value is .80.",
        " ",
        "  -of=str",
        "     OUtput file for final segmentation.",
        "     Default value is <fstem>_<mode>_<nprotos>_<minspan>_<scost>[_<ldim>_<simt>].seg",
        " ",
        "  -seed=int",
        "     Specifies the seed for the random number generator. Use current time: -seed=t.",
        "     Default value is 1.",
        " ",
        "  -w",
        "     Outputs the final segmentation of -of file: seg order + seg points + seg protos.",
        " ",
        "  -wftrs",
        "     Outputs the key features for each proto.",
        " ",
        "  -dbglvl=int",
        "     Specifies the level of debugging information to be displayed.",
        "     Default value is 0.",
        " ",
        "  -nthreads=int",
        "     Specifies the number of threads that should be used in parallel computations.",
        "     Default value is 8.",
        " ",
        "  -help",
        "     Prints this message.",
        " ",
};



/*************************************************************************/
/*! This is the entry point of the command-line argument parser */
/*************************************************************************/
params_t *getcmdline_params(int argc, char *argv[])
{
    gk_idx_t i, j, k;
    int type=0;
    int c, option_index;
    params_t *params;
    char fname[100];

    params = (params_t *)gk_malloc(sizeof(params_t), "cmdline_parse: params");
    memset(params, 0, sizeof(params_t));

    /* print the command line */
    for (i=0; i<argc; i++)
        printf("%s ", argv[i]);
    printf("\n");

    /* initialize the params data structure */
    params->nprotos    = -1;
    params->datafile   = NULL;
    params->dataset    = NULL;
    params->ofile      = NULL;

    params->minspan    = 5;
    params->scost      = .01;
    params->mintp      = .2;
    params->n2frac     = .8;
    params->niters     = 20;
    params->simtype    = SIMTYPE_SQE;
    params->rowmodel   = ROWMODEL_LOG;
    params->mode       = MODE_DIFF;
    params->diffmode   = DIFFMODE_TPW;
    params->simt       = 0.9;
    params->ldim       = 20;
    params->diffglobal = 0;
    params->svdmaxit   = 500;

    params->writepaths = 0;
    params->writeftrs  = 0;

    params->dbglvl     = 0;
    params->nthreads   = 8;


    /* Parse the command line arguments  */
    while ((c = gk_getopt_long_only(argc, argv, "", long_options, &option_index)) != -1) {
        switch (c) {
            case CMD_SIMTYPE:
                if (gk_optarg) {
                    if ((params->simtype = gk_GetStringID(simtype_options, gk_optarg)) == -1)
                        errexit("Invalid simtype of %s.\n", gk_optarg);
                }
                break;

            case CMD_ROWMODEL:
                if (gk_optarg) {
                    if ((params->rowmodel = gk_GetStringID(rowmodel_options, gk_optarg)) == -1)
                        errexit("Invalid rowmodel of %s.\n", gk_optarg);
                }
                break;

            case CMD_MODE:
                if (gk_optarg) {
                    if ((params->mode = gk_GetStringID(mode_options, gk_optarg)) == -1)
                        errexit("Invalid mode of %s.\n", gk_optarg);
                }
                break;

            case CMD_DIFFMODE:
                if (gk_optarg) {
                    if ((params->diffmode = gk_GetStringID(diffmode_options, gk_optarg)) == -1)
                        errexit("Invalid diffmode of %s.\n", gk_optarg);
                }
                break;

            case CMD_MINSPAN:
                if (gk_optarg) {
                    if ((params->minspan = atoi(gk_optarg)) < 1)
                        errexit("The -minspan must be greater than 0.\n");
                }
                break;

            case CMD_SCOST:
                if (gk_optarg) {
                    params->scost = atof(gk_optarg);
                    if (params->scost < 0 || params->scost > 1)
                        errexit("The -scost must be between [0...1].\n");
                }
                break;

            case CMD_DIFFSIMT:
                if (gk_optarg) {
                    params->simt = atof(gk_optarg);
                    if (params->simt < 0 || params->simt > 1)
                        errexit("The -simt must be between [0...1].\n");
                }
                break;

            case CMD_DIFFLDIM:
                if (gk_optarg) {
                    params->ldim = atoi(gk_optarg);
                    if (params->ldim < 1)
                        errexit("The -ldim must be greater than 0.\n");
                }
                break;

            case CMD_SVDMAXIT:
                if (gk_optarg) 
                    params->svdmaxit = atoi(gk_optarg);
                break;

            case CMD_MINTP:
                if (gk_optarg) {
                    params->mintp = atof(gk_optarg);
                    if (params->mintp < 0 || params->mintp > 1)
                        errexit("The -mintp must be between [0...1].\n");
                }
                break;

            case CMD_N2FRAC:
                if (gk_optarg) {
                    params->n2frac = atof(gk_optarg);
                    if (params->n2frac < 0 || params->n2frac > 1)
                        errexit("The -n2frac must be between [0...1].\n");
                }
                break;

            case CMD_NITERS:
                if (gk_optarg) {
                    if ((params->niters = atoi(gk_optarg)) < 1)
                        errexit("The -niters must be greater than 0.\n");
                }
                break;

            case CMD_DIFFGLOBAL:
                params->diffglobal = 1;
                break;

            case CMD_WRITEPATHS:
                params->writepaths = 1;
                break;

            case CMD_WRITEFTRS:
                params->writeftrs = 1;
                break;

            case CMD_SEED:
                if (gk_optarg) {
                    if(gk_optarg[0] == 't' || gk_optarg[0] == 'T')
                        params->seed = time(NULL);
                    else if ((params->seed = atoi(gk_optarg)) < 1)
                        errexit("The -seed must be greater than 0.\n");
                }
                break;

            case CMD_DBGLVL:
                if (gk_optarg) {
                    params->dbglvl = atoi(gk_optarg);
                    if (params->dbglvl < 0)
                        errexit("The -dbglvl must be non-negative.\n");
                }
                break;

            case CMD_OFILE:
                if(gk_optarg)
                    params->ofile = gk_strdup(gk_optarg);
                break;

            case CMD_NUMTHREADS:
                if (gk_optarg) {
                    if ((params->nthreads = atoi(gk_optarg)) < 1)
                        gk_errexit(SIGERR,  "The -nthreads must be greater than 0.\n");
                } else
                    gk_errexit(SIGERR, "-nthreads requires an argument.");
                break;

            case CMD_HELP:
                for (i=0; strlen(helpstr[i]) > 0; i++)
                    printf("%s\n", helpstr[i]);
                exit(EXIT_SUCCESS);
                break;

            default:
                printf("Illegal command-line option(s)\nUse %s -help for a summary of the options.\n", argv[0]);
                exit(EXIT_FAILURE);
        }
    }

    /* Get the operation to be performed */
    if (argc-gk_optind < 2) {
        printf("Missing required parameters.\n  Use %s -help for a summary of the options.\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    params->datafile = gk_strdup(argv[gk_optind++]);

    if ((params->nprotos = atoi(argv[gk_optind])) < 0)
        errexit("Invalid operation %s.\n", argv[gk_optind]);
    gk_optind++;

    params->dataset = get_dataset(params);
    if (params->ofile == NULL){
        /* build output file: <fstem>_<mode>_<nprotos>_<minspan>_<scost>[_<diff>_<ldim>_<simt>].seg */
        switch(params->mode){
            case MODE_ORION:
            case MODE_OTM:
                sprintf(fname, "%s_%s_%d_%d_%.3f.seg", params->dataset, modenames[params->mode],
                        params->nprotos, params->minspan, params->scost);
                break;
            case MODE_DIFF:
                sprintf(fname, "%s_%s_%d_%d_%.3f_%s_%d_%.3f.seg", params->dataset, modenames[params->mode],
                        params->nprotos, params->minspan, params->scost,
                        diffusionnames[params->diffmode], params->ldim, params->simt);
                break;
            default:
                gk_errexit(SIGERR, "Invalid mode chosen.");
        }
        params->ofile = gk_strdup(fname);
    }

    /* set number of threads */
    omp_set_num_threads(params->nthreads);

#pragma omp parallel
{
    if(omp_get_thread_num() == 0)
        params->nthreads   = omp_get_num_threads();
}
    /* set random seed */
    srand(params->seed);

    return params;
}


