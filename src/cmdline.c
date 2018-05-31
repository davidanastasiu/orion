/*!
\file  
\brief Parsing of command-line arguments
 
This file parses the command line arguments

\date   Started 6/16/14
\author George
\version\verbatim $Id: cmdline.c 17345 2014-06-21 19:48:20Z karypis $ \endverbatim
*/


#include "orion.h"

/*-------------------------------------------------------------------
 * Command-line options 
 *-------------------------------------------------------------------*/
static struct gk_option long_options[] = {
/* {"simtype",           1,      0,      CMD_SIMTYPE}, */

  {"ncliters",          1,      0,      CMD_NCLITERS},
  {"niters",            1,      0,      CMD_NITERS},
  {"rowmodel",          1,      0,      CMD_ROWMODEL},
  {"minspan",           1,      0,      CMD_MINSPAN},
  {"mintp",             1,      0,      CMD_MINTP},
  {"n2frac",            1,      0,      CMD_N2FRAC},
  {"scost",             1,      0,      CMD_SCOST},
  {"seed",              1,      0,      CMD_SEED},

  {"countsfile",        1,      0,      CMD_COUNTSFILE},
  {"clabelsfile",       1,      0,      CMD_CLABELSFILE},
  {"writepaths",        0,      0,      CMD_WRITEPATHS},
  {"pathsfile",         1,      0,      CMD_PATHSFILE},
  {"writeftrs",         0,      0,      CMD_WRITEFTRS},
  {"ftrsfile",          1,      0,      CMD_FTRSFILE},
  {"writetrans",        0,      0,      CMD_WRITETRANS},
  {"transfile",         1,      0,      CMD_TRANSFILE},
  {"translevel",        1,      0,      CMD_TRANSLEVEL},
  {"writep2pt",         0,      0,      CMD_WRITEP2PT},
  {"p2ptfile",          1,      0,      CMD_P2PTFILE},
  {"writeprotos",       0,      0,      CMD_WRITEPROTOS},
  {"protosfile",        1,      0,      CMD_PROTOSFILE},
  {"writetinfo",        0,      0,      CMD_WRITETINFO},
  {"tinfofile",         1,      0,      CMD_TINFOFILE},

  {"dbglvl",            1,      0,      CMD_DBGLVL},
  {"help",              0,      0,      CMD_HELP},
  {0,                   0,      0,      0}
};

/*-------------------------------------------------------------------
 * Mappings for the various parameter values
 *-------------------------------------------------------------------*/
gk_StringMap_t simtype_options[] = {
  {"cos",                SIMTYPE_COS},
  {"ejc",                SIMTYPE_EJC},
  {"sqe",                SIMTYPE_SQE},
  {NULL,                 0}
};


gk_StringMap_t rowmodel_options[] = {
  {"none",               ROWMODEL_NONE},
  {"log",                ROWMODEL_LOG},
  {"sqrt",               ROWMODEL_SQRT},
  {NULL,                 0}
};

gk_StringMap_t fmt_options[] = {
  {"clu",               GK_CSR_FMT_CLUTO},
  {"csr",               GK_CSR_FMT_CSR},
  {"ijv",               GK_CSR_FMT_IJV},
  {"binr",              GK_CSR_FMT_BINROW},
  {"binc",              GK_CSR_FMT_BINCOL},
  {"bijv",              GK_CSR_FMT_BIJV},
  {NULL,                 0}
};



/*-------------------------------------------------------------------
 * Mini help
 *-------------------------------------------------------------------*/
static char helpstr[][100] =
{
" ",
"Usage: orion [options] fstem nprotos",
" ",
" ",        
" Description",
" ",
"  Orion performs multivariate resource utilization time series evolution analysis.",
" ",
" ",
" Options",
" ",
"  -rowmodel=text",
"     Specifies how the values will be scaled.",
"     Possible values are:",
"        none    No scaling.",
"        log     Take the log of the raw values. [Default]",
"        sqrt    Take the square-root of the raw values.",
" ",
"  -minspan=int",
"     Specifies the minimum number of consecutive weeks that a proto must",
"     cover in a user's time-series.",
"     Default value is 5.",
" ",
"  -ncliters=int",
"     Specifies the maximum number of initial clustering iterations.",
"     Default value is 20.",
" ",
"  -niters=int",
"     Specifies the maximum number of segmentation refinement iterations.",
"     Default value is 20.",
" ",
"  -scost=float",
"     A per segment cost as a fraction of the error.",
"     Default value is .01.",
" ",
"  -mintp=float",
"     Specifies the minimum proto-to-proto transition probability to be analyzed.",
"     Default value is .20.",
" ",
"  -n2frac=float",
"     Specifies the maximum fraction of a vector to be analyzed for features.",
"     Default value is .80.",
" ",
"  -countsfile=text",
"     Required input file containing the number of rows in <fstem> for each sequence.",
"     Default value is <fstem>.counts.",
" ",
"  -clabelsfile=text",
"     Optional input file for the feature labels.",
"     Default value is <fstem>.clabels (labels will be generated if file does not exit).",
" ",
"  -writepaths",
"     Outputs the proto paths of each sequence.",
" ",
"  -pathsfile=text",
"     Output file for the proto paths of each sequence.",
"     Default value is <fstem>.paths.c<scost>.s<minspan>.p<nprotos>.",
" ",
"  -writeftrs",
"     Outputs the key features for each proto.",
" ",
"  -ftrsfile=text",
"     Output file for the key features for each proto.",
"     Default value is <fstem>.ftrs.c<scost>.s<minspan>.p<nprotos>.",
" ",
"  -writep2pt",
"     Outputs the global proto-to-proto transition probabilities.",
" ",
"  -p2ptfile=text",
"     Output file for the global proto-to-proto transition probabilities.",
"     Default value is <fstem>.p2pt.c<scost>.s<minspan>.p<nprotos>.",
" ",
"  -writetrans",
"     Outputs transition probabilities for a given transition level.",
" ",
"  -transfile=text",
"     Output file for the transition probabilities for a given transition level.",
"     Default value is <fstem>.trans.c<scost>.s<minspan>.p<nprotos>.l<translevel>.",
" ",
"  -translevel=int",
"     Transition level to output transitions for. Must be > 0.",
"     Default value is 1.",
" ",
"  -writetinfo",
"     Outputs frequent transitions and their key discriminative features.",
" ",
"  -tinfofile=text",
"     Output file for frequent transitions and their key discriminative features.",
"     Default value is <fstem>.tinfo.c<scost>.s<minspan>.p<nprotos>.",
" ",
"  -writeprotos",
"     Outputs the prototypical usage vectors.",
" ",
"  -protosfile=text",
"     Output file for the prototypical usage vectors.",
"     Default value is <fstem>.protos.c<scost>.s<minspan>.p<nprotos>.",
" ",
"  -seed=int",
"     Specifies the seed for clustering.",
"     Default value is 1.",
" ",
"  -dbglvl=int",
"     Specifies the level of debugging information to be displayed.",
"     Possible values are (2 and above can be combined in binary form):",
"         0  Disable all standard output.",
"         1  Print general progress and output information. [Default]",
"         2  Include the protos.",
"         4  Include per-sequence errors.",
"         8  Include dynamic programming algorithm progress information.",
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

    params = (params_t *)gk_malloc(sizeof(params_t), "cmdline_parse: params");
    memset(params, 0, sizeof(params_t));

    /* initialize the params data structure */
    params->nprotos     = -1;
    params->datafile    = NULL;
    params->countsfile  = NULL;
    params->clabelsfile = NULL;
    params->ftrsfile    = NULL;
    params->transfile   = NULL;
    params->pathsfile   = NULL;
    params->protosfile  = NULL;
    params->p2ptfile    = NULL;
    params->tinfofile   = NULL;

    params->minspan     = 5;
    params->scost       = .01;
    params->mintp       = .2;
    params->n2frac      = .8;
    params->ncliters    = 20;
    params->niters      = 20;
    params->translevel  = 1;
    params->simtype     = SIMTYPE_SQE;
    params->rowmodel    = ROWMODEL_LOG;

    params->writepaths  = 0;
    params->writeftrs   = 0;
    params->writetrans  = 0;
    params->writeprotos = 0;
    params->writep2pt   = 0;
    params->writetinfo  = 0;

    params->dbglvl      = 1;


    /* Parse the command line arguments  */
    while ((c = gk_getopt_long_only(argc, argv, "", long_options, &option_index)) != -1) {
        switch (c) {
            case CMD_SIMTYPE:
                if (gk_optarg) {
                    if ((params->simtype = gk_GetStringID(simtype_options, gk_optarg)) == -1)
                        errexit("Invalid simtype %s.\n", gk_optarg);
                }
                break;

            case CMD_ROWMODEL:
                if (gk_optarg) {
                    if ((params->rowmodel = gk_GetStringID(rowmodel_options, gk_optarg)) == -1)
                        errexit("Invalid rowmodel %s.\n", gk_optarg);
                }
                break;

            case CMD_MINSPAN:
                if (gk_optarg) {
                    if ((params->minspan = atoi(gk_optarg)) < 1)
                        errexit("The -minspan value must be greater than 0.\n");
                }
                break;

            case CMD_SCOST:
                if (gk_optarg) {
                    params->scost = atof(gk_optarg);
                    if (params->scost < 0 || params->scost > 1)
                        errexit("The -scost value must be between [0...1].\n");
                }
                break;

            case CMD_MINTP:
                if (gk_optarg) {
                    params->mintp = atof(gk_optarg);
                    if (params->mintp < 0 || params->mintp > 1)
                        errexit("The -mintp value must be between [0...1].\n");
                }
                break;

            case CMD_N2FRAC:
                if (gk_optarg) {
                    params->n2frac = atof(gk_optarg);
                    if (params->n2frac < 0 || params->n2frac > 1)
                        errexit("The -n2frac value must be between [0...1].\n");
                }
                break;

            case CMD_NCLITERS:
                if (gk_optarg) {
                    if ((params->ncliters = atoi(gk_optarg)) < 1)
                        errexit("The -ncliters value must be greater than 0.\n");
                }
                break;

            case CMD_NITERS:
                if (gk_optarg) {
                    if ((params->niters = atoi(gk_optarg)) < 1)
                        errexit("The -niters value must be greater than 0.\n");
                }
                break;

            case CMD_COUNTSFILE:
                params->countsfile = gk_strdup(gk_optarg);
                break;

            case CMD_CLABELSFILE:
                params->clabelsfile = gk_strdup(gk_optarg);
                break;

            case CMD_WRITEPATHS:
                params->writepaths = 1;
                break;

            case CMD_PATHSFILE:
                params->pathsfile = gk_strdup(gk_optarg);
                break;

            case CMD_WRITEFTRS:
                params->writeftrs = 1;
                break;

            case CMD_FTRSFILE:
                params->ftrsfile = gk_strdup(gk_optarg);
                break;

            case CMD_WRITETRANS:
                params->writetrans = 1;
                break;

            case CMD_TRANSFILE:
                params->transfile = gk_strdup(gk_optarg);
                break;

            case CMD_TRANSLEVEL:
                if (gk_optarg) {
                    if ((params->translevel = atoi(gk_optarg)) < 1)
                        errexit("The -translevel value must be greater than 0.\n");
                }
                break;

            case CMD_WRITEP2PT:
                params->writep2pt = 1;
                break;

            case CMD_P2PTFILE:
                params->p2ptfile = gk_strdup(gk_optarg);
                break;

            case CMD_WRITEPROTOS:
                params->writeprotos = 1;
                break;

            case CMD_PROTOSFILE:
                params->protosfile = gk_strdup(gk_optarg);
                break;

            case CMD_WRITETINFO:
                params->writetinfo = 1;
                break;

            case CMD_TINFOFILE:
                params->tinfofile = gk_strdup(gk_optarg);
                break;

            case CMD_SEED:
                if (gk_optarg) {
                    if ((params->seed = atoi(gk_optarg)) < 1)
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


    /* print the command line */
    if (params->dbglvl>0){
        for (i=0; i<argc; i++)
            printf("%s ", argv[i]);
        printf("\n");
    }

    /* Get the operation to be performed */
    if (argc-gk_optind < 2) {
        printf("Missing required parameters.\n  Use %s -help for a summary of the options.\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    params->datafile = gk_strdup(argv[gk_optind++]);

    if ((params->nprotos = atoi(argv[gk_optind])) < 0){
        printf("Invalid operation %s.\n\n", argv[gk_optind]);
        for (i=0; strlen(helpstr[i]) > 0; i++)
            printf("%s\n", helpstr[i]);
        exit(EXIT_FAILURE);
    }
    gk_optind++;

    return params;
}


