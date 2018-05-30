/*!
\file
\brief This file contains various constant definitions
\date Started 6/16/14
\author George, David
 */

#ifndef __DEF_H__
#define __DEF_H__


#define PROGRAM_NAME        "orion+"
#define VER_MAJOR           0
#define VER_MINOR           1
#define VER_SUBMINOR        0
#define VER_COMMENT         "dev"

//#define NOMP

#define MAXLINE         1024*128
#define MAX_STRLEN      1024*128

#define MFACTOR         10000.0
#define SPOWER          1.0

#define PRUNEFRACTION   0.75
#define PRUNEFRACTION2  0.75

#define NCOMMON         20

#define MINCOS          .40

#define NIMPROVE        10 /* stop iterating unless solution has improved in max NIMPROVE last iterations */

#define LSI_SIM_SCALE           1
#define LSI_SIM_ABS             2
#define LSI_SIM_TRUNC           3
#define LSI_SIM_MINMAX          4
#define LSI_SIM_NONE            5
#define LSI_SIM                 LSI_SIM_TRUNC

#define LSI_SHOW_ERROR          0


/* command-line options */
#define CMD_SIMTYPE             1
#define CMD_NITERS              2
#define CMD_ROWMODEL            3
#define CMD_MINSPAN             4
#define CMD_MINTP               5
#define CMD_N2FRAC              6
#define CMD_SCOST               7
#define CMD_MODE                30
#define CMD_DIFFMODE            31
#define CMD_DIFFSIMT            32
#define CMD_DIFFLDIM            33
#define CMD_DIFFGLOBAL          35
#define CMD_SVDMAXIT            39
#define CMD_OFILE               40
#define CMD_WRITEPATHS          50
#define CMD_WRITEFTRS           51
#define CMD_SEED                70
#define CMD_DBGLVL              100
#define CMD_NUMTHREADS          190
#define CMD_HELP                200

/* Execution modes */
#define MODE_ORION              1   /* Original Orion */
#define MODE_DIFF               2   /* Orion with diffused usage vectors */
#define MODE_OTM                3   /* Micro-proto topic modeling */

/* Diffusion modes */
#define DIFFMODE_PW             1   /* pairwise similarities, S = X*X^T */
#define DIFFMODE_TPW            2   /* thresholded pairwise similarities, S = L2AP(S, t) */
#define DIFFMODE_DR             3   /* Dimensionality reduction */
#define DIFFMODE_TDR            4   /* Dimensionality reduction + thresholded pairwise similarities */
#define DIFFMODE_PRW            5   /* Personalized random walk */

/* simtypes */
#define SIMTYPE_COS             1
#define SIMTYPE_EJC             2
#define SIMTYPE_SQE             3

/* rowmodels */
#define ROWMODEL_NONE           1
#define ROWMODEL_LOG            2
#define ROWMODEL_SQRT           3


/* The text labels for the different simtypes */
static char simtypenames[][10] = {"", "cos", "ejc", "sqe", ""};

/* The text labels for the different rowmodels */
static char rowmodelnames[][10] = {"", "none", "log", "sqrt", ""};

/* The text labels for the different simtypes */
static char modenames[][10] = {"", "orion", "diff", "otm", ""};

/* The text labels for the different simtypes */
static char diffusionnames[][10] = {"", "pw", "tpw", "dr", "tdr", "prw", ""};


#endif

