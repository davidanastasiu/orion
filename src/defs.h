/*!
\file
\brief This file contains various constant definitions
\date Started 6/16/14
\author George
\version\verbatim $Id: defs.h 17345 2014-06-21 19:48:20Z karypis $ \endverbatim
*/

#ifndef __DEF_H__
#define __DEF_H__


#define PROGRAM_NAME        "orion"
#define VER_MAJOR           1
#define VER_MINOR           0
#define VER_SUBMINOR        0
#define VER_COMMENT         "initial public release"

#define MAXLINE         1024*128
#define MAX_STRLEN      1024*128

#define MFACTOR         10000.0
#define SPOWER          1.0

#define PRUNEFRACTION   0.75
#define PRUNEFRACTION2  0.75

#define NCOMMON         20

#define MINCOS          .40


/* command-line options */
#define CMD_SIMTYPE             1
#define CMD_NCLITERS            2
#define CMD_NITERS              3
#define CMD_ROWMODEL            4
#define CMD_MINSPAN             5
#define CMD_MINTP               6
#define CMD_N2FRAC              7
#define CMD_SCOST               8
#define CMD_WRITEPATHS          52
#define CMD_WRITEFTRS           53
#define CMD_WRITETRANS          54
#define CMD_WRITEP2PT           55
#define CMD_WRITEPROTOS         56
#define CMD_WRITETINFO          57
#define CMD_COUNTSFILE          60
#define CMD_CLABELSFILE         61
#define CMD_PATHSFILE           62
#define CMD_FTRSFILE            63
#define CMD_TRANSFILE           64
#define CMD_P2PTFILE            65
#define CMD_PROTOSFILE          66
#define CMD_TINFOFILE           67
#define CMD_TRANSLEVEL          70
#define CMD_SEED                80
#define CMD_DBGLVL              100
#define CMD_HELP                200

/* simtypes */
#define SIMTYPE_COS             1
#define SIMTYPE_EJC             2
#define SIMTYPE_SQE             3

/* rowmodels */
#define ROWMODEL_NONE           1
#define ROWMODEL_LOG            2
#define ROWMODEL_SQRT           3


/* The text labels for the different simtypes */
static char simtypenames[][10] = 
                {"", "cos", "ejc", "sqe", ""}; 

/* The text labels for the different rowmodels */
static char rowmodelnames[][10] = 
                {"", "none", "log", "sqrt", ""}; 


extern gk_StringMap_t simtype_options[];
extern gk_StringMap_t rowmodel_options[];
extern gk_StringMap_t fmt_options[];


#endif

