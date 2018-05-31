/*!
\file
\brief This file contains function prototypes
\date Started 6/16/2014
\author George
*/

#ifndef _PROTO_H_
#define _PROTO_H_


/* main.c */
void free_memory(params_t *params, vault_t *vault, solution_t *sol);

/* io.c */
vault_t *loadData(params_t *params);
char* da_getStringKey(gk_StringMap_t *strmap, const char id);
int da_getStringID(gk_StringMap_t *strmap, char *key);
char da_getFileFormat(char *file, const char format);

/* cmdline.c */
params_t *getcmdline_params(int argc, char *argv[]);

/* protos.c */
solution_t *protos_Cluster(params_t *params, vault_t *vault);
void protos_PreprocessData(params_t *params, vault_t *vault);
solution_t *protos_FindInitial(params_t *params, vault_t *vault);
void protos_Refine(params_t *params, vault_t *vault, solution_t *sol);
void protos_DecodeSequence(params_t *params, vault_t *vault, solution_t *sol,
         int u);
void protos_PrintDistances(params_t *params, vault_t *vault, solution_t *sol);
void protos_ComputeProtos(params_t *params, vault_t *vault, solution_t *sol);
double protos_ObjValue(params_t *params, vault_t *vault, solution_t *sol);
void protos_Print(params_t *params, vault_t *vault, solution_t *sol);
void protos_Write(params_t *params, vault_t *vault, solution_t *sol);
void protos_PrintKeyFtrs(params_t *params, vault_t *vault, solution_t *sol);
void protos_WriteKeyFtrs(params_t *params, vault_t *vault, solution_t *sol);
void protos_PrintP2PTMat(params_t *params, vault_t *vault, solution_t *sol);
void protos_WriteP2PTMat(params_t *params, vault_t *vault, solution_t *sol);
void protos_PrintP2PTInfo(params_t *params, vault_t *vault, solution_t *sol);
void protos_WriteP2PTInfo(params_t *params, vault_t *vault, solution_t *sol);
void protos_WriteTransitions(params_t *params, vault_t *vault, solution_t *sol);
void protos_WritePaths(params_t *params, vault_t *vault, solution_t *sol);

#endif
