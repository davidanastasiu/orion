/*!
\file
\brief This file contains function prototypes
\date Started 6/16/2014
\author George, David
*/

#ifndef _PROTO_H_
#define _PROTO_H_

/* io.c */
vault_t *loadData(params_t *params);

/* cmdline.c */
params_t *getcmdline_params(int argc, char *argv[]);

/* orion.c */
solution_t *cluster_orion(params_t *params, vault_t *vault);

/* diff.c */
solution_t *cluster_diff(params_t *params, vault_t *vault);

/* otm.c */
solution_t *cluster_otm(params_t *params, vault_t *vault);

/* protos.c */
void protos_PreprocessData(params_t *params, vault_t *vault);
solution_t *protos_FindInitial(params_t *params, vault_t *vault);
void protos_PrintDistances(params_t *params, vault_t *vault, solution_t *sol);
void protos_ComputeProtos(params_t *params, vault_t *vault, solution_t *sol);
double protos_ObjValue(params_t *params, vault_t *vault, solution_t *sol);
void protos_Print(params_t *params, vault_t *vault, solution_t *sol);
void protos_PrintKeyFtrs(params_t *params, vault_t *vault, solution_t *sol);
void protos_WriteKeyFtrs(params_t *params, vault_t *vault, solution_t *sol);
void protos_PrintP2PTMat(params_t *params, vault_t *vault, solution_t *sol);
void protos_WritePaths(params_t *params, vault_t *vault, solution_t *sol);
void protos_PrintP2PTInfo(params_t *params, vault_t *vault, solution_t *sol);

/* util.c */
char* get_dataset(params_t* const params);
void write_dmat(char* filename, const double* const mat, const int nrows, const int ncols);
void write_fmat(char* filename, const float** const mat, const int nrows, const int ncols);
float** csr2fmat(gk_csr_t* mat, float** fmat);
double* csr2dmat(gk_csr_t* mat, double* dmat);
solution_t* copy_output_sol(solution_t *sa, solution_t* sb, int nrows, int ncols, int nprotos);
char compute_svd_dmat(const double* const A, const int nrows, const int ncols, 
        const int nSig, const int svdmaxit,
        double *Ut, double *S, double *Vt);
double compute_svd_error_dmat(const double* const A, const int nrows, const int ncols, const int nSig,
        double *Ut, double *S, double *Vt);
char compute_svd_csr(const gk_csr_t* const A, const int nSig, const int svdmaxit,
        double *Ut, double *S, double *Vt);
double compute_svd_error_csr(const gk_csr_t* const pmat, const int nSig,
        double *Ut, double *S, double *Vt);
double* pairwise_sims_dmat(const double* const Ut, const int nrows, const int nSig, const double* S,
        double* psim);

#endif
