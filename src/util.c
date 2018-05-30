/*!
\file
\brief A number of utility functions shared among code modules.
\date Started 10/04/2014
\author David
 */
#include "orion.h"



/**
 * Retrieve the dataset (file - extension) from the input filename
 */
char* get_dataset(params_t* const params)
{
    if(params->dataset)
        return params->dataset;
    int len;
    char *end, *start = strrchr(params->datafile, '/');
    if(start)
        start++;
    else
        start = params->datafile;
    end = strstr(start, ".");
    if(!end)
        end = params->datafile + strlen(params->datafile);
    len = end - start;
    params->dataset = gk_cmalloc(len+1, "get_dataset: dataset");
    strncpy(params->dataset, start, len);
    params->dataset[len] = '\0';
    return params->dataset;
}


/**
 * Transfer sparse matrix to dense one
 */
float** csr2fmat(gk_csr_t* mat, float** fmat)
{
    int nrows, ncols, i, j;

    assert(mat->rowptr || mat->colptr);
    nrows = mat->nrows;
    ncols = mat->ncols;

    if(fmat == NULL)
        fmat = gk_fAllocMatrix(nrows, ncols, 0.0, "csr2fmat: fmat");
    else
        gk_fSetMatrix(fmat, nrows, ncols, 0.0);

    if(mat->rowptr){
        for(i=0; i < nrows; ++i)
            for(j=mat->rowptr[i]; j < mat->rowptr[i+1]; ++j)
                fmat[i][mat->rowind[j]] = mat->rowval[j];
    } else if (mat->colptr){
        for(i=0; i < ncols; ++i)
            for(j=mat->colptr[i]; j < mat->colptr[i+1]; ++j)
                fmat[mat->colind[j]][i] = mat->colval[j];
    } else
        gk_errexit(SIGERR, "csr2fmat: no ptr data in input matrix.");

    return fmat;
}


/**
 * Transfer sparse matrix to dense one stored as 1-D array
 */
double* csr2dmat(gk_csr_t* mat, double* dmat)
{
    int nrows, ncols, i, j;

    assert(mat->rowptr || mat->colptr);
    nrows = mat->nrows;
    ncols = mat->ncols;

    if(dmat == NULL)
        dmat = gk_dsmalloc(nrows * ncols, 0.0, "csr2dmat: dmat");
    else
        gk_dset(nrows * ncols, 0.0, dmat);

    if(mat->rowptr){
        for(i=0; i < nrows; ++i)
            for(j=mat->rowptr[i]; j < mat->rowptr[i+1]; ++j)
                dmat[i*mat->ncols + mat->rowind[j]] = mat->rowval[j];
    } else if (mat->colptr){
        for(i=0; i < ncols; ++i)
            for(j=mat->colptr[i]; j < mat->colptr[i+1]; ++j)
                dmat[mat->colind[j]*mat->nrows + i] = mat->colval[j];
    } else
        gk_errexit(SIGERR, "csr2dmat: no ptr data in input matrix.");

    return dmat;
}


/**
 * Write a matrix (stored as a 1-D array of doubles) to a file
 */
void write_dmat(char* filename, const double* const mat, const int nrows, const int ncols)
{
    ssize_t i, j;
    FILE* fpout = NULL;
    if (filename)
        fpout = gk_fopen(filename, "w", "write_dmat: fpout");
    else
        fpout = stdout;

    for(i=0; i < nrows; ++i){
        for(j=0; j < ncols; ++j)
            fprintf(fpout, "%.17g ", mat[i*ncols+j]);
        fprintf(fpout, "\n");
    }

    if (filename)
        gk_fclose(fpout);
}


/**
 * Write a matrix (stored as a 2-D array of floats) to a file
 */
void write_fmat(char* filename, const float** const mat, const int nrows, const int ncols)
{
    ssize_t i, j;
    FILE* fpout = NULL;
    if (filename)
        fpout = gk_fopen(filename, "w", "write_fmat: fpout");
    else
        fpout = stdout;

    for(i=0; i < nrows; ++i){
        for(j=0; j < ncols; ++j)
            fprintf(fpout, "%.17g ", mat[i][j]);
        fprintf(fpout, "\n");
    }

    if (filename)
        gk_fclose(fpout);
}


/**
 * Copy only the parts of the solution that are needed for output
 * from sa to sb.
 */
solution_t* copy_output_sol(solution_t *sa, solution_t* sb, int nrows, int ncols, int nprotos)
{
    size_t i;
    if(!sb){
        /* allocate memory for the solution */
        sb = (solution_t *)gk_malloc(sizeof(solution_t), "protos_FindInitial: sol");
        memset(sb, 0, sizeof(solution_t));
        sb->protos  = gk_fAllocMatrix(nprotos, ncols, 0.0, "copy_output_col: sb->protos");
        sb->vpart   = gk_imalloc(nrows, "vpart");
        sb->pnorms  = gk_fmalloc(nprotos, "pnorms");
        sb->pcounts = gk_imalloc(nprotos, "pcounts");
    }
    sb->objval  = sa->objval;
    sb->vpart   = memcpy(sb->vpart, sa->vpart, nrows * sizeof(int));
    sb->pcounts = memcpy(sb->pcounts, sa->pcounts, nprotos * sizeof(int));
    sb->pnorms  = memcpy(sb->pnorms, sa->pnorms, nprotos * sizeof(float));
    for(i=0; i < nprotos; ++i)
        sb->protos[i]  = memcpy(sb->protos[i], sa->protos[i], ncols * sizeof(float));

    return sb;
}




/**
 * Compute the Singular Value Decomposition of dense matrix A using the
 * SVDLIBC library and store results in Ut, S, and V^t.
 * \param A is the matrix to be decomposed, say of size mxn, stored as row-wise 1-D array
 * \param nrows Number of rows in A (m)
 * \param ncols Number of cols in A (n)
 * \param nSig Number of singular values to return (dictates size of U and V^t as well)
 *      nSig in [1, min(A->nrows, A->ncols)]
 * \param svdmaxit Max number of iterations for SVD
 * \param Ut pointer to a 1-D array to hold the resulting Ut component.
 *      If pointing to NULL, Ut will not be included in the result. Size: nSig x m
 * \param S the singular values (diagonal values of Sigma).
 *      If pointing to NULL, S will not be included in the result. Length: nSig
 * \param Vt pointer to a dense matrix to hold the resulting V^t component.
 *      If pointing to NULL, V^t will not be included in the result. Size: nSig x n
 */
char compute_svd_dmat(const double* const A, const int nrows, const int ncols, 
        const int nSig, const int svdmaxit,
        double *Ut, double *S, double *Vt)
{
    SMat sA = NULL;
    SVDRec R = NULL;
    DMat dA;
    ssize_t i;
    double las2end[2] = {-1.0e-30, 1.0e-30};
    double kappa = 1e-6;
    double residual = 0.0;

    assert(Ut != NULL || Vt != NULL || S != NULL);
    assert(A && nrows > 0 && ncols > 0);
    dA = (DMat) gk_malloc(sizeof(struct dmat),
            "da_dmat_ComputeSVD: temp dmat representation of A");
    dA->value = (double **) gk_malloc(nrows * sizeof(double *), "dA->value");
    dA->rows = nrows;
    dA->cols = ncols;
    dA->value[0] = (double*) A;
    for (i = 1; i < nrows; i++)
        dA->value[i] = dA->value[i-1] + ncols;

    /* get sparse representation of the matrix for the LAS2 Algo */
    sA = svdConvertDtoS(dA);
    gk_free((void**)&dA->value, &dA, LTERM);

    SVDVerbosity = 0;

    /** execute the decomposition */
    if(!(R = svdLAS2(sA, nSig, svdmaxit, las2end, kappa))){
        svdFreeSMat(sA);
        svdFreeSVDRec(R);
        return 0;
    }
    /* transfer results */
    //    Size Ut: Ut->nrows == nSig && Ut->ncols == A->nrows
    if(Ut){
        if(R->d == nSig){
            gk_dcopy(nSig * nrows, R->Ut->value[0], Ut);
        } else {
            memset(Ut, 0, nSig * nrows * sizeof(double));
            for(i=0; i < R->d; ++i)
                gk_dcopy(nrows, R->Ut->value[i], Ut + i*nrows);
        }
    }
    if(S){
        gk_dcopy(R->d, R->S, S);
        for(i=R->d; i < nSig; ++i)
            S[i] = 0; /* in case less than desired dimensions are returned */
    }
    //    Size Vt: Vt->nrows == nSig && Vt->ncols == A->ncols
    if(Vt){
        if(R->d == nSig){
            gk_dcopy(nSig * ncols, R->Vt->value[0], Vt);
        } else {
            memset(Vt, 0, nSig * ncols * sizeof(double));
            for(i=0; i < R->d; ++i)
                gk_dcopy(ncols, R->Vt->value[i], Vt);
        }
    }

    svdFreeSMat(sA);
    svdFreeSVDRec(R);

    return 1;
}

/**
 * Compute the reconstruction error of the Singular Value Decomposition of dense
 * matrix A, i.e. A=U*S*V^t, Err = ||A-U*S*V^t||_F^2.
 * \param A is the matrix to be decomposed, say of size mxn, stored as row-wise 1-D array
 * \param nrows Number of rows in A (m)
 * \param ncols Number of cols in A (n)
 * \param nSig Number of singular values to return (dictates size of U and V^t as well)
 *      nSig in [1, min(A->nrows, A->ncols)]
 * \param Ut pointer to a 1-D array to hold the resulting Ut component.
 *      If pointing to NULL, Ut will not be included in the result. Size: nSig x m
 * \param S the singular values (diagonal values of Sigma).
 *      If pointing to NULL, S will not be included in the result. Length: nSig
 * \param Vt pointer to a dense matrix to hold the resulting V^t component.
 *      If pointing to NULL, V^t will not be included in the result. Size: nSig x n
 */
double compute_svd_error_dmat(const double* const A, const int nrows, const int ncols, const int nSig,
        double *Ut, double *S, double *Vt)
{
    ssize_t i, j, k, l;
    double sum, vX;

    for(sum=0.0, i=0; i < nrows; ++i){
        for(j=0; j < ncols; ++j){
            for(vX=0.0, k=0; k < nSig; ++k)
                vX += Ut[k*nrows+i] * S[k] * Vt[k*ncols+j];
            sum += (A[i*ncols + j] - vX) * (A[i*ncols + j] - vX);
        }
    }
    return sum;
}



/**
 * Compute the Singular Value Decomposition of sparse CSR matrix A using the
 * SVDLIBC library and store results in Ut, S, and V^t.
 * \param A is the matrix to be decomposed, say of size mxn, stored as CSR
 * \param nSig Number of singular values to return (dictates size of U and V^t as well)
 *      nSig in [1, min(A->nrows, A->ncols)]
 * \param svdmaxit Max number of iterations for SVD
 * \param Ut pointer to a 1-D array to hold the resulting Ut component.
 *      If pointing to NULL, Ut will not be included in the result. Size: nSig x m
 * \param S the singular values (diagonal values of Sigma).
 *      If pointing to NULL, S will not be included in the result. Length: nSig
 * \param Vt pointer to a dense matrix to hold the resulting V^t component.
 *      If pointing to NULL, V^t will not be included in the result. Size: nSig x n
 */
char compute_svd_csr(const gk_csr_t* const A, const int nSig, const int svdmaxit,
        double *Ut, double *S, double *Vt)
{
    SMat sA = NULL;
    SVDRec R = NULL;
    ssize_t i, nrows, ncols, nnz;
    double las2end[2] = {-1.0e-30, 1.0e-30};
    double kappa = 1e-6;
    double residual = 0.0;
    gk_csr_t* tmp;

    assert(Ut != NULL || Vt != NULL || S != NULL);
    assert(A && A->nrows > 0 && A->ncols > 0 && (A->rowptr || A->colptr));

    nrows = A->nrows;
    ncols = A->ncols;
    nnz = A->rowptr ? A->rowptr[A->nrows] : A->colptr[A->ncols];
    sA = svdNewSMat(nrows, ncols, nnz);
    if(A->colptr){
        for(i=0; i < nnz; ++i)
            sA->rowind[i] = A->colind[i];
        for(i=0; i < nnz; ++i)
            sA->value[i] = A->colval[i];
        for(i=0; i < ncols+1; ++i)
            sA->pointr[i] = A->colptr[i];
    } else if(A->rowptr){
        tmp = gk_csr_Dup((gk_csr_t*) A);
        gk_csr_CreateIndex(tmp, GK_CSR_COL);
        gk_csr_SortIndices(tmp, GK_CSR_COL);
        for(i=0; i < nnz; ++i)
            sA->rowind[i] = tmp->colind[i];
        for(i=0; i < nnz; ++i)
            sA->value[i] = tmp->colval[i];
        for(i=0; i < ncols+1; ++i)
            sA->pointr[i] = tmp->colptr[i];
        gk_csr_Free(&tmp);
    }

    SVDVerbosity = 0;

    /** execute the decomposition */
    if(!(R = svdLAS2(sA, nSig, svdmaxit, las2end, kappa))){
        svdFreeSMat(sA);
        svdFreeSVDRec(R);
        return 0;
    }

    /* transfer results */
    //    Size Ut: Ut->nrows == nSig && Ut->ncols == A->nrows
    if(Ut){
        if(R->d == nSig){
            gk_dcopy(nSig * nrows, R->Ut->value[0], Ut);
        } else {
            memset(Ut, 0, nSig * nrows * sizeof(double));
            for(i=0; i < R->d; ++i)
                gk_dcopy(nrows, R->Ut->value[i], Ut + i*nrows);
        }
    }
    if(S){
        gk_dcopy(R->d, R->S, S);
        for(i=R->d; i < nSig; ++i)
            S[i] = 0.0; /* in case less than desired dimensions are returned */
    }
    //    Size Vt: Vt->nrows == nSig && Vt->ncols == A->ncols
    if(Vt){
        if(R->d == nSig){
            gk_dcopy(nSig * ncols, R->Vt->value[0], Vt);
        } else {
            memset(Vt, 0, nSig * ncols * sizeof(double));
            for(i=0; i < R->d; ++i)
                gk_dcopy(ncols, R->Vt->value[i], Vt + i*ncols);
        }
    }

    svdFreeSMat(sA);
    svdFreeSVDRec(R);

    return 1;
}


/**
 * Compute the reconstruction error of the Singular Value Decomposition of sparse
 * CSR matrix A, i.e. A=U*S*V^t, Err = ||A-U*S*V^t||_F^2.
 * \param A is the matrix to be decomposed, say of size mxn, stored as CSR
 * \param nSig Number of singular values to return (dictates size of U and V^t as well)
 *      nSig in [1, min(A->nrows, A->ncols)]
 * \param Ut pointer to a 1-D array to hold the resulting Ut component.
 *      If pointing to NULL, Ut will not be included in the result. Size: nSig x m
 * \param S the singular values (diagonal values of Sigma).
 *      If pointing to NULL, S will not be included in the result. Length: nSig
 * \param Vt pointer to a dense matrix to hold the resulting V^t component.
 *      If pointing to NULL, V^t will not be included in the result. Size: nSig x n
 */
double compute_svd_error_csr(const gk_csr_t* const pmat, const int nSig,
        double *Ut, double *S, double *Vt)
{
    ssize_t i, j, k, l, e;
    double sum, vX;

    sum = 0.0;
    if(pmat->rowptr){
        for(i=0; i < pmat->nrows; ++i){
            for(e=0, l=pmat->rowptr[i]; l < pmat->rowptr[i+1]; ++l){
                j = pmat->rowind[l];
                while(e < j){
                    for(vX=0.0, k=0; k < nSig; ++k)
                        vX += Ut[k*pmat->nrows+i] * S[k] * Vt[k*pmat->ncols+e];
                    sum += vX*vX;
                    e++;
                }
                for(vX=0.0, k=0; k < nSig; ++k)
                    vX += Ut[k*pmat->nrows+i] * S[k] * Vt[k*pmat->ncols+j];
                sum += (pmat->rowval[l] - vX) * (pmat->rowval[l] - vX);
                e++;
            }

            while(e < pmat->ncols){
                for(vX=0.0, k=0; k < nSig; ++k)
                    vX += Ut[k*pmat->nrows+i] * S[k] * Vt[k*pmat->ncols+e];
                sum += vX*vX;
                e++;
            }
        }
    } else if(pmat->colptr) {
        for(j=0; j < pmat->ncols; ++j){
            for(e=0, l=pmat->colptr[j]; l < pmat->colptr[j+1]; ++l){
                i = pmat->colind[l];
                while(e < i){
                    for(vX=0.0, k=0; k < nSig; ++k)
                        vX += Ut[k*pmat->nrows+e] * S[k] * Vt[k*pmat->ncols+j];
                    sum += vX*vX;
                    e++;
                }
                for(vX=0.0, k=0; k < nSig; ++k)
                    vX += Ut[k*pmat->nrows+i] * S[k] * Vt[k*pmat->ncols+j];
                sum += (pmat->colval[l] - vX) * (pmat->colval[l] - vX);
                e++;
            }
            while(e < pmat->nrows){
                for(vX=0.0, k=0; k < nSig; ++k)
                    vX += Ut[k*pmat->nrows+e] * S[k] * Vt[k*pmat->ncols+j];
                sum += vX*vX;
                e++;
            }
        }
    } else
        gk_errexit(SIGERR, "compute_svd_error_csr: input matrix does not have data.");

    return sum;
}

/**
 * Compute pairwise cosine similarities between the lower-dimensional representation of rows of the matrix (U*S)
 * given results from the SVD decomposition.
 * \param Ut pointer to a 1-D array to hold the resulting Ut component. Size: nSig x nrows
 * \param nrows Number of rows in original matrix (also number of rows in U, but number of columns in Ut)
 * \param nSig Number of latent dimensions used in SVD (also number of cols in U, but number of rows in Ut)
 * \param S Singular values returned by SVD, length: nSig
 * \param psim matrix as 1-D array to hold pairwise similarities in (will be created if NULL)
 */
double* pairwise_sims_dmat(const double* const Ut, const int nrows, const int nSig, const double* S,
        double* psim)
{
    ssize_t i, j, k, n, nnz;
    double sum, v, min, max, mean, lmean, var;
    double *sim = NULL, *mat, *m;

    if(psim){
        sim = psim;
    } else {
        sim = gk_dmalloc(nrows * nrows, "dmat[p]");
    }
    /* zero-out the pairwise similarities */
    memset(sim, 0, nrows*nrows*sizeof(double));

    /* create row-based l2-normalized matrix to speed up cosine sim computations */
    mat = gk_dmalloc(nrows * nSig, "pairwise_sims_dmat: mat");
    for(i=0; i < nrows; ++i)
        for(j=0; j < nSig; ++j)
            mat[i*nSig+j] = Ut[j*nrows+i] * S[j];
    for(i=0; i < nrows; ++i){
        for(m=mat+(i*nSig), v=0.0, j=0; j < nSig; ++j)
            v += m[j] * m[j];
        if(v != 0)
            gk_dscale(nSig, 1.0/sqrt(v), mat + (i*nSig), 1);
    }

    for(min=DBL_MAX, max=-DBL_MAX, i=0; i < nrows; ++i){
        for(j=0; j < i; ++j){
            v = gk_ddot(nSig, mat+(i*nSig), 1, mat+(j*nSig), 1);
#if LSI_SIM == LSI_SIM_SCALE
            sim[i*nrows+j] = sim[j*nrows+i] = (1.0+v)/2;
#elif LSI_SIM == LSI_SIM_ABS
            sim[i*nrows+j] = sim[j*nrows+i] = v >= 0.0 ? v : -v;
#elif LSI_SIM == LSI_SIM_TRUNC
            sim[i*nrows+j] = sim[j*nrows+i] = v > 0.0 ? v : 0.0;
#elif LSI_SIM == LSI_SIM_MINMAX
            sim[i*nrows+j] = sim[j*nrows+i] = v;
            if(v > max)
                max = v;
            if(v < min)
                min = v;
#elif LSI_SIM == LSI_SIM_NONE
            sim[i*nrows+j] = sim[j*nrows+i] = v;
#else
            gk_errexit(SIGERR, "undefined LSI_SIM in defs.h");
#endif
        }
        sim[i*nrows+i] = 1.0;
        if(1.0 > max)
            max = 1.0;
        if(1.0 < min)
            min = 1.0;

    }

#if LSI_SIM == LSI_SIM_MINMAX   /* min-max normalize */
    for(v=max-min, i=0; i < nrows; ++i){
        for(j=0; j < i; ++j){
            sim[i*nrows+j] = sim[j*nrows+i] = (sim[j*nrows+i] - min)/v;
        }
        sim[i*nrows+i] = (1.0-min)/v;
    }
#endif

    /* compute statistics on similarities */
    for(mean=sim[0], var=0.0, min=DBL_MAX, max=-DBL_MAX, n=1, nnz=0, i=0; i < nrows; ++i){
        for(j=0; j < i; ++j){
            v = sim[i*nrows+j];
            if(v > max)
                max = v;
            if(v < min)
                min = v;
            if(v != 0.0)
                nnz+=2;
            lmean = mean;
            mean += (1.0/n++)*(v-mean);
            var += (v - lmean)*(v - mean);
            lmean = mean;
            mean += (1.0/n++)*(v-mean);
            var += (v - lmean)*(v - mean);
        }
        if(1.0 > max)
            max = 1.0;
        if(1.0 < min)
            min = 1.0;
        nnz++;
        lmean = mean;
        mean += (1.0/n++)*(1.0-mean);
        var += (v - lmean)*(1.0 - mean);
    }
    var /= n;

    printf("\tsim stats - nnz: %zu, min: %.3f, max: %.3f, mean: %.3f, stdev: %.5f\n",
            nnz, min, max, mean, sqrt(var));

    gk_free((void**)&mat, LTERM);

    return sim;
}
