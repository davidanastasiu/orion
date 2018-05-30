/*!
\file
\brief Contains the core routines for proto-based clustering of sequences
\date Started 6/16/2014
\author George, David
 */

#include "orion.h"




/*************************************************************************/
/*! Preprocesses the input data */
/*************************************************************************/
void protos_PreprocessData(params_t *params, vault_t *vault)
{
    ssize_t ir, ic, id, it, nnz1, nnz2, nnz3, nrows, nnz;
    ssize_t *rowptr;
    int *rowind;
    float *rowval;

    switch (params->rowmodel) {
        case ROWMODEL_NONE:
            break;
        case ROWMODEL_LOG:
            gk_csr_Scale(vault->mat, GK_CSR_LOG);
            break;
        case ROWMODEL_SQRT:
            gk_csr_Scale(vault->mat, GK_CSR_SQRT);
            break;
    }

    switch (params->simtype) {
        case SIMTYPE_COS:
            gk_csr_Normalize(vault->mat, GK_CSR_ROW, 2);
            break;
        case SIMTYPE_EJC:
            gk_csr_ComputeSums(vault->mat, GK_CSR_ROW);
            break;
        case SIMTYPE_SQE:
            gk_csr_ComputeSquaredNorms(vault->mat, GK_CSR_ROW);
            break;
    }

}


/*************************************************************************/
/*! Finds the initial set of protos */
/*************************************************************************/
solution_t *protos_FindInitial(params_t *params, vault_t *vault)
{
    ssize_t i, j, k, iter;
    int nrows, ncols, nprotos, nconsecutive;
    ssize_t *rowptr;
    int *rowind, *part;
    float *rowval, *rnorms, **protos, *pnorms;
    float minsqe, sqe, totsqe, scale;
    solution_t *sol;

    nprotos = params->nprotos;

    nrows  = vault->mat->nrows;
    ncols  = vault->mat->ncols;
    rowptr = vault->mat->rowptr;
    rowind = vault->mat->rowind;
    rowval = vault->mat->rowval;
    rnorms = vault->mat->rnorms;
    minsqe = FLT_MAX;

    /* allocate memory for the solution */
    sol = (solution_t *)gk_malloc(sizeof(solution_t), "protos_FindInitial: sol");
    memset(sol, 0, sizeof(solution_t));

    sol->protos  = gk_fAllocMatrix(nprotos, ncols, 0.0, "protos_FindInitial: sol->protos");
    sol->vpart   = gk_imalloc(nrows, "vpart");
    sol->pnorms  = gk_fmalloc(nprotos, "pnorms");
    sol->pcounts = gk_imalloc(nprotos, "pcounts");
    sol->sqes    = gk_fmalloc(vault->nseqs, "sqes");

    protos = sol->protos;
    part   = sol->vpart;
    pnorms = sol->pnorms;


    /* pick initial centers randomly */
    gk_irandArrayPermute(nrows, part, 10*nrows, 1);
    for (k=0; k<nprotos; k++) {
        gk_fset(ncols, 0.0, protos[k]);
        pnorms[k] = 0.0;

        i = part[k];
        for (j=rowptr[i]; j<rowptr[i+1]; j++) {
            protos[k][rowind[j]] = rowval[j];
            pnorms[k] += rowval[j]*rowval[j];
        }
    }


    /* perform a small number of k-means iterations */
    for (iter=0; iter<20; iter++) {
        /* determine the new partition of the rows */
        nconsecutive = 0;
        totsqe = 0.0;
        for (i=0; i<nrows; i++) {
            part[i] = -1;
            for (k=0; k<nprotos; k++) {
                sqe = rnorms[i]+pnorms[k];
                for (j=rowptr[i]; j<rowptr[i+1]; j++)
                    sqe -= 2*rowval[j]*protos[k][rowind[j]];
                if (part[i] == -1 || sqe < minsqe) {
                    minsqe  = sqe;
                    part[i] = k;
                }
            }
            totsqe += minsqe;
            if (i > 0 && part[i-1] == part[i])
                nconsecutive++;
        }
        printf("Init: iter: %2zd, totsqe: %.4le, nconsecutive: %6d/%6d\n",
                iter, totsqe/nrows, nconsecutive, nrows);
        fflush(stdout);

        /* compute new centroids/protos */
        protos_ComputeProtos(params, vault, sol);
    }

    sol->objval  = protos_ObjValue(params, vault, sol);

    return sol;

}




/*************************************************************************/
/*! Computes the protos given a partitioning */
/*************************************************************************/
void protos_ComputeProtos(params_t *params, vault_t *vault, solution_t *sol)
{
    ssize_t i, j, iP;
    int nrows, ncols, nprotos;
    ssize_t *rowptr;
    int *rowind, *part, *pcounts;
    float *rowval, **protos, *pnorms, *proto;
    float scale, pnorm;

    nprotos = params->nprotos;

    nrows  = vault->mat->nrows;
    ncols  = vault->mat->ncols;
    rowptr = vault->mat->rowptr;
    rowind = vault->mat->rowind;
    rowval = vault->mat->rowval;


    protos  = sol->protos;
    part    = sol->vpart;
    pnorms  = sol->pnorms;
    pcounts = sol->pcounts;


    /* compute new centroids/protos */
    gk_iset(nprotos, 0, pcounts);
    for (iP=0; iP<nprotos; iP++)
        gk_fset(ncols, 0.0, protos[iP]);

    for (i=0; i<nrows; i++) {
        pcounts[part[i]]++;
        for (j=rowptr[i]; j<rowptr[i+1]; j++)
            protos[part[i]][rowind[j]] += rowval[j];
    }

#ifndef NOMP
#pragma omp parallel private(iP, j, scale, proto, pnorm)
{

    #pragma omp for schedule(static)
    for (iP=0; iP<nprotos; iP++) {
        scale = (pcounts[iP] > 0.0 ? 1.0/pcounts[iP] : 0);
        proto = protos[iP];

        for (pnorm=0.0, j=0; j<ncols; j++) {
            proto[j] *= scale;
            pnorm += proto[j]*proto[j];
        }
        pnorms[iP] = pnorm;
    }

}
#else
    for (iP=0; iP<nprotos; iP++) {
        scale = (pcounts[iP] > 0.0 ? 1.0/pcounts[iP] : 0);
        pnorms[iP] = 0.0;
        for (j=0; j<ncols; j++) {
            protos[iP][j] *= scale;
            pnorms[iP] += protos[iP][j]*protos[iP][j];
        }
    }
#endif

}


/*************************************************************************/
/*! Computes the objective value of the solution and the per-sequence
    errors. */
/*************************************************************************/
double protos_ObjValue(params_t *params, vault_t *vault, solution_t *sol)
{
    ssize_t iS, iV, iP, iD;
    int nrows, nseqs;
    ssize_t *rowptr;
    int *rowind, *vpart, *sptr;
    float *rowval, *rnorms, **protos, *pnorms, *sqes;
    double totsqe, vnorms;

    nseqs  = vault->nseqs;
    sptr   = vault->sptr;

    nrows  = vault->mat->nrows;
    rowptr = vault->mat->rowptr;
    rowind = vault->mat->rowind;
    rowval = vault->mat->rowval;
    rnorms = vault->mat->rnorms;

    protos = sol->protos;
    pnorms = sol->pnorms;
    vpart  = sol->vpart;
    sqes   = sol->sqes;

    totsqe = 0.0;
    for (iS=0; iS<nseqs; iS++) {
        sqes[iS] = 0.0;
        vnorms   = 0.0;
        for (iV=sptr[iS]; iV<sptr[iS+1]; iV++) {
            iP = vpart[iV];
            vnorms += rnorms[iV];
            sqes[iS] += rnorms[iV] + pnorms[iP];
            for (iD=rowptr[iV]; iD<rowptr[iV+1]; iD++)
                sqes[iS] -= 2*rowval[iD]*protos[iP][rowind[iD]];
        }
        totsqe += sqes[iS];
        sqes[iS] = sqes[iS]/(sptr[iS+1]-sptr[iS]);
        vnorms   = vnorms/(sptr[iS+1]-sptr[iS]);
        if (params->dbglvl&4){
            printf("%4zd %3d %.4f %.4f\n", iS, sptr[iS+1]-sptr[iS], sol->sqes[iS], sol->sqes[iS]/vnorms);
            fflush(stdout);
        }
    }

    return (double) (totsqe/nrows);

}


/*************************************************************************/
/*! Prints the pairwise distances of the protos */
/*************************************************************************/
void protos_PrintDistances(params_t *params, vault_t *vault, solution_t *sol)
{
    ssize_t iP, jP, iD;
    int ndims, nprotos;
    float **protos;
    float sqe;

    nprotos = params->nprotos;
    ndims   = vault->mat->ncols;
    protos  = sol->protos;

    printf("\nProto similarity -------------\n");
    for (iP=0; iP<nprotos; iP++) {
        printf("%6.1f => ", sol->pnorms[iP]);
        for (jP=0; jP<nprotos; jP++) {
            sqe = 0.0;
            for (iD=0; iD<ndims; iD++)
                sqe += (protos[iP][iD]-protos[jP][iD])*(protos[iP][iD]-protos[jP][iD]);
            printf(" %6.1f", sqe);
        }
        printf("\n");
        fflush(stdout);
    }
    printf("End proto similarity --------------\n\n");
    fflush(stdout);
}


/*************************************************************************/
/*! Prints the protos */
/*************************************************************************/
void protos_Print(params_t *params, vault_t *vault, solution_t *sol)
{
    ssize_t iP, iC;

    printf("\nProtos -----\n");
    for (iC=0; iC<vault->mat->ncols; iC++) {
        for (iP=0; iP<params->nprotos; iP++)
            printf(" %.3f,", sol->protos[iP][iC]);
        if (vault->clabels)
            printf(" %s", vault->clabels[iC]);
        printf("\n");
        fflush(stdout);
    }
    printf("End protos -----\n\n");
    fflush(stdout);

}


/*************************************************************************/
/*! Prints key features from the protos (stdout) */
/*************************************************************************/
void protos_PrintKeyFtrs(params_t *params, vault_t *vault, solution_t *sol)
{
    ssize_t iS, iV, iP, jP, iD;
    int nvecs, ndims, nprotos, sum;
    ssize_t *rowptr;
    int *rowind, *pcounts;
    float *rowval, **protos, *pnorms, *center;
    float scale, sqe, rnrm, cnorm;
    gk_fkv_t *dims;

    GKASSERT(vault->clabels != NULL);

    nprotos = params->nprotos;

    ndims = vault->ndims;

    nvecs  = vault->mat->nrows;
    rowptr = vault->mat->rowptr;
    rowind = vault->mat->rowind;
    rowval = vault->mat->rowval;

    protos  = sol->protos;
    pnorms  = sol->pnorms;
    pcounts = sol->pcounts;

    center = gk_fsmalloc(ndims, 0.0, "center");
    dims   = gk_fkvmalloc(ndims, "dims");

    /* compute the global center and its norm */
    for (iV=0; iV<nvecs; iV++) {
        for (iD=rowptr[iV]; iD<rowptr[iV+1]; iD++)
            center[rowind[iD]] += rowval[iD];
    }
    scale = 1.0/nvecs;
    cnorm = 0.0;
    for (iD=0; iD<ndims; iD++) {
        center[iD] *= scale;
        cnorm += center[iD]*center[iD];
    }


    printf("\nProtos features -----\n");

    for (iP=0; iP<nprotos; iP++) {


        for (sqe=0.0, iD=0; iD<ndims; iD++) {
            dims[iD].val = iD;
            dims[iD].key = protos[iP][iD]*protos[iP][iD];
            sqe += (protos[iP][iD]-center[iD])*(protos[iP][iD]-center[iD]);
        }

        printf("Proto: %zd [%d %.2f %.2f]\n", iP, pcounts[iP], pnorms[iP], sqe);

        if(pcounts[iP] == 0){
            printf("\tNO ASSIGNMENTS!\n");
            continue;
        }

        gk_fkvsortd(ndims, dims);
        for (rnrm=0.0, iD=0; iD<ndims; iD++) {
            rnrm += dims[iD].key;

            printf("  %.3f %+.3f %.3f %s\n",
                    dims[iD].key/pnorms[iP],
                    (protos[iP][dims[iD].val]>center[dims[iD].val] ? +1 : -1)*(protos[iP][dims[iD].val]-center[dims[iD].val])*(protos[iP][dims[iD].val]-center[dims[iD].val])/sqe,
                    rnrm/pnorms[iP],
                    vault->clabels[dims[iD].val]);

            if (rnrm > pnorms[iP]*params->n2frac || dims[iD].key/pnorms[iP] < .005)
                break;
        }
        fflush(stdout);
    }

    printf("\nEnd protos features-----\n\n");
    fflush(stdout);

    gk_free((void **)&center, &dims, LTERM);

}


/*************************************************************************/
/*! Writes key features from the protos (file) */
/*************************************************************************/
void protos_WriteKeyFtrs(params_t *params, vault_t *vault, solution_t *sol)
{
    ssize_t iS, iV, iP, jP, iD;
    int nvecs, ndims, nprotos, sum;
    ssize_t *rowptr;
    int *rowind, *pcounts;
    float *rowval, **protos, *pnorms, *center;
    float scale, sqe, rnrm, cnorm;
    gk_fkv_t *dims;
    char filename[MAX_STRLEN];
    FILE *fpout;

    sprintf(filename, "%s.ftrs.c%.3f.s%d.p%d",
            params->datafile, params->scost, params->minspan, params->nprotos);
    fpout = gk_fopen(filename, "w", "ftrfile");


    GKASSERT(vault->clabels != NULL);

    nprotos = params->nprotos;

    ndims = vault->ndims;

    nvecs  = vault->mat->nrows;
    rowptr = vault->mat->rowptr;
    rowind = vault->mat->rowind;
    rowval = vault->mat->rowval;

    protos  = sol->protos;
    pnorms  = sol->pnorms;
    pcounts = sol->pcounts;

    center = gk_fsmalloc(ndims, 0.0, "center");
    dims   = gk_fkvmalloc(ndims, "dims");

    /* compute the global center and its norm */
    for (iV=0; iV<nvecs; iV++) {
        for (iD=rowptr[iV]; iD<rowptr[iV+1]; iD++)
            center[rowind[iD]] += rowval[iD];
    }
    scale = 1.0/nvecs;
    cnorm = 0.0;
    for (iD=0; iD<ndims; iD++) {
        center[iD] *= scale;
        cnorm += center[iD]*center[iD];
    }


    for (iP=0; iP<nprotos; iP++) {
        for (sqe=0.0, iD=0; iD<ndims; iD++) {
            dims[iD].val = iD;
            dims[iD].key = protos[iP][iD]*protos[iP][iD];
            sqe += (protos[iP][iD]-center[iD])*(protos[iP][iD]-center[iD]);
        }

        fprintf(fpout, "proto:%zd vcounts:%d norm2:%.2f sqe:%.2f ", iP, pcounts[iP], pnorms[iP], sqe);

        gk_fkvsortd(ndims, dims);
        for (rnrm=0.0, iD=0; iD<ndims; iD++) {
            rnrm += dims[iD].key;

            fprintf(fpout, " pwgt:%.3f cwgt:%.3f label:%s",
                    dims[iD].key/pnorms[iP],
                    center[dims[iD].val]*center[dims[iD].val]/cnorm,
                    vault->clabels[dims[iD].val]);

            if (rnrm > pnorms[iP]*params->n2frac || dims[iD].key/pnorms[iP] < .005)
                break;
        }
        fprintf(fpout, "\n");
    }
    gk_fclose(fpout);

    gk_free((void **)&center, &dims, LTERM);

}


/*************************************************************************/
/*! Computes the proto-to-proto transition matrix */
/*************************************************************************/
void protos_PrintP2PTMat(params_t *params, vault_t *vault, solution_t *sol)
{
    ssize_t iS, iV, iP, jP;
    int nseqs, nprotos, sum;
    int *vpart, *sptr, **p2pfreq;
    float scale;

    nprotos = params->nprotos;

    nseqs = vault->nseqs;
    sptr  = vault->sptr;
    vpart = sol->vpart;

    gk_AllocMatrix((void ***)&p2pfreq, sizeof(int), nprotos, nprotos);
    for (iP=0; iP<nprotos; iP++)
        gk_iset(nprotos, 0, p2pfreq[iP]);

    /* compute the frequency of transitions across all sequences */
    for (iS=0; iS<nseqs; iS++) {
        for (iV=sptr[iS]+1; iV<sptr[iS+1]; iV++) {
            if (vpart[iV-1] != vpart[iV])
                p2pfreq[vpart[iV-1]][vpart[iV]]++;
        }
    }

    printf("\nProto-to-Proto transition matrix -------------\n");
    for (iP=0; iP<nprotos; iP++) {
        for (sum=0, jP=0; jP<nprotos; jP++)
            sum += p2pfreq[iP][jP];

        printf("%5d => ", sum);
        for (scale=1.0/(1+sum), jP=0; jP<nprotos; jP++)
            printf(" %.3f", scale*p2pfreq[iP][jP]);
        printf("\n");
        fflush(stdout);
    }
    printf("End Proto-to-Proto transition matrix -------------\n\n");
    fflush(stdout);

    gk_FreeMatrix((void ***)&p2pfreq, nprotos, nprotos);
}


/*************************************************************************/
/*! Writes the paths taken by each sequence */
/*************************************************************************/
void protos_WritePaths(params_t *params, vault_t *vault, solution_t *sol)
{
    ssize_t iS, iV, j, n;
    int nseqs, slen;
    int *vpart, *sptr, **ptr = NULL;
    FILE *fpout;

    fpout = gk_fopen(params->ofile, "w", "pathfile");

    nseqs = vault->nseqs;
    sptr  = vault->sptr;
    vpart = sol->vpart;
    slen = nseqs > 0 ? 10 * (sptr[1] - sptr[0]) : 0;  /* initial max seq length */

    if(slen)
        ptr = gk_iAllocMatrix(2, slen, 0, "protos_WritePaths: ptr");

    /* output the paths */
    for (iS=0; iS<nseqs; iS++) {
        if(sptr[iS+1] - sptr[iS] > slen){
            gk_iFreeMatrix(&ptr, 2, slen);
            slen = 2 * sptr[iS+1] - sptr[iS];
            ptr = gk_iAllocMatrix(2, slen, 0, "protos_WritePaths: ptr");
        }
        ptr[1][0] = vpart[sptr[iS]];
        for (n=0, j=1, iV=sptr[iS]+1; iV<sptr[iS+1]; iV++, j++) {
            if (vpart[iV-1] != vpart[iV]){
                ptr[1][n+1] = vpart[iV];
                ptr[0][n++] = j;
            }
        }
        ptr[0][n++] = sptr[iS+1] - sptr[iS];

        fprintf(fpout, "%zu ", n);
        for (j=0; j<n; ++j)
            fprintf(fpout, "%d ", ptr[0][j]);
        for (j=0; j<n; ++j)
            fprintf(fpout, "%d ", ptr[1][j]);

        fprintf(fpout, "\n");
    }

    gk_fclose(fpout);
    if(ptr)
        gk_iFreeMatrix(&ptr, 2, slen);

}


/*************************************************************************/
/*! Prints information about the most frequent transitions */
/*************************************************************************/
void protos_PrintP2PTInfo(params_t *params, vault_t *vault, solution_t *sol)
{
    ssize_t iS, iV, iP, jP, iT, iD;
    int nseqs, ndims, nprotos, sum, ntrans;
    int *vpart, *sptr, **p2pfreq;
    float **protos, *pnorms;
    float scale, sqe, rsqe;
    gk_fkv_t *trans, *dims;

    GKASSERT(vault->clabels != NULL);

    nprotos = params->nprotos;

    ndims = vault->ndims;
    nseqs = vault->nseqs;
    sptr  = vault->sptr;

    protos = sol->protos;
    pnorms = sol->pnorms;
    vpart  = sol->vpart;

    gk_AllocMatrix((void ***)&p2pfreq, sizeof(int), nprotos, nprotos);
    for (iP=0; iP<nprotos; iP++)
        gk_iset(nprotos, 0, p2pfreq[iP]);

    trans = gk_fkvmalloc(nprotos*nprotos, "trans");
    dims  = gk_fkvmalloc(ndims, "dims");

    /* compute the frequency of transitions across all sequences */
    for (iS=0; iS<nseqs; iS++) {
        for (iV=sptr[iS]+1; iV<sptr[iS+1]; iV++) {
            if (vpart[iV-1] != vpart[iV])
                p2pfreq[vpart[iV-1]][vpart[iV]]++;
        }
    }

    /* identify the transitions that satisfy the minimum cut-offs */
    ntrans = 0;
    for (iP=0; iP<nprotos; iP++) {
        for (sum=0, jP=0; jP<nprotos; jP++)
            sum += p2pfreq[iP][jP];

        for (scale=1.0/(1+sum), jP=0; jP<nprotos; jP++) {
            if (scale*p2pfreq[iP][jP] > params->mintp) {
                trans[ntrans].val = iP*nprotos+jP;
                trans[ntrans].key = scale*p2pfreq[iP][jP];
                ntrans++;
            }
        }
    }
    gk_fkvsortd(ntrans, trans);

    printf("\nFrequent Proto-to-Proto transitions -------------\n");

    for (iT=0; iT<ntrans; iT++) {
        iP = trans[iT].val/nprotos;
        jP = trans[iT].val%nprotos;

        for (rsqe=0.0, sqe=0.0, iD=0; iD<ndims; iD++) {
            dims[iD].val = iD;
            dims[iD].key = (protos[iP][iD]-protos[jP][iD])*(protos[iP][iD]-protos[jP][iD]);
            sqe += dims[iD].key;
            rsqe += protos[iP][iD]*protos[jP][iD];
        }

        printf("%zd => %zd [%d %.3f][%.2f %.2f %.2f %.2f]\n",
                iP, jP, p2pfreq[iP][jP], trans[iT].key, pnorms[iP], pnorms[jP], sqe,
                pnorms[iP]+pnorms[jP]-2.0*rsqe);
        fflush(stdout);

        gk_fkvsortd(ndims, dims);
        for (rsqe=0.0, iD=0; iD<ndims; iD++) {
            rsqe += dims[iD].key;

            printf("  %+6.2f %+6.3f %.3f %.3f %.2f %s\n",
                    (protos[iP][dims[iD].val] > protos[jP][dims[iD].val] ? -1 : 1)*sqrt(dims[iD].key)/protos[iP][dims[iD].val],
                    log(protos[jP][dims[iD].val]/(.01+protos[iP][dims[iD].val])),
                    protos[iP][dims[iD].val]*protos[iP][dims[iD].val]/pnorms[iP],
                    protos[jP][dims[iD].val]*protos[jP][dims[iD].val]/pnorms[jP],
                    rsqe/sqe,
                    vault->clabels[dims[iD].val]);
            fflush(stdout);

            if (rsqe > sqe*params->n2frac)
                break;
        }
    }

    printf("\nEnd Frequent Proto-to-Proto transitions -------------\n");

    gk_FreeMatrix((void ***)&p2pfreq, nprotos, nprotos);
    gk_free((void **)&trans, &dims, LTERM);

}



