/*!
\file
\brief Contains the core routines for proto-based clustering of sequences
\date Started 6/16/2014
\author George, David
 */

#include "orion.h"

/** internal function prototypes **/
solution_t* diff_Refine(params_t *params, vault_t *vault, solution_t *sol);
double diff_ObjValue(params_t *params, vault_t *vault, solution_t *sol);
void diff_DecodeSequence(params_t *params, vault_t *vault, solution_t *sol, int u);
void diff_SetProtoAppCounts(params_t *params, vault_t *vault, solution_t *sol);
void diff_SetGlobalDiffusionMatrix(params_t *params, vault_t *vault, solution_t *sol);
void diff_UpdateDiffusionMatrix(params_t *params, vault_t *vault, solution_t *sol, int p);
void diff_ComputeProtos(params_t *params, vault_t *vault, solution_t *sol);


/*************************************************************************/
/*! Top-level proto-based clustering routine */
/*************************************************************************/
solution_t *cluster_diff(params_t *params, vault_t *vault)
{
    solution_t *csol;

    srand(params->seed);

    protos_PreprocessData(params, vault);

    /* need sorted indices if using sparse diffusion matrix */
    if(params->diffmode == DIFFMODE_TPW || params->diffmode == DIFFMODE_TDR)
        gk_csr_SortIndices(vault->mat, GK_CSR_ROW);

    csol = protos_FindInitial(params, vault);

    if (params->dbglvl&2)
        protos_Print(params, vault, csol);

    protos_PrintDistances(params, vault, csol);
    protos_PrintP2PTMat(params, vault, csol);

    csol = diff_Refine(params, vault, csol);

    if (params->dbglvl&2)
        protos_Print(params, vault, csol);

    protos_PrintDistances(params, vault, csol);
    protos_PrintP2PTMat(params, vault, csol);
    if (params->dbglvl&1){
        protos_PrintKeyFtrs(params, vault, csol);
        protos_PrintP2PTInfo(params, vault, csol);
    }

    if (params->writepaths)
        protos_WritePaths(params, vault, csol);

    if (params->writeftrs)
        protos_WriteKeyFtrs(params, vault, csol);

    return csol;

}


/*************************************************************************/
/*! Iteratively refines the segmentation */
/*************************************************************************/
solution_t* diff_Refine(params_t *params, vault_t *vault, solution_t *sol)
{
    ssize_t iS, iP, imp;
    int nrows, ncols, nseqs, napps, nprotos, iter;
    solution_t *bsol = NULL;  /* best solution to date */

    nrows   = vault->mat->nrows;
    ncols   = vault->mat->ncols;
    nseqs   = vault->nseqs;
    napps   = vault->mat->ncols;
    nprotos = params->nprotos;
    bsol    = NULL;

    /** allocate additional required memory */
    sol->apcnts = gk_iAllocMatrix(nprotos, napps, 0, "diff_Refine: vault->apcnts");
    sol->smat = (gk_csr_t**) gk_malloc(nprotos * sizeof(gk_csr_t*), "diff_Refine: vault->smat");
    sol->dmat = (double**) gk_malloc(nprotos * sizeof(double*), "diff_Refine: vault->dmat");
    memset(sol->smat, 0, nprotos * sizeof(gk_csr_t*));
    memset(sol->dmat, 0, nprotos * sizeof(double*));

    printf("Initial objval: %.4lf\n", sol->objval);

    /** Set initial proto app counts OR global diffusion matrix */
    if(params->diffglobal)
        diff_SetGlobalDiffusionMatrix(params, vault, sol);
    else
        diff_SetProtoAppCounts(params, vault, sol);

    for (imp=0, iter=0; iter < params->niters; iter++) {

        /* update the proto diffusion matrices X_p */
        if(!params->diffglobal){
            printf("\tupdate diffusion matrices...\n");
            if(params->diffmode == DIFFMODE_TPW || params->diffmode == DIFFMODE_TDR)
                printf("\t");
            fflush(stdout);

            if(params->diffmode != DIFFMODE_DR && params->diffmode != DIFFMODE_TDR){
#ifndef NOMP
#pragma omp parallel private(iP)
{
                #pragma omp for schedule(static)
                for(iP=0; iP < nprotos; ++iP)
                    diff_UpdateDiffusionMatrix(params, vault, sol, iP);
}
#else
                for(iP=0; iP < nprotos; ++iP)
                    diff_UpdateDiffusionMatrix(params, vault, sol, iP);
#endif
            } else {  /** svdlibc does not appear to be thread safe :-( **/
                for(iP=0; iP < nprotos; ++iP)
                    diff_UpdateDiffusionMatrix(params, vault, sol, iP);
            }
            if(params->diffmode == DIFFMODE_TPW || params->diffmode == DIFFMODE_TDR)
                printf("\n");
        }

        printf("\tdecode sequences...\n");
        fflush(stdout);
        for (iS=0; iS<nseqs; ++iS)
            diff_DecodeSequence(params, vault, sol, iS);

        printf("\tre-computing protos...\n");
        fflush(stdout);
        gk_iSetMatrix(sol->apcnts, nprotos, napps, 0); /* start new counts */
        diff_ComputeProtos(params, vault, sol);
        sol->objval = protos_ObjValue(params, vault, sol);

        if(bsol == NULL || sol->objval < bsol->objval){
            bsol = copy_output_sol(sol, bsol, nrows, ncols, nprotos);
            imp = 0;
        } else
            imp++;
        printf("  iter: %3d, objval: %.4lf, best objval: %.4lf\n", iter, sol->objval, bsol->objval);
        fflush(stdout);

        if(imp == NIMPROVE){
            printf("  No improvement in %d iterations. Stopping...\n", NIMPROVE);
            break;
        }
    }

    fflush(stdout);

    /** free additional required memory */
    if(params->diffglobal){
        if(sol->dmat)
            for(iP=1; iP<nprotos; ++iP)
                sol->dmat[iP] = NULL;
        if(sol->smat)
            for(iP=1; iP<nprotos; ++iP)
                sol->smat[iP] = NULL;
    }
    switch(params->diffmode){
        case DIFFMODE_PW:
            for(iP=0; iP<nprotos; ++iP)
                if(sol->dmat[iP])
                    free(sol->dmat[iP]);
            break;
        case DIFFMODE_TPW:
            for(iP=0; iP<nprotos; ++iP)
                if(sol->smat[iP])
                    l2ap_csr_free((l2ap_csr_t*)sol->smat[iP]);
            break;
        case DIFFMODE_DR:
            for(iP=0; iP<nprotos; ++iP)
                if(sol->dmat[iP])
                    gk_free((void**)&sol->dmat[iP], LTERM);
            break;
        case DIFFMODE_TDR:
            for(iP=0; iP<nprotos; ++iP){
                if(sol->dmat[iP])
                    gk_free((void**)&sol->dmat[iP], LTERM);
                if(sol->smat[iP])
                    gk_csr_Free(&sol->smat[iP]);
            }
            break;
        case DIFFMODE_PRW:
            gk_errexit(SIGERR, "Diff mode PRW not yet implemented.");
            break;
        default:
            gk_errexit(SIGERR, "Invalid diff mode.");
            break;
    }
    /* free solution matrix */
    gk_fFreeMatrix(&(sol->protos), nprotos, ncols);
    gk_iFreeMatrix(&(sol->apcnts), nprotos, napps);
    gk_free((void**)&sol->smat, &sol->dmat, &sol->pcounts, &sol->pnorms, &sol->sqes, &sol->vpart, &sol, LTERM);

    return bsol;

}


/*************************************************************************/
/*! Sets app proto counts                                                */
/*************************************************************************/
void diff_SetProtoAppCounts(params_t *params, vault_t *vault, solution_t *sol)
{
    ssize_t i, j;
    int nrows, ncols, nprotos;
    ssize_t *rowptr;
    int *rowind, *part, **apcnts;

    nprotos = params->nprotos;

    nrows  = vault->mat->nrows;
    ncols  = vault->mat->ncols;
    rowptr = vault->mat->rowptr;
    rowind = vault->mat->rowind;
    part    = sol->vpart;
    apcnts  = sol->apcnts;

    /* reset proto app counts */
    gk_iSetMatrix(apcnts, nprotos, ncols, 0);

    /* compute new proto app counts */
    for (i=0; i<nrows; i++)
        for (j=rowptr[i]; j<rowptr[i+1]; j++)
            apcnts[part[i]][rowind[j]]++;

}


/*************************************************************************/
/*! Computes the protos given a partitioning + updates app/proto counts  */
/*************************************************************************/
void diff_ComputeProtos(params_t *params, vault_t *vault, solution_t *sol)
{
    ssize_t i, j, iP;
    int nrows, ncols, nprotos;
    ssize_t *rowptr;
    int *rowind, *part, *pcounts, **apcnts;
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
    apcnts  = sol->apcnts;

    /* reset proto/app counts */
    gk_iSetMatrix(apcnts, nprotos, ncols, 0);

    /* compute new centroids/protos */
    gk_iset(nprotos, 0, pcounts);
    for (iP=0; iP<nprotos; iP++)
        gk_fset(ncols, 0.0, protos[iP]);

    for (i=0; i<nrows; i++) {
        pcounts[part[i]]++;
        for (j=rowptr[i]; j<rowptr[i+1]; j++){
            protos[part[i]][rowind[j]] += rowval[j];
            apcnts[part[i]][rowind[j]]++;
        }
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
/*! Updates the diffusion matrix for proto p                             */
/*************************************************************************/
void diff_UpdateDiffusionMatrix(params_t *params, vault_t *vault, solution_t *sol, int p)
{
    ssize_t iA, iV, i, j, k, nnz, nnz2;
    int nrows, ncols;
    ssize_t *rowptr, *ptr, *cnts;
    int *rowind, *ind, *part, *pcnts;
    float *rowval, *val;
    gk_csr_t *pmat, *smat;
    l2ap_output_t *psol = NULL;
    double *dmat, sum, *Ut=NULL, *S=NULL, *Vt=NULL;
    char fname[100];

    nrows   = vault->mat->nrows;
    ncols   = vault->mat->ncols;
    rowptr  = vault->mat->rowptr;
    rowind  = vault->mat->rowind;
    rowval  = vault->mat->rowval;
    part    = sol->vpart;
    pcnts   = sol->apcnts[p]; /* proto app counts for this proto */

    /* create sparse matrix containing app data for this proto */
    nnz = gk_isum(ncols, pcnts, 1);
    pmat = gk_csr_Create();
    ptr = gk_zmalloc(ncols+1, "diff_UpdateDiffusionMatrix: ptr");
    ind = gk_imalloc(nnz, "diff_UpdateDiffusionMatrix: ind");
    val = gk_fmalloc(nnz, "diff_UpdateDiffusionMatrix: val");

    for(ptr[0]=0, i=0; i < ncols; ++i)
        ptr[i+1] = ptr[i] + pcnts[i];

    for (nnz2=0, iV=0, i=0; i<nrows; ++i) {
        if(part[i] != p)
            continue;
        for (j=rowptr[i]; j<rowptr[i+1]; j++){
            iA = rowind[j]; /* app */
            ind[ptr[iA]] = iV;
            val[ptr[iA]++] = rowval[j];
            nnz2++;
        }
        iV++;
    }
    GKASSERT(nnz == nnz2);

    /** rest pointer array */
    for(ptr[0]=0, i=0; i < ncols; ++i)
        ptr[i+1] = ptr[i] + pcnts[i];

    pmat->nrows = ncols;
    pmat->ncols = iV;
    pmat->rowptr = ptr;
    pmat->rowind = ind;
    pmat->rowval = val;

    /* Normalize rows */
    gk_csr_Normalize(pmat, GK_CSR_ROW, 2);

    switch(params->diffmode){
        case DIFFMODE_PW:
            /* pairwise similarities stored in dense matrix as 1-D array */
            if(!sol->dmat[p])
                sol->dmat[p] = gk_dmalloc(ncols * ncols, "sol->dmat[p]");
            if(pmat->nrows > 0 && pmat->ncols > 0){
                psol = psim_cos_d((const l2ap_csr_t*) pmat, sol->dmat[p], 0.0, L2AP_MODE_IDXJOIN, 0);
                l2ap_dmat_set_eye(psol->dneighbors, ncols, ncols);
                dmat = psol->dneighbors;
                /* normalize rows w/ l-1 norm */
                for(i=0; i < ncols; ++i){
                    sum = 1.0/gk_dsum(ncols, dmat + (i*ncols), 1);
                    gk_dscale(ncols, sum, dmat + (i*ncols), 1);
                }
                sol->dmat[p] = dmat;
                psol->dneighbors = NULL;
                l2ap_out_free(psol);
            } else{
                memset(sol->dmat[p], 0, ncols * ncols);
                printf("\tempty input for proto %d\n", p);
            }

            break;

        case DIFFMODE_TPW:
            /* pairwise similarities stored in CSR matrix, including diagonal */
            if(sol->smat[p])
                l2ap_csr_free((l2ap_csr_t*)sol->smat[p]);
            psol = psim_cos_c((const l2ap_csr_t*) pmat, params->simt, L2AP_MODE_IDXJOIN);
            /* add upper triangular and Identity - also sorts columns */
            l2ap_csr_set_ut(psol->neighbors, 1, 0);
            /* normalize rows of the sim matrix */
            l2ap_csr_normalize_l1(psol->neighbors);
            sol->smat[p] = (gk_csr_t*) psol->neighbors;
            psol->neighbors = NULL;
            l2ap_out_free(psol);
            gk_csr_CreateIndex(sol->smat[p], GK_CSR_COL);

            printf(" p%d: %zu,", p, sol->smat[p]->rowptr[sol->smat[p]->nrows]);
            fflush(stdout);

            break;

        case DIFFMODE_DR:
        case DIFFMODE_TDR:

            if(!sol->dmat[p])
                sol->dmat[p] = gk_dmalloc(ncols * ncols, "sol->dmat[p]");
            if(pmat->nrows > 0 && pmat->ncols > 0){
                /* first execute dimensionality reduction on the input matrix */
//                l2ap_csr_print_info((l2ap_csr_t*)pmat, "\tComputing SVD on... ", "\n");
//                fflush(stdout);
                Ut = gk_dmalloc(params->ldim * pmat->nrows, "Ut");
                S  = gk_dmalloc(params->ldim, "S");
#if LSI_SHOW_ERROR
                Vt = gk_dmalloc(params->ldim * pmat->ncols, "Vt");
#endif
                compute_svd_csr(pmat, params->ldim, params->svdmaxit, Ut, S, NULL);
                /* figure out reconstruction error */
#if LSI_SHOW_ERROR
                printf("\tsvd error: %f\n", compute_svd_error_csr(pmat, params->ldim, Ut, S, Vt));
                fflush(stdout);
                gk_free((void**)&Vt, LTERM);
#endif
                sol->dmat[p] = pairwise_sims_dmat(Ut, pmat->nrows, params->ldim, S, sol->dmat[p]);
                gk_free((void**)&Ut, &S, LTERM);

                if(params->diffmode == DIFFMODE_DR){
                    /* normalize rows with l1 norm */
                    for(i=0; i < ncols; ++i){
                        sum = 1.0/gk_dsum(ncols, sol->dmat[p] + (i*ncols), 1);
                        gk_dscale(ncols, sum, sol->dmat[p] + (i*ncols), 1);
                    }
                } else {
                    if(sol->smat[p])
                        gk_csr_Free(&sol->smat[p]);
                    dmat = sol->dmat[p];
                    smat = gk_csr_Create();
                    smat->nrows = smat->ncols = ncols;
                    cnts = gk_zsmalloc(ncols+1, 0, "tptr");
                    /* count non-zeros */
                    for(i=0; i < ncols; ++i)
                        for(j=0; j < ncols; ++j)
                            if(dmat[i*ncols+j] >= params->simt)
                                cnts[j]++;
                    nnz = gk_zsum(ncols, cnts, 1);
                    smat->rowind = gk_imalloc(nnz, "smat->rowind");
                    smat->rowval = gk_fmalloc(nnz, "smat->rowval");
                    smat->rowptr = gk_zmalloc(ncols+1, "smat->rowptr");
                    /* transfer data */
                    smat->rowptr[0] = 0;
                    for(i=0; i < ncols; ++i)
                        smat->rowptr[i+1] = smat->rowptr[i] + cnts[i];
                    gk_zcopy(ncols+1, smat->rowptr, cnts);
                    for(i=0; i < ncols; ++i)
                        for(j=0; j < ncols; ++j)
                            if(dmat[i*ncols+j] >= params->simt){
                                smat->rowind[cnts[j]] = j;
                                smat->rowval[cnts[j]++] = dmat[i*ncols+j];
                            }
                    printf(" p%d: %zu,", p, smat->rowptr[smat->nrows]);
                    gk_csr_CreateIndex(smat, GK_CSR_COL);
                    sol->smat[p] = smat;
                    gk_free((void**)&cnts, LTERM);
                }

            } else {
                memset(sol->dmat[p], 0, ncols * ncols);
                printf("\tempty input for proto %d\n", p);
            }

            break;

        case DIFFMODE_PRW:
            gk_errexit(SIGERR, "Diff mode PRW not yet implemented.");
            break;
        default:
            gk_errexit(SIGERR, "Invalid diff mode.");
            break;
    }

    /** free temporary data matrix */
    gk_csr_Free(&pmat);

}


/*************************************************************************/
/*! Updates the diffusion matrix for proto p                             */
/*************************************************************************/
void diff_SetGlobalDiffusionMatrix(params_t *params, vault_t *vault, solution_t *sol)
{
    ssize_t iA, iV, i, j, k, nnz, nnz2;
    int nrows, ncols, nprotos;
    ssize_t *ptr, *cnts;
    int *ind, *pcnts;
    float *val;
    gk_csr_t *pmat, *smat;
    l2ap_output_t *psol = NULL;
    double *dmat, sum, *Ut=NULL, *S=NULL, *Vt=NULL;
    char fname[100];

    nprotos = params->nprotos;

    nrows   = vault->mat->nrows;
    ncols   = vault->mat->ncols;

    /* create sparse matrix containing app data for this proto */
    if(!vault->mat->colptr)
        gk_csr_CreateIndex(vault->mat, GK_CSR_COL);
    pmat = gk_csr_Create();
    ptr = vault->mat->colptr;
    ind = vault->mat->colind;
    val = vault->mat->colval;

    pmat->nrows = ncols;
    pmat->ncols = nrows;
    pmat->rowptr = ptr;
    pmat->rowind = ind;
    pmat->rowval = val;

    /* Normalize rows */
    gk_csr_Normalize(pmat, GK_CSR_ROW, 2);

    switch(params->diffmode){
        case DIFFMODE_PW:
            /* pairwise similarities stored in dense matrix as 1-D array */
            if(!sol->dmat[0])
                sol->dmat[0] = gk_dmalloc(ncols * ncols, "sol->dmat[p]");
            psol = psim_cos_d((const l2ap_csr_t*) pmat, sol->dmat[0], 0.0, L2AP_MODE_IDXJOIN, 0);
            l2ap_dmat_set_eye(psol->dneighbors, ncols, ncols);
            dmat = psol->dneighbors;
            /* normalize rows w/ l-1 norm */
            for(i=0; i < ncols; ++i){
                sum = 1.0/gk_dsum(ncols, dmat + (i*ncols), 1);
                gk_dscale(ncols, sum, dmat + (i*ncols), 1);
            }
            sol->dmat[0] = dmat;
            psol->dneighbors = NULL;
            l2ap_out_free(psol);

            /* all protos share matrix */
            for(i=1; i < nprotos; ++i)
                sol->dmat[i] = sol->dmat[0];

            break;

        case DIFFMODE_TPW:
            /* pairwise similarities stored in CSR matrix, including diagonal */
            if(sol->smat[0])
                l2ap_csr_free((l2ap_csr_t*)sol->smat[0]);
            psol = psim_cos_c((const l2ap_csr_t*) pmat, params->simt, L2AP_MODE_IDXJOIN);
            /* add upper triangular and Identity - also sorts columns */
            l2ap_csr_set_ut(psol->neighbors, 1, 0);
            /* normalize rows of the sim matrix */
            l2ap_csr_normalize_l1(psol->neighbors);
            sol->smat[0] = (gk_csr_t*) psol->neighbors;
            psol->neighbors = NULL;
            l2ap_out_free(psol);
            gk_csr_CreateIndex(sol->smat[0], GK_CSR_COL);

            printf(" p%d: %zu,", 0, sol->smat[0]->rowptr[sol->smat[0]->nrows]);
            fflush(stdout);

            /* all protos share matrix */
            for(i=1; i < nprotos; ++i)
                sol->smat[i] = sol->smat[0];

            break;

        case DIFFMODE_DR:
        case DIFFMODE_TDR:

            if(!sol->dmat[0])
                sol->dmat[0] = gk_dmalloc(ncols * ncols, "sol->dmat[p]");
            /* first execute dimensionality reduction on the input matrix */
            Ut = gk_dmalloc(params->ldim * pmat->nrows, "Ut");
            S  = gk_dmalloc(params->ldim, "S");
#if LSI_SHOW_ERROR
            Vt = gk_dmalloc(params->ldim * pmat->ncols, "Vt");
#endif
            compute_svd_csr(pmat, params->ldim, params->svdmaxit, Ut, S, NULL);
            /* figure out reconstruction error */
#if LSI_SHOW_ERROR
            printf("\tsvd error: %f\n", compute_svd_error_csr(pmat, params->ldim, Ut, S, Vt));
            fflush(stdout);
            gk_free((void**)&Vt, LTERM);
#endif
            sol->dmat[0] = pairwise_sims_dmat(Ut, pmat->nrows, params->ldim, S, sol->dmat[0]);
            gk_free((void**)&Ut, &S, LTERM);

            if(params->diffmode == DIFFMODE_DR){
                /* normalize rows with l1 norm */
                for(i=0; i < ncols; ++i){
                    sum = 1.0/gk_dsum(ncols, sol->dmat[0] + (i*ncols), 1);
                    gk_dscale(ncols, sum, sol->dmat[0] + (i*ncols), 1);
                }
                /* all protos share matrix */
                for(i=1; i < nprotos; ++i)
                    sol->dmat[i] = sol->dmat[0];
            } else {
                if(sol->smat[0])
                    gk_csr_Free(&sol->smat[0]);
                dmat = sol->dmat[0];
                smat = gk_csr_Create();
                smat->nrows = smat->ncols = ncols;
                cnts = gk_zsmalloc(ncols+1, 0, "tptr");
                /* count non-zeros */
                for(i=0; i < ncols; ++i)
                    for(j=0; j < ncols; ++j)
                        if(dmat[i*ncols+j] >= params->simt)
                            cnts[j]++;
                nnz = gk_zsum(ncols, cnts, 1);
                smat->rowind = gk_imalloc(nnz, "smat->rowind");
                smat->rowval = gk_fmalloc(nnz, "smat->rowval");
                smat->rowptr = gk_zmalloc(ncols+1, "smat->rowptr");
                /* transfer data */
                smat->rowptr[0] = 0;
                for(i=0; i < ncols; ++i)
                    smat->rowptr[i+1] = smat->rowptr[i] + cnts[i];
                gk_zcopy(ncols+1, smat->rowptr, cnts);
                for(i=0; i < ncols; ++i)
                    for(j=0; j < ncols; ++j)
                        if(dmat[i*ncols+j] >= params->simt){
                            smat->rowind[cnts[j]] = j;
                            smat->rowval[cnts[j]++] = dmat[i*ncols+j];
                        }
                printf(" p%d: %zu,", 0, smat->rowptr[smat->nrows]);
                gk_csr_CreateIndex(smat, GK_CSR_COL);
                sol->smat[0] = smat;
                gk_free((void**)&cnts, LTERM);
                /* all protos share matrix */
                for(i=1; i < nprotos; ++i)
                    sol->smat[i] = sol->smat[0];
            }

            break;

        case DIFFMODE_PRW:
            gk_errexit(SIGERR, "Diff mode PRW not yet implemented.");
            break;
        default:
            gk_errexit(SIGERR, "Invalid diff mode.");
            break;
    }

    /** free temporary data matrix */
    pmat->rowptr = NULL;
    pmat->rowind = NULL;
    pmat->rowval = NULL;
    gk_csr_Free(&pmat);

}


/*************************************************************************/
/*! Computes the objective value of the solution and the per-sequence
    errors. */
/*************************************************************************/
double diff_ObjValue(params_t *params, vault_t *vault, solution_t *sol)
{
    ssize_t iS, iV, iP, iD, j, l, k;
    int nrows, ncols, nseqs;
    ssize_t *rowptr, *ptr;
    int *rowind, *ind, *vpart, *sptr;
    float *rowval, *val, **protos, *proto, *pnorms, *sqes;
    double totsqe, vnorms, **dmat, *dmatP, dtmp, dp, vnorm, v;
    gk_csr_t **smat;

    nseqs  = vault->nseqs;
    sptr   = vault->sptr;

    nrows  = vault->mat->nrows;
    ncols  = vault->mat->ncols;
    rowptr = vault->mat->rowptr;
    rowind = vault->mat->rowind;
    rowval = vault->mat->rowval;
    protos = sol->protos;
    pnorms = sol->pnorms;
    vpart  = sol->vpart;
    sqes   = sol->sqes;
    dmat   = sol->dmat;             /* dense proto diffusion matrices */
    smat   = sol->smat;             /* sparse proto diffusion matrices */

    totsqe = 0.0;
//    for (iS=0; iS<nseqs; iS++) {
//        sqes[iS] = 0.0;
//        vnorms   = 0.0;
//        for (iV=sptr[iS]; iV<sptr[iS+1]; iV++) {
//            iP = vpart[iV];
//            vnorms += rnorms[iV];
//            sqes[iS] += rnorms[iV] + pnorms[iP];
//            for (iD=rowptr[iV]; iD<rowptr[iV+1]; iD++)
//                sqes[iS] -= 2*rowval[iD]*protos[iP][rowind[iD]];
//        }
//        totsqe += sqes[iS];
//        sqes[iS] = sqes[iS]/(sptr[iS+1]-sptr[iS]);
//        vnorms   = vnorms/(sptr[iS+1]-sptr[iS]);
//        if (params->dbglvl&4){
//            printf("%4zd %3d %.4f %.4f\n", iS, sptr[iS+1]-sptr[iS], sol->sqes[iS], sol->sqes[iS]/vnorms);
//            fflush(stdout);
//        }
//    }

    for (iS=0; iS<nseqs; iS++) {
        sqes[iS] = 0.0;
        vnorms   = 0.0;
        for (iV=sptr[iS]; iV<sptr[iS+1]; iV++) {
            iP = vpart[iV];
            /* compute norm of D_p*v */
            vnorm = 0.0;
            dtmp = 0.0;
            if(smat[iP]){
                assert(smat[iP]->colptr);
                for (proto=protos[iP], j=rowptr[iV]; j<rowptr[iV+1]; j++){
                    ptr = smat[iP]->colptr;
                    ind = smat[iP]->colind;
                    val = smat[iP]->colval;
                    for (l=ptr[rowind[j]]; l<ptr[rowind[j]+1]; ++l){
                        dtmp += val[l]*rowval[j] * proto[ind[l]];
                        vnorm += val[l]*val[l] * rowval[j]*rowval[j];
                    }
                }
            } else if(dmat[iP]){
                dmatP = dmat[iP];
                for (proto=protos[iP], j=rowptr[iV]; j<rowptr[iV+1]; j++){
                    for (k=0; k<ncols; ++k){
                        v = dmatP[k*ncols + rowind[j]];
                        if(v > 0.0){
                            dtmp += v*rowval[j] * proto[k];
                            vnorm += v*v * rowval[j]*rowval[j];
                        }
                    }
                }
            } else
                gk_errexit(SIGERR, "Missing proto diffusion matrix for proto %d.", iP);

            vnorms += vnorm;
            sqes[iS] += vnorm + pnorms[iP] - 2*dtmp;

//            vnorms += rnorms[iV];
//            sqes[iS] += rnorms[iV] + pnorms[iP];
//            for (iD=rowptr[iV]; iD<rowptr[iV+1]; iD++)
//                sqes[iS] -= 2*rowval[iD]*protos[iP][rowind[iD]];
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
/*! Computes an optimal proto-based decoding of a multi-variate sequence */
/*************************************************************************/
void diff_DecodeSequence(params_t *params, vault_t *vault, solution_t *sol, int u)
{
    ssize_t i, j, k, l, iV, jV, iP, iS;
    int nvecs, ncols, nprotos, minspan, minp, minj, minseg;
    ssize_t *rowptr, *ptr;
    int *rowind, *ind, *sptr, *vpart, *pos, *ids, *nsegs;
    float *rowval, *val, **protos, *proto, *pnorms, **dotps, **vnorms, *scores;
    float minsqe, sqe, totsqe, oldcost, newcost;
    double **dmat, *dmatP, dtmp, dp, vnorm, v;
    gk_csr_t **smat;
    char fname[100];

    nprotos = params->nprotos;
    minspan = params->minspan;

    sptr = vault->sptr;

    ncols  = vault->mat->ncols;
    rowptr = vault->mat->rowptr;
    rowind = vault->mat->rowind;
    rowval = vault->mat->rowval;

    protos = sol->protos;
    pnorms = sol->pnorms;
    vpart  = sol->vpart + sptr[u];  /* to simplify accessing vpart[] */
    dmat   = sol->dmat;             /* dense proto diffusion matrices */
    smat   = sol->smat;             /* sparse proto diffusion matrices */

    nvecs = sptr[u+1]-sptr[u];

    GKASSERT(nvecs >= 2*minspan);

    /* allocate memory for the solutions */
    scores = gk_fsmalloc(nvecs+1, 0.0, "scores");
    pos    = gk_ismalloc(nvecs+1, -1, "pos");
    ids    = gk_ismalloc(nvecs+1, -1, "ids");
    nsegs  = gk_ismalloc(nvecs+1, 0, "nsegs");
    scores++;
    pos++;
    ids++;
    nsegs++;


    /* compute the dot-products between weeks & protos and the norms of the diffused vectors */
    gk_AllocMatrix((void ***)&dotps, sizeof(float), nvecs, nprotos);
    gk_AllocMatrix((void ***)&vnorms, sizeof(float), nvecs, nprotos);


    for (iV=0; iV<nvecs; iV++) {
        i = sptr[u]+iV;
#ifndef NOMP
#pragma omp parallel private(iP, j, k, l, v, proto, dtmp, vnorm, ptr, ind, val, dmatP)
{
        #pragma omp for schedule(static)
        for (iP=0; iP<nprotos; iP++) {
            dtmp=0.0;
            vnorm=0.0;
            if(smat[iP]){
                assert(smat[iP]->colptr);
                for (proto=protos[iP], j=rowptr[i]; j<rowptr[i+1]; j++){
                    ptr = smat[iP]->colptr;
                    ind = smat[iP]->colind;
                    val = smat[iP]->colval;
                    for (l=ptr[rowind[j]]; l<ptr[rowind[j]+1]; ++l){
                        dtmp += val[l]*rowval[j] * proto[ind[l]];
                        vnorm += val[l]*val[l] * rowval[j]*rowval[j];
                    }
                }
            } else if(dmat[iP]){
                dmatP = dmat[iP];
                for (proto=protos[iP], j=rowptr[i]; j<rowptr[i+1]; j++){
                    for (k=0; k<ncols; ++k){
                        v = dmatP[k*ncols + rowind[j]];
                        if(v > 0.0){
                            dtmp += v*rowval[j] * proto[k];
                            vnorm += v*v * rowval[j]*rowval[j];
                        }
                    }
                }
            } else
                gk_errexit(SIGERR, "Missing proto diffusion matrix for proto %d.", iP);
            vnorms[iV][iP] = vnorm;
            dotps[iV][iP] = dtmp;
        }
}
#else
        for (iP=0; iP<nprotos; iP++) {
            dtmp=0.0;
            vnorm=0.0;
            if(smat[iP]){
                assert(smat[iP]->colptr);
                for (proto=protos[iP], j=rowptr[i]; j<rowptr[i+1]; j++){
                    ptr = smat[iP]->colptr;
                    ind = smat[iP]->colind;
                    val = smat[iP]->colval;
                    for (l=ptr[rowind[j]]; l<ptr[rowind[j]+1]; ++l){
                        dtmp += val[l]*rowval[j] * proto[ind[l]];
                        vnorm += val[l]*val[l] * rowval[j]*rowval[j];
                    }
                }
            } else if(dmat[iP]){
                dmatP = dmat[iP];
                for (proto=protos[iP], j=rowptr[i]; j<rowptr[i+1]; j++){
                    for (k=0; k<ncols; ++k){
                        v = dmatP[k*ncols + rowind[j]];
                        if(v > 0.0){
                            dtmp += v*rowval[j] * proto[k];
                            vnorm += v*v * rowval[j]*rowval[j];
                        }
                    }
                }
            } else
                gk_errexit(SIGERR, "Missing proto diffusion matrix for proto %d.", iP);
            vnorms[iV][iP] = vnorm;
            dotps[iV][iP] = dtmp;
        }
#endif
    }

    /* setup the initial conditions
     the scores[minspan-1 ... 2*minspan-2] are set based on a single best scoring proto */
    for (iP=0; iP<nprotos; iP++) {
        for (sqe=0.0, iV=0; iV<minspan; iV++)
            sqe += vnorms[iV][iP] + pnorms[iP] - 2.0*dotps[iV][iP];

        if (ids[iV-1] == -1 || scores[iV-1] > sqe) {
            scores[iV-1] = sqe;
            ids[iV-1]    = iP;
            pos[iV-1]    = 0;
            nsegs[iV-1]  = 1;
        }

        for (iV=minspan; iV<2*minspan-1; iV++) {
            sqe += vnorms[iV][iP] + pnorms[iP] - 2.0*dotps[iV][iP];

            if (ids[iV] == -1 || scores[iV] > sqe) {
                scores[iV] = sqe;
                ids[iV]    = iP;
                pos[iV]    = 0;
                nsegs[iV]  = 1;
            }
        }
    }


    /* now apply the DP algorithm to fill in the rest */
    for (iV=2*minspan-1; iV<nvecs; iV++) {
        minsqe = 0.0;
        minj   = -1;
        minp   = -1;
        for (iP=0; iP<nprotos; iP++) {
            /* find best [jV...iV] span for iP */
            for (sqe=0.0, jV=iV; jV>=0; jV--) {
                sqe += vnorms[iV][iP] + pnorms[iP] - 2.0*dotps[jV][iP];

                if (iV-jV < minspan-1)
                    continue; /* skip short segments */

                if (jV > 0 && jV < minspan)
                    continue; /* skip any positions within the first minspan excluding the start */

                if (minj == -1) {
                    minsqe = sqe + scores[jV-1];
                    minj   = jV;
                    minp   = iP;
                    minseg = nsegs[minj-1] + 1;
                    if (params->dbglvl&1)
                        printf("minj==-1 :: iV=%3zd, minj=jV=%3d, minp=iP=%3d, sqe=%.3f, minsqe=%.3f, minseg=%d\n",
                                iV, minj, minp, sqe, minsqe, minseg);
                }
                else {
                    oldcost = (1.0+minseg*params->scost)*minsqe;
                    newcost = (1.0+(nsegs[jV-1]+1)*params->scost)*(sqe + scores[jV-1]);
                    if (1.00001*oldcost > newcost) {
                        minsqe = sqe + scores[jV-1];
                        minseg = nsegs[jV-1] + 1;
                        minj   = jV;
                        minp   = iP;
                    }
                }
                if (params->dbglvl&1){
                    printf("jV=%3zd, iP=%3zd, sqe=%.3f, minsqe=%.3f, minj=%d, minp=%d\n",
                            jV, iP, sqe, minsqe, minj, minp);
                    fflush(stdout);
                }
            }
        }

        scores[iV] = minsqe;
        ids[iV]    = minp;
        pos[iV]    = minj;
        nsegs[iV]  = nsegs[pos[iV]-1] + 1;

        if (params->dbglvl&1){
            printf("iV=%3zd, scores[iV]=%.3f, ids[iV]=%d, pos[iV]=%d, nsegs[iV]=%d, ids[iV-1]=%d, scores[iV-1]=%.3f\n",
                    iV, scores[iV], ids[iV], pos[iV], nsegs[iV], ids[iV-1], scores[iV-1]);
            fflush(stdout);
        }
    }
    //printf("u:%d, iV=%3zd, scores[iV]=%.1f, ids[iV]=%d\n", u, iV-1, scores[iV-1], ids[iV-1]);

    for (iV=nvecs-1; iV>=0; ) {
        for (jV=iV; jV>=pos[iV]; jV--) {
            if (params->dbglvl&1){
                printf("jV=%3zd, old-vpart[jV]=%2d new-vpart[jV]=%2d\n", jV, vpart[jV], ids[iV]);
                fflush(stdout);
            }

            vpart[jV] = ids[iV];
        }
        iV = pos[iV]-1;
    }

    /* adjust the +1's */
    scores--;
    pos--;
    ids--;
    nsegs--;

    gk_FreeMatrix((void ***)&dotps, nvecs, nprotos);
    gk_FreeMatrix((void ***)&vnorms, nvecs, nprotos);
    gk_free((void **)&scores, &ids, &pos, &nsegs, LTERM);

}


