/*!
\file
\brief Contains the core routines for proto-based clustering of sequences
\date Started 6/16/2014
\author George
 */

#include "orion.h"


/*************************************************************************/
/*! Top-level proto-based clustering routine */
/*************************************************************************/
solution_t *protos_Cluster(params_t *params, vault_t *vault)
{
    solution_t *csol;

    srand(params->seed);

    protos_PreprocessData(params, vault);

    csol = protos_FindInitial(params, vault);

    if (params->dbglvl&2)
        protos_Print(params, vault, csol);


    if (params->dbglvl>0){
        protos_PrintDistances(params, vault, csol);
        protos_PrintP2PTMat(params, vault, csol);
    }

    protos_Refine(params, vault, csol);

    if (params->dbglvl&2)
        protos_Print(params, vault, csol);


    if (params->dbglvl>0){
        printf("\n\nAnalysis results:\n\n");
        protos_PrintDistances(params, vault, csol);
        protos_PrintP2PTMat(params, vault, csol);
        protos_PrintKeyFtrs(params, vault, csol);
        protos_PrintP2PTInfo(params, vault, csol);
    }

    if (params->writepaths)
        protos_WritePaths(params, vault, csol);

    if (params->writeftrs)
        protos_WriteKeyFtrs(params, vault, csol);

    if (params->writetrans)
        protos_WriteTransitions(params, vault, csol);

    if (params->writep2pt)
        protos_WriteP2PTMat(params, vault, csol);

    if (params->writeprotos)
        protos_Write(params, vault, csol);

    if (params->writetinfo)
        protos_WriteP2PTInfo(params, vault, csol);

    return csol;

}


/*************************************************************************/
/*! Pre-process the input data */
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
/*! Find the initial set of protos */
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
    minsqe  = FLT_MAX;

    nrows  = vault->mat->nrows;
    ncols  = vault->mat->ncols;
    rowptr = vault->mat->rowptr;
    rowind = vault->mat->rowind;
    rowval = vault->mat->rowval;
    rnorms = vault->mat->rnorms;

    /* allocate memory for the solution */
    sol = (solution_t *)gk_malloc(sizeof(solution_t), "protos_FindInitial: sol");
    memset(sol, 0, sizeof(solution_t));

    gk_AllocMatrix((void ***)&(sol->protos), sizeof(float), nprotos, ncols);
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

    if (params->dbglvl>0)
        printf("Initial clustering:\n");

    /* perform a small number of k-means iterations */
    for (iter=0; iter<params->ncliters; iter++) {
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
        if (params->dbglvl>0)
            printf("  iter: %2zd, totsqe: %.4le, nconsecutive: %6d/%6d\n",
                iter, totsqe/nrows, nconsecutive, nrows);

        /* compute new centroids/protos */
        protos_ComputeProtos(params, vault, sol);
    }

    sol->objval  = protos_ObjValue(params, vault, sol);

    return sol;

}


/*************************************************************************/
/*! Iteratively refine the segmentation */
/*************************************************************************/
void protos_Refine(params_t *params, vault_t *vault, solution_t *sol)
{
    ssize_t iS;
    int nseqs, iter;

    nseqs = vault->nseqs;
    if (params->dbglvl>0)
        printf("Initial objval: %.4lf\n", sol->objval);

    for (iter=0; iter<params->niters; iter++) {
        for (iS=0; iS<nseqs; iS++)
            protos_DecodeSequence(params, vault, sol, iS);

        protos_ComputeProtos(params, vault, sol);
        sol->objval = protos_ObjValue(params, vault, sol);
        if (params->dbglvl>0)
            printf("  iter: %3d, objval: %.4lf\n", iter, sol->objval);
    }
}


/*************************************************************************/
/*! Compute an optimal proto-based decoding of a multi-variate sequence */
/*************************************************************************/
void protos_DecodeSequence(params_t *params, vault_t *vault, solution_t *sol, 
        int u)
{
    ssize_t i, j, iV, jV, iP;
    int nvecs, ncols, nprotos, minspan, minp, minj, minseg;
    ssize_t *rowptr;
    int *rowind, *sptr, *vpart, *pos, *ids, *nsegs;
    float *rowval, *rnorms, **protos, *proto, *pnorms, **dotps, *scores;
    float minsqe, sqe, totsqe, ftmp, oldcost, newcost;


    nprotos = params->nprotos;
    minspan = params->minspan;

    sptr = vault->sptr;

    ncols  = vault->mat->ncols;
    rowptr = vault->mat->rowptr;
    rowind = vault->mat->rowind;
    rowval = vault->mat->rowval;
    rnorms = vault->mat->rnorms + sptr[u];  /* to simplify accessing rnorms[] */

    protos = sol->protos;
    pnorms = sol->pnorms;
    vpart  = sol->vpart + sptr[u];  /* to simplify accessing vpart[] */

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


    /* compute the dot-products between weeks & protos */
    gk_AllocMatrix((void ***)&dotps, sizeof(float), nvecs, nprotos);
    for (iV=0; iV<nvecs; iV++) {
        i = sptr[u]+iV;
        for (iP=0; iP<nprotos; iP++) {
            for (ftmp=0.0, proto=protos[iP], j=rowptr[i]; j<rowptr[i+1]; j++)
                ftmp += rowval[j]*proto[rowind[j]];
            dotps[iV][iP] = ftmp;
        }
    }

    /* setup the initial conditions
     the scores[minspan-1 ... 2*minspan-2] are set based on a single best scoring proto */
    for (iP=0; iP<nprotos; iP++) {
        for (sqe=0.0, iV=0; iV<minspan; iV++)
            sqe += rnorms[iV] + pnorms[iP] - 2.0*dotps[iV][iP];

        if (ids[iV-1] == -1 || scores[iV-1] > sqe) {
            scores[iV-1] = sqe;
            ids[iV-1]    = iP;
            pos[iV-1]    = 0;
            nsegs[iV-1]  = 1;
        }

        for (iV=minspan; iV<2*minspan-1; iV++) {
            sqe += rnorms[iV] + pnorms[iP] - 2.0*dotps[iV][iP];

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
                sqe += rnorms[jV] + pnorms[iP] - 2.0*dotps[jV][iP];

                if (iV-jV < minspan-1)
                    continue; /* skip short segments */

                if (jV > 0 && jV < minspan)
                    continue; /* skip any positions within the first minspan excluding the start */

                if (minj == -1) {
                    minsqe = sqe + scores[jV-1];
                    minj   = jV;
                    minp   = iP;
                    minseg = nsegs[minj-1] + 1;
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
                if (params->dbglvl&8)
                    printf("jV=%3zd, iP=%3zd, sqe=%.3f, minsqe=%.3f, minj=%d, minp=%d\n",
                            jV, iP, sqe, minsqe, minj, minp);
            }
        }

        scores[iV] = minsqe;
        ids[iV]    = minp;
        pos[iV]    = minj;
        nsegs[iV]  = nsegs[pos[iV]-1] + 1;

        if (params->dbglvl&8)
            printf("iV=%3zd, scores[iV]=%.3f, ids[iV]=%d, pos[iV]=%d, nsegs[iV]=%d, ids[iV-1]=%d, scores[iV-1]=%.3f\n",
                    iV, scores[iV], ids[iV], pos[iV], nsegs[iV], ids[iV-1], scores[iV-1]);
    }
    //printf("u:%d, iV=%3zd, scores[iV]=%.1f, ids[iV]=%d\n", u, iV-1, scores[iV-1], ids[iV-1]);

    for (iV=nvecs-1; iV>=0; ) {
        for (jV=iV; jV>=pos[iV]; jV--) {
            if (params->dbglvl&8)
                printf("jV=%3zd, old-vpart[jV]=%2d new-vpart[jV]=%2d\n", jV, vpart[jV], ids[iV]);

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
    gk_free((void **)&scores, &ids, &pos, &nsegs, LTERM);

}


/*************************************************************************/
/*! Compute the protos given a partitioning */
/*************************************************************************/
void protos_ComputeProtos(params_t *params, vault_t *vault, solution_t *sol)
{
    ssize_t i, j, iP;
    int nrows, ncols, nprotos;
    ssize_t *rowptr;
    int *rowind, *part, *pcounts;
    float *rowval, *rnorms, **protos, *pnorms;
    float scale;

    nprotos = params->nprotos;

    nrows  = vault->mat->nrows;
    ncols  = vault->mat->ncols;
    rowptr = vault->mat->rowptr;
    rowind = vault->mat->rowind;
    rowval = vault->mat->rowval;
    rnorms = vault->mat->rnorms;


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

    for (iP=0; iP<nprotos; iP++) {
        scale = (pcounts[iP] > 0.0 ? 1.0/pcounts[iP] : 0);
        pnorms[iP] = 0.0;
        for (j=0; j<ncols; j++) {
            protos[iP][j] *= scale;
            pnorms[iP] += protos[iP][j]*protos[iP][j];
        }
    }

}


/*************************************************************************/
/*! Compute the objective value of the solution and the per-sequence errors. */
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
        if (params->dbglvl&4)
            printf("%4zd %3d %.4f %.4f\n", iS, sptr[iS+1]-sptr[iS], sol->sqes[iS], sol->sqes[iS]/vnorms);
    }

    return (double) (totsqe/nrows);

}


/*************************************************************************/
/*! Print the pairwise distances of the protos */
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

    printf("\nProto distances -------------\n");
    for (iP=0; iP<nprotos; iP++) {
        printf("%6.1f => ", sol->pnorms[iP]);
        for (jP=0; jP<nprotos; jP++) {
            sqe = 0.0;
            for (iD=0; iD<ndims; iD++)
                sqe += (protos[iP][iD]-protos[jP][iD])*(protos[iP][iD]-protos[jP][iD]);
            printf(" %6.1f", sqe);
        }
        printf("\n");
    }
    printf("End proto distances --------------\n\n");

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
    }
    printf("End protos -----\n\n");

}


/*************************************************************************/
/*! Print the protos */
/*************************************************************************/
void protos_Write(params_t *params, vault_t *vault, solution_t *sol)
{
    ssize_t iP, iC;
    char filename[MAX_STRLEN];
    FILE *fpout;

    if(params->protosfile){
        fpout = gk_fopen(params->protosfile, "w", "protosfile");
    } else {
        sprintf(filename, "%s.protos.c%.3f.s%d.p%d",
                params->datafile, params->scost, params->minspan, params->nprotos);
        fpout = gk_fopen(filename, "w", "protosfile");
    }
    for (iC=0; iC<vault->mat->ncols; iC++) {
        for (iP=0; iP<params->nprotos; iP++)
            fprintf(fpout, " %.3f", sol->protos[iP][iC]);
        if (vault->clabels)
            fprintf(fpout, " %s", vault->clabels[iC]);
        fprintf(fpout, "\n");
    }

    gk_fclose(fpout);

}


/*************************************************************************/
/*! Print key features from the protos (stdout) */
/*************************************************************************/
void protos_PrintKeyFtrs(params_t *params, vault_t *vault, solution_t *sol)
{
    ssize_t iS, iV, iP, jP, iD;
    int nseqs, nvecs, ndims, nprotos, sum;
    ssize_t *rowptr;
    int *rowind, *vpart, *sptr, *pcounts;
    float *rowval, **protos, *pnorms, *center;
    float scale, sqe, rnrm, cnorm;
    gk_fkv_t *dims;

    GKASSERT(vault->clabels != NULL);

    nprotos = params->nprotos;

    ndims = vault->ndims;
    nseqs = vault->nseqs;
    sptr  = vault->sptr;

    nvecs  = vault->mat->nrows;
    rowptr = vault->mat->rowptr;
    rowind = vault->mat->rowind;
    rowval = vault->mat->rowval;

    protos  = sol->protos;
    pnorms  = sol->pnorms;
    pcounts = sol->pcounts;
    vpart   = sol->vpart;

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


    printf("\nProto features -----\n");

    for (iP=0; iP<nprotos; iP++) {
        for (sqe=0.0, iD=0; iD<ndims; iD++) {
            dims[iD].val = iD;
            dims[iD].key = protos[iP][iD]*protos[iP][iD];
            sqe += (protos[iP][iD]-center[iD])*(protos[iP][iD]-center[iD]);
        }

        printf("Proto: %zd [%d %.2f %.2f]\n", iP, pcounts[iP], pnorms[iP], sqe);

        if(pcounts[iP] == 0)
            continue;

        gk_fkvsortd(ndims, dims);
        for (rnrm=0.0, iD=0; iD<ndims; iD++) {
            rnrm += dims[iD].key;

            printf("  %.3f %.3f %+.3f %.3f %s\n",
                    dims[iD].key/pnorms[iP],
                    center[dims[iD].val]*center[dims[iD].val]/cnorm,
                    sqe > 0 ?
                        (protos[iP][dims[iD].val]>center[dims[iD].val] ? +1 : -1) *
                        (protos[iP][dims[iD].val]-center[dims[iD].val]) *
                        (protos[iP][dims[iD].val]-center[dims[iD].val]) / sqe
                        : 0,
                    rnrm/pnorms[iP],
                    vault->clabels[dims[iD].val]);

            if (rnrm > pnorms[iP]*params->n2frac || dims[iD].key/pnorms[iP] < .005)
                break;
        }
    }

    printf("End proto features-----\n\n");

    gk_free((void **)&center, &dims, LTERM);

}


/*************************************************************************/
/*! Write key features from the protos (file) */
/*************************************************************************/
void protos_WriteKeyFtrs(params_t *params, vault_t *vault, solution_t *sol)
{
    ssize_t iS, iV, iP, jP, iD;
    int nseqs, nvecs, ndims, nprotos, sum;
    ssize_t *rowptr;
    int *rowind, *vpart, *sptr, *pcounts;
    float *rowval, **protos, *pnorms, *center;
    float scale, sqe, rnrm, cnorm;
    gk_fkv_t *dims;
    char filename[MAX_STRLEN];
    FILE *fpout;

    if(params->ftrsfile){
        fpout = gk_fopen(params->ftrsfile, "w", "ftrfile");
    } else {
        sprintf(filename, "%s.ftrs.c%.3f.s%d.p%d",
                params->datafile, params->scost, params->minspan, params->nprotos);
        fpout = gk_fopen(filename, "w", "ftrfile");
    }

    GKASSERT(vault->clabels != NULL);

    nprotos = params->nprotos;

    ndims = vault->ndims;
    nseqs = vault->nseqs;
    sptr  = vault->sptr;

    nvecs  = vault->mat->nrows;
    rowptr = vault->mat->rowptr;
    rowind = vault->mat->rowind;
    rowval = vault->mat->rowval;

    protos  = sol->protos;
    pnorms  = sol->pnorms;
    pcounts = sol->pcounts;
    vpart   = sol->vpart;

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

        fprintf(fpout, "%zd %d %.2f %.2f\t", iP, pcounts[iP], pnorms[iP], sqe);

        if(pcounts[iP] == 0){
            fprintf(fpout, "\n");
            continue;
        }

        gk_fkvsortd(ndims, dims);
        for (rnrm=0.0, iD=0; iD<ndims; iD++) {
            rnrm += dims[iD].key;

            fprintf(fpout, "%.3f %.3f %+.3f %.3f %s ",
                dims[iD].key/pnorms[iP],
                center[dims[iD].val]*center[dims[iD].val]/cnorm,
                sqe > 0 ?
                    (protos[iP][dims[iD].val]>center[dims[iD].val] ? +1 : -1) *
                    (protos[iP][dims[iD].val]-center[dims[iD].val]) *
                    (protos[iP][dims[iD].val]-center[dims[iD].val]) / sqe
                    : 0,
                rnrm/pnorms[iP],
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
/*! Compute the proto-to-proto transition matrix */
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

    printf("\nProto-to-proto transition matrix -------------\n");
    for (iP=0; iP<nprotos; iP++) {
        for (sum=0, jP=0; jP<nprotos; jP++)
            sum += p2pfreq[iP][jP];

        printf("%5d => ", sum);
        for (scale=1.0/(1+sum), jP=0; jP<nprotos; jP++)
            printf(" %.3f", scale*p2pfreq[iP][jP]);
        printf("\n");
    }
    printf("End proto-to-proto transition matrix -------------\n\n");

    gk_FreeMatrix((void ***)&p2pfreq, nprotos, nprotos);
}


/*************************************************************************/
/*! Write the global proto-to-proto transition matrix */
/*************************************************************************/
void protos_WriteP2PTMat(params_t *params, vault_t *vault, solution_t *sol)
{
    ssize_t iS, iV, iP, jP;
    int nseqs, nprotos, sum;
    int *vpart, *sptr, **p2pfreq;
    float scale;
    char filename[MAX_STRLEN];
    FILE *fpout;

    if(params->p2ptfile){
        fpout = gk_fopen(params->p2ptfile, "w", "p2ptfile");
    } else {
        sprintf(filename, "%s.p2pt.c%.3f.s%d.p%d",
                params->datafile, params->scost, params->minspan, params->nprotos);
        fpout = gk_fopen(filename, "w", "p2ptfile");
    }

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

    for (iP=0; iP<nprotos; iP++) {
        for (sum=0, jP=0; jP<nprotos; jP++)
            sum += p2pfreq[iP][jP];

        for (scale=1.0/(1+sum), jP=0; jP<nprotos; jP++)
            fprintf(fpout, " %.3f", scale*p2pfreq[iP][jP]);
        fprintf(fpout, "\n");
    }

    gk_FreeMatrix((void ***)&p2pfreq, nprotos, nprotos);

    gk_fclose(fpout);
}


/*************************************************************************/
/*! Write proto-to-proto transition matrix for a certain level */
/*************************************************************************/
void protos_WriteTransitions(params_t *params, vault_t *vault, solution_t *sol)
{
    ssize_t iS, iV, iP, jP, l;
    int nseqs, nprotos, sum;
    int *vpart, *sptr, **p2pfreq, *sfreq;
    char filename[MAX_STRLEN];
    float scale, nfseqs;
    FILE *fpout;

    if(params->transfile){
        fpout = gk_fopen(params->transfile, "w", "transfile");
    } else {
        sprintf(filename, "%s.trans.c%.3f.s%d.p%d.l%d",
                params->datafile, params->scost, params->minspan, params->nprotos, params->translevel);
        fpout = gk_fopen(filename, "w", "transfile");
    }

    nprotos = params->nprotos;

    nseqs  = vault->nseqs;
    sptr   = vault->sptr;
    vpart  = sol->vpart;

    sfreq = gk_ismalloc(nprotos, 0, "start frequencies");
    gk_AllocMatrix((void ***)&p2pfreq, sizeof(int), nprotos, nprotos);
    for (iP=0; iP<nprotos; iP++)
        gk_iset(nprotos, 0, p2pfreq[iP]);

    /* compute the frequency of transitions in level <translevel> across all sequences */
    for (iS=0; iS<nseqs; iS++) {
        for (l=0, iV=sptr[iS]+1; iV<sptr[iS+1]; iV++) {
            if(l == params->translevel-1){
                sfreq[vpart[iV-1]]++;
                l++;
            }
            if (vpart[iV-1] != vpart[iV]){
                if(l == params->translevel){
                    p2pfreq[vpart[iV-1]][vpart[iV]]++;
                    break;
                }
                l++;
            }
        }
    }

    nfseqs = (float) gk_isum(nprotos, sfreq, 1);
    if(nfseqs == 0.0)
        gk_errexit(SIGERR, "\nWriteTransitions error: -translevel set too high. No transitions found on this level.");

    for (iP=0; iP<nprotos; iP++) {
        for (sum=0, jP=0; jP<nprotos; jP++)
            sum += p2pfreq[iP][jP];

        /* write out the start state probability for iP */
        fprintf(fpout, "S %zu %.3f\n", iP, sfreq[iP]/nfseqs);

        /* write proto-to-proto transition probabilities */
        for (scale=1.0/(float)sfreq[iP], jP=0; jP<nprotos; jP++)
            if(p2pfreq[iP][jP] > 0)
                fprintf(fpout, "%zu %zu %.3f\n", iP, jP, scale*p2pfreq[iP][jP]);

        /* write out the end state probability for iP */
        fprintf(fpout, "%zu E %.3f\n", iP, sfreq[iP] > 0 ? (sfreq[iP]-sum)/(float)sfreq[iP] : 0);

    }

    gk_FreeMatrix((void ***)&p2pfreq, nprotos, nprotos);
    gk_free((void**)&sfreq, LTERM);

    gk_fclose(fpout);
}



/*************************************************************************/
/*! Write the paths taken by each sequence */
/*************************************************************************/
void protos_WritePaths(params_t *params, vault_t *vault, solution_t *sol)
{
    ssize_t iS, iV;
    int nseqs, nprotos;
    int *vpart, *sptr;
    char filename[MAX_STRLEN];
    FILE *fpout;

    if(params->pathsfile){
        fpout = gk_fopen(params->pathsfile, "w", "pathfile");
    } else {
        sprintf(filename, "%s.paths.c%.3f.s%d.p%d",
                params->datafile, params->scost, params->minspan, params->nprotos);
        fpout = gk_fopen(filename, "w", "pathfile");
    }

    nprotos = params->nprotos;

    nseqs = vault->nseqs;
    sptr  = vault->sptr;
    vpart = sol->vpart;

    /* output the paths */
    for (iS=0; iS<nseqs; iS++) {
        fprintf(fpout, "%d ", vpart[sptr[iS]]);
        for (iV=sptr[iS]+1; iV<sptr[iS+1]; iV++) {
            if (vpart[iV-1] != vpart[iV])
                fprintf(fpout, "%d ", vpart[iV]);
        }
        fprintf(fpout, "\n");
    }

    gk_fclose(fpout);

}


/*************************************************************************/
/*! Print information about the most frequent transitions */
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

    printf("\nFrequent proto-to-proto transitions -------------\n");

    for (iT=0; iT<ntrans; iT++) {
        iP = trans[iT].val/nprotos;
        jP = trans[iT].val%nprotos;

        for (rsqe=0.0, sqe=0.0, iD=0; iD<ndims; iD++) {
            dims[iD].val = iD;
            dims[iD].key = (protos[iP][iD]-protos[jP][iD])*(protos[iP][iD]-protos[jP][iD]);
            sqe += dims[iD].key;
        }

        printf("%zd => %zd [%d %.3f][%.2f %.2f %.2f]\n",
                iP, jP, p2pfreq[iP][jP], trans[iT].key, pnorms[iP], pnorms[jP], sqe);

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

            if (rsqe > sqe*params->n2frac)
                break;
        }
    }

    printf("End frequent proto-to-proto transitions -------------\n");

    gk_FreeMatrix((void ***)&p2pfreq, nprotos, nprotos);
    gk_free((void **)&trans, &dims, LTERM);

}





/*************************************************************************/
/*! Print information about the most frequent transitions */
/*************************************************************************/
void protos_WriteP2PTInfo(params_t *params, vault_t *vault, solution_t *sol)
{
    ssize_t iS, iV, iP, jP, iT, iD;
    int nseqs, ndims, nprotos, sum, ntrans;
    int *vpart, *sptr, **p2pfreq;
    float **protos, *pnorms;
    float scale, sqe, rsqe;
    gk_fkv_t *trans, *dims;
    char filename[MAX_STRLEN];
    FILE *fpout;

    if(params->tinfofile){
        fpout = gk_fopen(params->tinfofile, "w", "tinfofile");
    } else {
        sprintf(filename, "%s.tinfo.c%.3f.s%d.p%d",
                params->datafile, params->scost, params->minspan, params->nprotos);
        fpout = gk_fopen(filename, "w", "tinfofile");
    }

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

    for (iT=0; iT<ntrans; iT++) {
        iP = trans[iT].val/nprotos;
        jP = trans[iT].val%nprotos;

        for (rsqe=0.0, sqe=0.0, iD=0; iD<ndims; iD++) {
            dims[iD].val = iD;
            dims[iD].key = (protos[iP][iD]-protos[jP][iD])*(protos[iP][iD]-protos[jP][iD]);
            sqe += dims[iD].key;
        }

        fprintf(fpout, "%zd %zd %d %.3f %.2f %.2f %.2f\t",
                iP, jP, p2pfreq[iP][jP], trans[iT].key, pnorms[iP], pnorms[jP], sqe);

        gk_fkvsortd(ndims, dims);
        for (rsqe=0.0, iD=0; iD<ndims; iD++) {
            rsqe += dims[iD].key;

            fprintf(fpout, "%+6.2f %+6.3f %.3f %.3f %.2f %s ",
                    (protos[iP][dims[iD].val] > protos[jP][dims[iD].val] ? -1 : 1)*sqrt(dims[iD].key)/protos[iP][dims[iD].val],
                    log(protos[jP][dims[iD].val]/(.01+protos[iP][dims[iD].val])),
                    protos[iP][dims[iD].val]*protos[iP][dims[iD].val]/pnorms[iP],
                    protos[jP][dims[iD].val]*protos[jP][dims[iD].val]/pnorms[jP],
                    rsqe/sqe,
                    vault->clabels[dims[iD].val]);

            if (rsqe > sqe*params->n2frac)
                break;
        }
        fprintf(fpout, "\n");
    }

    gk_FreeMatrix((void ***)&p2pfreq, nprotos, nprotos);
    gk_free((void **)&trans, &dims, LTERM);

    gk_fclose(fpout);
}

