/*!
\file
\brief Contains the core routines for proto-based clustering of sequences
\date Started 6/16/2014
\author George, David
 */

#include "orion.h"

/** internal function prototypes **/
void orion_Refine(params_t *params, vault_t *vault, solution_t *sol);
void orion_DecodeSequence(params_t *params, vault_t *vault, solution_t *sol,
         int u);


/*************************************************************************/
/*! Top-level proto-based clustering routine */
/*************************************************************************/
solution_t *cluster_orion(params_t *params, vault_t *vault)
{
    solution_t *csol;

    protos_PreprocessData(params, vault);

    csol = protos_FindInitial(params, vault);

    if (params->dbglvl&2)
        protos_Print(params, vault, csol);


    protos_PrintDistances(params, vault, csol);
    protos_PrintP2PTMat(params, vault, csol);

    //protos_DecodeSequence(params, vault, csol, 79);
    //exit(0);

    orion_Refine(params, vault, csol);

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
void orion_Refine(params_t *params, vault_t *vault, solution_t *sol)
{
    ssize_t iS;
    int nseqs, iter, imp;
    double objval;

    nseqs = vault->nseqs;

    printf("Initial objval: %.4lf\n", sol->objval);
    objval = sol->objval;

    for (imp=0, iter=0; iter<params->niters; iter++) {
        for (iS=0; iS<nseqs; iS++)
            orion_DecodeSequence(params, vault, sol, iS);

        protos_ComputeProtos(params, vault, sol);
        sol->objval = protos_ObjValue(params, vault, sol);

        printf("  iter: %3d, objval: %.4lf\n", iter, sol->objval);
        fflush(stdout);
        if(objval < sol->objval || gk_abs(sol->objval - objval) < 10e-3)
            imp++;
        else
            imp = 0;
        if(imp == NIMPROVE){
            printf("  No improvement in %d iterations. Stopping...\n", NIMPROVE);
            break;
        }

    }

    fflush(stdout);
}


/*************************************************************************/
/*! Computes an optimal proto-based decoding of a multi-variate sequence */
/*************************************************************************/
void orion_DecodeSequence(params_t *params, vault_t *vault, solution_t *sol, 
        int u)
{
    ssize_t i, j, iV, jV, iP;
    int nvecs, nprotos, minspan, minp, minj, minseg;
    ssize_t *rowptr;
    int *rowind, *sptr, *vpart, *pos, *ids, *nsegs;
    float *rowval, *rnorms, **protos, *proto, *pnorms, **dotps, *scores;
    float minsqe, sqe, totsqe, ftmp, oldcost, newcost;


    nprotos = params->nprotos;
    minspan = params->minspan;

    sptr = vault->sptr;

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
    gk_free((void **)&scores, &ids, &pos, &nsegs, LTERM);

}


