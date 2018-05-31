/*!
\file
\brief The entry point of the proto-based decoding 
\date Started 6/16/2014
\author George
 */

#include "orion.h"

/*************************************************************************
 * The entry point
 **************************************************************************/
int main(int argc, char *argv[])
{
    ssize_t i;
    float *rcmds;
    params_t *params;
    vault_t *vault;
    solution_t *sol;

    params = getcmdline_params(argc, argv);

    vault = loadData(params);

    if (params->dbglvl>0){
        printf("\n-----------------\n");
        printf("#seqs: %d, #vecs: %d, #dims: %d, #protos: %d\n",
                vault->nseqs, vault->mat->nrows, vault->ndims, params->nprotos);
        printf("rowmodel: %s, minspan: %d, niters: %d\n",
                rowmodelnames[params->rowmodel], params->minspan, params->niters);
        printf("scost: %.3f, mintp: %.3f, n2frac: %.3f\n",
                params->scost, params->mintp, params->n2frac);
        printf("\n");
    }

    sol = protos_Cluster(params, vault);

    if (params->dbglvl>0)
        printf("\n-----------------\n");

    /** free memory */
    free_memory(params, vault, sol);
}



void free_memory(params_t *params, vault_t *vault, solution_t *sol)
{
    ssize_t i, nprotos;
    nprotos = params->nprotos;

    if(sol){
        if(sol->protos){
            for(i=0; i < nprotos; ++i)
                gk_free((void**)&sol->protos[i], LTERM);
            gk_free((void**)&sol->protos, LTERM);
        }
        gk_free((void**)&sol->pcounts, &sol->pnorms, &sol->sqes, &sol->vpart, &sol, LTERM);
    }
    if(vault){
        if(vault->clabels){
            for(i=0; i < vault->ndims; ++i)
                gk_free((void**)&vault->clabels[i], LTERM);
            gk_free((void**)&vault->clabels, LTERM);
        }
        if(vault->mat)
            gk_csr_Free(&vault->mat);
        gk_free((void**)&vault->sptr, &vault, LTERM);
    }
    gk_free((void**)&params->datafile, &params->transfile, &params->p2ptfile,
            &params->pathsfile, &params->ftrsfile, &params->protosfile, &params->tinfofile,
            &params->countsfile, &params->clabelsfile,
            &params, LTERM);
}


