/*!
\file
\brief The entry point of the proto-based decoding 
\date Started 6/16/2014
\author George, David
 */

#include "orion.h"

void free_memory(params_t *params, vault_t *vault, solution_t *sol);

/*************************************************************************
 * The entry point
 **************************************************************************/
int main(int argc, char *argv[])
{
    ssize_t i;
    params_t *params = NULL;
    vault_t *vault = NULL;
    solution_t *sol = NULL;

    params = getcmdline_params(argc, argv);

    vault = loadData(params);

    printf("\n-----------------\n");
    printf("mode: %s, diffmode: %s, diffglobal: %d, nthreads: %d\n",
            modenames[params->mode], params->mode == MODE_DIFF ? diffusionnames[params->diffmode] : "N/A",
            params->diffglobal, params->nthreads);
    printf("#seqs: %d, #vecs: %d, #dims: %d, #protos: %d\n",
            vault->nseqs, vault->mat->nrows, vault->ndims, params->nprotos);
    printf("simtype: %s, rowmodel: %s, minspan: %d, niters: %d\n",
            simtypenames[params->simtype], rowmodelnames[params->rowmodel],
            params->minspan, params->niters);
    printf("scost: %.3f, mintp: %.3f, n2frac: %.3f, simt: %.3f\n",
            params->scost, params->mintp, params->n2frac, params->mode == MODE_DIFF ? params->simt : -1);
    printf("\n");

    switch(params->mode){
        case MODE_ORION:
            sol = cluster_orion(params, vault);
            break;
        case MODE_DIFF:
            sol = cluster_diff(params, vault);
            break;
        case MODE_OTM:
            sol = cluster_otm(params, vault);
            break;
        default:
            gk_errexit(SIGERR, "Invalid mode selected %d.", params->mode);
            break;
    }

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
        if(sol->apcnts){
            for(i=0; i < nprotos; ++i)
                gk_free((void**)&sol->apcnts[i], LTERM);
            gk_free((void**)&sol->apcnts, LTERM);
        }
        if(sol->dmat){
            for(i=0; i < nprotos; ++i)
                gk_free((void**)&sol->dmat[i], LTERM);
            gk_free((void**)&sol->dmat, LTERM);
        }
        if(sol->smat){
            for(i=0; i < nprotos; ++i)
                gk_csr_Free((gk_csr_t**)&sol->smat[i]);
            gk_free((void**)&sol->smat, LTERM);
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
    gk_free((void**)&params->datafile, &params->dataset, &params->ofile, &params, LTERM);
}

