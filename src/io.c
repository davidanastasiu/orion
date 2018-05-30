/*!
\file
\brief This file contains routines for reading in the data set
\date Started 6/16/2014
\author George, David
 */

#include "orion.h"


/**************************************************************************/
/*! Reads the input data */
/**************************************************************************/
vault_t *loadData(params_t *params)
{
    ssize_t i, j;
    size_t nlines;
    vault_t *vault;
    char filename[MAX_STRLEN];
    int32_t *counts;

    vault = (vault_t *)gk_malloc(sizeof(vault_t), "loadData: vault");
    memset(vault, 0, sizeof(vault_t));

    /* read the matrix */
    printf("Reading matrix %s...\n", params->datafile);
    GKASSERT(gk_fexists(params->datafile));
    vault->mat = gk_csr_Read(params->datafile, GK_CSR_FMT_CSR, 1, 1);

    /* read the counts */
    printf("Reading the week counts of the users...\n");
    sprintf(filename, "%s.counts", params->datafile);
    GKASSERT(gk_fexists(filename));
    counts = gk_i32readfile(filename, &nlines);

    /* setup vault->sptr from the counts */
    vault->sptr = gk_i32malloc(nlines+1, "sptr");
    gk_i32copy(nlines, counts, vault->sptr);
    MAKECSR(i, (ssize_t) nlines, vault->sptr);
    gk_free((void **)&counts, LTERM);

    /* set dataset size info and do a sanity check */
    vault->nseqs = nlines;
    vault->ndims = vault->mat->ncols;
    ASSERT(vault->mat->nrows == vault->sptr[vault->nseqs]);


    /* read the clabels */
    printf("Reading the feature labels...\n");
    sprintf(filename, "%s.clabels", params->datafile);
    if (gk_fexists(filename)){
        vault->clabels = gk_readfile(filename, &nlines);
        GKASSERT((ssize_t)nlines == vault->ndims);
    } else {
        vault->clabels = (char**) gk_malloc(vault->ndims * sizeof(char*), "vault->clabels");
        for(i=0; i < vault->ndims; ++i){
            sprintf(filename, "l%zu", i+1);
            vault->clabels[i] = gk_strdup(filename);
        }
    }


    return vault;

}

