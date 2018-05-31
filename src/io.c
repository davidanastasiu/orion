/*!
\file
\brief This file contains routines for reading in the data set
\date Started 6/16/2014
\author George
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
  if (params->dbglvl>0)
      printf("Reading matrix %s...\n", params->datafile);
  GKASSERTP((params->datafile && gk_fexists(params->datafile)),
          ("Missing <fstem> file %s.", params->datafile));
  vault->mat = gk_csr_Read(params->datafile, da_getFileFormat(params->datafile, 0), 1, 1);

  /* read the counts */
  if (params->dbglvl>0)
      printf("Reading user sequence counts...\n");
  if (params->countsfile && gk_fexists(params->countsfile)){
      counts = gk_i32readfile(params->countsfile, &nlines);
  } else {
      sprintf(filename, "%s.counts", params->datafile);
      GKASSERTP(gk_fexists(filename),
              ("Missing counts file %s.counts.", params->datafile));
      counts = gk_i32readfile(filename, &nlines);
  }

  /* setup vault->sptr from the counts */
  vault->sptr = gk_i32malloc(nlines+1, "sptr");
  gk_i32copy(nlines, counts, vault->sptr);
  MAKECSR(i, nlines, vault->sptr);
  gk_free((void **)&counts, LTERM);

  /* set dataset size info and do a sanity check */
  vault->nseqs = nlines;
  vault->ndims = vault->mat->ncols;
  ASSERT(vault->mat->nrows == vault->sptr[vault->nseqs]);


  /* read the clabels */
  sprintf(filename, "%s.clabels", params->datafile);
  if (params->clabelsfile && gk_fexists(params->clabelsfile)){
      vault->clabels = gk_readfile(params->clabelsfile, &nlines);
      GKASSERTP((ssize_t)nlines == vault->ndims,
              ("The number of lines in the feature file (%zu) does not match the number of features (%d).",
                      nlines, vault->ndims));
  } else if (gk_fexists(filename)){
      if (params->dbglvl>0)
          printf("Reading feature labels...\n");
      vault->clabels = gk_readfile(filename, &nlines);
      GKASSERTP((ssize_t)nlines == vault->ndims,
              ("The number of lines in the feature file (%zu) does not match the number of features (%d).",
                      nlines, vault->ndims));
  } else {
      if (params->dbglvl>0)
          printf("Generating feature labels...\n");
      vault->clabels = (char**) gk_malloc(vault->ndims * sizeof(char*), "vault->clabels");
      for(i=0; i < vault->ndims; ++i){
          sprintf(filename, "l%zu", i+1);
          vault->clabels[i] = gk_strdup(filename);
      }
  }
  return vault;

}



/*************************************************************************
* This function returns the key of a particular StringMap ID
**************************************************************************/
char* da_getStringKey(gk_StringMap_t *strmap, const char id)
{
  int i;

  for (i=0; strmap[i].name; i++) {
    if (strmap[i].id == id)
      return strmap[i].name;
  }

  return NULL;
}


/*************************************************************************
* This function returns the ID of a particular string based on the
* supplied StringMap array
**************************************************************************/
int da_getStringID(gk_StringMap_t *strmap, char *key)
{
  int i;

  for (i=0; strmap[i].name; i++) {
    if (gk_strcasecmp(key, strmap[i].name))
      return strmap[i].id;
  }

  return -1;
}


/*************************************************************************/
/*! If format not specifically given (> 0), check if a text (non-binary) text file
 *  containing a csr is in CLUTO or CSR format.
    \param file is the matrix file to be checked.
    \return the CSR format: DA_FMT_CLUTO or DA_FMT_CSR
 */
/*************************************************************************/
char da_getFileFormat(char *file, const char format)
{
    if(format > 0) return format;
    size_t nnz;
    char fmt;
    char *ext, *p;

    ext = strrchr(file, '.');
    if(ext){
        ext++;
        //make lowercase
        for (p=ext ; *p; ++p) *p = tolower(*p);
        if ((fmt = da_getStringID(fmt_options, ext)) > -1)
            return fmt;
    } else if(gk_fexists(file)){ // assume some sort of CSR. Can we guess?
        gk_getfilestats(file, NULL, &nnz, NULL, NULL);
        return (nnz%2 == 1) ? GK_CSR_FMT_CLUTO : GK_CSR_FMT_CSR;
    }
    return -1;
}
