/*!
\file
\brief Data structures used in the program
\date Started 6/16/2014
\author George, David
 */

#ifndef _STRUCT_ORION_H_
#define _STRUCT_ORION_H_




/*************************************************************************
 * the data vault
 **************************************************************************/
typedef struct {
    int nseqs;          /* # of sequences */
    int ndims;          /* # of features */

    int32_t *sptr;      /* the ptr structure for identifying the sequences in mat */
    gk_csr_t *mat;      /* the csr structure storing the vectors of the sequences */

    char **clabels;     /* the labels of the columns */

} vault_t; 


/*************************************************************************
 * the solution
 **************************************************************************/
typedef struct {
    float **protos;     /* the current set of protos */
    float *pnorms;      /* the squared norms of the protos */
    int *pcounts;       /* the occupancy counts for each proto */
    int *vpart;         /* the decoding of the vectors to the protos */
    float *sqes;        /* the errors associated with each sequence decoding */
    double objval;      /* the objective value of the decoding */

    int **apcnts;       /* value counts for each app/proto combination */
    gk_csr_t **smat;    /* csr structures storing sparse diffusion matrices for each proto */
    double **dmat;     /* LT-diag portions of pairwise sim matrices storing dense diffusion matrices */
                        /* a proto will either have a dense or a sparse diffusion matrix, depending on load */

} solution_t; 


/*************************************************************************
 * run parameters
 **************************************************************************/
typedef struct {
    int mode;             /* The algorithm to apply - orion, dsim, otp */
    int diffmode;         /* The type diffusion to apply - how to compute the pairwise similarity matrix in DSIM mode */
    int nprotos;          /* The number of protos */
    int niters;           /* The number of clustering iterations */
    int simtype;          /* The similarity type */
    int rowmodel;         /* The row model */
    int minspan;          /* The minimum number of weeks to be covered by a proto */
    int nthreads;         /* number of threads */
    int diffglobal;       /* whether to use a gloabl diffusion matrix */
    int ldim;             /* Number of latent dimensions for dimensionality reduction */
    int svdmaxit;         /* maximum number of SVD iterations */

    float mintp;          /* The minimum P2P transition probability */
    float n2frac;         /* The maximum fraction of a vector's L2^2 for analysis */
    float scost;          /* A per-segment cost */
    float simt;           /* Minimum similarity for thresholded pairwise similarity computation in DIFF mode */

    int writepaths;       /* Outputs the paths for each sequence */
    int writeftrs;        /* Outputs the key features of each proto */

    int seed;             /* Seed for the random number generator */
    int dbglvl;           /* Debugging information */

    char *datafile;       /* The file storing the CSR format of the input data */
    char *dataset;        /* The datafile - extension */
    char *ofile;          /* The file storing the output segmentation and chosen protos */
} params_t;



#endif 
