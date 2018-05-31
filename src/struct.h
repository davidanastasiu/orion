/*!
\file
\brief Data structures used in the program
\date Started 6/16/2014
\author George
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

} solution_t; 


/*************************************************************************
 * run parameters
 **************************************************************************/
typedef struct {
    int nprotos;          /* The number of protos */
    int ncliters;         /* The number of clustering iterations */
    int niters;           /* The number of sequence refinement iterations */
    int simtype;          /* The similarity type */
    int rowmodel;         /* The row model */
    int minspan;          /* The minimum number of weeks to be covered by a proto */
    int translevel;       /* The transitions level to output probabilities for */
    float mintp;          /* The minimum P2P transition probability */
    float n2frac;         /* The maximum fraction of a vector's L2^2 for analysis */
    float scost;          /* A per-segment cost */

    int writepaths;       /* Outputs the paths for each sequence */
    int writeftrs;        /* Outputs the key features of each proto */
    int writetrans;       /* Outputs the transition probabilities for a given transition level */
    int writetinfo;       /* Outputs frequent transitions and key discriminative features for those transitions */
    int writep2pt;        /* Outputs the global proto-to-proto transition probabilities matrix */
    int writeprotos;      /* Outputs the proto vectors */

    int seed;             /* Seed for the random number generator */
    int dbglvl;           /* Debugging information */

    char *datafile;       /* The file storing the CSR format of the input data */
    char *countsfile;     /* The file storing the row counts for each sequence in the datafile */
    char *clabelsfile;    /* The file storing the feature labels */
    char *pathsfile;      /* The output file to store the proto transition paths for each of the sequences */
    char *ftrsfile;       /* The output file to store the feature analysis data */
    char *transfile;      /* The output file to store transition probabilities for a given trasition level */
    char *tinfofile;      /* The output file to store frequent transitions and key discriminative features for those transitions */
    char *p2ptfile;       /* The output file to store the global proto-to-proto transition probabilities */
    char *protosfile;     /* The output file to store the final protos */
} params_t;



#endif 
