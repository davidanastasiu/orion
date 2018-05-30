/*!
\file
\brief Contains the core routines for proto-based clustering of sequences
\date Started 6/16/2014
\author George, David
 */

#include "orion.h"


/*************************************************************************/
/*! Top-level proto-based clustering routine */
/*************************************************************************/
solution_t *cluster_otm(params_t *params, vault_t *vault)
{
    solution_t *csol;

    srand(params->seed);

    gk_errexit(SIGERR, "OTM mode not yet implemented.");

    protos_PreprocessData(params, vault);

    csol = protos_FindInitial(params, vault);

    if (params->dbglvl&2)
        protos_Print(params, vault, csol);


    protos_PrintDistances(params, vault, csol);
    protos_PrintP2PTMat(params, vault, csol);

//    otm_Refine(params, vault, csol);

    if (params->dbglvl&2)
        protos_Print(params, vault, csol);

    protos_PrintDistances(params, vault, csol);
    protos_PrintP2PTMat(params, vault, csol);
    protos_PrintKeyFtrs(params, vault, csol);
    protos_PrintP2PTInfo(params, vault, csol);

    if (params->writepaths)
        protos_WritePaths(params, vault, csol);

    if (params->writeftrs)
        protos_WriteKeyFtrs(params, vault, csol);

    return csol;

}
