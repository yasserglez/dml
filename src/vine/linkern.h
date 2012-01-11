/****************************************************************************/
/*                                                                          */
/*  (c) Copyright 1995--1999 by David Applegate, Robert Bixby,              */
/*  Vasek Chvatal, and William Cook                                         */
/*                                                                          */
/*  Permission is granted for academic research use.  For other uses,       */
/*  contact the authors for licensing options.                              */
/*                                                                          */
/*  Use at your own risk.  We make no guarantees about the                  */
/*  correctness or usefulness of this code.                                 */
/*                                                                          */
/****************************************************************************/

#ifndef LINKERN_H_
#define LINKERN_H_

typedef struct CCrandstate {
    int a;
    int b;
    int arr[55];
} CCrandstate;

typedef struct CCdata_user {
    double *x;
    double *y;
} CCdata_user;

typedef struct CCdata_rhvector {
    int dist_00;
    int dist_01;
    int dist_02;
    int dist_12;
    int dist_22;
    double p;
    int rhlength;
    char *space;
    char **vectors;
} CCdata_rhvector;

typedef struct CCdatagroup {
    int (*edgelen)(int i, int j, struct CCdatagroup *dat);
    double *x;
    double *y;
    double *z;
    int **adj;
    int *adjspace;
    int **len;
    int *lenspace;
    int *degree;
    int norm;
    int dsjrand_param;
    int default_len;
    int sparse_ecount;
    double gridsize;
    double dsjrand_factor;
    CCdata_rhvector rhdat;
    CCdata_user userdat;
    int ndepot;
    int orig_ncount;
    int *depotcost;
    int *orig_names;
} CCdatagroup;

void CCutil_sprand(int seed, CCrandstate *r);
int CCutil_graph2dat_matrix(int ncount, int ecount, int *elist, int *elen, int defaultlen, CCdatagroup *dat);
void CCutil_freedatagroup(CCdatagroup *dat);
int CClinkern_tour(int ncount, CCdatagroup *dat, int ecount, int *elist, int stallcount, int repeatcount, int *incycle, int *outcycle, double *val, CCrandstate *rstate);

#endif /* LINKERN_H_ */
