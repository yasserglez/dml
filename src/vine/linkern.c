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

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "linkern.h"

#define CC_PRANDMAX 1000000007
#define MAXDEPTH 25
#define KICK_MAXDEPTH 50
#define IMPROVE_SWITCH -1
#define LATE_DEPTH 10
#define MARK_LEVEL 10
#define BACKTRACK 4
#define MAX_BACK 12
#define BIGINT 2000000000
#define GROUPSIZE_FACTOR 0.50
#define SEGMENT_SPLIT_CUTOFF 0.30

#define SAME_SEGMENT(a, b)                                                   \
     (a->parent == b->parent &&                                              \
      ((!((F->reversed)^(a->parent->rev)) && a->id <= b->id) ||              \
       (((F->reversed)^(a->parent->rev)) && a->id >= b->id)))

#define CC_PTRWORLD_LISTFREE_ROUTINE(type, ptr_listfree_r, ptr_free_r)       \
                                                                             \
static void ptr_listfree_r (CCptrworld *world, type *p)                      \
{                                                                            \
    type *next;                                                              \
                                                                             \
    while (p != (type *) NULL) {                                             \
        next = p->next;                                                      \
        ptr_free_r (world, p);                                               \
        p = next;                                                            \
    }                                                                        \
}

#define CC_PTRWORLD_ALLOC_ROUTINE(type, ptr_alloc_r, ptr_bulkalloc_r)        \
                                                                             \
static int ptr_bulkalloc_r (CCptrworld *world, int nalloc)                   \
{                                                                            \
    CCbigchunkptr *bp;                                                       \
    int i;                                                                   \
    int count = CC_BIGCHUNK / sizeof ( type );                               \
    type *p;                                                                 \
                                                                             \
    while (nalloc > 0) {                                                     \
        bp = CCutil_bigchunkalloc ();                                        \
        if (bp == (CCbigchunkptr *) NULL) {                                  \
            fprintf (stderr, "ptr alloc failed\n");                          \
            return 1;                                                        \
        }                                                                    \
        bp->next = world->chunklist ;                                        \
        world->chunklist = bp;                                               \
                                                                             \
        p = ( type * ) bp->this_one;                                         \
        for (i=count-2; i>=0; i--) {                                         \
            p[i].next = &p[i+1];                                             \
        }                                                                    \
        p[count - 1].next = (type *) world->freelist;                        \
        world->freelist = (void *) p;                                        \
        nalloc -= count;                                                     \
    }                                                                        \
    return 0;                                                                \
}                                                                            \
                                                                             \
static type *ptr_alloc_r (CCptrworld *world)                                 \
{                                                                            \
    type *p;                                                                 \
                                                                             \
    if (world->freelist == (void *) NULL) {                                  \
        if (ptr_bulkalloc_r (world, 1)) {                                    \
            fprintf (stderr, "ptr alloc failed\n");                          \
            return ( type * ) NULL;                                          \
        }                                                                    \
    }                                                                        \
    p = (type *) world->freelist ;                                           \
    world->freelist = (void *) p->next;                                      \
                                                                             \
    return p;                                                                \
}

#define CC_PTRWORLD_FREE_ROUTINE(type, ptr_free_r)                           \
                                                                             \
static void ptr_free_r (CCptrworld *world, type *p)                          \
{                                                                            \
    p->next = (type *) world->freelist ;                                     \
    world->freelist = (void *) p;                                            \
}

#define CC_PTRWORLD_ROUTINES(type, ptr_alloc_r, ptr_bulkalloc_r, ptr_free_r) \
CC_PTRWORLD_ALLOC_ROUTINE (type, ptr_alloc_r, ptr_bulkalloc_r)               \
CC_PTRWORLD_FREE_ROUTINE (type, ptr_free_r)

#define CC_SWAP(a,b,t) (((t)=(a)),((a)=(b)),((b)=(t)))

#define CC_BIGCHUNK ((int) ((1<<16) - sizeof (CCbigchunkptr) - 16))

#define CC_SAFE_MALLOC(nnum,type)                                            \
    (type *) CCutil_allocrus (((size_t) (nnum)) * sizeof (type))

#define CC_FREE(object,type) {                                               \
    CCutil_freerus ((void *) (object));                                      \
    object = (type *) NULL;                                                  \
}

#define CC_IFFREE(object,type) {                                             \
    if ((object)) CC_FREE ((object),type);                                   \
}

#define Edgelen(n1, n2, D)  dist (n1, n2, D)

#define FLIP(aprev, a, b, bnext, f, x) {                                     \
    CClinkern_flipper_flip ((x),(a), (b));                                   \
    (f)->stack[(f)->counter].first = (a);                                    \
    (f)->stack[(f)->counter++].last = (b);                                   \
}

#define UNFLIP(aprev, a, b, bnext, f, x) {                                   \
    CClinkern_flipper_flip ((x), (b), (a));                                  \
    (f)->counter--;                                                          \
}

#define MARK(xn, xQ, xF, xD, xG, xW)  turn ((xn), (xQ), (xW))

#define markedge_add(n1, n2, E)    E->add_edges[n1 ^ n2] = 1
#define markedge_del(n1, n2, E)    E->del_edges[n1 ^ n2] = 1
#define unmarkedge_add(n1, n2, E)  E->add_edges[n1 ^ n2] = 0
#define unmarkedge_del(n1, n2, E)  E->del_edges[n1 ^ n2] = 0
#define is_it_added(n1, n2, E)     E->add_edges[n1 ^ n2]
#define is_it_deleted(n1, n2, E)   E->del_edges[n1 ^ n2]

struct CCbigchunk;

typedef struct CCbigchunkptr {
    void *this_one;
    struct CCbigchunk *this_chunk;
    struct CCbigchunkptr *next;
} CCbigchunkptr;

typedef struct CCptrworld {
    int refcount;
    void *freelist;
    CCbigchunkptr *chunklist;
} CCptrworld;

typedef struct CCbigchunk {
    char space[CC_BIGCHUNK];
    CCbigchunkptr ptr;
} CCbigchunk;

typedef struct CCdheap {
    double *key;
    int *entry;
    int *loc;
    int total_space;
    int size;
} CCdheap;

typedef struct edge {
    int other;
    int weight;
} edge;

typedef struct edgelook {
    struct edgelook *next;
    int other;
    int diff;
    int over;
    int seq;
    int side;
    int mm;
} edgelook;

typedef struct intptr {
    int this;
    struct intptr *next;
} intptr;

typedef struct flippair {
    int firstprev;
    int first;
    int last;
    int lastnext;
} flippair;

typedef struct flipstack {
    flippair *stack;
    int counter;
    int max;
} flipstack;

typedef struct graph {
    edge **goodlist;
    edge *edgespace;
    int *degree;
    int *weirdmark;
    int weirdmagic;
    int ncount;
    CCrandstate *rstate;
} graph;

typedef struct distobj {
    CCdatagroup *dat;
    int *cacheval;
    int *cacheind;
    int cacheM;
} distobj;

typedef struct adddel {
    char *add_edges;
    char *del_edges;
} adddel;

typedef struct aqueue {
    char *active;
    intptr *active_queue;
    intptr *bottom_active_queue;
    CCdheap *h;
} aqueue;

typedef struct CClk_parentnode {
    struct CClk_parentnode *adj[2];
    struct CClk_childnode *ends[2];
    int size;
    int id;
    int rev;
} CClk_parentnode;

typedef struct CClk_childnode {
    struct CClk_parentnode *parent;
    struct CClk_childnode *adj[2];
    int id;
    int name;
} CClk_childnode;

typedef struct CClk_flipper {
    CClk_parentnode *parents;
    CClk_childnode *children;
    int reversed;
    int nsegments;
    int groupsize;
    int split_cutoff;
} CClk_flipper;

typedef struct CCkdnode {
    double cutval;
    struct CCkdnode *loson;
    struct CCkdnode *hison;
    struct CCkdnode *father;
    struct CCkdnode *next;
    struct CCkdbnds *bnds;
    int lopt;
    int hipt;
    char bucket;
    char empty;
    char cutdim;
} CCkdnode;

typedef struct CCkdtree {
    CCkdnode *root;
    CCkdnode **bucketptr;
    int *perm;
    CCptrworld kdnode_world;
    CCptrworld kdbnds_world;
} CCkdtree;

static const int backtrack_count[BACKTRACK] = { 4, 3, 3, 2 };
static const int weird_backtrack_count[3] = { 4, 3, 3 };

static int CCutil_lprand(CCrandstate *r);
static void free_rhdata(CCdata_rhvector *rhdat);
static void free_userdat(CCdata_user *userdat);
static int matrix_edgelen(int i, int j, CCdatagroup *dat);
static int edgelen_nonorm(int i, int j, CCdatagroup *dat);
static void init_userdat(CCdata_user *userdat);
static void init_rhdata(CCdata_rhvector *rhdat);
static void CCutil_init_datagroup(CCdatagroup *dat);
static void *CCutil_allocrus(size_t size);
static void CCutil_freerus(void *p);
static CCbigchunkptr *CCutil_bigchunkalloc(void);
static void CCutil_bigchunkfree(CCbigchunkptr *bp);
static void CCptrworld_init(CCptrworld *world);
static void CCptrworld_delete(CCptrworld *world);
static int repeated_lin_kernighan(graph *G, distobj *D, int *cyc,
        int stallcount, int count, double *val, CCptrworld *intptr_world,
        CCptrworld *edgelook_world, CCrandstate *rstate);
static void lin_kernighan(graph *G, distobj *D, adddel *E, aqueue *Q,
        CClk_flipper *F, double *val, int *win_cycle, flipstack *win,
        flipstack *fstack, CCptrworld *intptr_world, CCptrworld *edgelook_world);
static double improve_tour(graph *G, distobj *D, adddel *E, aqueue *Q,
        CClk_flipper *F, int t1, flipstack *fstack, CCptrworld *intptr_world,
        CCptrworld *edgelook_world);
static int step(graph *G, distobj *D, adddel *E, aqueue *Q, CClk_flipper *F,
        int level, int gain, int *Gstar, int first, int last, flipstack *fstack,
        CCptrworld *intptr_world, CCptrworld *edgelook_world);
static int step_noback(graph *G, distobj *D, adddel *E, aqueue *Q,
        CClk_flipper *F, int level, int gain, int *Gstar, int first, int last,
        flipstack *fstack, CCptrworld *intptr_world);
static double kick_improve(graph *G, distobj *D, adddel *E, aqueue *Q,
        CClk_flipper *F, flipstack *win, flipstack *fstack,
        CCptrworld *intptr_world);
static int kick_step_noback(graph *G, distobj *D, adddel *E, aqueue *Q,
        CClk_flipper *F, int level, int gain, int *Gstar, int first, int last,
        flipstack *win, flipstack *fstack, CCptrworld *intptr_world);
static int weird_second_step(graph *G, distobj *D, adddel *E, aqueue *Q,
        CClk_flipper *F, int len_t1_t2, int t1, int t2, flipstack *fstack,
        CCptrworld *intptr_world, CCptrworld *edgelook_world);
static edgelook *look_ahead(graph *G, distobj *D, adddel *E, CClk_flipper *F,
        int first, int last, int gain, int level, CCptrworld *edgelook_world);
static void look_ahead_noback(graph *G, distobj *D, adddel *E, CClk_flipper *F,
        int first, int last, int gain, edgelook *winner);
static edgelook *weird_look_ahead(graph *G, distobj *D, CClk_flipper *F,
        int gain, int t1, int t2, CCptrworld *edgelook_world);
static edgelook *weird_look_ahead2(graph *G, distobj *D, CClk_flipper *F,
        int gain, int t2, int t3, int t4, CCptrworld *edgelook_world);
static edgelook *weird_look_ahead3(graph *G, distobj *D, CClk_flipper *F,
        int gain, int t2, int t3, int t6, CCptrworld *edgelook_world);
static double cycle_length(int ncount, int *cyc, distobj *D);
static int random_four_swap(graph *G, distobj *D, aqueue *Q, CClk_flipper *F,
        CCkdtree *kdt, int *delta, flipstack *win, flipstack *fstack,
        CCptrworld *intptr_world, CCrandstate *rstate);
static void first_kicker(graph *G, distobj *D, CClk_flipper *F, int *t1,
        int *t2);
static void find_walk_four(graph *G, distobj *D, CClk_flipper *F, int *t1,
        int *t2, int *t3, int *t4, int *t5, int *t6, int *t7, int *t8);
static void turn(int n, aqueue *Q, CCptrworld *intptr_world);
static void kickturn(int n, aqueue *Q, distobj *D, graph *G, CClk_flipper *F,
        CCptrworld *intptr_world);
static void bigturn(graph *G, int n, int tonext, aqueue *Q, CClk_flipper *F,
        distobj *D, CCptrworld *intptr_world);
static void randcycle(int ncount, int *cyc, CCrandstate *rstate);
static void initgraph(graph *G);
static void freegraph(graph *G);
static int buildgraph(graph *G, int ncount, int ecount, int *elist, distobj *D);
static void insertedge(graph *G, int n1, int n2, int w);
static void linkern_free_world(CCptrworld *intptr_world,
        CCptrworld *edgelook_world);
static int init_flipstack(flipstack *f, int total, int single);
static void free_flipstack(flipstack *f);
static void init_adddel(adddel *E);
static void free_adddel(adddel *E);
static int build_adddel(adddel *E, int ncount);
static void init_aqueue(aqueue *Q);
static void free_aqueue(aqueue *Q, CCptrworld *intptr_world);
static int build_aqueue(aqueue *Q, int ncount, CCptrworld *intptr_world);
static void add_to_active_queue(int n, aqueue *Q, CCptrworld *intptr_world);
static int pop_from_active_queue(aqueue *Q, CCptrworld *intptr_world);
static void init_distobj(distobj *D);
static void free_distobj(distobj *D);
static int build_distobj(distobj *D, int ncount, CCdatagroup *dat);
static int CCutil_dat_edgelen(int i, int j, CCdatagroup *dat);
static int dist(int i, int j, distobj *D);
static int CClinkern_flipper_init(CClk_flipper *F, int ncount, int *cyc);
static void CClinkern_flipper_cycle(CClk_flipper *F, int *x);
static void CClinkern_flipper_finish(CClk_flipper *F);
static int CClinkern_flipper_next(CClk_flipper *F, int x);
static int CClinkern_flipper_prev(CClk_flipper *F, int x);
static void CClinkern_flipper_flip(CClk_flipper *F, int x, int y);
static void same_segment_flip(CClk_flipper *F, CClk_childnode *a,
        CClk_childnode *b);
static void consecutive_segment_flip(CClk_flipper *F, CClk_parentnode *a,
        CClk_parentnode *b);
static void segment_split(CClk_flipper *F, CClk_parentnode *p,
        CClk_childnode *aprev, CClk_childnode *a, int left_or_right);
static int CClinkern_flipper_sequence(CClk_flipper *F, int x, int y, int z);
static void init_flipper(CClk_flipper *Fl);
static void free_flipper(CClk_flipper *Fl);
static int build_flipper(CClk_flipper *Fl, int ncount);
static void CCutil_dheap_free(CCdheap *h);

CC_PTRWORLD_ROUTINES(intptr, intptralloc, intptr_bulkalloc, intptrfree)
CC_PTRWORLD_LISTFREE_ROUTINE(intptr, intptr_listfree, intptrfree)
CC_PTRWORLD_ROUTINES(edgelook, edgelookalloc, edgelook_bulkalloc, edgelookfree)
CC_PTRWORLD_LISTFREE_ROUTINE(edgelook, edgelook_listfree, edgelookfree)

void CCutil_sprand(int seed, CCrandstate *r)
{
    int i, ii;
    int last, next;
    int *arr = r->arr;

    seed %= CC_PRANDMAX;
    if (seed < 0)
        seed += CC_PRANDMAX;

    arr[0] = last = seed;
    next = 1;
    for (i = 1; i < 55; i++) {
        ii = (21 * i) % 55;
        arr[ii] = next;
        next = last - next;
        if (next < 0)
            next += CC_PRANDMAX;
        last = arr[ii];
    }
    r->a = 0;
    r->b = 24;
    for (i = 0; i < 165; i++)
        last = CCutil_lprand(r);
}

int CCutil_graph2dat_matrix(int ncount, int ecount, int *elist, int *elen,
        int defaultlen, CCdatagroup *dat)
{
    int i, j, k;
    int rval;

    CCutil_init_datagroup(dat);
    dat->adj = CC_SAFE_MALLOC(ncount, int *);
    dat->adjspace = CC_SAFE_MALLOC(ncount * (ncount+1) / 2, int);
    if (dat->adj == (int **) NULL || dat->adjspace == (int *) NULL) {
        fprintf(stderr, "Our of memory in CCutil_graph2dat\n");
        rval = 1;
        goto CLEANUP;
    }

    for (i = 0, j = 0; i < ncount; i++) {
        dat->adj[i] = dat->adjspace + j;
        j += (i + 1);
    }
    for (i = 0, k = 0; i < ncount; i++) {
        for (j = 0; j < i; j++) {
            dat->adj[i][j] = defaultlen;
        }
        dat->adj[i][i] = 0;
    }
    for (i = 0; i < ecount; i++) {
        j = elist[2 * i];
        k = elist[2 * i + 1];
        if (j < k)
            dat->adj[k][j] = elen[i];
        else
            dat->adj[j][k] = elen[i];
    }
    dat->edgelen = matrix_edgelen;
    rval = 0;

CLEANUP:

    if (rval) {
        CCutil_freedatagroup(dat);
    }
    return rval;
}

void CCutil_freedatagroup(CCdatagroup *dat)
{
    CC_IFFREE(dat->x, double);
    CC_IFFREE(dat->y, double);
    CC_IFFREE(dat->z, double);
    CC_IFFREE(dat->adj, int *);
    CC_IFFREE(dat->adjspace, int);
    CC_IFFREE(dat->len, int *);
    CC_IFFREE(dat->lenspace, int);
    CC_IFFREE(dat->degree, int);
    free_userdat(&dat->userdat);
    free_rhdata(&dat->rhdat);
    CC_IFFREE(dat->depotcost, int);
    CC_IFFREE(dat->orig_names, int);
}

int CClinkern_tour(int ncount, CCdatagroup *dat, int ecount, int *elist,
        int stallcount, int repeatcount, int *incycle, int *outcycle,
        double *val, CCrandstate *rstate)
{
    int rval = 0;
    int i;
    int *tcyc = (int *) NULL;
    graph G;
    distobj D;
    CCptrworld intptr_world;
    CCptrworld edgelook_world;

    initgraph(&G);
    init_distobj(&D);
    CCptrworld_init(&intptr_world);
    CCptrworld_init(&edgelook_world);
    G.rstate = rstate;

    if (ncount < 10 && repeatcount > 0) {
        repeatcount = 0;
    }

    rval = intptr_bulkalloc(&intptr_world, ncount);
    if (rval) {
        fprintf(stderr, "Unable to allocate initial intptrs\n");
        goto CLEANUP;
    }

    rval = edgelook_bulkalloc(&edgelook_world, MAX_BACK * (BACKTRACK + 3));
    if (rval) {
        fprintf(stderr, "Unable to allocate initial edgelooks\n");
        goto CLEANUP;
    }

    tcyc = CC_SAFE_MALLOC (ncount, int);
    if (tcyc == (int *) NULL) {
        fprintf(stderr, "out of memory in linkern\n");
        rval = 1;
        goto CLEANUP;
    }

    rval = build_distobj(&D, ncount, dat);
    if (rval)
        goto CLEANUP;

    rval = buildgraph(&G, ncount, ecount, elist, &D);
    if (rval) {
        fprintf(stderr, "buildgraph failed\n");
        goto CLEANUP;
    }

    if (incycle) {
        for (i = 0; i < ncount; i++)
            tcyc[i] = incycle[i];
    } else {
        randcycle(ncount, tcyc, G.rstate);
    }
    *val = cycle_length(ncount, tcyc, &D);

    rval = repeated_lin_kernighan(&G, &D, tcyc, stallcount, repeatcount, val,
            &intptr_world, &edgelook_world, rstate);
    if (rval) {
        fprintf(stderr, "repeated_lin_kernighan failed\n");
        goto CLEANUP;
    }

    if (outcycle) {
        for (i = 0; i < ncount; i++)
            outcycle[i] = tcyc[i];
    }

CLEANUP:

    CC_IFFREE(tcyc, int);
    freegraph(&G);
    free_distobj(&D);
    linkern_free_world(&intptr_world, &edgelook_world);

    return rval;
}

static int CCutil_lprand(CCrandstate *r)
{
    int t;

    if (r->a-- == 0)
        r->a = 54;
    if (r->b-- == 0)
        r->b = 54;

    t = r->arr[r->a] - r->arr[r->b];

    if (t < 0)
        t += CC_PRANDMAX;

    r->arr[r->a] = t;

    return t;
}

static void free_rhdata(CCdata_rhvector *rhdat)
{
    CC_IFFREE(rhdat->space, char);
    CC_IFFREE(rhdat->vectors, char *);
    rhdat->rhlength = 0;
}

static void free_userdat(CCdata_user *userdat)
{
    CC_IFFREE(userdat->x, double);
    CC_IFFREE(userdat->y, double);
}

static int matrix_edgelen(int i, int j, CCdatagroup *dat)
{
    if (i > j)
        return (dat->adj[i])[j];
    else
        return (dat->adj[j])[i];
}

static int edgelen_nonorm(int i, int j, CCdatagroup *dat)
{
    return -1;
}

static void init_userdat(CCdata_user *userdat)
{
    userdat->x = (double *) NULL;
    userdat->y = (double *) NULL;
}

static void init_rhdata(CCdata_rhvector *rhdat)
{
    rhdat->space = (char *) NULL;
    rhdat->vectors = (char **) NULL;
    rhdat->rhlength = 0;
    rhdat->dist_00 = 0;
    rhdat->dist_01 = 0;
    rhdat->dist_02 = 0;
    rhdat->dist_22 = 0;
    rhdat->p = 0.0;
}

static void CCutil_init_datagroup(CCdatagroup *dat)
{
    dat->x = (double *) NULL;
    dat->y = (double *) NULL;
    dat->z = (double *) NULL;
    dat->adj = (int **) NULL;
    dat->adjspace = (int *) NULL;
    dat->len = (int **) NULL;
    dat->lenspace = (int *) NULL;
    dat->degree = (int *) NULL;
    dat->norm = 0;
    dat->dsjrand_param = 1;
    dat->dsjrand_factor = 1.0;
    dat->default_len = 100000;
    dat->sparse_ecount = 0;
    dat->edgelen = edgelen_nonorm;
    init_userdat(&dat->userdat);
    init_rhdata(&dat->rhdat);
    dat->ndepot = 0;
    dat->orig_ncount = 0;
    dat->depotcost = (int *) NULL;
    dat->orig_names = (int *) NULL;
}

static void *CCutil_allocrus(size_t size)
{
    void *mem = (void *) NULL;

    if (size == 0) {
        fprintf(stderr, "Warning: 0 bytes allocated\n");
    }

    mem = (void *) malloc(size);
    if (mem == (void *) NULL) {
        fprintf(stderr, "Out of memory. Asked for %d bytes\n", (int) size);
    }
    return mem;
}

static void CCutil_freerus(void *p)
{
    if (!p) {
        fprintf(stderr, "Warning: null pointer freed\n");
        return;
    }

    free(p);
}

static CCbigchunkptr *CCutil_bigchunkalloc(void)
{
    CCbigchunk *p = CC_SAFE_MALLOC (1, CCbigchunk);

    if (p == (CCbigchunk *) NULL) {
        fprintf(stderr, "Out of memory in CCutil_bigchunkalloc\n");
        return (CCbigchunkptr *) NULL;
    }
    p->ptr.this_chunk = p;
    p->ptr.this_one = (void *) p->space;
    return &(p->ptr);
}

static void CCutil_bigchunkfree(CCbigchunkptr *bp)
{
    CCbigchunk *p = bp->this_chunk;
    CC_FREE(p, CCbigchunk);
}

static void CCptrworld_init(CCptrworld *world)
{
    world->refcount = 1;
    world->freelist = (void *) NULL;
    world->chunklist = (CCbigchunkptr *) NULL;
}

static void CCptrworld_delete(CCptrworld *world)
{
    world->refcount--;
    if (world->refcount <= 0) {
        CCbigchunkptr *bp, *bpnext;

        for (bp = world->chunklist; bp; bp = bpnext) {
            bpnext = bp->next;
            CCutil_bigchunkfree(bp);
        }
        world->chunklist = (CCbigchunkptr *) NULL;
        world->freelist = (void *) NULL;
        world->refcount = 0;
    }
}

static int repeated_lin_kernighan(graph *G, distobj *D, int *cyc,
        int stallcount, int count, double *val, CCptrworld *intptr_world,
        CCptrworld *edgelook_world, CCrandstate *rstate)
{
    int rval = 0;
    int round = 0;
    int quitcount, hit, delta;
    int *win_cycle = (int *) NULL;
    CCkdtree kdt;
    flipstack winstack, fstack;
    double t, best = *val;
    int ncount = G->ncount;
    adddel E;
    CClk_flipper F;
    aqueue Q;
    int *tcyc = (int *) NULL;
    int i;

    init_aqueue(&Q);
    init_adddel(&E);
    rval = build_aqueue(&Q, ncount, intptr_world);
    if (rval) {
        fprintf(stderr, "build_aqueue failed\n");
        goto CLEANUP;
    }
    rval = build_adddel(&E, ncount);
    if (rval) {
        fprintf(stderr, "build_adddel failed\n");
        goto CLEANUP;
    }

    hit = 2 * (MAXDEPTH + 7 + KICK_MAXDEPTH);
    rval = init_flipstack(&fstack, hit, 0);
    if (rval) {
        fprintf(stderr, "init_flipstack failed\n");
        goto CLEANUP;
    }
    rval = init_flipstack(&winstack, 500 + ncount / 50, hit);
    if (rval) {
        fprintf(stderr, "init_flipstack failed\n");
        goto CLEANUP;
    }

    win_cycle = CC_SAFE_MALLOC (ncount, int);
    if (win_cycle == (int *) NULL) {
        fprintf(stderr, "out of memory in repeated_lin_kernighan\n");
        rval = 1;
        goto CLEANUP;
    }
    win_cycle[0] = -1;

    quitcount = stallcount;
    if (quitcount > count)
        quitcount = count;

    CClinkern_flipper_init(&F, ncount, cyc);
    fstack.counter = 0;
    winstack.counter = 0;
    win_cycle[0] = -1;

    tcyc = CC_SAFE_MALLOC (ncount, int);
    if (tcyc == (int *) NULL) {
        fprintf(stderr, "out of memory in repeated_lin_kernighan\n");
        rval = 1;
        goto CLEANUP;
    }
    /* init active_queue with random order */
    randcycle(ncount, tcyc, G->rstate);
    for (i = 0; i < ncount; i++) {
        add_to_active_queue(tcyc[i], &Q, intptr_world);
    }CC_IFFREE(tcyc, int);

    lin_kernighan(G, D, &E, &Q, &F, &best, win_cycle, &winstack, &fstack,
            intptr_world, edgelook_world);

    winstack.counter = 0;
    win_cycle[0] = -1;

    while (round < quitcount) {
        hit = 0;
        fstack.counter = 0;

        if (IMPROVE_SWITCH == -1 || round < IMPROVE_SWITCH) {
            rval = random_four_swap(G, D, &Q, &F, &kdt, &delta, &winstack,
                    &fstack, intptr_world, rstate);
            if (rval) {
                fprintf(stderr, "random_four_swap failed\n");
                goto CLEANUP;
            }
        } else {
            delta = kick_improve(G, D, &E, &Q, &F, &winstack, &fstack,
                    intptr_world);
        }

        fstack.counter = 0;
        t = best + delta;
        lin_kernighan(G, D, &E, &Q, &F, &t, win_cycle, &winstack, &fstack,
                intptr_world, edgelook_world);

        if (t <= best) {
            winstack.counter = 0;
            win_cycle[0] = -1;
            if (t < best) {
                best = t;
                quitcount = round + stallcount;
                if (quitcount > count)
                    quitcount = count;
                hit++;
            }
        } else {
            if (win_cycle[0] == -1) {
                while (winstack.counter) {
                    winstack.counter--;
                    CClinkern_flipper_flip(&F,
                            winstack.stack[winstack.counter].last,
                            winstack.stack[winstack.counter].first);
                }
            } else {
                CClinkern_flipper_finish(&F);
                CClinkern_flipper_init(&F, ncount, win_cycle);
                while (winstack.counter) {
                    winstack.counter--;
                    CClinkern_flipper_flip(&F,
                            winstack.stack[winstack.counter].last,
                            winstack.stack[winstack.counter].first);
                }
                win_cycle[0] = -1;
            }
        }

        round++;
    }

    CClinkern_flipper_cycle(&F, cyc);
    CClinkern_flipper_finish(&F);

    t = cycle_length(ncount, cyc, D);
    if (t != best) {
        best = t;
    }
    *val = best;

CLEANUP:

    free_aqueue(&Q, intptr_world);
    free_adddel(&E);
    free_flipstack(&fstack);
    free_flipstack(&winstack);
    CC_IFFREE(win_cycle, int);
    return rval;
}

static void lin_kernighan(graph *G, distobj *D, adddel *E, aqueue *Q,
        CClk_flipper *F, double *val, int *win_cycle, flipstack *win,
        flipstack *fstack, CCptrworld *intptr_world, CCptrworld *edgelook_world)
{
    int start, i;
    double delta, totalwin = 0.0;

    while (1) {
        start = pop_from_active_queue(Q, intptr_world);
        if (start == -1)
            break;

        delta = improve_tour(G, D, E, Q, F, start, fstack, intptr_world,
                edgelook_world);
        if (delta > 0.0) {
            totalwin += delta;
            if (win->counter < win->max) {
                for (i = 0; i < fstack->counter; i++) {
                    win->stack[win->counter].first = fstack->stack[i].first;
                    win->stack[win->counter].last = fstack->stack[i].last;
                    win->stack[win->counter].firstprev =
                            fstack->stack[i].firstprev;
                    win->stack[win->counter].lastnext =
                            fstack->stack[i].lastnext;
                    win->counter++;
                }
            } else if (win_cycle[0] == -1) {
                for (i = 0; i < fstack->counter; i++) {
                    win->stack[win->counter].first = fstack->stack[i].first;
                    win->stack[win->counter].last = fstack->stack[i].last;
                    win->counter++;
                }
                CClinkern_flipper_cycle(F, win_cycle);
            }
            fstack->counter = 0;
        }
    }

    if (win_cycle[0] == -1) {
        for (i = 0; i < fstack->counter; i++) {
            win->stack[win->counter].first = fstack->stack[i].first;
            win->stack[win->counter].last = fstack->stack[i].last;
            win->counter++;
        }
    }
    (*val) -= totalwin;
}

static double improve_tour(graph *G, distobj *D, adddel *E, aqueue *Q,
        CClk_flipper *F, int t1, flipstack *fstack, CCptrworld *intptr_world,
        CCptrworld *edgelook_world)
{
    int t2 = CClinkern_flipper_next(F, t1);
    int gain, Gstar = 0;

    gain = Edgelen(t1, t2, D);
    markedge_del(t1, t2, E);

    if (step(G, D, E, Q, F, 0, gain, &Gstar, t1, t2, fstack, intptr_world,
            edgelook_world) == 0) {
        Gstar = weird_second_step(G, D, E, Q, F, gain, t1, t2, fstack,
                intptr_world, edgelook_world);
    }unmarkedge_del(t1, t2, E);

    if (Gstar) {
        MARK(t1, Q, F, D, G, intptr_world);
        MARK(t2, Q, F, D, G, intptr_world);
    }
    return (double) Gstar;
}

static int step(graph *G, distobj *D, adddel *E, aqueue *Q, CClk_flipper *F,
        int level, int gain, int *Gstar, int first, int last, flipstack *fstack,
        CCptrworld *intptr_world, CCptrworld *edgelook_world)
{
    int val, this, newlast, hit = 0, oldG = gain;
    edgelook *list, *e;

    if (level >= BACKTRACK) {
        return step_noback(G, D, E, Q, F, level, gain, Gstar, first, last,
                fstack, intptr_world);
    }

    list = look_ahead(G, D, E, F, first, last, gain, level, edgelook_world);
    for (e = list; e; e = e->next) {
        this = e->other;
        newlast = e->over;

        gain = oldG - e->diff;
        val = gain - Edgelen (newlast, first, D);
        if (val > *Gstar) {
            *Gstar = val;
            hit++;
        }

        FLIP(first, last, newlast, this, fstack, F);

        if (level < MAXDEPTH) {
            markedge_add(last, this, E);
            markedge_del(this, newlast, E);
            hit += step(G, D, E, Q, F, level + 1, gain, Gstar, first, newlast,
                    fstack, intptr_world, edgelook_world);
            unmarkedge_add(last, this, E);
            unmarkedge_del(this, newlast, E);
        }

        if (!hit) {
            UNFLIP(first, last, newlast, this, fstack, F);
        } else {
            MARK(this, Q, F, D, G, intptr_world);
            MARK(newlast, Q, F, D, G, intptr_world);
            edgelook_listfree(edgelook_world, list);
            return 1;
        }
    }
    edgelook_listfree(edgelook_world, list);
    return 0;
}

static int step_noback(graph *G, distobj *D, adddel *E, aqueue *Q,
        CClk_flipper *F, int level, int gain, int *Gstar, int first, int last,
        flipstack *fstack, CCptrworld *intptr_world)
{
    edgelook e; e.other = -1; e.over = -1;

    look_ahead_noback(G, D, E, F, first, last, gain - *Gstar - level, &e);

    if (e.diff < BIGINT) {
        {
            if (e.mm) {
                int hit = 0;
                int this = e.other;
                int newfirst = e.over;
                int val;

                gain -= e.diff;
                val = gain - Edgelen (newfirst, last, D);
                if (val > *Gstar) {
                    *Gstar = val;
                    hit++;
                }
                FLIP(this, newfirst, first, last, fstack, F);

                if (level < MAXDEPTH) {
                    markedge_add(first, this, E);
                    markedge_del(this, newfirst, E);
                    hit += step_noback(G, D, E, Q, F, level + 1, gain, Gstar,
                            newfirst, last, fstack, intptr_world);
                    unmarkedge_add(first, this, E);
                    unmarkedge_del(this, newfirst, E);
                }

                if (!hit) {
                    UNFLIP(this, newfirst, first, last, fstack, F);
                    return 0;
                } else {
                    MARK(this, Q, F, D, G, intptr_world);
                    MARK(newfirst, Q, F, D, G, intptr_world);
                    return 1;
                }
            } else {
                int hit = 0;
                int this = e.other;
                int newlast = e.over;
                int val;

                gain -= e.diff;
                val = gain - Edgelen (newlast, first, D);
                if (val > *Gstar) {
                    *Gstar = val;
                    hit++;
                }

                FLIP(first, last, newlast, this, fstack, F);

                if (level < MAXDEPTH) {
                    markedge_add(last, this, E);
                    markedge_del(this, newlast, E);
                    hit += step_noback(G, D, E, Q, F, level + 1, gain, Gstar,
                            first, newlast, fstack, intptr_world);
                    unmarkedge_add(last, this, E);
                    unmarkedge_del(this, newlast, E);
                }

                if (!hit) {
                    UNFLIP(first, last, newlast, this, fstack, F);
                    return 0;
                } else {
                    MARK(this, Q, F, D, G, intptr_world);
                    MARK(newlast, Q, F, D, G, intptr_world);
                    return 1;
                }
            }
        }
    } else {
        return 0;
    }
}

static double kick_improve(graph *G, distobj *D, adddel *E, aqueue *Q,
        CClk_flipper *F, flipstack *win, flipstack *fstack,
        CCptrworld *intptr_world)
{
    int t1, t2;
    int gain, Gstar = 0;
    int hit = 0;

    do {
        first_kicker(G, D, F, &t1, &t2);
        gain = Edgelen (t1, t2, D);
        markedge_del(t1, t2, E);
        hit = kick_step_noback(G, D, E, Q, F, 0, gain, &Gstar, t1, t2, win,
                fstack, intptr_world);
        unmarkedge_del(t1, t2, E);
    } while (!hit);

    kickturn(t1, Q, D, G, F, intptr_world);
    kickturn(t2, Q, D, G, F, intptr_world);

    return (double) -Gstar;
}

#define G_MULT 1.5

static int kick_step_noback(graph *G, distobj *D, adddel *E, aqueue *Q,
        CClk_flipper *F, int level, int gain, int *Gstar, int first, int last,
        flipstack *win, flipstack *fstack, CCptrworld *intptr_world)
{
    edgelook winner;
    int val;
    int this, prev, newlast;
    int lastnext = CClinkern_flipper_next(F, last);
    int i;
    int cutoff = (int) (G_MULT * (double) gain);
    edge **goodlist = G->goodlist;

    winner.diff = BIGINT;
    for (i = 0; goodlist[last][i].weight < cutoff; i++) {
        this = goodlist[last][i].other;
        if (!is_it_deleted (last, this, E) && this != first
                && this != lastnext) {
            prev = CClinkern_flipper_prev(F, this);
            if (!is_it_added (this, prev, E)) {
                val = goodlist[last][i].weight - Edgelen (this, prev, D);
                if (val < winner.diff) {
                    winner.diff = val;
                    winner.other = this;
                    winner.over = prev;
                }
            }
        }
    }

    if (winner.diff < BIGINT) {
        this = winner.other;
        newlast = winner.over;
        gain -= winner.diff;
        *Gstar = gain - Edgelen (newlast, first, D);

        FLIP(first, last, newlast, this, fstack, F);
        kickturn(this, Q, D, G, F, intptr_world);
        kickturn(newlast, Q, D, G, F, intptr_world);
        if (win->counter < win->max) {
            win->stack[win->counter].first = last;
            win->stack[win->counter].last = newlast;
            win->counter++;
        }

        if (level < KICK_MAXDEPTH) {
            markedge_add(last, this, E);
            markedge_del(this, newlast, E);
            kick_step_noback(G, D, E, Q, F, level + 1, gain, Gstar, first,
                    newlast, win, fstack, intptr_world);
            unmarkedge_add(last, this, E);
            unmarkedge_del(this, newlast, E);
        }
        return 1;
    } else {
        return 0;
    }
}

static int weird_second_step(graph *G, distobj *D, adddel *E, aqueue *Q,
        CClk_flipper *F, int len_t1_t2, int t1, int t2, flipstack *fstack,
        CCptrworld *intptr_world, CCptrworld *edgelook_world)
{
    int t3, t4, t5, t6, t7, t8;
    int oldG, gain, tG, Gstar = 0, val, hit;
    int t4next;
    edgelook *e, *f, *h, *list, *list2, *list3;

    list = weird_look_ahead(G, D, F, len_t1_t2, t1, t2, edgelook_world);
    for (h = list; h; h = h->next) {
        t3 = h->other;
        t4 = h->over;

        oldG = len_t1_t2 - h->diff;

        t4next = CClinkern_flipper_next(F, t4);

        markedge_add(t2, t3, E);
        markedge_del(t3, t4, E);
        G->weirdmagic++;
        G->weirdmark[t1] = G->weirdmagic;
        G->weirdmark[t2] = G->weirdmagic;
        G->weirdmark[t3] = G->weirdmagic;
        G->weirdmark[t4next] = G->weirdmagic;

        list2 = weird_look_ahead2(G, D, F, oldG, t2, t3, t4, edgelook_world);
        for (e = list2; e; e = e->next) {
            t5 = e->other;
            t6 = e->over;

            markedge_add(t4, t5, E);
            if (e->seq) {
                if (!e->side) {
                    gain = oldG - e->diff;
                    val = gain - Edgelen (t6, t1, D);
                    if (val > Gstar)
                        Gstar = val;
                    FLIP(t1, t2, t6, t5, fstack, F);
                    FLIP(t2, t5, t3, t4, fstack, F);

                    markedge_del(t5, t6, E);
                    hit = step(G, D, E, Q, F, 2, gain, &Gstar, t1, t6, fstack,
                            intptr_world, edgelook_world);
                    unmarkedge_del(t5, t6, E);

                    if (!hit && Gstar)
                        hit = 1;

                    if (!hit) {
                        UNFLIP(t2, t5, t3, t4, fstack, F);
                        UNFLIP(t1, t2, t6, t5, fstack, F);
                    } else {
                        unmarkedge_add(t2, t3, E);
                        unmarkedge_del(t3, t4, E);
                        unmarkedge_add(t4, t5, E);
                        MARK(t3, Q, F, D, G, intptr_world);
                        MARK(t4, Q, F, D, G, intptr_world);
                        MARK(t5, Q, F, D, G, intptr_world);
                        MARK(t6, Q, F, D, G, intptr_world);
                        edgelook_listfree(edgelook_world, list);
                        edgelook_listfree(edgelook_world, list2);
                        return Gstar;
                    }
                } else {
                    gain = oldG - e->diff;
                    val = gain - Edgelen (t6, t1, D);
                    if (val > Gstar)
                        Gstar = val;
                    FLIP(t1, t2, t3, t4, fstack, F);
                    FLIP(t6, t5, t2, t4, fstack, F);
                    FLIP(t1, t3, t6, t2, fstack, F);

                    markedge_del(t5, t6, E);
                    hit = step(G, D, E, Q, F, 2, gain, &Gstar, t1, t6, fstack,
                            intptr_world, edgelook_world);
                    unmarkedge_del(t5, t6, E);

                    if (!hit && Gstar)
                        hit = 1;

                    if (!hit) {
                        UNFLIP(t1, t3, t6, t2, fstack, F);
                        UNFLIP(t6, t5, t2, t4, fstack, F);
                        UNFLIP(t1, t2, t3, t4, fstack, F);
                    } else {
                        unmarkedge_add(t2, t3, E);
                        unmarkedge_del(t3, t4, E);
                        unmarkedge_add(t4, t5, E);
                        MARK(t3, Q, F, D, G, intptr_world);
                        MARK(t4, Q, F, D, G, intptr_world);
                        MARK(t5, Q, F, D, G, intptr_world);
                        MARK(t6, Q, F, D, G, intptr_world);
                        edgelook_listfree(edgelook_world, list);
                        edgelook_listfree(edgelook_world, list2);
                        return Gstar;
                    }
                }
            } else {
                tG = oldG - e->diff;
                markedge_del(t5, t6, E);
                list3 = weird_look_ahead3(G, D, F, tG, t2, t3, t6,
                        edgelook_world);
                for (f = list3; f; f = f->next) {
                    t7 = f->other;
                    t8 = f->over;
                    gain = tG - f->diff;
                    if (!f->side) {
                        val = gain - Edgelen (t8, t1, D);
                        if (val > Gstar)
                            Gstar = val;
                        FLIP(t1, t2, t8, t7, fstack, F);
                        FLIP(t2, t7, t3, t4, fstack, F);
                        FLIP(t7, t4, t6, t5, fstack, F);

                        markedge_add(t6, t7, E);
                        markedge_del(t7, t8, E);
                        hit = step(G, D, E, Q, F, 3, gain, &Gstar, t1, t8,
                                fstack, intptr_world, edgelook_world);
                        unmarkedge_del(t6, t7, E);
                        unmarkedge_del(t7, t8, E);

                        if (!hit && Gstar)
                            hit = 1;

                        if (!hit) {
                            UNFLIP(t7, t4, t6, t5, fstack, F);
                            UNFLIP(t2, t7, t3, t4, fstack, F);
                            UNFLIP(t1, t2, t8, t7, fstack, F);
                        } else {
                            unmarkedge_add(t2, t3, E);
                            unmarkedge_del(t3, t4, E);
                            unmarkedge_add(t4, t5, E);
                            unmarkedge_del(t5, t6, E);
                            MARK(t3, Q, F, D, G, intptr_world);
                            MARK(t4, Q, F, D, G, intptr_world);
                            MARK(t5, Q, F, D, G, intptr_world);
                            MARK(t6, Q, F, D, G, intptr_world);
                            MARK(t7, Q, F, D, G, intptr_world);
                            MARK(t8, Q, F, D, G, intptr_world);
                            edgelook_listfree(edgelook_world, list);
                            edgelook_listfree(edgelook_world, list2);
                            edgelook_listfree(edgelook_world, list3);
                            return Gstar;
                        }
                    } else {
                        val = gain - Edgelen (t8, t1, D);
                        if (val > Gstar)
                            Gstar = val;
                        FLIP(t1, t2, t6, t5, fstack, F);
                        FLIP(t1, t6, t8, t7, fstack, F);
                        FLIP(t3, t4, t2, t5, fstack, F);

                        markedge_add(t6, t7, E);
                        markedge_del(t7, t8, E);
                        hit = step(G, D, E, Q, F, 3, gain, &Gstar, t1, t8,
                                fstack, intptr_world, edgelook_world);
                        unmarkedge_add(t6, t7, E);
                        unmarkedge_del(t7, t8, E);

                        if (!hit && Gstar)
                            hit = 1;

                        if (!hit) {
                            UNFLIP(t3, t4, t2, t5, fstack, F);
                            UNFLIP(t1, t6, t8, t7, fstack, F);
                            UNFLIP(t1, t2, t6, t5, fstack, F);
                        } else {
                            unmarkedge_add(t2, t3, E);
                            unmarkedge_del(t3, t4, E);
                            unmarkedge_add(t4, t5, E);
                            unmarkedge_del(t5, t6, E);
                            MARK(t3, Q, F, D, G, intptr_world);
                            MARK(t4, Q, F, D, G, intptr_world);
                            MARK(t5, Q, F, D, G, intptr_world);
                            MARK(t6, Q, F, D, G, intptr_world);
                            MARK(t7, Q, F, D, G, intptr_world);
                            MARK(t8, Q, F, D, G, intptr_world);
                            edgelook_listfree(edgelook_world, list);
                            edgelook_listfree(edgelook_world, list2);
                            edgelook_listfree(edgelook_world, list3);
                            return Gstar;
                        }
                    }
                }
                edgelook_listfree(edgelook_world, list3);
                unmarkedge_del(t5, t6, E);
            }
            unmarkedge_add(t4, t5, E);
        }
        edgelook_listfree(edgelook_world, list2);
        unmarkedge_add(t2, t3, E);
        unmarkedge_del(t3, t4, E);
    }
    edgelook_listfree(edgelook_world, list);
    return 0;
}

static edgelook *look_ahead(graph *G, distobj *D, adddel *E, CClk_flipper *F,
        int first, int last, int gain, int level, CCptrworld *edgelook_world)
{
    edgelook *list = (edgelook *) NULL, *el;
    int i, val;
    int this, prev;
    int lastnext = CClinkern_flipper_next(F, last);
    int other[MAX_BACK], save[MAX_BACK];
    int value[MAX_BACK + 1];int
    k, ahead = backtrack_count[level];
    edge **goodlist = G->goodlist;

    for (i = 0; i < ahead; i++) {
        value[i] = BIGINT;
    }
    value[ahead] = -BIGINT;

    for (i = 0; goodlist[last][i].weight <= gain; i++) {
        this = goodlist[last][i].other;
        if (!is_it_deleted (last, this, E) && this != first
                && this != lastnext) {
            prev = CClinkern_flipper_prev(F, this);
            if (!is_it_added (this, prev, E)) {
                val = goodlist[last][i].weight - Edgelen (this, prev, D);
                if (val < value[0]) {
                    for (k = 0; value[k + 1] > val; k++) {
                        value[k] = value[k + 1];
                        other[k] = other[k + 1];
                        save[k] = save[k + 1];
                    }
                    value[k] = val;
                    other[k] = this;
                    save[k] = prev;
                }
            }
        }
    }

    for (i = 0; i < ahead; i++) {
        if (value[i] < BIGINT) {
            el = edgelookalloc(edgelook_world);
            el->diff = value[i];
            el->other = other[i];
            el->over = save[i];
            el->next = list;
            list = el;
        }
    }

    return list;
}

static void look_ahead_noback(graph *G, distobj *D, adddel *E, CClk_flipper *F,
        int first, int last, int gain, edgelook *winner)
{
    int val;
    int this, prev;
    int lastnext = CClinkern_flipper_next(F, last);
    int i;
    int next;
    edge **goodlist = G->goodlist;
    int firstprev;

    winner->diff = BIGINT;
    for (i = 0; goodlist[last][i].weight < gain; i++) {
        this = goodlist[last][i].other;
        if (!is_it_deleted (last, this, E) && this != first
                && this != lastnext) {
            prev = CClinkern_flipper_prev(F, this);
            if (!is_it_added (this, prev, E)) {
                val = goodlist[last][i].weight - Edgelen (this, prev, D);
                if (val < winner->diff) {
                    winner->diff = val;
                    winner->other = this;
                    winner->over = prev;
                    winner->mm = 0;
                }
            }
        }
    }
    firstprev = CClinkern_flipper_prev(F, first);

    for (i = 0; goodlist[first][i].weight < gain; i++) {
        this = goodlist[first][i].other;
        if (!is_it_deleted (first, this, E) && this != last
                && this != firstprev) {
            next = CClinkern_flipper_next(F, this);
            if (!is_it_added (this, next, E)) {
                val = goodlist[first][i].weight - Edgelen (this, next, D);
                if (val < winner->diff) {
                    winner->diff = val;
                    winner->other = this;
                    winner->over = next;
                    winner->mm = 1;
                }
            }
        }
    }
}

static edgelook *weird_look_ahead(graph *G, distobj *D, CClk_flipper *F,
        int gain, int t1, int t2, CCptrworld *edgelook_world)
{
    edgelook *list, *el;
    int i, this, next;
    int other[MAX_BACK], save[MAX_BACK];
    int value[MAX_BACK + 1];int
    k, val, ahead;
    edge **goodlist = G->goodlist;

    list = (edgelook *) NULL;
    ahead = weird_backtrack_count[0];
    for (i = 0; i < ahead; i++)
        value[i] = BIGINT;
    value[ahead] = -BIGINT;

    for (i = 0; goodlist[t2][i].weight <= gain; i++) {
        this = goodlist[t2][i].other;
        if (this != t1) {
            next = CClinkern_flipper_next(F, this);
            val = goodlist[t2][i].weight - Edgelen (this, next, D);
            if (val < value[0]) {
                for (k = 0; value[k + 1] > val; k++) {
                    value[k] = value[k + 1];
                    other[k] = other[k + 1];
                    save[k] = save[k + 1];
                }
                value[k] = val;
                other[k] = this;
                save[k] = next;
            }
        }
    }
    for (i = 0; i < ahead; i++) {
        if (value[i] < BIGINT) {
            el = edgelookalloc(edgelook_world);
            el->diff = value[i];
            el->other = other[i];
            el->over = save[i];
            el->next = list;
            list = el;
        }
    }
    return list;
}

static edgelook *weird_look_ahead2(graph *G, distobj *D, CClk_flipper *F,
        int gain, int t2, int t3, int t4, CCptrworld *edgelook_world)
{
    edgelook *list = (edgelook *) NULL;
    edgelook *el;
    int i, t5, t6;
    int other[MAX_BACK], save[MAX_BACK], seq[MAX_BACK], side[MAX_BACK];
    int value[MAX_BACK + 1];int
    k, val;
    int ahead = weird_backtrack_count[1];
    edge **goodlist = G->goodlist;
    int *weirdmark = G->weirdmark;
    int weirdmagic = G->weirdmagic;

    for (i = 0; i < ahead; i++)
        value[i] = BIGINT;
    value[ahead] = -BIGINT;

    for (i = 0; goodlist[t4][i].weight <= gain; i++) {
        t5 = goodlist[t4][i].other;
        if (weirdmark[t5] != weirdmagic) {
            if (CClinkern_flipper_sequence(F, t2, t5, t3)) {
                t6 = CClinkern_flipper_prev(F, t5);
                val = goodlist[t4][i].weight - Edgelen (t5, t6, D);
                if (val < value[0]) {
                    for (k = 0; value[k + 1] > val; k++) {
                        value[k] = value[k + 1];
                        other[k] = other[k + 1];
                        save[k] = save[k + 1];
                        seq[k] = seq[k + 1];
                        side[k] = side[k + 1];
                    }
                    value[k] = val;
                    other[k] = t5;
                    save[k] = t6;
                    seq[k] = 1;
                    side[k] = 0;
                }
                t6 = CClinkern_flipper_next(F, t5);
                val = goodlist[t4][i].weight - Edgelen (t5, t6, D);
                if (val < value[0]) {
                    for (k = 0; value[k + 1] > val; k++) {
                        value[k] = value[k + 1];
                        other[k] = other[k + 1];
                        save[k] = save[k + 1];
                        seq[k] = seq[k + 1];
                        side[k] = side[k + 1];
                    }
                    value[k] = val;
                    other[k] = t5;
                    save[k] = t6;
                    seq[k] = 1;
                    side[k] = 1;
                }
            } else {
                t6 = CClinkern_flipper_prev(F, t5);
                val = goodlist[t4][i].weight - Edgelen (t5, t6, D);
                if (val < value[0]) {
                    for (k = 0; value[k + 1] > val; k++) {
                        value[k] = value[k + 1];
                        other[k] = other[k + 1];
                        save[k] = save[k + 1];
                        seq[k] = seq[k + 1];
                        side[k] = side[k + 1];
                    }
                    value[k] = val;
                    other[k] = t5;
                    save[k] = t6;
                    seq[k] = 0;
                    side[k] = 0;
                }
            }
        }
    }

    for (i = 0; i < ahead; i++) {
        if (value[i] < BIGINT) {
            el = edgelookalloc(edgelook_world);
            el->diff = value[i];
            el->other = other[i];
            el->over = save[i];
            el->seq = seq[i];
            el->side = side[i];
            el->next = list;
            list = el;
        }
    }
    return list;
}

static edgelook *weird_look_ahead3(graph *G, distobj *D, CClk_flipper *F,
        int gain, int t2, int t3, int t6, CCptrworld *edgelook_world)
{
    edgelook *list = (edgelook *) NULL;
    edgelook *el;
    int i, t7, t8;
    int other[MAX_BACK], save[MAX_BACK], side[MAX_BACK];
    int value[MAX_BACK + 1];int
    k, val;
    int ahead = weird_backtrack_count[2];
    edge **goodlist = G->goodlist;
    int *weirdmark = G->weirdmark;
    int weirdmagic = G->weirdmagic;

    for (i = 0; i < ahead; i++)
        value[i] = BIGINT;
    value[ahead] = -BIGINT;

    for (i = 0; goodlist[t6][i].weight <= gain; i++) {
        t7 = goodlist[t6][i].other; /* Need t7 != t2, t3, t2next, t3prev */
        if (weirdmark[t7] != weirdmagic
                && CClinkern_flipper_sequence(F, t2, t7, t3)) {
            t8 = CClinkern_flipper_prev(F, t7);
            val = goodlist[t6][i].weight - Edgelen (t7, t8, D);
            if (val < value[0]) {
                for (k = 0; value[k + 1] > val; k++) {
                    value[k] = value[k + 1];
                    other[k] = other[k + 1];
                    save[k] = save[k + 1];
                    side[k] = side[k + 1];
                }
                value[k] = val;
                other[k] = t7;
                save[k] = t8;
                side[k] = 0;
            }
            t8 = CClinkern_flipper_next(F, t7);
            val = goodlist[t6][i].weight - Edgelen (t7, t8, D);
            if (val < value[0]) {
                for (k = 0; value[k + 1] > val; k++) {
                    value[k] = value[k + 1];
                    other[k] = other[k + 1];
                    save[k] = save[k + 1];
                    side[k] = side[k + 1];
                }
                value[k] = val;
                other[k] = t7;
                save[k] = t8;
                side[k] = 1;
            }
        }
    }

    for (i = 0; i < ahead; i++) {
        if (value[i] < BIGINT) {
            el = edgelookalloc(edgelook_world);
            el->diff = value[i];
            el->other = other[i];
            el->over = save[i];
            el->side = side[i];
            el->next = list;
            list = el;
        }
    }
    return list;
}

static double cycle_length(int ncount, int *cyc, distobj *D)
{
    int i;
    double val = 0.0;

    for (i = 1; i < ncount; i++) {
        val += (double) Edgelen (cyc[i - 1], cyc[i], D);
    }
    val += (double) Edgelen (cyc[0], cyc[ncount - 1], D);

    return val;
}

static int random_four_swap(graph *G, distobj *D, aqueue *Q, CClk_flipper *F,
        CCkdtree *kdt, int *delta, flipstack *win, flipstack *fstack,
        CCptrworld *intptr_world, CCrandstate *rstate)
{
    int t1, t2, t3, t4, t5, t6, t7, t8, temp;

    find_walk_four(G, D, F, &t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8);

    if (!CClinkern_flipper_sequence(F, t1, t3, t5)) {
        CC_SWAP(t3, t5, temp);
        CC_SWAP(t4, t6, temp);
    }
    if (!CClinkern_flipper_sequence(F, t1, t5, t7)) {
        CC_SWAP(t5, t7, temp);
        CC_SWAP(t6, t8, temp);
        if (!CClinkern_flipper_sequence(F, t1, t3, t5)) {
            CC_SWAP(t3, t5, temp);
            CC_SWAP(t4, t6, temp);
        }
    }
    FLIP(t1, t2, t5, t6, fstack, F);
    FLIP(t4, t3, t7, t8, fstack, F);
    FLIP(t1, t5, t6, t2, fstack, F);

    if (win->counter < win->max) {
        win->stack[win->counter].first = t2;
        win->stack[win->counter].last = t5;
        win->counter++;
    }
    if (win->counter < win->max) {
        win->stack[win->counter].first = t3;
        win->stack[win->counter].last = t7;
        win->counter++;
    }
    if (win->counter < win->max) {
        win->stack[win->counter].first = t5;
        win->stack[win->counter].last = t6;
        win->counter++;
    }

    bigturn(G, t1, 0, Q, F, D, intptr_world);
    bigturn(G, t2, 1, Q, F, D, intptr_world);
    bigturn(G, t3, 0, Q, F, D, intptr_world);
    bigturn(G, t4, 1, Q, F, D, intptr_world);
    bigturn(G, t5, 0, Q, F, D, intptr_world);
    bigturn(G, t6, 1, Q, F, D, intptr_world);
    bigturn(G, t7, 0, Q, F, D, intptr_world);
    bigturn(G, t8, 1, Q, F, D, intptr_world);

    *delta = Edgelen (t1, t6, D) + Edgelen (t2, t5, D) + Edgelen (t3, t8, D)
            + Edgelen (t4, t7, D) - Edgelen (t1, t2, D) - Edgelen (t3, t4, D)
            - Edgelen (t5, t6, D) - Edgelen (t7, t8, D);
    return 0;
}

#define HUNT_PORTION_LONG 0.001

static void first_kicker(graph *G, distobj *D, CClk_flipper *F, int *t1,
        int *t2)
{
    int longcount = (int) ((double) G->ncount * HUNT_PORTION_LONG) + 10;
    int i, best, try1, len, next, prev, nextl, prevl;
    int ncount = G->ncount;
    edge **goodlist = G->goodlist;

    try1 = CCutil_lprand(G->rstate) % ncount;
    next = CClinkern_flipper_next(F, try1);
    prev = CClinkern_flipper_prev(F, try1);
    nextl = Edgelen (try1, next, D);
    prevl = Edgelen (try1, prev, D);
    if (nextl >= prevl) {
        *t1 = try1;
        *t2 = next;
        best = nextl - goodlist[*t1][0].weight;
    } else {
        *t1 = prev;
        *t2 = try1;
        best = prevl - goodlist[*t1][0].weight;
    }

    for (i = 0; i < longcount; i++) {
        try1 = CCutil_lprand(G->rstate) % ncount;
        next = CClinkern_flipper_next(F, try1);
        prev = CClinkern_flipper_prev(F, try1);
        nextl = Edgelen (try1, next, D);
        prevl = Edgelen (try1, prev, D);
        if (nextl >= prevl) {
            len = nextl - goodlist[try1][0].weight;
            if (len > best) {
                *t1 = try1;
                *t2 = next;
            }
        } else {
            len = prevl - goodlist[try1][0].weight;
            if (len > best) {
                *t1 = prev;
                *t2 = try1;
            }
        }
    }
}

#define HUNT_PORTION 0.03
#define RAND_TRYS 6
#define GEO_FACTOR 50
#define GEO_MAX 250
#define WALK_STEPS 50

static void find_walk_four(graph *G, distobj *D, CClk_flipper *F, int *t1,
        int *t2, int *t3, int *t4, int *t5, int *t6, int *t7, int *t8)
{
    int s1, s2, s3, s4, s5, s6, s7, s8;
    int old, n, i, j;

    first_kicker(G, D, F, &s1, &s2);

    do {
        old = -1;
        n = s2;

        for (i = 0; i < WALK_STEPS; i++) {
            j = CCutil_lprand(G->rstate) % (G->degree[n]);
            if (old != G->goodlist[n][j].other) {
                old = n;
                n = G->goodlist[n][j].other;
            }
        }
        s3 = n;
        s4 = CClinkern_flipper_next(F, s3);

        n = s4;
        for (i = 0; i < WALK_STEPS; i++) {
            j = CCutil_lprand(G->rstate) % (G->degree[n]);
            if (old != G->goodlist[n][j].other) {
                old = n;
                n = G->goodlist[n][j].other;
            }
        }
        s5 = n;
        s6 = CClinkern_flipper_next(F, s5);

        n = s6;
        for (i = 0; i < WALK_STEPS; i++) {
            j = CCutil_lprand(G->rstate) % (G->degree[n]);
            if (old != G->goodlist[n][j].other) {
                old = n;
                n = G->goodlist[n][j].other;
            }
        }
        s7 = n;
        s8 = CClinkern_flipper_next(F, s7);
    } while (s1 == s3 || s1 == s4 || s1 == s5 || s1 == s6 || s1 == s7
            || s1 == s8 || s2 == s3 || s2 == s4 || s2 == s5 || s2 == s6
            || s2 == s7 || s2 == s8 || s3 == s5 || s3 == s6 || s3 == s7
            || s3 == s8 || s4 == s5 || s4 == s6 || s4 == s7 || s4 == s8
            || s5 == s7 || s5 == s8 || s6 == s7 || s6 == s8);

    *t1 = s1;
    *t2 = s2;
    *t3 = s3;
    *t4 = s4;
    *t5 = s5;
    *t6 = s6;
    *t7 = s7;
    *t8 = s8;
}

static void turn(int n, aqueue *Q, CCptrworld *intptr_world)
{
    add_to_active_queue(n, Q, intptr_world);
}

static void kickturn(int n, aqueue *Q, distobj *D, graph *G, CClk_flipper *F,
        CCptrworld *intptr_world)
{
    add_to_active_queue(n, Q, intptr_world);
    {
        int k;
        k = CClinkern_flipper_next(F, n);
        add_to_active_queue(k, Q, intptr_world);
        k = CClinkern_flipper_next(F, k);
        add_to_active_queue(k, Q, intptr_world);
        k = CClinkern_flipper_prev(F, n);
        add_to_active_queue(k, Q, intptr_world);
        k = CClinkern_flipper_prev(F, k);
        add_to_active_queue(k, Q, intptr_world);
    }
}

static void bigturn(graph *G, int n, int tonext, aqueue *Q, CClk_flipper *F,
        distobj *D, CCptrworld *intptr_world)
{
    int i, k;

    add_to_active_queue(n, Q, intptr_world);
    if (tonext) {
        for (i = 0, k = n; i < MARK_LEVEL; i++) {
            k = CClinkern_flipper_next(F, k);
            add_to_active_queue(k, Q, intptr_world);
        }
    } else {
        for (i = 0, k = n; i < MARK_LEVEL; i++) {
            k = CClinkern_flipper_prev(F, k);
            add_to_active_queue(k, Q, intptr_world);
        }
    }

    for (i = 0; i < G->degree[n]; i++) {
        add_to_active_queue(G->goodlist[n][i].other, Q, intptr_world);
    }
}

static void randcycle(int ncount, int *cyc, CCrandstate *rstate)
{
    int i, k, temp;

    for (i = 0; i < ncount; i++)
        cyc[i] = i;
    for (i = ncount; i > 1; i--) {
        k = CCutil_lprand(rstate) % i;
        CC_SWAP(cyc[i - 1], cyc[k], temp);
    }
}

static void initgraph(graph *G)
{
    G->goodlist = (edge **) NULL;
    G->edgespace = (edge *) NULL;
    G->degree = (int *) NULL;
    G->weirdmark = (int *) NULL;
    G->weirdmagic = 0;
    G->ncount = 0;
}

static void freegraph(graph *G)
{
    if (G) {
        CC_IFFREE(G->goodlist, edge *);
        CC_IFFREE(G->edgespace, edge);
        CC_IFFREE(G->degree, int);
        CC_IFFREE(G->weirdmark, int);
        G->weirdmagic = 0;
        G->ncount = 0;
    }
}

static int buildgraph(graph *G, int ncount, int ecount, int *elist, distobj *D)
{
    int rval = 0;
    int n1, n2, w, i;
    edge *p;

    G->goodlist = CC_SAFE_MALLOC (ncount, edge *);
    G->degree = CC_SAFE_MALLOC (ncount, int);
    G->weirdmark = CC_SAFE_MALLOC (ncount, int);
    G->edgespace = CC_SAFE_MALLOC ((2 * ecount) + ncount, edge);
    if (G->goodlist == (edge **) NULL || G->degree == (int *) NULL
            || G->edgespace == (edge *) NULL) {
        fprintf(stderr, "out of memory in buildgraph\n");
        rval = 1;
        goto CLEANUP;
    }

    for (i = 0; i < ncount; i++) {
        G->degree[i] = 1;
        G->weirdmark[i] = 0;
    }
    for (i = ecount - 1; i >= 0; i--) {
        G->degree[elist[2 * i]]++;
        G->degree[elist[(2 * i) + 1]]++;
    }

    for (i = 0, p = G->edgespace; i < ncount; i++) {
        G->goodlist[i] = p;
        p += (G->degree[i]);
        G->goodlist[i][G->degree[i] - 1].weight = BIGINT;
        G->degree[i] = 0;
    }

    for (i = ecount - 1; i >= 0; i--) {
        n1 = elist[2 * i];
        n2 = elist[(2 * i) + 1];
        w = Edgelen (n1, n2, D);
        insertedge(G, n1, n2, w);
        insertedge(G, n2, n1, w);
    }
    G->ncount = ncount;
    G->weirdmagic = 0;

CLEANUP:

    if (rval)
        freegraph(G);
    return rval;
}

static void insertedge(graph *G, int n1, int n2, int w)
{
    int i;
    edge *e = G->goodlist[n1];

    for (i = G->degree[n1] - 1; i >= 0 && e[i].weight >= w; i--) {
        e[i + 1].weight = e[i].weight;
        e[i + 1].other = e[i].other;
    }
    e[i + 1].weight = w;
    e[i + 1].other = n2;
    G->degree[n1]++;
}

static void linkern_free_world(CCptrworld *intptr_world,
        CCptrworld *edgelook_world)
{
    CCptrworld_delete(intptr_world);
    CCptrworld_delete(edgelook_world);
}

static int init_flipstack(flipstack *f, int total, int single)
{
    f->counter = 0;
    f->max = 0;
    f->stack = (flippair *) NULL;

    f->stack = CC_SAFE_MALLOC (total + single, flippair);
    if (f->stack == (flippair *) NULL) {
        fprintf(stderr, "out of memory in init_flipstack\n");
        return 1;
    }
    f->max = total;

    return 0;
}

static void free_flipstack(flipstack *f)
{
    f->counter = 0;
    f->max = 0;
    CC_IFFREE(f->stack, flippair);
}

static void init_adddel(adddel *E)
{
    E->add_edges = (char *) NULL;
    E->del_edges = (char *) NULL;
}

static void free_adddel(adddel *E)
{
    if (E) {
        CC_IFFREE(E->add_edges, char);
        CC_IFFREE(E->del_edges, char);
    }
}

static int build_adddel(adddel *E, int ncount)
{
    int rval = 0;
    int i, M;

    i = 0;
    while ((1 << i) < ncount)
        i++;
    M = (1 << i);

    E->add_edges = CC_SAFE_MALLOC (M, char);
    E->del_edges = CC_SAFE_MALLOC (M, char);
    if (E->add_edges == (char *) NULL || E->del_edges == (char *) NULL) {
        fprintf(stderr, "out of memory in build_adddel\n");
        rval = 1;
        goto CLEANUP;
    }
    for (i = 0; i < M; i++) {
        E->add_edges[i] = 0;
        E->del_edges[i] = 0;
    }

CLEANUP:

    if (rval) {
        free_adddel(E);
    }
    return rval;
}

static void init_aqueue(aqueue *Q)
{
    Q->active = (char *) NULL;
    Q->active_queue = (intptr *) NULL;
    Q->bottom_active_queue = (intptr *) NULL;
    Q->h = (CCdheap *) NULL;
}

static void free_aqueue(aqueue *Q, CCptrworld *intptr_world)
{
    if (Q) {
        CC_IFFREE(Q->active, char);
        intptr_listfree(intptr_world, Q->active_queue);
        Q->active_queue = (intptr *) NULL;
        Q->bottom_active_queue = (intptr *) NULL;
        if (Q->h) {
            CCutil_dheap_free(Q->h);
            Q->h = (CCdheap *) NULL;
        }
    }
}

static int build_aqueue(aqueue *Q, int ncount, CCptrworld *intptr_world)
{
    int rval = 0;
    int i;

    init_aqueue(Q);

    Q->active = CC_SAFE_MALLOC (ncount, char);
    if (Q->active == (char *) NULL) {
        fprintf(stderr, "out of memory in build_aqueue\n");
        rval = 1;
        goto CLEANUP;
    }
    for (i = 0; i < ncount; i++)
        Q->active[i] = 0;

CLEANUP:

    if (rval) {
        free_aqueue(Q, intptr_world);
    }
    return rval;
}

static void add_to_active_queue(int n, aqueue *Q, CCptrworld *intptr_world)
{
    intptr *ip;

    if (Q->active[n] == 0) {
        Q->active[n] = 1;
        ip = intptralloc(intptr_world);
        ip->this = n;
        ip->next = (intptr *) NULL;
        if (Q->bottom_active_queue) {
            Q->bottom_active_queue->next = ip;
        } else {
            Q->active_queue = ip;
        }
        Q->bottom_active_queue = ip;
    }
}

static int pop_from_active_queue(aqueue *Q, CCptrworld *intptr_world)
{
    intptr *ip;
    int n = -1;

    if (Q->active_queue != (intptr *) NULL) {
        ip = Q->active_queue;
        n = ip->this;
        Q->active_queue = ip->next;
        if (ip == Q->bottom_active_queue) {
            Q->bottom_active_queue = (intptr *) NULL;
        }
        intptrfree(intptr_world, ip);
        Q->active[n] = 0;
    }
    return n;
}

static void init_distobj(distobj *D)
{
    D->dat = (CCdatagroup *) NULL;
    D->cacheind = (int *) NULL;
    D->cacheval = (int *) NULL;
    D->cacheM = 0;
}

static void free_distobj(distobj *D)
{
    if (D) {
        D->dat = (CCdatagroup *) NULL;
        CC_IFFREE(D->cacheind, int);
        CC_IFFREE(D->cacheval, int);
        D->cacheM = 0;
    }
}

static int build_distobj(distobj *D, int ncount, CCdatagroup *dat)
{
    int rval = 0;
    int i;

    init_distobj(D);
    D->dat = dat;

    i = 0;
    while ((1 << i) < (ncount << 2))
        i++;
    D->cacheM = (1 << i);
    D->cacheind = CC_SAFE_MALLOC (D->cacheM, int);
    D->cacheval = CC_SAFE_MALLOC (D->cacheM, int);
    if (D->cacheind == (int *) NULL || D->cacheval == (int *) NULL) {
        fprintf(stderr, "out of memory in build_distobj\n");
        rval = 1;
        goto CLEANUP;
    }
    for (i = 0; i < D->cacheM; i++) {
        D->cacheind[i] = -1;
    }
    D->cacheM--;

CLEANUP:

    if (rval) {
        free_distobj(D);
    }
    return rval;
}

static int CCutil_dat_edgelen(int i, int j, CCdatagroup *dat)
{
    if (dat->ndepot) {
        if (i >= dat->orig_ncount) {
            return dat->depotcost[j];
        } else if (j >= dat->orig_ncount) {
            return dat->depotcost[i];
        }
    }
    return (dat->edgelen)(i, j, dat);
}

static int dist(int i, int j, distobj *D)
{
    int ind;

    if (i > j) {
        int temp;
        CC_SWAP(i, j, temp);
    }
    ind = (((i << 8) + i + j) & (D->cacheM));

    if (D->cacheind[ind] != i) {
        D->cacheind[ind] = i;
        D->cacheval[ind] = CCutil_dat_edgelen(i, j, D->dat);
    }
    return D->cacheval[ind];
}

static int CClinkern_flipper_init(CClk_flipper *F, int ncount, int *cyc)
{
    int i, j, cind, remain;
    int rval = 0;
    CClk_childnode *c, *cprev;
    CClk_parentnode *p;

    init_flipper(F);
    rval = build_flipper(F, ncount);
    if (rval) {
        fprintf(stderr, "build_flipper failed\n");
        goto CLEANUP;
    }

    remain = ncount;
    i = 0;
    j = 2 * F->groupsize;
    while (remain >= j) {
        F->parents[i].size = F->groupsize;
        remain -= F->groupsize;
        i++;
    }
    if (remain > F->groupsize) {
        F->parents[i].size = remain / 2;
        remain -= (remain / 2);
        i++;
    }
    F->parents[i].size = remain;
    i++;

    if (i != F->nsegments) {
        fprintf(stderr, "seg count is wrong\n");
        rval = 1;
        goto CLEANUP;
    }

    c = &(F->children[cyc[ncount - 1]]);
    for (i = 0, p = F->parents, cind = 0; i < F->nsegments; p++, i++) {
        p->id = i;
        p->rev = 0;
        p->ends[0] = &(F->children[cyc[cind]]);
        for (j = p->size; j > 0; j--) {
            cprev = c;
            c = &(F->children[cyc[cind]]);
            c->id = cind;
            c->name = cyc[cind];
            c->parent = p;
            c->adj[0] = cprev;
            cprev->adj[1] = c;
            cind++;
        }
        p->ends[1] = c;
        p->adj[0] = p - 1;
        p->adj[1] = p + 1;
    }
    F->parents[0].adj[0] = &(F->parents[F->nsegments - 1]);
    F->parents[F->nsegments - 1].adj[1] = &(F->parents[0]);

CLEANUP:

    if (rval) {
        free_flipper(F);
    }
    return rval;
}

static void CClinkern_flipper_cycle(CClk_flipper *F, int *x)
{
    CClk_childnode *c, *start;
    int k = 0;

    start = &(F->children[0]);
    c = start->adj[!((F->reversed) ^ (start->parent->rev))];

    x[k++] = start->name;
    while (c != start) {
        x[k++] = c->name;
        c = c->adj[!((F->reversed) ^ (c->parent->rev))];
    }
}

static void CClinkern_flipper_finish(CClk_flipper *F)
{
    free_flipper(F);
}

static int CClinkern_flipper_next(CClk_flipper *F, int x)
{
    return F->children[x].adj[!((F->reversed) ^ (F->children[x].parent->rev))]->name;
}

static int CClinkern_flipper_prev(CClk_flipper *F, int x)
{
    return F->children[x].adj[(F->reversed) ^ (F->children[x].parent->rev)]->name;
}

static void CClinkern_flipper_flip(CClk_flipper *F, int x, int y)
{
    CClk_childnode *xc = &(F->children[x]);
    CClk_childnode *yc = &(F->children[y]);

    if (SAME_SEGMENT (xc, yc)) {
        if (xc != yc) {
            same_segment_flip(F, xc, yc);
        }
    } else {
        int xdir = ((F->reversed) ^ (xc->parent->rev));
        int ydir = ((F->reversed) ^ (yc->parent->rev));
        CClk_childnode *xprev = xc->adj[xdir];
        CClk_childnode *ynext = yc->adj[!ydir];
        if (SAME_SEGMENT (ynext, xprev)) {
            if (ynext != xprev) {
                same_segment_flip(F, ynext, xprev);
            }
            (F->reversed) ^= 1;
        } else {
            int side;
            if (xc->parent->ends[xdir] == xc && yc->parent->ends[!ydir] == yc) {
                if (F->reversed)
                    side = xc->parent->id - yc->parent->id;
                else
                    side = yc->parent->id - xc->parent->id;
                if (side < 0)
                    side += F->nsegments;
                if (side < F->nsegments / 2) {
                    consecutive_segment_flip(F, xc->parent, yc->parent);
                } else {
                    consecutive_segment_flip(F, yc->parent->adj[!F->reversed],
                            xc->parent->adj[F->reversed]);
                    (F->reversed) ^= 1;
                }
            } else {
                if (xprev->parent == xc->parent) {
                    segment_split(F, xc->parent, xprev, xc, 0);
                    if (SAME_SEGMENT (xc, yc)) {
                        if (xc != yc)
                            same_segment_flip(F, xc, yc);
                        return;
                    } else if (SAME_SEGMENT (ynext, xprev)) {
                        if (ynext != xprev) {
                            same_segment_flip(F, ynext, xprev);
                        }
                        (F->reversed) ^= 1;
                        return;
                    }
                }
                if (ynext->parent == yc->parent) {
                    segment_split(F, yc->parent, yc, ynext, 0);
                    if (SAME_SEGMENT (xc, yc)) {
                        if (xc != yc)
                            same_segment_flip(F, xc, yc);
                        return;
                    } else if (SAME_SEGMENT (ynext, xprev)) {
                        if (ynext != xprev) {
                            same_segment_flip(F, ynext, xprev);
                        }
                        (F->reversed) ^= 1;
                        return;
                    }
                }
                if (F->reversed)
                    side = xc->parent->id - yc->parent->id;
                else
                    side = yc->parent->id - xc->parent->id;
                if (side < 0)
                    side += F->nsegments;
                if (side < F->nsegments / 2) {
                    consecutive_segment_flip(F, xc->parent, yc->parent);
                } else {
                    consecutive_segment_flip(F, yc->parent->adj[!F->reversed],
                            xc->parent->adj[F->reversed]);
                    (F->reversed) ^= 1;
                }

            }
        }
    }
}

static void same_segment_flip(CClk_flipper *F, CClk_childnode *a,
        CClk_childnode *b)
{
    CClk_parentnode *parent = a->parent;
    int dir = ((F->reversed) ^ (parent->rev));
    CClk_childnode *aprev = a->adj[dir];
    CClk_childnode *bnext = b->adj[!dir];
    CClk_childnode *c, *cnext;

    if ((dir && a->id - b->id > F->split_cutoff)
            || (!dir && b->id - a->id > F->split_cutoff)) {
        if (aprev->parent == parent)
            segment_split(F, parent, aprev, a, 1);
        if (bnext->parent == parent)
            segment_split(F, parent, b, bnext, 2);
        aprev->adj[!((F->reversed) ^ (aprev->parent->rev))] = b;
        bnext->adj[(F->reversed) ^ (bnext->parent->rev)] = a;
        a->adj[dir] = bnext;
        b->adj[!dir] = aprev;
        parent->rev ^= 1;
        return;
    }

    if (dir) {
        int id = a->id;
        aprev->adj[!((F->reversed) ^ (aprev->parent->rev))] = b;
        bnext->adj[(F->reversed) ^ (bnext->parent->rev)] = a;
        cnext = b->adj[1];
        b->adj[1] = aprev;
        b->adj[0] = cnext;
        b->id = id--;
        c = cnext;
        while (c != a) {
            cnext = c->adj[1];
            c->adj[1] = c->adj[0];
            c->adj[0] = cnext;
            c->id = id--;
            c = cnext;
        }
        a->adj[1] = a->adj[0];
        a->adj[0] = bnext;
        a->id = id;
        if (parent->ends[1] == a)
            parent->ends[1] = b;
        if (parent->ends[0] == b)
            parent->ends[0] = a;
    } else {
        int id = a->id;
        aprev->adj[!((F->reversed) ^ (aprev->parent->rev))] = b;
        bnext->adj[(F->reversed) ^ (bnext->parent->rev)] = a;
        c = b->adj[0];
        b->adj[0] = aprev;
        b->adj[1] = c;
        b->id = id++;
        while (c != a) {
            cnext = c->adj[0];
            c->adj[0] = c->adj[1];
            c->adj[1] = cnext;
            c->id = id++;
            c = cnext;
        }
        a->adj[0] = a->adj[1];
        a->adj[1] = bnext;
        a->id = id;
        if (parent->ends[0] == a)
            parent->ends[0] = b;
        if (parent->ends[1] == b)
            parent->ends[1] = a;
    }
}

static void consecutive_segment_flip(CClk_flipper *F, CClk_parentnode *a,
        CClk_parentnode *b)
{
    CClk_parentnode *aprev = a->adj[F->reversed];
    CClk_parentnode *bnext = b->adj[!F->reversed];
    CClk_parentnode *c, *cnext;
    CClk_childnode *achild = a->ends[(F->reversed) ^ (a->rev)];
    CClk_childnode *bchild = b->ends[!((F->reversed) ^ (b->rev))];
    CClk_childnode *childprev, *childnext;
    int id = a->id;

    if (F->reversed) {
        childprev = achild->adj[!a->rev];
        childnext = bchild->adj[b->rev];
        childprev->adj[childprev->parent->rev] = bchild;
        childnext->adj[!childnext->parent->rev] = achild;
        bchild->adj[b->rev] = childprev;
        achild->adj[!a->rev] = childnext;

        aprev->adj[0] = b;
        bnext->adj[1] = a;
        c = b->adj[1];
        b->adj[1] = aprev;
        b->adj[0] = c;
        b->id = id--;
        b->rev ^= 1;
        while (c != a) {
            cnext = c->adj[1];
            c->adj[1] = c->adj[0];
            c->adj[0] = cnext;
            c->id = id--;
            c->rev ^= 1;
            c = cnext;
        }
        a->adj[1] = a->adj[0];
        a->adj[0] = bnext;
        a->id = id;
        a->rev ^= 1;
    } else {
        childprev = achild->adj[a->rev];
        childnext = bchild->adj[!b->rev];
        childprev->adj[!childprev->parent->rev] = bchild;
        childnext->adj[childnext->parent->rev] = achild;
        bchild->adj[!b->rev] = childprev;
        achild->adj[a->rev] = childnext;

        aprev->adj[1] = b;
        bnext->adj[0] = a;
        c = b->adj[0];
        b->adj[0] = aprev;
        b->adj[1] = c;
        b->id = id++;
        b->rev ^= 1;
        while (c != a) {
            cnext = c->adj[0];
            c->adj[0] = c->adj[1];
            c->adj[1] = cnext;
            c->id = id++;
            c->rev ^= 1;
            c = cnext;
        }
        a->adj[0] = a->adj[1];
        a->adj[1] = bnext;
        a->id = id;
        a->rev ^= 1;
    }
}

static void segment_split(CClk_flipper *F, CClk_parentnode *p,
        CClk_childnode *aprev, CClk_childnode *a, int left_or_right)
{
    int side;
    int dir = ((F->reversed) ^ (p->rev));
    int id;
    CClk_parentnode *pnext;
    CClk_childnode *b, *bnext;

    if (dir)
        side = p->ends[1]->id - aprev->id + 1;
    else
        side = aprev->id - p->ends[0]->id + 1;

    if ((left_or_right == 0 && side <= p->size / 2) || left_or_right == 1) {
        pnext = p->adj[F->reversed];
        pnext->size += side;
        p->size -= side;
        if (pnext->rev == p->rev) {
            b = pnext->ends[!dir];
            id = b->id;
            if (dir) {
                do {
                    b = b->adj[0];
                    b->id = --id;
                    b->parent = pnext;
                } while (b != aprev);
            } else {
                do {
                    b = b->adj[1];
                    b->id = ++id;
                    b->parent = pnext;
                } while (b != aprev);
            }
            pnext->ends[!dir] = aprev;
            p->ends[dir] = a;
        } else {
            b = pnext->ends[dir];
            id = b->id;
            if (!dir) {
                bnext = b->adj[0];
                do {
                    b = bnext;
                    b->id = --id;
                    b->parent = pnext;
                    bnext = b->adj[1];
                    b->adj[1] = b->adj[0];
                    b->adj[0] = bnext;
                } while (b != aprev);
            } else {
                bnext = b->adj[1];
                do {
                    b = bnext;
                    b->id = ++id;
                    b->parent = pnext;
                    bnext = b->adj[0];
                    b->adj[0] = b->adj[1];
                    b->adj[1] = bnext;
                } while (b != aprev);
            }
            pnext->ends[dir] = aprev;
            p->ends[dir] = a;
        }
    } else {
        pnext = p->adj[!F->reversed];
        pnext->size += (p->size - side);
        p->size = side;
        if (pnext->rev == p->rev) {
            b = pnext->ends[dir];
            id = b->id;
            if (dir) {
                do {
                    b = b->adj[1];
                    b->id = ++id;
                    b->parent = pnext;
                } while (b != a);
            } else {
                do {
                    b = b->adj[0];
                    b->id = --id;
                    b->parent = pnext;
                } while (b != a);
            }
            pnext->ends[dir] = a;
            p->ends[!dir] = aprev;
        } else {
            b = pnext->ends[!dir];
            id = b->id;
            if (!dir) {
                bnext = b->adj[1];
                do {
                    b = bnext;
                    b->id = ++id;
                    b->parent = pnext;
                    bnext = b->adj[0];
                    b->adj[0] = b->adj[1];
                    b->adj[1] = bnext;
                } while (b != a);
            } else {
                bnext = b->adj[0];
                do {
                    b = bnext;
                    b->id = --id;
                    b->parent = pnext;
                    bnext = b->adj[1];
                    b->adj[1] = b->adj[0];
                    b->adj[0] = bnext;
                } while (b != a);
            }
            pnext->ends[!dir] = a;
            p->ends[!dir] = aprev;
        }
    }
}

static int CClinkern_flipper_sequence(CClk_flipper *F, int x, int y, int z)
{
    CClk_childnode *a = &(F->children[x]);
    CClk_childnode *b = &(F->children[y]);
    CClk_childnode *c = &(F->children[z]);
    CClk_parentnode *pa = a->parent;
    CClk_parentnode *pb = b->parent;
    CClk_parentnode *pc = c->parent;

    if (pa == pb) {
        if (pa == pc) {
            if ((F->reversed) ^ (pa->rev)) {
                if (a->id >= b->id) {
                    return (b->id >= c->id || c->id >= a->id);
                } else {
                    return (b->id >= c->id && c->id >= a->id);
                }
            } else {
                if (a->id <= b->id) {
                    return (b->id <= c->id || c->id <= a->id);
                } else {
                    return (b->id <= c->id && c->id <= a->id);
                }
            }
        } else {
            if ((F->reversed) ^ (pa->rev)) {
                return (a->id >= b->id);
            } else {
                return (a->id <= b->id);
            }
        }
    } else if (pa == pc) {
        if ((F->reversed) ^ (pa->rev)) {
            return (a->id <= c->id);
        } else {
            return (a->id >= c->id);
        }
    } else if (pb == pc) {
        if ((F->reversed) ^ (pb->rev)) {
            return (b->id >= c->id);
        } else {
            return (b->id <= c->id);
        }
    } else {
        if (F->reversed) {
            if (pa->id >= pb->id) {
                return (pb->id >= pc->id || pc->id >= pa->id);
            } else {
                return (pb->id >= pc->id && pc->id >= pa->id);
            }
        } else {
            if (pa->id <= pb->id) {
                return (pb->id <= pc->id || pc->id <= pa->id);
            } else {
                return (pb->id <= pc->id && pc->id <= pa->id);
            }
        }
    }
}

static void init_flipper(CClk_flipper *Fl)
{
    Fl->parents = (CClk_parentnode *) NULL;
    Fl->children = (CClk_childnode *) NULL;
    Fl->reversed = 0;
    Fl->nsegments = 0;
    Fl->groupsize = 100;
    Fl->split_cutoff = 100;
}

static void free_flipper(CClk_flipper *Fl)
{
    if (Fl) {
        CC_IFFREE(Fl->parents, CClk_parentnode);
        CC_IFFREE(Fl->children, CClk_childnode);
        Fl->reversed = 0;
        Fl->nsegments = 0;
        Fl->groupsize = 0;
        Fl->split_cutoff = 0;
    }
}

static int build_flipper(CClk_flipper *Fl, int ncount)
{
    int rval = 0;

    Fl->reversed = 0;
    Fl->groupsize = (int) (sqrt((double) ncount) * GROUPSIZE_FACTOR);
    Fl->nsegments = (ncount + Fl->groupsize - 1) / Fl->groupsize;
    Fl->split_cutoff = Fl->groupsize * SEGMENT_SPLIT_CUTOFF;

    Fl->parents = CC_SAFE_MALLOC (Fl->nsegments, CClk_parentnode);
    Fl->children = CC_SAFE_MALLOC (ncount + 1, CClk_childnode);
    if (Fl->parents == (CClk_parentnode *) NULL
            || Fl->children == (CClk_childnode *) NULL) {
        fprintf(stderr, "out of memory in build_flipper\n");
        rval = 1;
        goto CLEANUP;
    }

CLEANUP:

    if (rval) {
        free_flipper(Fl);
    }
    return rval;
}

static void CCutil_dheap_free(CCdheap *h)
{
    CC_IFFREE(h->entry, int);
    CC_IFFREE(h->loc, int);
    CC_IFFREE(h->key, double);
}
