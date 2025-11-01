// approx_core.c — WL(1) (no hashing) + greedy mapping (exact vs additions),
// plus utilities: findGraphSize and spectral_distance_directed_normalized.

#include "graph.h"
#include "matrix.h"
#include "approx_task.h"

#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>
#include <math.h>

// ---------------------------------------------------------------------
// 0) Graph size (exact; V + total directed multiplicities)
// ---------------------------------------------------------------------

int findGraphSize(const MultiGraph* graph)
{
    if (!graph || !graph->adj) return 0;
    int n = graph->V, e = 0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            e += graph->adj[i][j];
    return n + e;
}

// ---------------------------------------------------------------------
// Small helpers
// ---------------------------------------------------------------------

static inline int indeg(const MultiGraph* g, int v){
    int s=0; for(int u=0; u<g->V; ++u) s += g->adj[u][v]; return s;
}
static inline int outdeg(const MultiGraph* g, int u){
    int s=0; for(int v=0; v<g->V; ++v) s += g->adj[u][v]; return s;
}

// ---------------------------------------------------------------------
// 1-WL color refinement WITHOUT hashing (canonical signatures via sorting)
// Directed multigraphs; produces compact integer colors 0..K-1.
// ---------------------------------------------------------------------

typedef struct {
    int *data;   // pointer into a contiguous buffer
    int len;     // length of signature
    int vertex;  // vertex id
} Sig;

static int cmp_int_lex(const int *a, int la, const int *b, int lb) {
    int m = (la < lb) ? la : lb;
    for (int i = 0; i < m; ++i) {
        if (a[i] != b[i]) return (a[i] < b[i]) ? -1 : 1;
    }
    if (la != lb) return (la < lb) ? -1 : 1;
    return 0;
}

static int cmp_sig(const void *A, const void *B) {
    const Sig *a = (const Sig*)A, *b = (const Sig*)B;
    int r = cmp_int_lex(a->data, a->len, b->data, b->len);
    if (r != 0) return r;
    return (a->vertex < b->vertex) ? -1 : (a->vertex > b->vertex);
}

static void wl_color_refinement_nohash(const MultiGraph *g, int rounds, int *out_color) {
    int n = g->V;
    if (n <= 0) return;

    // Round 0: initial signatures (IN, OUT, LOOP) → canonical colors
    Sig *sigs = (Sig*)malloc(sizeof(Sig)*(size_t)n);
    int *buf  = (int*)malloc(sizeof(int)*3*(size_t)n);
    int *cur  = (int*)malloc(sizeof(int)*n);
    if (!sigs || !buf || !cur) {
        for (int i=0;i<n;++i) out_color[i]=0;
        free(sigs); free(buf); free(cur);
        return;
    }
    for (int u = 0; u < n; ++u) {
        buf[3*u+0] = 0; for (int x=0; x<n; ++x) buf[3*u+0] += g->adj[x][u]; // IN
        buf[3*u+1] = 0; for (int v=0; v<n; ++v) buf[3*u+1] += g->adj[u][v]; // OUT
        buf[3*u+2] = g->adj[u][u]; // LOOP
        sigs[u].data = &buf[3*u];
        sigs[u].len = 3;
        sigs[u].vertex = u;
    }
    qsort(sigs, (size_t)n, sizeof(Sig), cmp_sig);
    int color_id = 0;
    cur[sigs[0].vertex] = color_id;
    for (int i = 1; i < n; ++i) {
        if (cmp_int_lex(sigs[i-1].data, sigs[i-1].len, sigs[i].data, sigs[i].len) != 0) ++color_id;
        cur[sigs[i].vertex] = color_id;
    }

    for (int r = 0; r < rounds; ++r) {
        // Signature: [cur[u], -1, (OUT color,count pairs sorted), -2, (IN color,count pairs sorted)]
        int *tmp_color = (int*)malloc(sizeof(int)*(size_t)n);
        int *tmp_count = (int*)malloc(sizeof(int)*(size_t)n);
        int sig_cap = 2 + 4*n;
        int *sig_storage = (int*)malloc(sizeof(int)*(size_t)sig_cap*(size_t)n);
        if (!tmp_color || !tmp_count || !sig_storage) {
            free(tmp_color); free(tmp_count); free(sig_storage);
            break;
        }

        for (int u = 0; u < n; ++u) {
            int *dst = &sig_storage[u*sig_cap];
            int k = 0;

            dst[k++] = cur[u];
            dst[k++] = -1; // OUT marker

            int m = 0;
            // OUT multiset by neighbor color, counting multiplicities
            for (int v=0; v<n; ++v) {
                int mult = g->adj[u][v];
                if (mult <= 0) continue;
                int c = cur[v];
                int pos = -1;
                for (int t=0; t<m; ++t) if (tmp_color[t]==c) { pos=t; break; }
                if (pos<0) { tmp_color[m]=c; tmp_count[m]=mult; ++m; }
                else       { tmp_count[pos]+=mult; }
            }
            // sort by color (simple selection sort for clarity)
            for (int a=0;a<m;++a){
                int best=a;
                for(int b=a+1;b<m;++b) if(tmp_color[b]<tmp_color[best]) best=b;
                if(best!=a){
                    int tc=tmp_color[a]; tmp_color[a]=tmp_color[best]; tmp_color[best]=tc;
                    int tn=tmp_count[a]; tmp_count[a]=tmp_count[best]; tmp_count[best]=tn;
                }
            }
            for (int t=0;t<m;++t){ dst[k++]=tmp_color[t]; dst[k++]=tmp_count[t]; }

            dst[k++] = -2; // IN marker
            m = 0;
            for (int x=0; x<n; ++x) {
                int mult = g->adj[x][u];
                if (mult <= 0) continue;
                int c = cur[x];
                int pos=-1;
                for (int t=0; t<m; ++t) if (tmp_color[t]==c) { pos=t; break; }
                if (pos<0) { tmp_color[m]=c; tmp_count[m]=mult; ++m; }
                else       { tmp_count[pos]+=mult; }
            }
            for (int a=0;a<m;++a){
                int best=a;
                for(int b=a+1;b<m;++b) if(tmp_color[b]<tmp_color[best]) best=b;
                if(best!=a){
                    int tc=tmp_color[a]; tmp_color[a]=tmp_color[best]; tmp_color[best]=tc;
                    int tn=tmp_count[a]; tmp_count[a]=tmp_count[best]; tmp_count[best]=tn;
                }
            }
            for (int t=0;t<m;++t){ dst[k++]=tmp_color[t]; dst[k++]=tmp_count[t]; }

            sigs[u].data = dst;
            sigs[u].len = k;
            sigs[u].vertex = u;
        }

        qsort(sigs, (size_t)n, sizeof(Sig), cmp_sig);
        int *newc = (int*)malloc(sizeof(int)*n);
        if (!newc) { free(tmp_color); free(tmp_count); free(sig_storage); break; }
        int next_color = 0;
        newc[sigs[0].vertex] = next_color;
        for (int i=1;i<n;++i) {
            if (cmp_int_lex(sigs[i-1].data, sigs[i-1].len, sigs[i].data, sigs[i].len) != 0) ++next_color;
            newc[sigs[i].vertex] = next_color;
        }
        memcpy(cur, newc, sizeof(int)*n);
        free(newc);
        free(tmp_color); free(tmp_count); free(sig_storage);
    }

    for (int i=0;i<n;++i) out_color[i]=cur[i];
    free(sigs); free(buf); free(cur);
}

// ---------------------------------------------------------------------
// Ordering for G1: rarer WL color first, then higher total degree
// ---------------------------------------------------------------------

typedef struct { int wl; int cnt; } WLCount;
typedef struct { int u; int wl; int degSum; int rarity; } OrderKey;

static int find_wl_idx(const WLCount* bag, int bagN, int wl){
    for(int i=0;i<bagN;++i) if(bag[i].wl==wl) return i;
    return -1;
}
static int cmp_order_key(const void* A, const void* B){
    const OrderKey *a=(const OrderKey*)A, *b=(const OrderKey*)B;
    if (a->rarity != b->rarity) return a->rarity - b->rarity;      // rarer first
    if (a->degSum != b->degSum) return b->degSum - a->degSum;      // then higher degree
    if (a->wl != b->wl) return (a->wl < b->wl) ? -1 : 1;           // stable
    return a->u - b->u;
}
static void order_vertices_by_rarity_then_degree(const MultiGraph* g1, const int *wl1, int *order){
    int n1 = g1->V;
    WLCount *bag = (WLCount*)malloc(sizeof(WLCount)*(size_t)n1);
    int bagN=0;
    for (int i=0;i<n1;++i){
        int k = find_wl_idx(bag, bagN, wl1[i]);
        if(k<0){ bag[bagN].wl=wl1[i]; bag[bagN].cnt=1; bagN++; }
        else    bag[k].cnt++;
    }
    OrderKey *keys=(OrderKey*)malloc(sizeof(OrderKey)*(size_t)n1);
    for (int i=0;i<n1;++i){
        keys[i].u = i;
        keys[i].wl = wl1[i];
        keys[i].degSum = indeg(g1,i)+outdeg(g1,i);
        int k = find_wl_idx(bag, bagN, wl1[i]);
        keys[i].rarity = (k>=0)? bag[k].cnt : INT_MAX;
    }
    qsort(keys, (size_t)n1, sizeof(OrderKey), cmp_order_key);
    for (int i=0;i<n1;++i) order[i]=keys[i].u;
    free(keys); free(bag);
}

// ---------------------------------------------------------------------
// Consistency & incremental deficit helpers
// ---------------------------------------------------------------------

static int consistent_with_partial(const MultiGraph* g1, const MultiGraph* g2,
                                   int u, int v, const int *map, int n1)
{
    for (int w = 0; w < n1; ++w) {
        int x = map[w];
        if (x < 0) continue; // not mapped yet
        if (g1->adj[u][w] > g2->adj[v][x]) return 0;
        if (g1->adj[w][u] > g2->adj[x][v]) return 0;
    }
    return 1;
}

static int incremental_deficit_if(const MultiGraph* g1, const MultiGraph* g2,
                                  int u, int v, const int *map, int n1)
{
    int add = 0;
    for (int w = 0; w < n1; ++w) {
        int x = map[w];
        if (x < 0) continue;
        int need_fw = g1->adj[u][w] - g2->adj[v][x];
        int need_bw = g1->adj[w][u] - g2->adj[x][v];
        if (need_fw > 0) add += need_fw;
        if (need_bw > 0) add += need_bw;
    }
    return add;
}

// ---------------------------------------------------------------------
// Two greedy matchers:
//  - exact: for graphIsSubgraphApprox (enforce degree feasibility + partial consistency)
//  - additions: for minEdgeAdditionsToEmbedApprox / getMinimalEdgeAdditionsApprox
//               (no feasibility filters; always assign; choose min incremental deficit)
// ---------------------------------------------------------------------

static void greedy_wl_match_exact(const MultiGraph* g1, const MultiGraph* g2,
                                  const int *wl1, const int *wl2,
                                  const int *order, int *map_out)
{
    int n1 = g1->V, n2 = g2->V;
    char *used = (char*)calloc((size_t)n2, 1);
    for (int i = 0; i < n1; ++i) map_out[i] = -1;

    for (int t = 0; t < n1; ++t) {
        int u = order[t];
        int want = wl1[u];

        int in_u = 0, out_u = 0;
        for (int x = 0; x < n1; ++x) in_u  += g1->adj[x][u];
        for (int y = 0; y < n1; ++y) out_u += g1->adj[u][y];

        int best_v = -1;
        double best_d = INFINITY;

        // prefer same-color, enforce degree-feasibility and partial-consistency
        for (int v = 0; v < n2; ++v) {
            if (used[v]) continue;
            if (wl2[v] != want) continue;

            int in_v = 0, out_v = 0;
            for (int x = 0; x < n2; ++x) in_v  += g2->adj[x][v];
            for (int y = 0; y < n2; ++y) out_v += g2->adj[v][y];

            if (out_v < out_u || in_v < in_u) continue;
            if (!consistent_with_partial(g1, g2, u, v, map_out, n1)) continue;

            double d = fabs((double)out_u - (double)out_v)
                     + fabs((double) in_u - (double) in_v)
                     + 0.3 * fabs((double)g1->adj[u][u] - (double)g2->adj[v][v]);
            if (d < best_d) { best_d = d; best_v = v; }
        }

        // fallback: allow different color if still feasible/consistent
        if (best_v < 0) {
            for (int v = 0; v < n2; ++v) {
                if (used[v]) continue;

                int in_v = 0, out_v = 0;
                for (int x = 0; x < n2; ++x) in_v  += g2->adj[x][v];
                for (int y = 0; y < n2; ++y) out_v += g2->adj[v][y];

                if (out_v < out_u || in_v < in_u) continue;
                if (!consistent_with_partial(g1, g2, u, v, map_out, n1)) continue;

                double d = fabs((double)out_u - (double)out_v)
                         + fabs((double) in_u - (double) in_v)
                         + 0.3 * fabs((double)g1->adj[u][u] - (double)g2->adj[v][v]);
                if (d < best_d) { best_d = d; best_v = v; }
            }
        }

        if (best_v >= 0) { map_out[u] = best_v; used[best_v] = 1; }
        // else remain -1 (exact verification will fail → safe)
    }
    free(used);
}

static void greedy_wl_match_additions(const MultiGraph* g1, const MultiGraph* g2,
                                      const int *wl1, const int *wl2,
                                      const int *order, int *map_out)
{
    int n1 = g1->V, n2 = g2->V;
    char *used = (char*)calloc((size_t)n2, 1);
    for (int i = 0; i < n1; ++i) map_out[i] = -1;

    for (int t = 0; t < n1; ++t) {
        int u = order[t];
        int want = wl1[u];

        int best_v = -1;
        int best_inc = INT_MAX;

        // try same-color unused targets first; choose min incremental deficit
        for (int v = 0; v < n2; ++v) {
            if (used[v]) continue;
            if (wl2[v] != want) continue;
            int inc = incremental_deficit_if(g1, g2, u, v, map_out, n1);
            if (inc < best_inc) { best_inc = inc; best_v = v; }
        }
        // fallback: any unused target, pick min incremental deficit
        if (best_v < 0) {
            for (int v = 0; v < n2; ++v) {
                if (used[v]) continue;
                int inc = incremental_deficit_if(g1, g2, u, v, map_out, n1);
                if (inc < best_inc) { best_inc = inc; best_v = v; }
            }
        }
        // must always assign (since |V1| <= |V2|)
        if (best_v < 0) {
            for (int v = 0; v < n2; ++v) if (!used[v]) { best_v = v; break; }
        }
        map_out[u] = best_v;
        if (best_v >= 0) used[best_v] = 1;
    }
    free(used);
}

// ---------------------------------------------------------------------
// Deficit cost and deficit graph
// ---------------------------------------------------------------------

static long long total_deficit_cost(const MultiGraph* g1, const MultiGraph* g2, const int *map){
    int n1=g1->V;
    long long cost=0;
    for (int i=0;i<n1;++i){
        int v = map[i]; if(v<0) continue;
        for (int j=0;j<n1;++j){
            int w = map[j]; if(w<0) continue;
            int need = g1->adj[i][j] - g2->adj[v][w];
            if (need > 0) cost += need;
        }
    }
    return cost;
}

static MultiGraph* build_deficit_graph(const MultiGraph* g1, const MultiGraph* g2, const int *map){
    MultiGraph *R = createGraph(g2->V);
    if (!R) return NULL;
    int n1=g1->V;
    for (int i=0;i<n1;++i){
        int v = map[i]; if(v<0) continue;
        for (int j=0;j<n1;++j){
            int w = map[j]; if(w<0) continue;
            int need = g1->adj[i][j] - g2->adj[v][w];
            if (need > 0) R->adj[v][w] += need;
        }
    }
    return R;
}

// ---------------------------------------------------------------------
// Public approx API
// ---------------------------------------------------------------------

// Returns true only if the WL-guided mapping (with feasibility+consistency) verifies exactly.
// May be false-negative; never false-positive. Non-induced semantics (extras in g2 allowed).
bool graphIsSubgraphApprox(const MultiGraph* g1, const MultiGraph* g2, int wl_rounds)
{
    if (!g1 || !g2 || g1->V > g2->V) return false;
    if (wl_rounds < 0) wl_rounds = 0;

    int n1=g1->V;

    int *wl1=(int*)malloc(sizeof(int)*n1);
    int *wl2=(int*)malloc(sizeof(int)*g2->V);
    int *order=(int*)malloc(sizeof(int)*n1);
    int *map  =(int*)malloc(sizeof(int)*n1);
    if (!wl1 || !wl2 || !order || !map){ free(wl1); free(wl2); free(order); free(map); return false; }

    wl_color_refinement_nohash(g1, wl_rounds, wl1);
    wl_color_refinement_nohash(g2, wl_rounds, wl2);
    order_vertices_by_rarity_then_degree(g1, wl1, order);
    greedy_wl_match_exact(g1, g2, wl1, wl2, order, map);

    // verify exactly
    bool ok = true;
    for (int i=0;i<n1 && ok;++i){
        int v = map[i]; if (v<0) { ok=false; break; }
        for (int j=0;j<n1;++j){
            int w = map[j]; if (w<0) { ok=false; break; }
            if (g1->adj[i][j] > g2->adj[v][w]) { ok=false; break; }
        }
    }

    free(wl1); free(wl2); free(order); free(map);
    return ok;
}

// Approximate deficit (minimum edges to add) under additions-matcher mapping.
// Returns -1 on failure.
int minEdgeAdditionsToEmbedApprox(const MultiGraph* g1, const MultiGraph* g2, int wl_rounds)
{
    if (!g1 || !g2) return -1;
    if (g1->V > g2->V) return -1;
    if (wl_rounds < 0) wl_rounds = 0;

    int n1=g1->V;

    int *wl1=(int*)malloc(sizeof(int)*n1);
    int *wl2=(int*)malloc(sizeof(int)*g2->V);
    int *order=(int*)malloc(sizeof(int)*n1);
    int *map  =(int*)malloc(sizeof(int)*n1);
    if (!wl1 || !wl2 || !order || !map){ free(wl1); free(wl2); free(order); free(map); return -1; }

    wl_color_refinement_nohash(g1, wl_rounds, wl1);
    wl_color_refinement_nohash(g2, wl_rounds, wl2);
    order_vertices_by_rarity_then_degree(g1, wl1, order);
    greedy_wl_match_additions(g1, g2, wl1, wl2, order, map);

    long long cost = total_deficit_cost(g1, g2, map);

    free(wl1); free(wl2); free(order); free(map);
    if (cost > INT_MAX) return INT_MAX;
    return (int)cost;
}

// Returns a graph with only the edges to add to g2 (under additions-matcher mapping).
// Returns NULL if nothing to add or on failure.
MultiGraph* getMinimalEdgeAdditionsApprox(const MultiGraph* g1, const MultiGraph* g2, int wl_rounds)
{
    if (!g1 || !g2) return NULL;
    if (g1->V > g2->V) return NULL;
    if (wl_rounds < 0) wl_rounds = 0;

    int n1=g1->V;

    int *wl1=(int*)malloc(sizeof(int)*n1);
    int *wl2=(int*)malloc(sizeof(int)*g2->V);
    int *order=(int*)malloc(sizeof(int)*n1);
    int *map  =(int*)malloc(sizeof(int)*n1);
    if (!wl1 || !wl2 || !order || !map){ free(wl1); free(wl2); free(order); free(map); return NULL; }

    wl_color_refinement_nohash(g1, wl_rounds, wl1);
    wl_color_refinement_nohash(g2, wl_rounds, wl2);
    order_vertices_by_rarity_then_degree(g1, wl1, order);
    greedy_wl_match_additions(g1, g2, wl1, wl2, order, map);

    MultiGraph *R = build_deficit_graph(g1, g2, map);

    // return NULL if R has no edges to add
    if (R) {
        int all_zero = 1;
        for (int i=0;i<g2->V && all_zero;++i)
            for (int j=0;j<g2->V;++j)
                if (R->adj[i][j] != 0) { all_zero = 0; break; }
        if (all_zero) { /* freeGraph(R); */ R = NULL; }
    }

    free(wl1); free(wl2); free(order); free(map);
    return R;
}

// // ---------------------------------------------------------------------
// // Spectral distance (directed) using symmetric lift H = [0 A; A^T 0]
// // ---------------------------------------------------------------------

// #ifndef SPECTRAL_EPS
// #define SPECTRAL_EPS 1e-10
// #endif

// static int cmp_desc_double(const void *a, const void *b) {
//     double da = *(const double*)a, db = *(const double*)b;
//     return (db > da) - (db < da);
// }

// static double** build_hermitian_lift(const MultiGraph *g) {
//     if (!g || !g->adj || g->V <= 0) return NULL;
//     int n = g->V, N = 2 * n;
//     double** H = alloc_matrix_nxn_double(N);
//     if (!H) return NULL;
//     for (int i = 0; i < N; ++i)
//         for (int j = 0; j < N; ++j)
//             H[i][j] = 0.0;
//     for (int i = 0; i < n; ++i)
//         for (int j = 0; j < n; ++j) {
//             double aij = (double)g->adj[i][j];
//             H[i][n + j] = aij;
//             H[n + j][i] = aij;
//         }
//     return H;
// }

// static int jacobi_symmetric_eigvals_inplace(double **A, int N, int max_sweeps) {
//     if (!A || N <= 0) return -1;
//     for (int sweep = 0; sweep < max_sweeps; ++sweep) {
//         int p = 0, q = 1;
//         double max_off = 0.0;
//         for (int i = 0; i < N; ++i)
//             for (int j = i + 1; j < N; ++j) {
//                 double aij = fabs(A[i][j]);
//                 if (aij > max_off) { max_off = aij; p = i; q = j; }
//             }
//         if (max_off < SPECTRAL_EPS) break;

//         double app = A[p][p], aqq = A[q][q], apq = A[p][q];
//         double theta = (aqq - app) / (2.0 * apq);
//         double t = (theta >= 0.0)
//             ? 1.0 / (fabs(theta) + sqrt(theta*theta + 1.0))
//             : -1.0 / (fabs(theta) + sqrt(theta*theta + 1.0));
//         double c = 1.0 / sqrt(1.0 + t*t);
//         double s = t * c;

//         for (int k = 0; k < N; ++k) {
//             if (k == p || k == q) continue;
//             double aip = A[p][k], aiq = A[q][k];
//             double vip = c * aip - s * aiq;
//             double viq = s * aip + c * aiq;
//             A[p][k] = vip; A[k][p] = vip;
//             A[q][k] = viq; A[k][q] = viq;
//         }
//         double app_new = c*c*app - 2.0*c*s*apq + s*s*aqq;
//         double aqq_new = s*s*app + 2.0*c*s*apq + c*c*aqq;
//         A[p][p] = app_new;
//         A[q][q] = aqq_new;
//         A[p][q] = 0.0;
//         A[q][p] = 0.0;
//     }
//     return 0;
// }

// static double* positive_eigenvalues_sorted(double **H, int N) {
//     if (!H || N <= 0) return NULL;
//     double **A = alloc_matrix_nxn_double(N);
//     if (!A) return NULL;
//     for (int i=0;i<N;++i)
//         for (int j=0;j<N;++j)
//             A[i][j] = H[i][j];

//     const int max_sweeps = 10 * N * N;
//     if (jacobi_symmetric_eigvals_inplace(A, N, max_sweeps) != 0) {
//         free_matrix_nxn((void**)A, N);
//         return NULL;
//     }

//     double *diag = (double*)malloc(sizeof(double)*(size_t)N);
//     if (!diag) { free_matrix_nxn((void**)A, N); return NULL; }
//     int m = 0;
//     for (int i=0;i<N;++i){
//         double val = A[i][i];
//         if (val >= -SPECTRAL_EPS) diag[m++] = (val < 0.0) ? 0.0 : val;
//     }
//     qsort(diag, (size_t)m, sizeof(double), cmp_desc_double);

//     int expect = N/2;               // positive spectrum size of H
//     if (m > expect) m = expect;

//     double *vals = (double*)malloc(sizeof(double)*(size_t)m);
//     if (!vals) { free(diag); free_matrix_nxn((void**)A, N); return NULL; }
//     for (int i=0;i<m;++i) vals[i]=diag[i];

//     free(diag);
//     free_matrix_nxn((void**)A, N);
//     return vals;
// }

// double spectral_distance_directed_normalized(const MultiGraph *g1, const MultiGraph *g2)
// {
//     if (!g1 || !g2) return -1.0;

//     int N1 = 2 * g1->V;
//     int N2 = 2 * g2->V;
//     double **H1 = build_hermitian_lift(g1);
//     double **H2 = build_hermitian_lift(g2);
//     if (!H1 || !H2) {
//         if (H1) free_matrix_nxn((void**)H1, N1);
//         if (H2) free_matrix_nxn((void**)H2, N2);
//         return -1.0;
//     }

//     double *S1 = positive_eigenvalues_sorted(H1, N1);
//     double *S2 = positive_eigenvalues_sorted(H2, N2);
//     if (!S1 || !S2) {
//         if (S1) free(S1);
//         if (S2) free(S2);
//         free_matrix_nxn((void**)H1, N1);
//         free_matrix_nxn((void**)H2, N2);
//         return -1.0;
//     }

//     int m1 = g1->V, m2 = g2->V;
//     int m  = (m1 > m2 ? m1 : m2);
//     double num2 = 0.0, den2 = 0.0;
//     for (int i = 0; i < m; ++i) {
//         double a = (i < m1) ? S1[i] : 0.0;
//         double b = (i < m2) ? S2[i] : 0.0;
//         double d = a - b;
//         num2 += d*d;
//         den2 += a*a + b*b;
//     }

//     free(S1); free(S2);
//     free_matrix_nxn((void**)H1, N1);
//     free_matrix_nxn((void**)H2, N2);

//     if (den2 <= SPECTRAL_EPS) return 0.0;
//     return sqrt(num2 / den2);
// }
