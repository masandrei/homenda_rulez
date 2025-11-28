#include "graph.h"
#include "matrix.h"
#include "stdlib.h"
#include "stdbool.h"
#include "limits.h"
#include "math.h"
#include "string.h"
#include "errno.h"
#include "exact_task.h"

#ifndef SPECTRAL_EPS
#define SPECTRAL_EPS 1e-10
#endif

// ========================= Degrees & ordering helpers =========================

static int computeOutDegree(const MultiGraph* graph, int v)
{
    int degree = 0;
    for (int j = 0; j < graph->V; j++)
        degree += graph->adj[v][j];
    return degree;
}

static void sortVerticesByDegreeDesc(const MultiGraph* graph, int* order)
{
    for (int i = 0; i < graph->V; i++) order[i] = i;
    for (int i = 0; i < graph->V; i++) {
        int best = i, bestDeg = computeOutDegree(graph, order[i]);
        for (int j = i + 1; j < graph->V; j++) {
            int dj = computeOutDegree(graph, order[j]);
            if (dj > bestDeg || (dj == bestDeg && order[j] < order[best])) {
                best = j; bestDeg = dj;
            }
        }
        if (best != i) { int t = order[i]; order[i] = order[best]; order[best] = t; }
    }
}

// ========================= Tarjan SCC (directed) =========================

typedef struct { int *stack; char *onstack; int *disc; int *low; int time; int top; } TarjanWS;

static void tarjan_dfs(const MultiGraph *g, int u, TarjanWS *ws, int *comp_id, int *comp_count) {
    ws->disc[u] = ws->low[u] = ++ws->time;
    ws->stack[++ws->top] = u;
    ws->onstack[u] = 1;

    for (int v = 0; v < g->V; ++v) {
        if (g->adj[u][v] <= 0) continue;
        if (ws->disc[v] == 0) {
            tarjan_dfs(g, v, ws, comp_id, comp_count);
            if (ws->low[v] < ws->low[u]) ws->low[u] = ws->low[v];
        } else if (ws->onstack[v]) {
            if (ws->disc[v] < ws->low[u]) ws->low[u] = ws->disc[v];
        }
    }

    if (ws->low[u] == ws->disc[u]) {
        for (;;) {
            int w = ws->stack[ws->top--];
            ws->onstack[w] = 0;
            comp_id[w] = *comp_count;
            if (w == u) break;
        }
        (*comp_count)++;
    }
}

static int tarjan_scc_labels(const MultiGraph *g, int **out_comp_id) {
    if (!g || g->V <= 0) { *out_comp_id = NULL; return 0; }

    int n = g->V;
    int *comp_id = (int*)malloc(sizeof(int)*n);
    TarjanWS ws = {
        .stack = (int*)malloc(sizeof(int)*n),
        .onstack = (char*)calloc((size_t)n, 1),
        .disc = (int*)calloc((size_t)n, sizeof(int)),
        .low = (int*)calloc((size_t)n, sizeof(int)),
        .time = 0, .top = -1
    };
    if (!comp_id || !ws.stack || !ws.onstack || !ws.disc || !ws.low) {
        free(comp_id); free(ws.stack); free(ws.onstack); free(ws.disc); free(ws.low);
        *out_comp_id = NULL; return 0;
    }

    int comp_count = 0;
    for (int u = 0; u < n; ++u)
        if (ws.disc[u] == 0)
            tarjan_dfs(g, u, &ws, comp_id, &comp_count);

    free(ws.stack); free(ws.onstack); free(ws.disc); free(ws.low);
    *out_comp_id = comp_id;
    return comp_count;
}

static int* scc_sizes_from_labels(int n, const int *cid, int num_comps) {
    int *sz = (int*)calloc((size_t)num_comps, sizeof(int));
    if (!sz) return NULL;
    for (int i = 0; i < n; ++i) sz[cid[i]]++;
    return sz;
}

static void sortVerticesBySCCThenDegreeDesc(const MultiGraph* graph, const int *comp_id, const int *comp_size, int *order) {
    for (int i = 0; i < graph->V; ++i) order[i] = i;
    for (int i = 0; i < graph->V; ++i) {
        int best = i;
        for (int j = i + 1; j < graph->V; ++j) {
            int ci = comp_id[order[best]], cj = comp_id[order[j]];
            int si = comp_size[ci], sj = comp_size[cj];
            if (sj > si) { best = j; continue; }
            if (sj == si) {
                int di = computeOutDegree(graph, order[best]);
                int dj = computeOutDegree(graph, order[j]);
                if (dj > di || (dj == di && order[j] < order[best])) best = j;
            }
        }
        if (best != i) { int t = order[i]; order[i] = order[best]; order[best] = t; }
    }
}

// ========================= Exact subgraph isomorphism (non-induced) =========================

static bool isConsistentPartialMapping_SCC(
    const MultiGraph* g1, const MultiGraph* g2,
    const int* order, int depth, const int* map, int candidateV2,
    const int* scc1, const int* scc2
){
    int u1 = order[depth];

    // Degree monotonicity (safe for exact non-induced check):
    if (computeOutDegree(g1, u1) > computeOutDegree(g2, candidateV2)) return false;

    for (int k = 0; k < depth; ++k) {
        int w1 = order[k];
        int mappedW2 = map[w1];

        // SCC: vertices from same SCC in g1 must land in same SCC in g2
        if (scc1[u1] == scc1[w1] && scc2[candidateV2] != scc2[mappedW2]) return false;

        // Multiplicity must be sufficient (non-induced)
        if (g1->adj[u1][w1] > g2->adj[candidateV2][mappedW2]) return false;
        if (g1->adj[w1][u1] > g2->adj[mappedW2][candidateV2]) return false;
    }
    return true;
}

static bool backtrackSubgraphIso_SCC(
    const MultiGraph* g1, const MultiGraph* g2,
    int* order, int depth, int* map, char* used,
    const int* scc1, const int* scc2
){
    if (depth == g1->V) return true;

    for (int v2 = 0; v2 < g2->V; v2++) {
        if (used[v2]) continue;
        if (!isConsistentPartialMapping_SCC(g1, g2, order, depth, map, v2, scc1, scc2)) continue;

        int u1 = order[depth];
        map[u1] = v2;
        used[v2] = 1;
        if (backtrackSubgraphIso_SCC(g1, g2, order, depth + 1, map, used, scc1, scc2)) return true;
        used[v2] = 0;
        map[u1] = -1;
    }
    return false;
}

bool graphIsSubgraph(const MultiGraph* graph1, const MultiGraph* graph2)
{
    if (!graph1 || !graph2 || graph1->V > graph2->V) return false;

    int *scc1 = NULL, *scc2 = NULL;
    int c1 = tarjan_scc_labels(graph1, &scc1);
    int c2 = tarjan_scc_labels(graph2, &scc2);
    if (c1 <= 0 || c2 <= 0 || !scc1 || !scc2) { free(scc1); free(scc2); return false; }

    int *sz1 = scc_sizes_from_labels(graph1->V, scc1, c1);
    int *sz2 = scc_sizes_from_labels(graph2->V, scc2, c2);
    if (!sz1 || !sz2) { free(scc1); free(scc2); free(sz1); free(sz2); return false; }

    int n1 = graph1->V;
    int* order = (int*)malloc(sizeof(int) * n1);
    int* map   = (int*)malloc(sizeof(int) * n1);
    char* used = (char*)calloc((size_t)graph2->V, 1);
    if (!order || !map || !used) {
        free(order); free(map); free(used);
        free(scc1); free(scc2); free(sz1); free(sz2);
        return false;
    }

    sortVerticesBySCCThenDegreeDesc(graph1, scc1, sz1, order);
    for (int i = 0; i < n1; ++i) map[i] = -1;

    bool ok = backtrackSubgraphIso_SCC(graph1, graph2, order, 0, map, used, scc1, scc2);

    free(order); free(map); free(used);
    free(scc1); free(scc2); free(sz1); free(sz2);
    return ok;
}

// ========================= Min-edge additions (add edges only; non-induced) =========================
//
// IMPORTANT: Do NOT prune using g2's current degrees/SCCs — additions will change them.
// We only count deficits between mapped pairs. Optional induced mode rejects mappings
// where g2 already has extra edges between mapped vertices.

static int computeIncrementalDeficit_minadd(
    const MultiGraph* g1, const MultiGraph* g2,
    int* order, int depth, int* map, int candidateV2
){
    int u1 = order[depth];
    int add = 0;

    for (int k = 0; k < depth; ++k) {
        int w1 = order[k];
        int v2 = candidateV2;
        int w2 = map[w1];

        int g1_fw = g1->adj[u1][w1], g2_fw = g2->adj[v2][w2];
        int g1_bw = g1->adj[w1][u1], g2_bw = g2->adj[w2][v2];
        if (g1_fw > g2_fw) add += (g1_fw - g2_fw);
        if (g1_bw > g2_bw) add += (g1_bw - g2_bw);
    }
    return add;
}

static int optimisticLowerBoundNext_minadd(
    const MultiGraph* g1, const MultiGraph* g2, int* order, int depth, const char* used
){
    (void)g1; (void)g2; (void)order; (void)depth; (void)used;
    return 0; // safe lower bound; degrees/SCCs can increase via additions
}

static void backtrackMinAdditions_core(
    const MultiGraph* g1, const MultiGraph* g2, int* order, int depth,
    int* map, char* used, int currentCost, int* bestCost, int* bestMap
){
    if (currentCost >= *bestCost) return;

    if (depth == g1->V) {
        if (currentCost < *bestCost) {
            *bestCost = currentCost;
            for (int i = 0; i < g1->V; ++i) bestMap[i] = map[i];
        }
        return;
    }

    int lb = optimisticLowerBoundNext_minadd(g1, g2, order, depth, used);
    if (currentCost + lb >= *bestCost) return;

    int u1 = order[depth];
    for (int v2 = 0; v2 < g2->V; ++v2) {
        if (used[v2]) continue; // injective

        int inc = computeIncrementalDeficit_minadd(g1, g2, order, depth, map, v2);
        if (inc == INT_MAX) continue; // impossible under induced semantics

        int nextCost = currentCost + inc;
        if (nextCost >= *bestCost) continue;

        map[u1] = v2;
        used[v2] = 1;
        backtrackMinAdditions_core(g1, g2, order, depth + 1, map, used, nextCost, bestCost, bestMap);
        used[v2] = 0;
        map[u1] = -1;
    }
}

int minEdgeAdditionsToEmbed(const MultiGraph* graph1, const MultiGraph* graph2)
{
    if (!graph1 || !graph2) return -1;
    if (graph1->V > graph2->V) return -1;

    const int n1 = graph1->V;
    int *order = (int*)malloc(sizeof(int)*n1);
    int *map   = (int*)malloc(sizeof(int)*n1);
    int *best  = (int*)malloc(sizeof(int)*n1);
    char *used = (char*)calloc((size_t)graph2->V, 1);
    if (!order || !map || !best || !used) { free(order); free(map); free(best); free(used); return -1; }

    // Ordering: use SCC-of-g1 (best) or degree (fallback)
    int *scc1 = NULL;
    int c1 = tarjan_scc_labels(graph1, &scc1);
    if (c1 > 0 && scc1) {
        int *sz1 = scc_sizes_from_labels(n1, scc1, c1);
        if (sz1) {
            sortVerticesBySCCThenDegreeDesc(graph1, scc1, sz1, order);
            free(sz1);
        } else {
            sortVerticesByDegreeDesc(graph1, order);
        }
    } else {
        sortVerticesByDegreeDesc(graph1, order);
    }
    free(scc1);

    for (int i = 0; i < n1; ++i) map[i] = -1;

    int bestCost = INT_MAX;
    backtrackMinAdditions_core(graph1, graph2, order, 0, map, used, 0, &bestCost, best);

    free(order); free(map); free(best); free(used);
    return bestCost == INT_MAX ? -1 : bestCost;
}

MultiGraph* getMinimalEdgeAdditions(const MultiGraph* graph1, const MultiGraph* graph2)
{
    if (!graph1 || !graph2) return NULL;
    if (graph1->V > graph2->V) return NULL;

    const int n1 = graph1->V;
    int *order = (int*)malloc(sizeof(int)*n1);
    int *map   = (int*)malloc(sizeof(int)*n1);
    int *best  = (int*)malloc(sizeof(int)*n1);
    char *used = (char*)calloc((size_t)graph2->V, 1);
    if (!order || !map || !best || !used) { free(order); free(map); free(best); free(used); return NULL; }

    // Same ordering strategy as above
    int *scc1 = NULL;
    int c1 = tarjan_scc_labels(graph1, &scc1);
    if (c1 > 0 && scc1) {
        int *sz1 = scc_sizes_from_labels(n1, scc1, c1);
        if (sz1) { sortVerticesBySCCThenDegreeDesc(graph1, scc1, sz1, order); free(sz1); }
        else { sortVerticesByDegreeDesc(graph1, order); }
    } else {
        sortVerticesByDegreeDesc(graph1, order);
    }
    free(scc1);

    for (int i = 0; i < n1; ++i) map[i] = -1;

    int bestCost = INT_MAX;
    backtrackMinAdditions_core(graph1, graph2, order, 0, map, used, 0, &bestCost, best);

    if (bestCost == INT_MAX || bestCost == 0) {
        free(order); free(map); free(best); free(used);
        return NULL; // none needed or impossible (induced mode)
    }

    MultiGraph* result = createGraph(graph2->V);
    if (!result) { free(order); free(map); free(best); free(used); return NULL; }

    // Build only the missing edges for the best mapping
    for (int i = 0; i < n1; ++i) {
        int u1 = order[i];
        int v2 = best[u1];
        for (int j = 0; j < n1; ++j) {
            int w1 = order[j];
            int w2 = best[w1];
#ifdef MINADD_INDUCED
            // Any excess would have been pruned already
#endif
            int need_fw = graph1->adj[u1][w1] - graph2->adj[v2][w2];
            int need_bw = graph1->adj[w1][u1] - graph2->adj[w2][v2];
            if (need_fw > 0) result->adj[v2][w2] = need_fw;
            if (need_bw > 0) result->adj[w2][v2] = need_bw;
        }
    }

    free(order); free(map); free(best); free(used);
    return result;
}

// ========================= Spectral distance (safe) =========================

static int cmp_desc_double(const void *a, const void *b) {
    double da = *(const double*)a, db = *(const double*)b;
    return (db > da) - (db < da); // descending
}

static double** build_hermitian_lift(const MultiGraph *g) {
    if (!g || !g->adj || g->V <= 0) return NULL;
    int n = g->V, N = 2 * n;
    double** H = alloc_matrix_nxn_double(N);
    if (!H) return NULL;
    // H = [0 A; A^T 0]
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            double aij = (double)g->adj[i][j];
            H[i][n + j] = aij;
            H[n + j][i] = aij;
        }
    return H;
}

static int jacobi_symmetric_eigvals_inplace(double **A, int N, int max_sweeps) {
    if (!A || N <= 0) return -1;

    for (int sweep = 0; sweep < max_sweeps; ++sweep) {
        int p = 0, q = 1;
        double max_off = 0.0;
        for (int i = 0; i < N; ++i)
            for (int j = i + 1; j < N; ++j) {
                double aij = fabs(A[i][j]);
                if (aij > max_off) { max_off = aij; p = i; q = j; }
            }
        if (max_off < SPECTRAL_EPS) break;

        double app = A[p][p], aqq = A[q][q], apq = A[p][q];
        double theta = (aqq - app) / (2.0 * apq);
        double t = (theta >= 0.0)
            ? 1.0 / (fabs(theta) + sqrt(theta*theta + 1.0))
            : -1.0 / (fabs(theta) + sqrt(theta*theta + 1.0));
        double c = 1.0 / sqrt(1.0 + t*t);
        double s = t * c;

        for (int k = 0; k < N; ++k) {
            if (k == p || k == q) continue;
            double aip = A[p][k], aiq = A[q][k];
            double vip = c * aip - s * aiq;
            double viq = s * aip + c * aiq;
            A[p][k] = vip; A[k][p] = vip;
            A[q][k] = viq; A[k][q] = viq;
        }

        double app_new = c*c*app - 2.0*c*s*apq + s*s*aqq;
        double aqq_new = s*s*app + 2.0*c*s*apq + c*c*aqq;
        A[p][p] = app_new;
        A[q][q] = aqq_new;
        A[p][q] = 0.0;
        A[q][p] = 0.0;
    }
    return 0;
}

static double* positive_eigenvalues_sorted(double **H, int N) {
    if (!H || N <= 0) return NULL;

    double **A = alloc_matrix_nxn_double(N);
    if (!A) return NULL;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            A[i][j] = H[i][j];

    const int max_sweeps = 10 * N * N;
    if (jacobi_symmetric_eigvals_inplace(A, N, max_sweeps) != 0) {
        free_matrix_nxn((void**)A, N);
        return NULL;
    }

    double *diag = (double*)malloc((size_t)N * sizeof(double));
    if (!diag) { free_matrix_nxn((void**)A, N); return NULL; }
    int m = 0;
    for (int i = 0; i < N; ++i) {
        double val = A[i][i];
        if (val >= -SPECTRAL_EPS) diag[m++] = (val < 0.0) ? 0.0 : val;
    }

    qsort(diag, (size_t)m, sizeof(double), cmp_desc_double);
    int expect = N / 2;
    if (m > expect) m = expect;

    double *vals = (double*)malloc((size_t)m * sizeof(double));
    if (!vals) { free(diag); free_matrix_nxn((void**)A, N); return NULL; }
    for (int i = 0; i < m; ++i) vals[i] = diag[i];

    free(diag);
    free_matrix_nxn((void**)A, N);
    return vals;
}

double spectral_distance_directed_normalized(const MultiGraph *g1, const MultiGraph *g2)
{
    if (!g1 || !g2) return -1.0;

    int N1 = 2 * g1->V, N2 = 2 * g2->V;
    double **H1 = build_hermitian_lift(g1);
    double **H2 = build_hermitian_lift(g2);
    if (!H1 || !H2) {
        if (H1) free_matrix_nxn((void**)H1, N1);
        if (H2) free_matrix_nxn((void**)H2, N2);
        return -1.0;
    }

    double *S1 = positive_eigenvalues_sorted(H1, N1);
    if (!S1) { free_matrix_nxn((void**)H1, N1); free_matrix_nxn((void**)H2, N2); return -1.0; }
    double *S2 = positive_eigenvalues_sorted(H2, N2);
    if (!S2) { free_matrix_nxn((void**)H1, N1); free_matrix_nxn((void**)H2, N2); free(S1); return -1.0; }

    double num2 = 0.0, den2 = 0.0;
    int m = (g1->V > g2->V ? g1->V : g2->V);
    for (int i = 0; i < m; ++i) {
        double a = (i < g1->V) ? S1[i] : 0.0;
        double b = (i < g2->V) ? S2[i] : 0.0;
        double d = a - b;
        num2 += d * d;
        den2 += a * a + b * b;
    }

    free_matrix_nxn((void**)H1, N1);
    free_matrix_nxn((void**)H2, N2);
    free(S1);
    free(S2);

    if (den2 <= SPECTRAL_EPS) return 0.0;
    return sqrt(num2 / den2);
}


int findGraphSizeExact(const MultiGraph* graph)
{
    if (!graph) return 0;
    int n = graph->V, e = 0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            e += graph->adj[i][j];
    return n + e;
}