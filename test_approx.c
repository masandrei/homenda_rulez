// test_approx_core.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "graph.h"
#include "approx_task.h"  // contains prototypes for the approx functions under test

#ifndef TEST_EPS
#define TEST_EPS 1e-9
#endif

// ---------- helpers ----------

static MultiGraph* make_graph_from_matrix(int V, const int *flat) {
    MultiGraph *g = createGraph(V);
    assert(g && g->adj);
    for (int i = 0; i < V; ++i) {
        for (int j = 0; j < V; ++j) {
            g->adj[i][j] = flat[i*V + j];
        }
    }
    return g;
}

static int matrix_sum(const MultiGraph *g) {
    int s = 0;
    for (int i = 0; i < g->V; ++i)
        for (int j = 0; j < g->V; ++j)
            s += g->adj[i][j];
    return s;
}

static int graph_is_all_zero(const MultiGraph *g) {
    for (int i = 0; i < g->V; ++i)
        for (int j = 0; j < g->V; ++j)
            if (g->adj[i][j] != 0) return 0;
    return 1;
}

static void require_int_equal(const char *msg, int a, int b) {
    if (a != b) {
        fprintf(stderr, "[FAIL] %s: got %d, expected %d\n", msg, a, b);
        exit(1);
    } else {
        // printf("[OK] %s\n", msg);
    }
}

static void require_true(const char *msg, int cond) {
    if (!cond) {
        fprintf(stderr, "[FAIL] %s: expected true\n", msg);
        exit(1);
    }
}

static void require_false(const char *msg, int cond) {
    if (cond) {
        fprintf(stderr, "[FAIL] %s: expected false\n", msg);
        exit(1);
    }
}

static void require_double_close(const char *msg, double a, double b, double eps) {
    if (isnan(a) || isnan(b) || fabs(a - b) > eps) {
        fprintf(stderr, "[FAIL] %s: got %.12g, expected %.12g (eps=%.1e)\n", msg, a, b, eps);
        exit(1);
    }
}

// ---------- tests for findGraphSize ----------

static void test_findGraphSize_basic(void) {
    // 3x3 with edges: sum(adj)=3, so size = V(3) + E(3) = 6
    const int V = 3;
    const int A[] = {
        0,1,0,
        0,0,1,
        1,0,0
    };
    MultiGraph *g = make_graph_from_matrix(V, A);
    int size = findGraphSize(g);
    require_int_equal("findGraphSize basic", size, 3 + 3);
    freeGraph(g);
}

static void test_findGraphSize_multiplicities(void) {
    // multiplicities count
    const int V = 2;
    const int A[] = {
        2,0,
        0,3
    };
    MultiGraph *g = make_graph_from_matrix(V, A);
    // E = 2 + 3 = 5
    require_int_equal("findGraphSize multiplicities", findGraphSize(g), 2 + 5);
    freeGraph(g);
}

// ---------- tests for graphIsSubgraphApprox ----------

static void test_isSubgraphApprox_trivial_empty(void) {
    // g1 empty, g2 arbitrary → true
    const int V = 3;
    const int A1[] = {
        0,0,0,
        0,0,0,
        0,0,0
    };
    const int A2[] = {
        0,0,1,
        0,1,0,
        2,0,0
    };
    MultiGraph *g1 = make_graph_from_matrix(V, A1);
    MultiGraph *g2 = make_graph_from_matrix(V, A2);
    require_true("isSubgraphApprox empty in arbitrary", graphIsSubgraphApprox(g1, g2, 2));
    freeGraph(g1); freeGraph(g2);
}

static void test_isSubgraphApprox_equal_graphs(void) {
    const int V = 3;
    const int A[] = {
        0,1,0,
        0,0,2,
        0,0,0
    };
    MultiGraph *g1 = make_graph_from_matrix(V, A);
    MultiGraph *g2 = make_graph_from_matrix(V, A);
    require_true("isSubgraphApprox equal graphs", graphIsSubgraphApprox(g1, g2, 2));
    freeGraph(g1); freeGraph(g2);
}

static void test_isSubgraphApprox_vertex_count_fail(void) {
    // g1 has more vertices than g2 → false immediately
    const int A1[] = { 0,1,0, 0,0,1, 0,0,0 }; // 3x3
    const int A2[] = { 0,1, 0,0 };           // 2x2
    MultiGraph *g1 = make_graph_from_matrix(3, A1);
    MultiGraph *g2 = make_graph_from_matrix(2, A2);
    require_false("isSubgraphApprox vertex count", graphIsSubgraphApprox(g1, g2, 2));
    freeGraph(g1); freeGraph(g2);
}

static void test_isSubgraphApprox_simple_true(void) {
    // g1: single edge 0->1 ; g2: supergraph with that edge present
    const int A1[] = {
        0,1,
        0,0
    };
    const int A2[] = {
        1,1,
        0,0
    };
    MultiGraph *g1 = make_graph_from_matrix(2, A1);
    MultiGraph *g2 = make_graph_from_matrix(2, A2);
    require_true("isSubgraphApprox simple true", graphIsSubgraphApprox(g1, g2, 2));
    freeGraph(g1); freeGraph(g2);
}

// ---------- tests for minEdgeAdditionsToEmbedApprox ----------

static void test_minAdditions_empty_target(void) {
    // g1 has 3 edges; g2 is empty with same V → deficit equals sum of g1 edges
    const int V = 3;
    const int A1[] = {
        0,1,0,
        0,0,1,
        2,0,0
    }; // total edges = 1 + 1 + 2 = 4
    const int A2[] = {
        0,0,0,
        0,0,0,
        0,0,0
    };
    MultiGraph *g1 = make_graph_from_matrix(V, A1);
    MultiGraph *g2 = make_graph_from_matrix(V, A2);
    int k = minEdgeAdditionsToEmbedApprox(g1, g2, 2);
    require_int_equal("minAdditions empty target equals |E1|", k, 4);
    freeGraph(g1); freeGraph(g2);
}

static void test_minAdditions_partial_overlap(void) {
    // g1 has edges 0->1 and 1->0; g2 has only 0->1, need 1->0 added → deficit = 1
    const int A1[] = {
        0,1,
        1,0
    };
    const int A2[] = {
        0,1,
        0,0
    };
    MultiGraph *g1 = make_graph_from_matrix(2, A1);
    MultiGraph *g2 = make_graph_from_matrix(2, A2);
    int k = minEdgeAdditionsToEmbedApprox(g1, g2, 2);
    require_int_equal("minAdditions partial overlap", k, 1);
    freeGraph(g1); freeGraph(g2);
}

// ---------- tests for getMinimalEdgeAdditionsApprox ----------

static void test_getMinimalAdditions_empty_target_simple(void) {
    // Deterministic mapping on empty g2: 0->0, 1->1 (WL & greedy tie-break)
    // g1 has a single edge 0->1 → deficit graph must have exactly that at (0,1)
    const int A1[] = {
        0,1,
        0,0
    };
    const int A2[] = {
        0,0,
        0,0
    };
    MultiGraph *g1 = make_graph_from_matrix(2, A1);
    MultiGraph *g2 = make_graph_from_matrix(2, A2);
    MultiGraph *R = getMinimalEdgeAdditionsApprox(g1, g2, 2);
    require_true("getMinimalEdgeAdditions returned non-NULL", R != NULL);
    require_int_equal("R same V as g2", R->V, g2->V);
    require_int_equal("R[0][1]==1", R->adj[0][1], 1);
    // all other entries zero
    int sum = matrix_sum(R);
    require_int_equal("R single edge only", sum, 1);
    freeGraph(R); freeGraph(g1); freeGraph(g2);
}

static void test_getMinimalAdditions_nothing_to_add_returns_null(void) {
    // When g2 already contains g1 among mapped vertices, function returns NULL by convention
    const int A[] = {
        0,1,
        0,0
    };
    MultiGraph *g1 = make_graph_from_matrix(2, A);
    MultiGraph *g2 = make_graph_from_matrix(2, A);
    MultiGraph *R = getMinimalEdgeAdditionsApprox(g1, g2, 2);
    require_true("getMinimalEdgeAdditions NULL when nothing to add", R == NULL);
    freeGraph(g1); freeGraph(g2);
}

// ---------- tests for spectral_distance_directed_normalized ----------

static void test_spectral_distance_equal_zero(void) {
    const int A[] = {
        0,1,0,
        0,0,1,
        0,0,0
    };
    MultiGraph *g1 = make_graph_from_matrix(3, A);
    MultiGraph *g2 = make_graph_from_matrix(3, A);
    double d = spectral_distance_directed_normalized(g1, g2);
    require_double_close("spectral distance equal graphs", d, 0.0, TEST_EPS);
    freeGraph(g1); freeGraph(g2);
}

static void test_spectral_distance_nonzero(void) {
    const int A1[] = {
        0,0,
        0,0
    };
    const int A2[] = {
        0,1,
        0,0
    };
    MultiGraph *g1 = make_graph_from_matrix(2, A1);
    MultiGraph *g2 = make_graph_from_matrix(2, A2);
    double d = spectral_distance_directed_normalized(g1, g2);
    require_true("spectral distance > 0 for different graphs", d > 0.0);
    freeGraph(g1); freeGraph(g2);
}

// ---------- main ----------

int main(void) {
    // size
    test_findGraphSize_basic();
    test_findGraphSize_multiplicities();

    // isSubgraphApprox
    test_isSubgraphApprox_trivial_empty();
    test_isSubgraphApprox_equal_graphs();
    test_isSubgraphApprox_vertex_count_fail();
    test_isSubgraphApprox_simple_true();

    // minEdgeAdditionsToEmbedApprox
    test_minAdditions_empty_target();
    test_minAdditions_partial_overlap();

    // getMinimalEdgeAdditionsApprox
    test_getMinimalAdditions_empty_target_simple();
    test_getMinimalAdditions_nothing_to_add_returns_null();

    // spectral distance
    test_spectral_distance_equal_zero();
    test_spectral_distance_nonzero();

    printf("All tests passed.\n");
    return 0;
}
