#include <stdio.h>
#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include "graph.h"
#include "exact_task.h"
// from fileio.c
MultiGraph* read_file(const char* name);
void write_file(const char* name, MultiGraph* graph);

// Test fixtures
static MultiGraph* g1 = NULL;
static MultiGraph* g2 = NULL;
static MultiGraph* g3 = NULL;

static void setup(void)
{
	g1 = read_file("data/graph1.txt");
	g2 = read_file("data/graph2.txt");
	g3 = read_file("data/graph3.txt");
	assert(g1 && g2 && g3);
}

static void teardown(void)
{
	if (g1) { freeGraph(g1); g1 = NULL; }
	if (g2) { freeGraph(g2); g2 = NULL; }
	if (g3) { freeGraph(g3); g3 = NULL; }
}

// Test: findGraphSize
static void test_findGraphSize(void)
{
	int s1 = findGraphSize(g1);
	int s2 = findGraphSize(g2);
	int s3 = findGraphSize(g3);
	printf("sizes: g1=%d g2=%d g3=%d\n", s1, s2, s3);
	assert(s3 == 3); // 3 vertices + 0 edges
	assert(s1 == 5);
	assert(s2 == 8);
}

// Test: graphIsSubgraph
static void test_graphIsSubgraph(void)
{
	assert(graphIsSubgraph(g1, g2) == true);
	assert(graphIsSubgraph(g1, g3) == false);
	assert(graphIsSubgraph(g3, g1) == true); // empty edges subgraph
}

// Test: getMinimalEdgeAdditions
static void test_getMinimalEdgeAdditions(void)
{
	MultiGraph* addGraph = getMinimalEdgeAdditions(g1, g3);
	// Should return NULL since g1->g3 mapping fails
	
	// Test with g1->g2 which should return NULL (no additions needed)
	MultiGraph* addGraph2 = getMinimalEdgeAdditions(g1, g2);
    assert(graphIsSubgraph(g1, g2));
	assert(addGraph2 == NULL); // No additions needed since g1 is subgraph of g2
	
	// Test with g3->g1 which should work (empty graph into non-empty)
	MultiGraph* addGraph3 = getMinimalEdgeAdditions(g1, g3);
    assert(!graphIsSubgraph(g1, g3));
	assert(addGraph3 != NULL);
	for(int i = 0; i < addGraph3->V; i++)
	{
		for(int j = 0; j < addGraph3->V; j++)
		{
			printf("%d ", addGraph3->adj[i][j]);
		}
		printf("\n");
	}
	// This might return NULL if the algorithm can't find a mapping
}

int main(void)
{
	setup();
	test_findGraphSize();
	test_graphIsSubgraph();
	test_getMinimalEdgeAdditions(); // Skip for now
	teardown();
	printf("All tests passed.\n");
	return 0;
}

