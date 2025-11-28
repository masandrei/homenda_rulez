#include <stdio.h>
#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <windows.h>
#include "graph.h"
#include "exact_task.h"

MultiGraph* read_file(const char* name);
void write_file(const char* name, MultiGraph* graph);

static void clean(MultiGraph* graph)
{
	if (graph) { freeGraph(graph); graph = NULL; }
}

void test_findGraphSize(MultiGraph* g1, MultiGraph* g2)
{
    int s1 = findGraphSize(g1);
    int s2 = findGraphSize(g2);
    printf("sizes: g1=%d g2=%d\n", s1, s2);
}

void test_graphIsSubgraph(MultiGraph* g1, MultiGraph* g2)
{
    bool result1 = graphIsSubgraph(g1, g2);
    bool result2 = graphIsSubgraph(g2, g1);
    printf("g1 is subgraph of g2? %s\n", result1 ? "true" : "false");
    printf("g2 is subgraph of g1? %s\n", result2 ? "true" : "false");
}

void test_getMinimalEdgeAdditions(MultiGraph* g1, MultiGraph* g2)
{
    MultiGraph* addGraph = getMinimalEdgeAdditions(g1, g2);
    if(addGraph)
    {
        printf("Minimal edge additions from g1 to g2:\n");
        for(int i = 0; i < addGraph->V; i++)
        {
            for(int j = 0; j < addGraph->V; j++)
                printf("%d ", addGraph->adj[i][j]);
            printf("\n");
        }
        freeGraph(addGraph);
    }
    else
    {
        printf("No edge additions needed or mapping failed.\n");
    }
}

int main(int argc, char* argv[])
{
    if(argc != 3)
    {
        fprintf(stderr, "Usage: %s <graph_file1> <graph_file2>\n", argv[0]);
        return 1;
    }

    MultiGraph* g1 = read_file(argv[1]);
    MultiGraph* g2 = read_file(argv[2]);
    if(!g1 || !g2)
    {
        fprintf(stderr, "Error reading graphs.\n");
        return 1;
    }

    test_findGraphSize(g1, g2);
    test_graphIsSubgraph(g1, g2);
    test_getMinimalEdgeAdditions(g1, g2);

    clean(g1);
    clean(g2);

    return 0;
}