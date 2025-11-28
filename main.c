#include <stdio.h>
#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <windows.h>
#include "graph.h"
#include "exact_task.h"
#include "approx_task.h"

MultiGraph* read_file(const char* name);
void write_file(const char* name, MultiGraph* graph);
void clean(MultiGraph* graph);
void test_findGraphSize(MultiGraph* g1, MultiGraph* g2);
void test_graphIsSubgraphExact(MultiGraph* g1, MultiGraph* g2);
void test_getMinimalEdgeAdditionsExact(MultiGraph* g1, MultiGraph* g2);
void test_graphIsSubgraphApprox(MultiGraph* g1, MultiGraph* g2, int wl);
void test_getMinimalEdgeAdditionsApprox(MultiGraph* g1, MultiGraph* g2, int wl);


int main(int argc, char* argv[])
{
    if(argc < 4 || argc > 5)
    {
        fprintf(stderr, "Usage: %s <graph_file1> <graph_file2> <option> [wl_number]\n", argv[0]);
        fprintf(stderr, "Options: exact, approx, exact_min, approx_min, all\n");
        fprintf(stderr, "wl_number: optional for approximate algorithms (default=1)\n");
        return 1;
    }

    const char* file1 = argv[1];
    const char* file2 = argv[2];
    const char* option = argv[3];
    int wl = 1;

    if(argc == 5)
        wl = atoi(argv[4]);

    MultiGraph* g1 = read_file(file1);
    MultiGraph* g2 = read_file(file2);
    if(!g1 || !g2)
    {
        fprintf(stderr, "Error reading graphs.\n");
        return 1;
    }

    int code = -1;
    if(strcmp(option, "exact") == 0) code = 0;
    else if(strcmp(option, "approx") == 0) code = 1;
    else if(strcmp(option, "exact_min_edge_addition") == 0) code = 2;
    else if(strcmp(option, "approx_min_edge_addition") == 0) code = 3;
    else if(strcmp(option, "all") == 0) code = 4;

    switch(code)
    {
        case 0:
            test_graphIsSubgraphExact(g1, g2);
            break;
        case 1:
            test_graphIsSubgraphApprox(g1, g2, wl);
            break;
        case 2:
            test_getMinimalEdgeAdditionsExact(g1, g2);
            break;
        case 3:
            test_getMinimalEdgeAdditionsApprox(g1, g2, wl);
            break;
        case 4:
            test_findGraphSize(g1, g2);
            test_graphIsSubgraphExact(g1, g2);
            test_getMinimalEdgeAdditionsExact(g1, g2);
            test_graphIsSubgraphApprox(g1, g2, wl);
            test_getMinimalEdgeAdditionsApprox(g1, g2, wl);
            break;
        default:
            fprintf(stderr, "Unknown option: %s\n", option);
            clean(g1);
            clean(g2);
            return 1;
    }

    clean(g1);
    clean(g2);
    return 0;
}


void clean(MultiGraph* graph)
{
	if (graph) { freeGraph(graph); graph = NULL; }
}

void test_findGraphSize(MultiGraph* g1, MultiGraph* g2)
{
    int s1 = findGraphSizeExact(g1);
    int s2 = findGraphSizeExact(g2);
    printf("sizes: g1=%d g2=%d\n", s1, s2);
}

void test_graphIsSubgraphExact(MultiGraph* g1, MultiGraph* g2)
{
    bool result1 = graphIsSubgraph(g1, g2);
    bool result2 = graphIsSubgraph(g2, g1);
    printf("g1 is subgraph of g2? %s\n", result1 ? "true" : "false");
    printf("g2 is subgraph of g1? %s\n", result2 ? "true" : "false");
}

void test_getMinimalEdgeAdditionsExact(MultiGraph* g1, MultiGraph* g2)
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

void test_graphIsSubgraphApprox(MultiGraph* g1, MultiGraph* g2, int wl)
{
    bool result1 = graphIsSubgraphApprox(g1, g2, wl);
    bool result2 = graphIsSubgraphApprox(g2, g1, wl);
    printf("g1 is subgraph of g2? %s\n", result1 ? "true" : "false");
    printf("g2 is subgraph of g1? %s\n", result2 ? "true" : "false");
}

void test_getMinimalEdgeAdditionsApprox(MultiGraph* g1, MultiGraph* g2, int wl)
{
    MultiGraph* addGraph =  getMinimalEdgeAdditionsApprox(g1, g2, wl);
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