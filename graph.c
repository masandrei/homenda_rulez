#include <stdio.h>
#include <stdlib.h>
#include "graph.h"
#include "matrix.h"

// Initializes the graph with the given number of vertices.
// In case of failure returns NULL
MultiGraph* createGraph(int V)
{
    MultiGraph* graph = (MultiGraph*)malloc(sizeof(MultiGraph));
    if(graph == NULL)
    {
        perror("Could not create a graph");
        return NULL;
    }
    graph->V = V;
    graph->adj = alloc_matrix_nxn_int(graph->V);
    return graph;
}

// Adds edge to the given graph.
void addEdge(MultiGraph* graph, int src, int dst)
{
    if(src >= 0 && dst >= 0 && src < graph->V && dst < graph->V)
    {
        graph->adj[src][dst]++;
    }
}

// Frees the memory from the graph
void freeGraph(MultiGraph* graph)
{
    free_matrix_nxn((void**)graph->adj, graph->V);
    free(graph);
}