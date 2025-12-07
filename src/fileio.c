#include <stdio.h>
#include "graph.h"

MultiGraph* read_file(const char* name)
{
    FILE* file = fopen(name, "r");
    if (!file) {
        perror("Error opening file");
        return NULL;
    }

    int V;
    if (fscanf(file, "%d", &V) != 1) {
        fprintf(stderr, "Invalid file format.\n");
        fclose(file);
        return NULL;
    }

    MultiGraph* graph = createGraph(V);

    for (int i = 0; i < V; i++) {
        for (int j = 0; j < V; j++) {
            if (fscanf(file, "%d", &graph->adj[i][j]) != 1) {
                fprintf(stderr, "Invalid adjacency matrix format.\n");
                freeGraph(graph);
                fclose(file);
                return NULL;
            }
        }
    }

    fclose(file);
    return graph;
}

void write_file(const char* name, MultiGraph* graph)
{
    FILE* file = fopen(name, "w");
    if (!file) {
        perror("Error opening file");
        return;
    }

    fprintf(file, "%d\n", graph->V);

    for (int i = 0; i < graph->V; i++) {
        for (int j = 0; j < graph->V; j++) {
            fprintf(file, "%d%s", graph->adj[i][j], (j + 1 < graph->V) ? " " : "");
        }
        fprintf(file, "\n");
    }

    fclose(file);
}