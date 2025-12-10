#include <stdio.h>
#include "graph.h"

int read_two_graphs_file(const char *filename, MultiGraph **g1, MultiGraph **g2)
{
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        return 0;
    }

    int V1, V2;

    /* ---- Read first graph ---- */
    if (fscanf(file, "%d", &V1) != 1) {
        fprintf(stderr, "Invalid file format (V1).\n");
        fclose(file);
        return 0;
    }

    *g1 = createGraph(V1);
    if (!*g1) {
        fclose(file);
        return 0;
    }

    for (int i = 0; i < V1; i++) {
        for (int j = 0; j < V1; j++) {
            if (fscanf(file, "%d", &((*g1)->adj[i][j])) != 1) {
                fprintf(stderr, "Invalid adjacency matrix for graph 1.\n");
                freeGraph(*g1);
                fclose(file);
                return 0;
            }
        }
    }

    /* ---- Read second graph ---- */
    if (fscanf(file, "%d", &V2) != 1) {
        fprintf(stderr, "Invalid file format (V2).\n");
        freeGraph(*g1);
        fclose(file);
        return 0;
    }

    *g2 = createGraph(V2);
    if (!*g2) {
        freeGraph(*g1);
        fclose(file);
        return 0;
    }

    for (int i = 0; i < V2; i++) {
        for (int j = 0; j < V2; j++) {
            if (fscanf(file, "%d", &((*g2)->adj[i][j])) != 1) {
                fprintf(stderr, "Invalid adjacency matrix for graph 2.\n");
                freeGraph(*g1);
                freeGraph(*g2);
                fclose(file);
                return 0;
            }
        }
    }

    fclose(file);
    return 1;
}


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