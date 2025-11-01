#pragma once

typedef struct {
	int V;
	int** adj;
} MultiGraph;

MultiGraph* createGraph(int V);
void addEdge(MultiGraph* graph, int src, int dst);
void freeGraph(MultiGraph* graph);

// Optional: expose minimal additions helpers
int minEdgeAdditionsToEmbed(const MultiGraph* graph1, const MultiGraph* graph2);
MultiGraph* getMinimalEdgeAdditions(const MultiGraph* graph1, const MultiGraph* graph2);
