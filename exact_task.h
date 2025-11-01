#include "graph.h"

int findGraphSize(const MultiGraph* graph);
bool graphIsSubgraph(const MultiGraph* graph1, const MultiGraph* graph2);
double spectral_distance_directed_normalized(const MultiGraph *g1, const MultiGraph *g2);