#include <stddef.h>
#include "graph.h"
#include "stdbool.h"

int minEdgeAdditionsToEmbedApprox(const MultiGraph* g1, const MultiGraph* g2, int wlRounds);

// Convenience wrapper using APPROX_* defaults.
MultiGraph* getMinimalEdgeAdditionsApprox(const MultiGraph* g1, const MultiGraph* g2, int wlRounds);

bool graphIsSubgraphApprox(const MultiGraph* g1, const MultiGraph* g2, int wlRounds);

double spectral_distance_directed_normalized(const MultiGraph* g1, const MultiGraph* g2);

int findGraphSizeApprox(const MultiGraph* graph);
