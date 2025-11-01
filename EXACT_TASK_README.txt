Subgraph check completion consists of 3 steps:
1. Preprocessing
    - Find the strongly-connected components of the graph using Tarjan's algorithm
    - Find the size of each component
    - Sort the vertices, firstly, by the sizes of corresponding SCCs, secondly, by the out-degree of the vertex. The results are placed into order array.
2. Backtracking
    - The core of the algorithm is a recursive Depth-First Search (DFS) over the space of all possible injective mappings.
    - The recursion tracks the current depth (which vertex u1 from the order array is being mapped), the current partial mapping (map), and which vertices in G2 are already used (used).
    - If the depth equals the number of vertices in G1, a complete valid mapping has been found. We're done, we can go home with a 5.0 and do nothing :D
    - For the current vertex u1 = order[depth], iterate through every available (unused) vertex v2 in G2 as a candidate for the mapping.
3. Check for consistency of partial mapping
    - Say v_1 is the vertex of graph G_1, we're trying to map, and v_2 is the candidate vertex in graph G_2, then mapping is consistent iff d_out(v1) <= d_out(v2)
    - If the current vertex u_1 and a previously mapped vertex w_1 belong to the same SCC in the source graph G_1, then their images, the candidate v_2 and the 
      already mapped vertex w_2, must also belong to the same SCC in the target graph G_2. If they do not, the mapping is invalid and this branch is pruned.
    - For every already mapped vertex w1 (mapped to w2), the number of edges between u1 and w1 in G1 must not exceed the number of edges between their images v2 and w2 in G2.
        G1.adj[u1][w1] <= G2.adj[v2][w2]
        G1.adj[w1][u1] <= G2.adj[w2][v2]

Minimal additions task is basically an upgraded version of subgraph check. If it produces a graph with no edges, it means that G1 is already a subgraph of G2.
However, subgraph check algorithm stops on success, whereas minimal additions won't stop until all possible mappings have been assessed, to make sure we found
minimal set of edges to add. Algorithm is the following set of steps:
1. Preprocessing
    - Vertices of the source graph (G1) are sorted into an order array (using SCC size and degree heuristics) to optimize the search sequence, similar to the
     standard Subgraph Isomorphism.
2. Backtracking
    - The search follows a Depth-First Search (DFS) pattern over the vertices of G1 (via the order array).
    - It tracks the currentCost (total deficit for the partial map) and a global bestCost (the minimum deficit found so far).
    - If all vertices of G1 are mapped, update bestCost and bestMap if the current total cost is lower.
3. Cost calculation and branching - when mapping any new vertex from G1 calculate the differences in the number of edges between vertex u1 from G1 and vartex w2 from G2.
If at some point current_cost(number of all edges to be added up to this moment) + difference for newly added vertex, if it exceeds the best_cost stop the computation.