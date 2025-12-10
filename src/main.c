#include <stdio.h>
#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <windows.h>
#include "graph.h"
#include "exact_task.h"
#include "approx_task.h"

#define MAX_FILE_PATH 256

MultiGraph* read_file(const char* name);
int read_two_graphs_file(const char *filename, MultiGraph **g1, MultiGraph **g2);
void write_file(const char* name, MultiGraph* graph);
void clean(MultiGraph* graph);
void print_graph(MultiGraph* graph);
int read_wl_option();
void test_findGraphSize(MultiGraph* g1, MultiGraph* g2);
void test_graphIsSubgraphExact(MultiGraph* g1, MultiGraph* g2);
void test_getMinimalEdgeAdditionsExact(MultiGraph* g1, MultiGraph* g2);
void test_graphIsSubgraphApprox(MultiGraph* g1, MultiGraph* g2, int wl);
void test_getMinimalEdgeAdditionsApprox(MultiGraph* g1, MultiGraph* g2, int wl);
void run_menu(MultiGraph** g1, MultiGraph** g2, char graphs_file[]);


int main(int argc, char* argv[])
{
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <graphs_file>\n", argv[0]);
        return 1;
    }

    MultiGraph *g1 = NULL, *g2 = NULL;

    char graphs_file[MAX_FILE_PATH];
    strncpy(graphs_file, argv[1], 255);

    if (!read_two_graphs_file(graphs_file, &g1, &g2)) {
        fprintf(stderr, "Error reading graphs.\n");
        return 1;
    }
    print_graph(g1);
    print_graph(g2);

    run_menu(&g1, &g2, graphs_file);

    clean(g1);
    clean(g2);

    return 0;
}



void run_menu(MultiGraph** g1, MultiGraph** g2, char graphs_file[]) {
    char command[64];
    char filename[MAX_FILE_PATH];
    int wl = 1;

    while (1) {
        printf("\n--- Menu ---\n");
        printf("1: Run exact algorithm\n");
        printf("2: Run approximate algorithm\n");
        printf("3: Change graphs filename\n");
        printf("4: Run exact minimal edge addition\n");
        printf("5: Run approximate minimal edge addition\n");
        printf("0: Exit\n");
        printf("Enter option: ");

        int option;

        if(!fgets(command, sizeof(command), stdin)) continue;
        if(sscanf(command, "%d", &option) != 1) continue;

        printf("\n");

        switch (option) {
            case 0:
                return;

            case 1:
                test_graphIsSubgraphExact(*g1, *g2);
                break;

            case 2:
                printf("Enter wl (default 1): ");
                wl = read_wl_option();
                test_graphIsSubgraphApprox(*g1, *g2, wl);
                break;

            case 3:
                printf("Enter new file path: ");
                if(fgets(filename, sizeof(char) * MAX_FILE_PATH, stdin)) {
                    filename[strcspn(filename, "\n")] = 0; 
                    if (filename[0] == '\0') {
                        printf("Filename unchanged.\n");
                        break;  
                    }
                    clean(*g1);
                    clean(*g2);

                    if (!read_two_graphs_file(filename, g1, g2)) {
                        fprintf(stderr, "Error reading graphs.\n");
                        break;
                    }
                }
                break;

            case 4: {
                MultiGraph* result = getMinimalEdgeAdditions(*g1, *g2);
                if (result) 
                {
                    printf("The result can be saved to file or printed otherwise.\n");
                    printf("Save result to file? (y/n): ");
                    char ans[8];
                    if(fgets(ans, sizeof(ans), stdin)) 
                    {
                        if(ans[0]=='y' || ans[0]=='Y') 
                        {
                            printf("Enter filename (default: result_exact.txt): ");
                            if (fgets(filename, sizeof(char) * MAX_FILE_PATH, stdin)) 
                            {
                                if (filename[0]=='\n') strcpy(filename, "result_exact.txt");
                                filename[strcspn(filename, "\n")] = 0;
                            } else 
                            {
                                strcpy(filename, "result_exact.txt");
                            }
                            write_file(filename, result);
                            printf("Saved to %s\n", filename);
                        }
                        else
                        {
                            print_graph(result);
                        }
                    }
                    clean(result);
                } 
                else printf("\nNo additions needed or failed\n");
                break;
            }

            case 5: {
                printf("Enter WL number (default 1): ");
                wl = read_wl_option();
                MultiGraph* result = getMinimalEdgeAdditionsApprox(*g1, *g2, wl);
                if (result) 
                {
                    printf("The result can be saved to file or printed otherwise.\n");
                    printf("Save to file? (y/n): ");
                    char ans[8];
                    if (fgets(ans, sizeof(ans), stdin)) 
                    {
                        if(ans[0]=='y' || ans[0]=='Y') 
                        {
                            printf("Enter filename (default: result_approx.txt): ");
                            if (fgets(filename, sizeof(filename), stdin)) 
                            {
                                if (filename[0]=='\n') strcpy(filename, "result_approx.txt");
                                filename[strcspn(filename, "\n")] = 0;
                            } else 
                            {
                                strcpy(filename, "result_approx.txt");
                            }
                            write_file(filename, result);
                            printf("Saved to %s\n", filename);
                        }
                        else
                        {
                            print_graph(result);
                        }
                    }
                    clean(result);
                } 
                else printf("\nNo additions needed or failed\n");
                break;
            }

            default:
                printf("Unknown option\n");
        }
    }
}

void clean(MultiGraph* graph)
{
	if (graph) { freeGraph(graph); graph = NULL; }
}

void print_graph(MultiGraph* g) 
{
    for (int i = 0; i< g->V; i++)
    {
        for (int j = 0; j < g->V; j++)
        {
            printf("%d ", g->adj[i][j]);
         }
        printf("\n");
    }
}

int read_wl_option()
{
    char line[64];
    int wl = 1;

    if(fgets(line, sizeof(line), stdin)) {
        if(sscanf(line, "%d", &wl) != 1) wl = 1;
    }
    return wl;
}

void test_findGraphSize(MultiGraph* g1, MultiGraph* g2)
{
    int s1 = findGraphSizeExact(g1);
    int s2 = findGraphSizeApprox(g2);
    printf("sizes: g1=%d g2=%d\n", s1, s2);
}

void test_graphIsSubgraphExact(MultiGraph* g1, MultiGraph* g2)
{
    bool result1 = graphIsSubgraph(g1, g2);
    bool result2 = graphIsSubgraph(g2, g1);
    printf("\n--- Exact Algorithm Result ---\n");
    printf("g1 is subgraph of g2? %s\n", result1 ? "true" : "false");
    printf("g2 is subgraph of g1? %s\n", result2 ? "true" : "false");
}

void test_getMinimalEdgeAdditionsExact(MultiGraph* g1, MultiGraph* g2)
{
    MultiGraph* addGraph = getMinimalEdgeAdditions(g1, g2);
    if(addGraph)
    {
        printf("\n--- Minimal exact edge additions from g1 to g2: ---\n");
        print_graph(addGraph);
        clean(addGraph);
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
    printf("\n--- Approx Algorithm Result ---\n");
    printf("g1 is subgraph of g2? %s\n", result1 ? "true" : "false");
    printf("g2 is subgraph of g1? %s\n", result2 ? "true" : "false");
}

void test_getMinimalEdgeAdditionsApprox(MultiGraph* g1, MultiGraph* g2, int wl)
{
    MultiGraph* addGraph =  getMinimalEdgeAdditionsApprox(g1, g2, wl);
    if(addGraph)
    {
        printf("\n--- Minimal approx edge additions from g1 to g2: ---\n");
        print_graph(addGraph);
        clean(addGraph);
    }
    else
    {
        printf("No edge additions needed or mapping failed.\n");
    }
}