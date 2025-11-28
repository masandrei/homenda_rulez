

command:
gcc main.c graph.c exact_task.c approx_task.c fileio.c matrix.c -o main


run:
 .\main.exe data\graph1.txt data\graph2.txt

Usage: %s <graph_file1> <graph_file2> <option> [wl_number]
Options: exact, approx, exact_min_edge_addition, approx_min_edge_addition, all
wl_number: optional for approximate algorithms (default=1)
