#pragma once
#include "stddef.h"

// Creates a square matrix of size N with each element having size elem_size
void** alloc_matrix_nxn(int rowCols, size_t elem_size);

void free_matrix_nxn(void** matrixToFree, int size);

#define alloc_matrix_nxn_int(rowCols) (int**)alloc_matrix_nxn(rowCols, sizeof(int))

#define alloc_matrix_nxn_double(rowCols) (double**)alloc_matrix_nxn(rowCols, sizeof(double))