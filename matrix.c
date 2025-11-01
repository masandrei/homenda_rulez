#include "matrix.h"
#include "stdlib.h"
#include <stdio.h>

void** alloc_matrix_nxn(int N, size_t size)
{
    void** matrix = (void**)malloc(N * sizeof(void*));
    if(matrix == NULL)
    {
        perror("Could not allocate matrix");
        return NULL;
    }
    for(int i = 0 ; i < N; i++)
    {
        matrix[i] = (void*)calloc(N, size);
        if(matrix[i] == NULL)
        {
            while(i-- > 0)
            {
                free(matrix[i]);
            }
            perror("Could not allocate rows");
            return NULL;
        }
    }
    return matrix;
}

void free_matrix_nxn(void **matrixToFree, int size)
{
    if(matrixToFree == NULL)
    {
        return;
    }
    for(int i = 0; i < size; i++)
    {
        free(matrixToFree[i]);
    }
    free(matrixToFree);
}