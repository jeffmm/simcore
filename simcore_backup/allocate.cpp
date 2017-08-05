#include "auxiliary.h"
#include "allocate.h"

/* Allocate 1-dimensional array. */
void *allocate_1d_array(size_t n, size_t size)
{
    void *ptr;

    ptr = gcalloc(n, size);

    return ptr;
}

/* Allocate 2-dimensional array. */
void **allocate_2d_array(size_t n1, size_t n2, size_t size)
{
    /* int i; */
    /* void **ptr; */

    /* ptr = gcalloc(n1, sizeof(void *)); */
    /* for (i = 0; i < n1; ++i) */
    /*     ptr[i] = gcalloc(n2, size); */

    /* return ptr; */
    unsigned int i;
    void *data_ptr;
    void **ptr;
    
    data_ptr = gcalloc(n1*n2, size);
    ptr = (void**) gcalloc(n1, sizeof(void*));
    for (i = 0; i < n1; ++i)
        ptr[i] = (char*) data_ptr + i * n2 * size;
    
    return ptr;

}


/* Allocate 3-dimensional array. */
void ***allocate_3d_array(size_t n1, size_t n2, size_t n3, size_t size)
{
    unsigned int i, j;
    void ***ptr;
    void *data_ptr;

    data_ptr = gcalloc(n1*n2*n3, size);

    ptr = (void ***) allocate_2d_array(n1, n2, sizeof(void*));
    for (i = 0; i < n1; ++i) {
        for (j = 0; j < n2; ++j) {
            ptr[i][j] = (char*) data_ptr + (i * n2 * n3   + j * n3 ) * size;
        }
    }
    return ptr;
}

/* Allocate 4-dimensional array. */
void ****allocate_4d_array(size_t n1, size_t n2, size_t n3, size_t n4, size_t size)
{
    unsigned int i, j, k;
    void ****ptr;

    ptr = (void ****) gcalloc(n1, sizeof(void ***));
    for (i = 0; i < n1; ++i) {
        ptr[i] = (void ***) gcalloc(n2, sizeof(void **));
        for (j = 0; j < n2; ++j) {
            ptr[i][j] = (void **) gcalloc(n3, sizeof(void *));
            for (k = 0; k < n3; ++k)
                ptr[i][j][k] = (void *) gcalloc(n4, size);
        }
    }

    return ptr;
}

/* Free 1-dimensional array. */
void free_1d_array(void *ptr)
{
    free(ptr);

    return;
}

/* Free 2-dimensional array. */
void free_2d_array(void **ptr, size_t n1)
{
    /* for (i = 0; i < n1; ++i) */
    /*     free(ptr[i]);; */
    /* free(ptr); */
    free(ptr[0]);
    free(ptr);

    return;
}

/* Free 3-dimensional array. */
void free_3d_array(void ***ptr, size_t n1, size_t n2)
{
    unsigned int i, j;

    for (i = 0; i < n1; ++i) {
        for (j = 0; j < n2; ++j)
            free(ptr[i][j]);
        free(ptr[i]);
    }
    free(ptr);

    return;
}

/* Free 4-dimensional array. */
void free_4d_array(void ****ptr, size_t n1, size_t n2, size_t n3)
{
    unsigned int i, j, k;

    for (i = 0; i < n1; ++i) {
        for (j = 0; j < n2; ++j) {
            for (k = 0; k < n3; ++k)
                free(ptr[i][j][k]);
            free(ptr[i][j]);
        }
        free(ptr[i]);
    }
    free(ptr);

    return;
}

/* Graceful malloc routine. */
void *gmalloc(size_t size)
{
    void *ptr;

    if ((ptr = malloc(size)) == NULL)
        error_exit("Cannot allocate memory in gmalloc");

    return ptr;
}

/* Graceful calloc routine. */
void *gcalloc(size_t n, size_t size)
{
    void *ptr = NULL;

    if ((ptr = calloc(n, size)) == NULL)
        error_exit("Cannot allocate memory in gcalloc");

    return ptr;
}

/* Graceful realloc routine. */
void *grealloc(void *ptr, size_t size)
{
    if (ptr == NULL)
        ptr = malloc(size);
    else
        ptr = realloc(ptr, size);
    if (ptr == NULL)
        error_exit("Cannot allocate memory in grealloc");

    return ptr;
}

/* Graceful realloc 1d array routine */
void *realloc_1d_array(void *ptr, size_t n1, size_t size)
{
    ptr = grealloc(ptr, size * n1);

    return ptr;
}

void **realloc_2d_array(void **ptr, size_t n1_old, size_t n1, size_t n2, size_t size)
{
    size_t i;

    /* If shrinking array, delete extra memory */
    if (n1 < n1_old)
        for (i = n1; i < n1_old; ++i)
            free(ptr[i]);

    ptr = (void **) grealloc(ptr, n1 * sizeof(void *));
    if (ptr == NULL)
        error_exit("Cannot allocate memory in realloc_2d_array\n");

    /* fix me */
    /* Play it safe and realloc everything */
    for (i = n1_old; i < n1; ++i)
        ptr[i] = (void *) gcalloc(n2, size);

    return ptr;
}

void ***realloc_3d_array(void ***ptr, size_t n1_old, size_t n2_old, size_t n1, size_t n2, size_t n3,
                         size_t size)
{
    size_t i, j;

    /* Delete trailing memory if for some reason you're shrinking the total number of n1 objects */
    for (i = n1; i < n1_old; ++i) {
        for (j = 0; j < n2_old; ++j)
            free(ptr[i][j]);
        free(ptr[i]);
    }

    /* Delete all references to memory if n2 is lower than originally */
    for (i = 0; i < n1_old; ++i)
        for (j = n2; j < n2_old; ++j)
            free(ptr[i][j]);

    ptr = (void ***) grealloc(ptr, n1 * sizeof(void **));

    /* Fix me */
    /* Play it safe and realloc all memory to make sure there are no holes. */
    for (i = n1_old; i < n1; ++i) {
        ptr[i] = (void **) gmalloc(n2 * sizeof(void *));
        for (j = 0; j < n2; ++j)
            ptr[i][j] = (void *) gcalloc(n3, size);
    }
    return ptr;
}
