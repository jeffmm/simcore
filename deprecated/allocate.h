#ifndef _ALLOCATE_H_
#define _ALLOCATE_H_

#include <stddef.h>

void *allocate_1d_array(size_t n, size_t size);

void **allocate_2d_array(size_t n1, size_t n2, size_t size);

void ***allocate_3d_array(size_t n1, size_t n2, size_t n3, size_t size);

void ****allocate_4d_array(size_t n1, size_t n2, size_t n3, size_t n4, size_t size);

void free_1d_array(void *ptr);

void free_2d_array(void **ptr, size_t n1);

void free_3d_array(void ***ptr, size_t n1, size_t n2); 

void free_4d_array(void ****ptr, size_t n1, size_t n2, size_t n3);

void *gmalloc(size_t size);

void *gcalloc(size_t n, size_t size);

void *grealloc(void *ptr, size_t size);

void *realloc_1d_array(void *ptr, size_t n1, size_t size);

void **realloc_2d_array(void **ptr, size_t n1_old, size_t n1, size_t n2, size_t size);

void ***realloc_3d_array(void ***ptr, size_t n1_old, size_t n2_old, size_t n1, size_t n2, size_t n3,
                         size_t size);

#endif
