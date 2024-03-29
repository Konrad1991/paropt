/* ----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * ----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ----------------------------------------------------------------
 * SUNDIALS SYCL memory helper header file.
 * ----------------------------------------------------------------*/

#ifndef _SUNDIALS_SYCLMEMORY_H
#define _SUNDIALS_SYCLMEMORY_H

#include <CL/sycl.hpp>
#include <sundials/sundials_memory.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/* Implementation specific functions */

SUNDIALS_EXPORT
SUNMemoryHelper SUNMemoryHelper_Sycl(SUNContext sunctx);

/* SUNMemoryHelper functions */

SUNDIALS_EXPORT
int SUNMemoryHelper_Alloc_Sycl(SUNMemoryHelper helper, SUNMemory* memptr,
                               size_t mem_size, SUNMemoryType mem_type,
                               void* queue);

SUNDIALS_EXPORT
SUNMemoryHelper SUNMemoryHelper_Clone_Sycl(SUNMemoryHelper helper);

SUNDIALS_EXPORT
int SUNMemoryHelper_Dealloc_Sycl(SUNMemoryHelper helper, SUNMemory mem,
                                 void* queue);

SUNDIALS_EXPORT
int SUNMemoryHelper_Copy_Sycl(SUNMemoryHelper helper, SUNMemory dst,
                              SUNMemory src, size_t memory_size, void* queue);

SUNDIALS_EXPORT
int SUNMemoryHelper_CopyAsync_Sycl(SUNMemoryHelper helper, SUNMemory dst,
                                   SUNMemory src, size_t memory_size,
                                   void* queue);

SUNDIALS_EXPORT
int SUNMemoryHelper_Destroy_Sycl(SUNMemoryHelper helper);

SUNDIALS_EXPORT
int SUNMemoryHelper_GetAllocStats_Sycl(SUNMemoryHelper helper, SUNMemoryType mem_type, unsigned long* num_allocations,
                                       unsigned long* num_deallocations, size_t* bytes_allocated,
                                       size_t* bytes_high_watermark);

#ifdef __cplusplus
}
#endif

#endif
