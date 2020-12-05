//==================================================================================
// Module       : gcblas.hpp
// Date         : 12/3/20 (DLR)
// Description  : Namespace for basic C-BLAS types, defs for CUDA & CPU interfaces.
//                These do not include 'GBLAS' types or calls
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GCBLAS_HPP)
#define _GCBLAS_HPP

#include "gtvector.hpp"

#if defined(USE_CBLAS)
#include "cblas.h"
#endif
#if defined(USE_CUBLAS)
  #include "cuda_runtime_api.h"
  #include "cublas_v2.h"
#endif

namespace GCBLAS
{

enum GBLAS_ORDER     {CblasRowMajor=101, CblasColMajor=102};
enum GBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113};

// Types, defs:
#if defined(USE_CUBLAS)
  typedef cudaStream_t      GCuStream;
  typedef cublasHandle_t    GBlasHandle;
  typedef cublasStatus_t    GBlasStatus;
  typedef cublasOperation_t GBlasOp;
  typedef cudaError_t       GBlasError;
#else
  typedef int               GCuStream;
  typedef int               GBlasHandle;
  typedef int               GBlasStatus;
  typedef int               GBlasOp;
  typedef int               GBlasError;
#endif

struct cuMatBlockDat {
  GSIZET               nbatch;  // total batched ('num elements')
  GTVector      <int>   ibblk;  // for each stream, starting blk index
  GTVector      <int>   ieblk;  // for each stream, ending blk index
  GBlasHandle         hcublas;  // handle for cuBLAS (1-deriv)
  GBlasHandle   hbatch_cublas;  // handle for batched-cuBLAS (2-deriv)
  GTVector<GCuStream> pStream;  // stream pointers
};

// Methods:
void handle_create (GBlasHandle &);
void handle_destroy(GBlasHandle &);
void stream_create (GCuStream &);
void stream_destroy(GCuStream &);

template<typename T>
void gemm(GBlasHandle h, 
          const enum GBLAS_ORDER Order,
          const enum GBLAS_TRANSPOSE TransA,
          const enum GBLAS_TRANSPOSE TransB,
          const int M, const int N, const int K,
          const T alpha, const T  *A, const int lda,
          const T *B, const int ldb, const T beta,
          T *C, const int ldc);

template<typename T>
void batched_gemm( cuMatBlockDat &cudat,
          const enum GBLAS_ORDER Order,
          const enum GBLAS_TRANSPOSE TransA,
          const enum GBLAS_TRANSPOSE TransB,
          const int M, const int N, const int K,
          const T alpha, const T  *A, const int lda,
          const T *B, const int ldb, const T beta,
          T *C, const int ldc);


} // end, namespace GCBLAS

#include "gcblas.ipp"

#endif // !defined(_GCBLAS_HPP)

