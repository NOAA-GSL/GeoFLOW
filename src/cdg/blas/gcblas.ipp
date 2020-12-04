//==================================================================================
// Module       : gcblas.ipp
// Date         : 12/3/20 (DLR)
// Description  : GCBLAS namespace template definitions
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#include "gcblas.hpp"

namespace GCBLAS
{

//**********************************************************************************
//**********************************************************************************
// METHOD : gemm
// DESC   : 
// ARGS   : 
// RETURNS: 
//**********************************************************************************
template<typename T>
void gemm(GBlasHandle h,  
          const enum CBLAS_ORDER Order,
          const enum CBLAS_TRANSPOSE TransA,
          const enum CBLAS_TRANSPOSE TransB,
          const int M, const int N, const int K,
          const T alpha, const T  *A, const int lda,
          const T *B, const int ldb, const T beta,
          T *C, const int ldc)
{
#if defined(_G_USE_CUDA)
  if      ( std::is_same<T,GFLOAT>::value ) {
    cublasSgemm( h, Order, TransA, TransB, M, N, K,
                 M, N, K, 
                 alpha, A, lda,
                 B, ldb, beta,
                 C, ldc);
  }
  else if ( std::is_same<T,GDOUBLE>::value ) {
    cublasDgemm( h, Order, TransA, TransB, M, N, K,
                 M, N, K, 
                 alpha, A, lda,
                 B, ldb, beta,
                 C, ldc);
  }
  else {
    assert(FALSE);
  }
#else
  if      ( std::is_same<T,GFLOAT>::value ) {
    cblas_fgemm( Order, TransA, TransB, M, N, K,
                 M, N, K, 
                 alpha, A, lda,
                 B, ldb, beta,
                 C, ldc);
  }
  else if ( std::is_same<T,GDOUBLE>::value ) {
    cblas_dgemm( Order, TransA, TransB, M, N, K,
                 M, N, K, 
                 alpha, A, lda,
                 B, ldb, beta,
                 C, ldc);
  }
  else {
    assert(FALSE);
  }
#endif

} // end, gemm


//**********************************************************************************
//**********************************************************************************
// METHOD : batched_gemm
// DESC   : 
// ARGS   : 
// RETURNS: 
//**********************************************************************************
template<typename T>
void batched_gemm(GBlasHandle h,
                  const enum CBLAS_ORDER Order,
                  const enum CBLAS_TRANSPOSE TransA,
                  const enum CBLAS_TRANSPOSE TransB,
                  const int M, const int N, const int K,
                  const T alpha, const T  *A, const int lda,
                  const T *B, const int ldb, const T beta,
                  T *C, const int ldc)
{
#if defined(_G_USE_CUDA)
  if      ( std::is_same<T,GFLOAT>::value ) {
    cublasSgemmBatched( h, Order, TransA, TransB, M, N, K,
                        M, N, K, 
                        alpha, A, lda,
                        B, ldb, beta,
                        C, ldc);
  }
  else if ( std::is_same<T,GDOUBLE>::value ) {
    cublasDgemmBatched( h, Order, TransA, TransB, M, N, K,
                        M, N, K, 
                        alpha, A, lda,
                        B, ldb, beta,
                        C, ldc);
  }
  else {
    assert(FALSE);
  }
#else
  if      ( std::is_same<T,GFLOAT>::value ) {
    cblas_fgemm( Order, TransA, TransB, M, N, K,
                 M, N, K, 
                 alpha, A, lda,
                 B, ldb, beta,
                 C, ldc);
  }
  else if ( std::is_same<T,GDOUBLE>::value ) {
    cblas_dgemm( Order, TransA, TransB, M, N, K,
                 M, N, K, 
                 alpha, A, lda,
                 B, ldb, beta,
                 C, ldc);
  }
  else {
    assert(FALSE);
  }
#endif

} // end, bacthed_gemm


} // end, namespace GCUDA


