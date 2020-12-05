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

#if defined(USE_CBLAS)
	  if ( std::is_same<T,GFLOAT>::value ) {
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

#elif defined(USE_CUBLAS)

  if ( std::is_same<T,GFLOAT>::value ) {
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

#else
  static_assert("Unrecognized BLAS Type");
#endif

} // end, gemm


//**********************************************************************************
//**********************************************************************************
// METHOD : batched_gemm
// DESC   : Args should be the matrix specs for a singe 'element';
//          A, C are of same dimension and layout; B is constant,
//         'small' 1d matrix
// ARGS   : 
// RETURNS: 
//**********************************************************************************
template<typename T>
void batched_gemm(cuMatBlockDat &cudat,
                  const enum CBLAS_ORDER Order,
                  const enum CBLAS_TRANSPOSE TransA,
                  const enum CBLAS_TRANSPOSE TransB,
                  const int M, const int N, const int K,
                  const T alpha, const T  *A, const int lda,
                  const T *B, const int ldb, const T beta,
                  T *C, const int ldc)
{
  GINT   nstreams=cudat.pStream.size();
  GSIZET iastart, iastride, ibstride, szblk;
 
#if defined(USE_CBLAS)

  szblk= M*K;
  for ( auto j=0; j<cudat.nbatch; j++ ) {
    GCBLAS::gemm<T>( cudat.hcublas, Order, TransA, TransB, M, N, K,
                     M, N, K,
                     alpha, A+j*szblk, lda,
                     B, ldb, beta,
                     C+j*szblk, ldc);
  }

#elif defined(USE_CUBLAS)

  for ( auto j=0; j<nstreams; j++ ) {
    cublasSetStream(cudat.hbatch_cublas, cudat.pStream[j]);
  }
  ibstride = 0;
  if      ( std::is_same<T,GFLOAT>::value ) {
    for ( auto j=0; j<nstreams; j++ ) {
      iastart  = ibblk[j]*M*K*sizeof(T);
      iastride = (ieblk[j] - ibblk[j] + 1)*M*K*sizeof(T);
      cublasSgemmStridedBatched( 
                   cudat.hbatch_cublas, 
                   Order, TransA, TransB, M, N, K,
                   M, N, K, 
                   alpha, A+iastart, lda, iastride,
                   B, ldb, ibstride, beta,
                   C+iastart, ldc, iastride);
    }
  }
  else if ( std::is_same<T,GDOUBLE>::value ) {
    for ( auto j=0; j<nstreams; j++ ) {
      iastart  = ibblk[j]*M*K*sizeof(T);
      iastride = (ieblk[j] - ibblk[j] + 1)*M*K*sizeof(T);
      cublasDgemmStridedBatched( 
                   cudat.hbatch_cublas, 
                   Order, TransA, TransB, M, N, K,
                   M, N, K, 
                   alpha, A+iastart, lda, iastride,
                   B, ldb, ibstride, beta,
                   C+iastart, ldc, iastride);
    }
  }
  else {
    assert(FALSE);
  }

#else
  static_assert("Unrecognized Batched BLAS Type");
#endif

} // end, bacthed_gemm


} // end, namespace GCUDA


