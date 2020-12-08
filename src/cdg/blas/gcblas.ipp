//==================================================================================
// Module       : gcblas.ipp
// Date         : 12/3/20 (DLR)
// Description  : GCBLAS namespace template definitions
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================

#include "tbox/pio.hpp"
#include "tbox/tracer.hpp"

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
          const GBLAS_ORDER Order,
          const GBLAS_TRANSPOSE TransA,
          const GBLAS_TRANSPOSE TransB,
          const int M, const int N, const int K,
          const T alpha, const T  *A, const int lda,
          const T *B, const int ldb, const T beta,
          T *C, const int ldc)
{
	GEOFLOW_TRACE();

#if defined(USE_CBLAS)
	  if constexpr ( std::is_same<T,float>::value ) {
		GEOFLOW_TRACE_MSG("cblas_sgemm(...)");
	    cblas_sgemm( (CBLAS_ORDER)Order, (CBLAS_TRANSPOSE)TransA, (CBLAS_TRANSPOSE)TransB,
	                 M, N, K,
	                 alpha, A, lda,
	                 B, ldb, beta,
	                 C, ldc);
	  }
	  else if constexpr ( std::is_same<T,double>::value ) {
		GEOFLOW_TRACE_MSG("cblas_dgemm(...)");
        cblas_dgemm( (CBLAS_ORDER)Order, (CBLAS_TRANSPOSE)TransA, (CBLAS_TRANSPOSE)TransB,
        		     M, N, K,
					 alpha, A, lda,
					 B, ldb, beta,
					 C, ldc);
	  }

#elif defined(USE_CUBLAS)

  static const GBlasOp cuBlasOps[] = { CUBLAS_OP_N, CUBLAS_OP_T, CUBLAS_OP_C };

  if constexpr ( std::is_same<T,GFLOAT>::value ) {
	GEOFLOW_TRACE_MSG("cublasSgemm(...)");
    cublasSgemm( h, cuBlasOps[TransA-CblasNoTrans], cuBlasOps[TransB-CblasNoTrans],
                 M, N, K, 
                 (const float*)(&alpha), (const float*)A, lda,
                 (const float*)B, ldb, (const float*)(&beta),
                 (float*)C, ldc);
  }
  else if constexpr ( std::is_same<T,GDOUBLE>::value ) {
	GEOFLOW_TRACE_MSG("cublasDgemm(...)");
    cublasDgemm( h, cuBlasOps[TransA-CblasNoTrans], cuBlasOps[TransB-CblasNoTrans],
                 M, N, K, 
                 (const double*)(&alpha), (const double*)A, lda,
                 (const double*)B, ldb, (const double*)(&beta),
                 (double*)C, ldc);
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
                  const GBLAS_ORDER Order,
                  const GBLAS_TRANSPOSE TransA,
                  const GBLAS_TRANSPOSE TransB,
                  const int M, const int N, const int K,
                  const T alpha, const T  *A, const int lda,
                  const T *B, const int ldb, const T beta,
                  T *C, const int ldc)
{
  GEOFLOW_TRACE();
  GINT   nstreams=cudat.pStream.size();
  GSIZET iastart, iastride, ibstride, szblk;
 
#if defined(USE_CBLAS)

  szblk= M*K;
  for ( auto j=0; j<cudat.nbatch; j++ ) {
    GCBLAS::gemm<T>( cudat.hcublas, Order, TransA, TransB,
                     M, N, K,
                     alpha, A+j*szblk, lda,
                     B, ldb, beta,
                     C+j*szblk, ldc);

  }

#elif defined(USE_CUBLAS)
  static const GBlasOp cuBlasOps[] = { CUBLAS_OP_N, CUBLAS_OP_T, CUBLAS_OP_C };

  for ( auto j=0; j<nstreams; j++ ) {
    cublasSetStream(cudat.hbatch_cublas, cudat.pStream[j]);
  }
  ibstride = 0;
  if      ( std::is_same<T,GFLOAT>::value ) {
    for ( auto j=0; j<nstreams; j++ ) {
      iastart  = cudat.ibblk[j]*M*K*sizeof(T);
      iastride = (cudat.ieblk[j] - cudat.ibblk[j] + 1)*M*K*sizeof(T);
      cublasSgemmStridedBatched( 
                   cudat.hbatch_cublas, 
                   cuBlasOps[TransA-CblasNoTrans], cuBlasOps[TransB-CblasNoTrans],
                   M, N, K, 
                   (const float*)&alpha, (const float*)(A+iastart), lda, iastride,
                   (const float*)B, ldb, ibstride, (const float*)&beta,
                   (float*)(C+iastart), ldc, iastride, nstreams);
    }
  }
  else if ( std::is_same<T,GDOUBLE>::value ) {
    for ( auto j=0; j<nstreams; j++ ) {
      iastart  = cudat.ibblk[j]*M*K*sizeof(T);
      iastride = (cudat.ieblk[j] - cudat.ibblk[j] + 1)*M*K*sizeof(T);
      cublasDgemmStridedBatched( 
                   cudat.hbatch_cublas, 
                   cuBlasOps[TransA-CblasNoTrans], cuBlasOps[TransB-CblasNoTrans],
                   M, N, K, 
                   (const double*)&alpha, (const double*)(A+iastart), lda, iastride,
                   (const double*)B, ldb, ibstride, (const double*)&beta,
                   (double*)(C+iastart), ldc, iastride, nstreams);
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


