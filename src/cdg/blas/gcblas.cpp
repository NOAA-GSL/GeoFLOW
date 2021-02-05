//==================================================================================
// Module       : gcblas.cpp
// Date         : 12/3/20 (DLR)
// Description  : GCBLASnamespace definitions
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#include "gcblas.hpp"

namespace GCBLAS
{

//**********************************************************************************
//**********************************************************************************
// METHOD : handle_create
// DESC   : 
// ARGS   : 
// RETURNS: GBlasHandle
//**********************************************************************************
void handle_create(GBlasHandle &h)
{
  GBlasStatus stat;
  #if defined(GEOFLOW_USE_CUBLAS)
    stat = cublasCreate(&h);
    assert(stat == CUBLAS_STATUS_SUCCESS);
  #else
    h = -1;
  #endif

} // end, handle_create


//**********************************************************************************
//**********************************************************************************
// METHOD : handle_destroy
// DESC   : 
// ARGS   : 
// RETURNS: none.
//**********************************************************************************
void handle_destroy(GBlasHandle &h)
{
  #if defined(GEOFLOW_USE_CUBLAS)
    cublasDestroy(h);
  #endif
} // end, handle_destroy


//**********************************************************************************
//**********************************************************************************
// METHOD : stream_create
// DESC   : 
// ARGS   : 
// RETURNS: none.
//**********************************************************************************
void stream_create(GCuStream &pstream)
{
  GBlasError ret;
  #if defined(GEOFLOW_USE_CUBLAS)
    ret = cudaStreamCreate(&pstream);
    assert(ret == cudaSuccess);
  #endif
} // end, stream_create


//**********************************************************************************
//**********************************************************************************
// METHOD : stream_destroy
// DESC   : 
// ARGS   : 
// RETURNS: none.
//**********************************************************************************
void stream_destroy(GCuStream &pstream)
{
  GBlasError ret;
  #if defined(GEOFLOW_USE_CUBLAS)
    ret = cudaStreamDestroy(pstream);
    assert(ret == cudaSuccess);
  #endif
} // end, stream_destroy



} // end, namespace GCUDA


