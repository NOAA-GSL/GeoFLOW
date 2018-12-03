//==================================================================================
// Module       : gmm.cpp
// Date         : 2/25/18 (DLR)
// Description  : Namespace to hold GeoFLOW memory manager/pool  methods.
//                  new, delete: Override new and delete operators by getting
//                               pointers to static memory block
//                  
// Copyright    : Copyright 2018. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================

#include "gmm.h"
using namespace GMM;

#if !defined(GMM_GLOBAL_DATA)
#define GMM_GLOBAL_DATA
GSIZET noccregions_ =0;                    // No. occupied regions
GUCHAR memblk_[G_MAX_MEMBLOCK_BYTES]= {};  // Memory block
GSIZET iocc_  [G_MAX_MEMBLOCK_BYTES]
                         = {G_MEMLOCNULL}; // Location of start of occupied mem region
GSIZET nocc_  [G_MAX_MEMBLOCK_BYTES]= {0}; // Number bytes occupied at each region start
#endif


//**********************************************************************************
//**********************************************************************************
// METHOD : New operator
// DESC   : Override std::new
// ARGS   : GSIZET n : no. bytes to allocate on stack
// RETURNS: Valid pointer on success; else NULLPTR. 
//**********************************************************************************
void *operator new(GSIZET n) 
{
  return ::operator new[](n);

} // end, new[] method


//**********************************************************************************
//**********************************************************************************
// METHOD : New operator[]
// DESC   : Override std::new
// ARGS   : GSIZET n : no. bytes to allocate on stack
// RETURNS: Valid pointer on success; else NULLPTR. 
//**********************************************************************************
void *operator new[](GSIZET n) // throw(std::bad_alloc)
{
  if ( n <= 0 ) 
  {
     std::cout << "GMM::new: cannot allocate block of size " << n << " bytes." << std::endl;
     throw std::bad_alloc();
  }

  // Find next valid contiguous pointer location. 
  //
  // Recall that arrays containing region data are
  // stored from smallest region index to largest:
  // Start at end of list of occupied regions, and work 
  // backwards to start of memblk, through occupied regions, 
  // until one is found that has more elements than required
  // for that region, and can accommodate request:
  GUCHAR *pret=NULLPTR;
  GSIZET j=noccregions_-1;              // Start at region furthest from beg
  GSIZET itst, ndel, iret;
  GSIZET ireg = G_MAX_MEMBLOCK_BYTES-1; // Boundary index of region
  while ( j>=0 && pret==NULLPTR )       // Cycle over existing regions
  {
     itst = iocc_[j] + nocc_[j];                  // Test index
     ndel = ireg - itst;                          // Max avail # bytes at test location 
     pret = ndel >= n ? &memblk_[ndel] : NULLPTR; // Set pointer, if enough space
     iret = pret != NULLPTR ? itst : G_MEMLOCNULL;// Set ret index if enough space
     ireg = iocc_[j];                             // Move region boundary
     j--;
  }
 
  if ( pret == NULLPTR ) 
  {
     std::cout << "GMM::new: cannot allocate block of size " << n << " bytes." << std::endl;
     throw std::bad_alloc();
  }

  iocc_[noccregions_] = iret;
  nocc_[noccregions_] = n;
  noccregions_++; 

  GMM::regsort();

  return pret;
 
} // end, new operator

//**********************************************************************************
//**********************************************************************************
// METHOD : delete[] operator
// DESC   : Override std::delete[]
// ARGS   : void *p : pointer to deallocate
// RETURNS: none.
//**********************************************************************************
void operator delete[](void *p) throw()
{
  ::operator delete(p);
} // end, operator delete[]

//**********************************************************************************
//**********************************************************************************
// METHOD : delete operator
// DESC   : Override std::delete
// ARGS   : void *p : pointer to deallocate
// RETURNS: none.
//**********************************************************************************
void operator delete(void *p) throw()
{

  // First, find index of p in memblk:
  GSIZET idx = (GUCHAR*)p - &memblk_[0];
  if ( idx <= 0 || idx >= G_MAX_MEMBLOCK_BYTES )
  {
     std::cout << "GMM::delete: invalid pointer: " << (GUCHAR*)p << std::endl;
     throw  std::unexpected;
  }

  // Find region that p occupies:
  GSIZET ireg = G_MEMLOCNULL;
  GSIZET j = 0;
  while ( ireg < 0 && j < noccregions_ )
  {
    ireg = iocc_[j] == idx ? j : G_MEMLOCNULL;
    j++;
  }
  
  if ( ireg == G_MEMLOCNULL )
  {
     std::cout << "GMM::delete: invalid region index: " << ireg  << std::endl;
     throw  std::unexpected;
  }
  iocc_[ireg] = G_MEMLOCNULL;
  nocc_[ireg] = 0;
  noccregions_--; 
  regsort();

} // end, delete operator
  
//**********************************************************************************
//**********************************************************************************
// METHOD : nfree
// DESC   : Total number of bytes remaining. Says nothing about the
//          amount of fragmenting of memory block.
// ARGS   : none.
// RETURNS: GSIZET number.
//**********************************************************************************
GSIZET GMM::nfree()
{
  GSIZET jsum=0.0;
  for ( GSIZET j=0; j<G_MAX_MEMBLOCK_BYTES; j++ )
  {
     jsum += nocc_[j]; 
  }
  
  return jsum;
} // end, method nfree


//**********************************************************************************
//**********************************************************************************
// METHOD : nused
// DESC   : Total number of bytes used. Says nothing about the
//          amount of fragmenting of memory block.
// ARGS   : none.
// RETURNS: GSIZET total.
//**********************************************************************************
GSIZET GMM::nused()
{
  return G_MAX_MEMBLOCK_BYTES-nfree();
} // end, method nused


//**********************************************************************************
//**********************************************************************************
// METHOD : fragfrac
// DESC   : Estimate of fragmentation fraction, computed as
//           (avg # unallocated bytes between data regions) /
//           (avg size of data regions)
//
// ARGS   : none.
// RETURNS: GDOUBLE fraction.
//**********************************************************************************
GDOUBLE GMM::fragfrac()
{
  GDOUBLE xunalloc=0.0;
  GDOUBLE xszreg  =0.0;

  // Don't count the final region, which is, 
  // by definition, not fragmented. Recall
  // that this 'final' region is represented by
  // region data at region index noccregions_-1, since this
  // data is sorted smallest region index to largest:
  GSIZET ireg; // largest bounding index for region
  for ( GSIZET j=0; j<noccregions_-1; j++ )
  {
    ireg      = iocc_[j+1];
    xszreg   += (GDOUBLE)nocc_[j];
    xunalloc += (GDOUBLE)( ireg - iocc_[j] - nocc_[j] );
  }
  return xszreg == 0.0 ? 0.0 : xunalloc/xszreg;
} // end, method fragfrac


//**********************************************************************************
//**********************************************************************************
// METHOD : nregn
// DESC   : Get number of contiguous memory regions
// ARGS   : none.
// RETURNS: GSIZET number regions
//**********************************************************************************
GSIZET GMM::nregn()
{
  return noccregions_;
} // end, method nregn


//**********************************************************************************
//**********************************************************************************
// METHOD : regsort
// DESC   : Sort pointer/location occupation region arrays from largest 
//          memblk index to smallest
// ARGS   : none.
// RETURNS: none.
//**********************************************************************************
void GMM::regsort()
{
  GSIZET idx;
 
  // Cycle through unsorted subarray.
  // NOTE: we look at noccregions_+1, in the case of 
  //       the delete operator, where a previously sorted 
  //       (when number of regions was noccregions_+1) region 
  //       index is set to G_MEMLOCNULL (<0):
  GSIZET nmax = MIN(G_MAX_MEMBLOCK_BYTES,noccregions_+1);
  GSIZET itmp;
  for ( GSIZET i=0;  i<nmax; i++)
  {
    // Find the maximum element in unsorted array
    idx = i;
    for (GSIZET j=i+1; j<nmax; j++)
    {
      idx = iocc_[j] >=0 && iocc_[j] < iocc_[idx] ? j : idx;
    }

    // Swap the newly found max element with the ith element:
    itmp           = iocc_[i]; 
    iocc_      [i] = iocc_[idx];
    iocc_    [idx] = itmp;

    itmp           = nocc_[i]; 
    nocc_      [i] = nocc_[idx];
    nocc_    [idx] = itmp;

  }
  
} // end, method regsort
  
