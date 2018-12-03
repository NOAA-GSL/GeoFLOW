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
#include <iostream>
#include <fstream>
#include <stdio.h>
#include "gmm.h"
using namespace GMM;

#if !defined(GMM_GLOBAL_DATA)
#define GMM_GLOBAL_DATA
GSIZET noccregions_ =0;                    // No. currently occupied regions (no. 'news')
GSIZET nnew_        =0;                    // Total no. 'new'
GSIZET nmaxused_    =0;                    // Max no. bytes used
GSIZET ndelete_     =0;                    // Total no. 'deletes'
GBOOL  bverbose_    =0;                    // Print verbose output to stdout
GUCHAR memblk_[G_MAX_MEMBLOCK_BYTES]= {};  // Memory block
GLONG  iocc_  [G_MAX_MEMBLOCK_BYTES]
                         = {G_MEMLOCNULL}; // Location of start of occupied mem region
                                           // Want to compare with G_MEMLOCNULL, 
                                           // so this is of GLONG type
                                              
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
void *operator new[](GSIZET n) throw(std::bad_alloc)
{
  if ( n <= 0 ) 
  {
     std::cout << "GMM::new: cannot allocate block of size " << n << " bytes." << std::endl;
     throw std::bad_alloc();
  }

  // Find next valid contiguous pointer location. 
  //
  // Recall that arrays containing region data are
  // stored from largest region index to smallest:
  // Start at end of list of occupied regions, and work 
  // backwards to start of memblk, through occupied regions, 
  // until one is found that has more elements than required
  // for that region, and can accommodate request:
  GSIZET j=0;                           // Start at region furthest from beg
  GSIZET ndel;
  GLONG  ireg = G_MAX_MEMBLOCK_BYTES-1; // Boundary index of region
  GLONG  itst, iret=0;                        
  GUCHAR *pret=noccregions_==0 && ireg>n ? memblk_ : NULLPTR;
  while ( j<noccregions_ && pret==NULLPTR )       // Cycle over existing regions
  {
     itst = iocc_[j] + nocc_[j];                  // Test index
     ndel = (GSIZET)ireg - (GSIZET)itst;          // Max avail # bytes at test location 
     pret = ndel >= n ? &memblk_[itst] : NULLPTR; // Set pointer, if enough space
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
  nnew_++;

  nmaxused_ = MAX(nmaxused_,bytesused());

  GMM::regsort();

  if ( bverbose_ > 0 ) 
  {
    std::cout << std::endl;
    for ( GSIZET j=0; j<noccregions_; j++ ) {
      std::cout << "GMM:new: iocc[" << j << "]=" << iocc_[j] << " nocc[" << j << "]=" << nocc_[j] << std::endl;
    }
    std::cout << std::endl;
  }

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
  GSIZET idx = (GSIZET)((GUCHAR*)p - memblk_);
  if ( idx <= 0 || idx >= G_MAX_MEMBLOCK_BYTES )
  {
     std::cout << "GMM::delete: invalid pointer: " << p << " idx: " << idx  << std::endl;
     throw  std::unexpected;
  }

  // Find region that p occupies:
  GLONG  ireg = G_MEMLOCNULL;
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
  ndelete_++;
  regsort();

  if ( bverbose_ > 0 ) 
  {
    std::cout << std::endl;
    for ( GSIZET j=0; j<noccregions_; j++ ) {
      std::cout << "GMM:delete: iocc[" << j << "]=" <<  iocc_[j] << " nocc[" << j << "]="<<  nocc_[j] << std::endl;
    }
    std::cout << std::endl;
  }

} // end, delete operator
  
//**********************************************************************************
//**********************************************************************************
// METHOD : bytesused
// DESC   : Current number of bytes remaining. Says nothing about the
//          amount of fragmenting of memory block.
// ARGS   : none.
// RETURNS: GSIZET number.
//**********************************************************************************
GSIZET GMM::bytesused()
{
  GSIZET jsum=0.0;
  for ( GSIZET j=0; j<noccregions_; j++ )
  {
     jsum += nocc_[j]; 
  }
  
  return jsum;
} // end, method bytesused


//**********************************************************************************
//**********************************************************************************
// METHOD : bytesfree
// DESC   : Current number of bytes available. Says nothing about the
//          amount of fragmenting of memory block.
// ARGS   : none.
// RETURNS: GSIZET total.
//**********************************************************************************
GSIZET GMM::bytesfree()
{
  return G_MAX_MEMBLOCK_BYTES-bytesused();
} // end, method bytesfree


//**********************************************************************************
//**********************************************************************************
// METHOD : maxbytesused
// DESC   : Max number of bytes used allocated as contiguous memory. 
// ARGS   : none.
// RETURNS: GSIZET total.
//**********************************************************************************
GSIZET GMM::maxbytesused()
{
  return nmaxused_;
} // end, method maxbytesused


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
  GLONG ireg; // largest bounding index for region
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
// METHOD : nnew
// DESC   : Get total number of new commands called
// ARGS   : none.
// RETURNS: GSIZET number 'news'
//**********************************************************************************
GSIZET GMM::nnew()
{
  return nnew_;
} // end, method nnew


//**********************************************************************************
//**********************************************************************************
// METHOD : ndelete
// DESC   : Get total number of 'deletes' called
// ARGS   : none.
// RETURNS: GSIZET number deletes
//**********************************************************************************
GSIZET GMM::ndelete()
{
  return ndelete_;
} // end, method ndelete


//**********************************************************************************
//**********************************************************************************
// METHOD : bverbose
// DESC   : Set verbose flag
// ARGS   : 0 (not verbose) or 1 (verbose)
// RETURNS: none.
//**********************************************************************************
void GMM::bverbose(GBOOL bv)
{
  bverbose_ = bv;
} // end, method bverbose


//**********************************************************************************
//**********************************************************************************
// METHOD : prstats
// DESC   : Prints GMM stats to file or stdout.
// ARGS   : const char *sfile: if NULLPTR, print to stdout; else to specified file.
// RETURNS: none.
//**********************************************************************************
void GMM::prstats(const char *sfile)
{
   std::ofstream os;

   if ( sfile == NULLPTR )
   {
     std::cout << "GMM Stats: " << std::endl;
     std::cout << "  Num current regions = " << GMM::nregn       () << std::endl;  
     std::cout << "  Num allocated bytes = " << GMM::bytesused   () << std::endl;
     std::cout << "  Num new calls       = " << GMM::nnew        () << std::endl;  
     std::cout << "  Num delete calls    = " << GMM::ndelete     () << std::endl;
     std::cout << "  Num free bytes      = " << GMM::bytesfree   () << std::endl;
     std::cout << "  Max num bytes used  = " << GMM::maxbytesused() << std::endl;
     std::cout << "  Fragmentaion measure= " << GMM::fragfrac    () << std::endl;
     return;
   }


/*
   os.open(sfile, std::ofstream::trunc | std::ofstream::out);

   os << "GMM Stats: " << std::endl 
      << "  Num current regions = " << GMM::nregn       () <<  std::endl 
      << "  Num allocated bytes = " << GMM::bytesused   () <<  std::endl
      << "  Num new calls       = " << GMM::nnew        () <<  std::endl   
      << "  Num delete calls    = " << GMM::ndelete     () <<  std::endl 
      << "  Num free bytes      = " << GMM::bytesfree   () <<  std::endl 
      << "  Max num bytes used  = " << GMM::maxbytesused() <<  std::endl 
      << "  Fragmentaion measure= " << GMM::fragfrac    () <<  std::endl;
   os.close();
*/
   FILE *fp;
   fp = fopen(sfile, "w");
   fprintf(fp,"GMM stats:\n");
   fprintf(fp,"  Num current regions = %zd\n",GMM::nregn       ()); 
   fprintf(fp,"  Num allocated bytes = %zd\n",GMM::bytesused   ());
   fprintf(fp,"  Num new calls       = %zd\n",GMM::nnew        ());
   fprintf(fp,"  Num delete calls    = %zd\n",GMM::ndelete     ());
   fprintf(fp,"  Num free bytes      = %zd\n",GMM::bytesfree   ());
   fprintf(fp,"  Max num bytes used  = %zd\n",GMM::maxbytesused());
   fprintf(fp,"  Fragmentaion measure= %zd\n",GMM::fragfrac    ());
   fclose(fp);

} // end, method prstats


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
  GSIZET ntmp;
  GLONG  itmp;
  for ( GSIZET i=0;  i<nmax; i++)
  {
    // Sort from max iocc to min iocc, so that iocc < 0 are placed
    // after valid (iocc>0) indices

    // Find the extremum element in unsorted array
    idx = i;
    for (GSIZET j=i+1; j<nmax; j++)
    {
      idx = iocc_[j] > iocc_[idx] ? j : idx;
    }

    // Swap the newly found extremum element with the ith element:
    itmp           = iocc_[i]; 
    iocc_      [i] = iocc_[idx];
    iocc_    [idx] = itmp;

    ntmp           = nocc_[i]; 
    nocc_      [i] = nocc_[idx];
    nocc_    [idx] = ntmp;

  }
  
} // end, method regsort
  
