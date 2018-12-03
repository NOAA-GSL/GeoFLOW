//==================================================================================
// Module       : gmm.h
// Date         : 2/25/18 (DLR)
// Description  : Namespace to hold GeoFLOW memory manager/pool methods.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================

#include "gtypes.h"

#if !defined(GMM_H)
#define GMM_H

extern GSIZET noccregions_;
extern GUCHAR memblk_[]; 
extern GSIZET iocc_  []; 
extern GSIZET nocc_  [];

void  *operator new(GSIZET n);
void  *operator new[](GSIZET n) throw(std::bad_alloc); 
void   operator delete(void *) throw();
void   operator delete[](void *) throw();
namespace GMM
{
  GSIZET  nregn();
  GSIZET  nfree();
  GSIZET  nused();
  void    regsort ();
  GDOUBLE fragfrac();
//void defrag();
} // end, namespace GMM

#endif
