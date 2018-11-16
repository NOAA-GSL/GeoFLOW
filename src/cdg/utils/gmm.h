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
extern GSIZET nnew_;
extern GSIZET ndelete_;
extern GUCHAR memblk_[]; 
extern GBOOL  bverbose_;
extern GLONG  iocc_  []; 
extern GSIZET nocc_  [];

void  *operator new(GSIZET n);
void  *operator new[](GSIZET n) throw(std::bad_alloc); 
void   operator delete(void *) throw();
void   operator delete[](void *) throw();
namespace GMM
{
  GSIZET  nregn();
  GSIZET  bytesfree();
  GSIZET  bytesused();
  GSIZET  maxbytesused();
  GSIZET  nnew();
  GSIZET  ndelete();
  void    regsort ();
  GDOUBLE fragfrac();
  void    bverbose(const GBOOL);
  void    prstats(const char *);
//void defrag();
} // end, namespace GMM

#endif
