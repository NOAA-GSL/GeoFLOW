//==================================================================================
// Module      : gmorton_keygen
// Date        : 7/1/2018 (DLR)
// Description : Encapsulates the access methods and data associated with
//               defining a Morton-type key-generator, as used in GASpAR.
// Copyright   : Copyright 2018. Colorado State University. All rights reserved
// Derived From: GKeyGen
//==================================================================================
#include <string>
#include <limits>
#include "gmorton_keygen.hpp"
#include "tbox/tracer.hpp"

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor Method (1)
// DESC   :
// ARGS   : none.
// RETURNS: none.
//**********************************************************************************
template<typename TK, typename TF> 
GMorton_KeyGen<TK,TF>::GMorton_KeyGen()
: GKeyGen<TK,TF>(),
itype_       (GMORTON_INTERLEAVE),
btakelog_    (FALSE),
bintlenset_  (FALSE),
delmax_      (std::numeric_limits<TF>::min()),
idelmax_     (1.0),
idel_        (1.0),
ttiny_       (100.0*std::numeric_limits<TF>::min()),
bb_          (NULLPTR)
{
	GEOFLOW_TRACE();
  idX_.x1  = idX_.x2 = idX_.x3 = ttiny_;
  logFact_ = log10(1.0/ttiny_);
  bb_       = new GBitBlock(BITSPERBYTE*sizeof(TK));
} // end of constructor method (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor
// DESC   :
// ARGS   : none.
// RETURNS: none.
//**********************************************************************************
template<typename TK, typename TF> 
GMorton_KeyGen<TK,TF>::~GMorton_KeyGen()
{
  if ( bb_ != NULLPTR ) delete bb_;
}

//**********************************************************************************
//**********************************************************************************
// METHOD     : setType
// DESCRIPTION: Morton-ordering method
// ARGUMENTS  : intype: type.
// RETURNS    : none.
//**********************************************************************************
template<typename TK, typename TF>
void GMorton_KeyGen<TK,TF>::setType(GMORTON_TYPE intype)
{
  itype_ = intype;
} // end of method setType


#if 0
//**********************************************************************************
//**********************************************************************************
// METHOD     : setOrigin
// DESCRIPTION: Set geometric origin for Morton-ordering
// ARGUMENTS  : inP0: GFPoint type for coordinate origin (Cartesian)
// RETURNS    : none.
//**********************************************************************************
template<typename TK, typename TF>
void GMorton_KeyGen<TK,TF>::setOrigin(GTPoint<TF> &inP0)
{
  P0_ = inP0;
} // end of method setOrigin


//**********************************************************************************
//**********************************************************************************
// METHOD     : setBox
// DESCRIPTION: Set geometric extent
// ARGUMENTS  : inP0: GFPoint type for coordinate origin (Cartesian)
//              inP1: Point type for coordinate diagonal
// RETURNS    : none.
//**********************************************************************************
template<typename TK, typename TF>
void GMorton_KeyGen<TK,TF>::setBox(GTPoint<TF> &inP0, GTPoint<TF> &inP1)
{
	GEOFLOW_TRACE();
  GINT i;

  P0_ = inP0;
  P1_ = inP1;
  for ( i=0,delmax_=0.0; i<P0_.dim(); i++ ) {
    dX_    [i] =       fabs(P1_[i] - P0_[i] + ttiny_);
    idX_   [i] = 1.0 / dX_[i];
    delmax_    = MAX(delmax_,dX_[i]);
  }
  logFact_ = log10(1.0 + delmax_/ttiny_);    
  idelmax_ = 1.0/delmax_;
  bintlenset_ = FALSE;
 
} // end of method setBox
#endif


//**********************************************************************************
//**********************************************************************************
// METHOD     : setDoLog
// DESCRIPTION: Set flag indicating whether or not to take log of difference before
//              integralizing.
// ARGUMENTS  : bDoLog: TRUE or FALSE flag
// RETURNS    : none.
//**********************************************************************************
template<typename TK, typename TF>
void GMorton_KeyGen<TK,TF>::setDoLog(GBOOL bDoLog)
{
  btakelog_ = bDoLog;
} // end of method setDoLog


//************************************************************************************
//************************************************************************************
// METHOD     : setIntegralLen
// DESCRIPTION: Set integral length for Morton-ordering. This is the smallest
//              length in each Cartesian direction such that it integralizes by
//              unity a quantity P-P0, where P is an arbitrary point, and P0 is
//              the coord origin, set in setBox
// ARGUMENTS  : inP0: GFPoint type for coordinate origin (Cartesian)
//              indX: GFPoint type for integral length(s) (Cartesian)
// RETURNS    : none.
//************************************************************************************
template<typename TK, typename TF>
void GMorton_KeyGen<TK,TF>::setIntegralLen(GTPoint<TF> &inP0, GTPoint<TF> &indX)
{
	GEOFLOW_TRACE();
  TF dmin=std::numeric_limits<TF>::max();

  delmax_ = std::numeric_limits<TF>::min();

  P0_.resize(inP0.dim());
  dX_.resize(inP0.dim());
  P0_ = inP0;
  dX_ = indX;
  for ( GSIZET i=0; i<P0_.dim(); i++ ) {
    if ( dX_[i] <= 0 ) {
      cout << "GMorton_KeyGen<TK,TF>::setIntegralLen: dX=" << dX_ << endl;
    }
    assert(dX_[i] > 0.0 && "Invalid integral length scale");
    idX_   [i] = 1.0 / dX_[i];
    delmax_    = MAX(delmax_,dX_[i]);
    dmin       = MIN(dmin,indX[i]);
  }
  logFact_ = log10(1.0 + delmax_/ttiny_);    
  idelmax_ = 1.0/delmax_;
  idel_    = 1.0/dmin;
  bintlenset_ = TRUE;

} // end of method setIntegralLen


//**********************************************************************************
//**********************************************************************************
// METHOD     : key (1)
// DESCRIPTION: Computes Morton-ordered key. GMORTON_TYPE is defined as follows:
//              Let X = x7 x6 x5 x4 x3 x2 x1 x0  and
//                  Y = y7 y6 y5 y4 y3 y2 y1 y0 
//              be the X,Y values of a point, where x0, ... x7 represent the bits 
//              0-7 of an 8-bit float, and the same for Y.
//              Then GMORTON_INTERLEAVE means that the computed key will be
//                  KEY = y7x7 y6x6 ... y2x2 y1x1 y0x0
//              while for GMORTON_STACKED the key is s.t.:
//                  KEY = y7y6y5...y2y1y0x7x6...x2x1x0.
//              with an obvious extension to 3D
// ARGUMENTS  : 
//              id   : Array of length 'n' containing the resultant keys. This
//                     must be allocated by caller and be of length 'n'.
//              idsz : size of 'id' array elements in bytes
//              point: Array of length 'n' of Cartesian points for which to 
//                     generate a Morton-type key. All points must have the same
//                     dimension; no checking done
//              n    : Number of points for which to generate keys.
//
// RETURNS    : none.
//**********************************************************************************
template<typename TK, typename TF>
void GMorton_KeyGen<TK,TF>::key(TK id[], GTPoint<TF> point[], GINT  n)
{
	GEOFLOW_TRACE();
  GString   serr = " GMorton_KeyGen<TK,TF>::key(1): ";
  GINT      i, idsz, j, k, ib, ishft, nbits, nbpc, tbits;
  GINT      gdim = point[0].dim();
  GUINT     ix[3]={0,0,0}, ixt;
  GUSHORT   ui;
  GDOUBLE   del;

  idsz = sizeof(TK);

  nbits = BITSPERBYTE * idsz / gdim; // - 1; // no bits per coord. direction
//nbpc  = BITSPERBYTE * sizeof(ix[0])   - 1; // no. bits per integral coordinate
  tbits = gdim * nbits;                      // total no. bits for Morton integers


  if ( !bintlenset_ ) {
    idel_ = pow(2.0,nbits) * idelmax_;
    bintlenset_ = TRUE;
  }

  logFact_ = fabs( pow(2.0,nbits) / log10(ttiny_) );

  // interleave integer bits into key:
  if ( itype_ == GMORTON_INTERLEAVE ) { // full-interleaving
    for ( i=0; i<n; i++ ) {
      bb_->reset();
      for ( k=0; k<gdim; k++ ) { // integralize position in each dir
        del   = fabs(point[i][k]-P0_[k]);
        ix[k] = btakelog_ ? (GINT)(logFact_*log10(del*idelmax_)) 
                          : (GINT)(               del*idel_+0.5);
      }
      for ( j=0,ib=0; j<nbits; j++ ) {
        for ( k=0; k<gdim; k++,ib++ ) {
          ixt   = ix[k];
#if 0
          ishft = ib % nbpc; 
#else
          ishft = j;
#endif
          ui    = ( (ixt >> ishft) & ~(~0 << 1) );     // shift & mask with 0000001
          bb_->setBits(ib, ui);
        }
      }
      bb_->transferBits(((GBYTE*)id+i*idsz), 1, idsz);

#if defined(_GBYTESWAP)
      GCUtils::swap(((GBYTE*)id+i*idsz),idsz);
#endif
    } // end, i loop
  } // end, GMORTON_INTERLEAVE
  else  {                                 // partial interleaving
    for ( i=0; i<n; i++ ) {
      bb_->reset();
      for ( k=0; k<gdim; k++ ) { // integralize position in each dir
        del   = fabs(point[i][k]-P0_[k]);
        ix[k] = btakelog_ ? (GINT)(logFact_*log10(del*idelmax_)) 
                          : (GINT)(               del*idel_); //+0.5);
      }
      for ( k=0,ib=0; k<gdim; k++ ) {
        ixt = ix[k];
        for ( j=0; j<nbits; j++,ib++ ) {
#if 0
          ishft = ib % nbpc; 
#else
          ishft = j;
#endif
          ui    = ( (ixt >> ishft) & ~(~0 << 1) );     // shift & mask with 0000001
          bb_->setBits(ib, ui);
        }
      }
      bb_->transferBits(((GBYTE*)id+i*idsz), 1, idsz);
#if defined(_GBYTESWAP)
      GCUtils::swap(((GBYTE*)id+i*idsz),idsz);
#endif
    } // end, i loop
  } // end, GMORTON_STACKED

} // end of method key (1)


//**********************************************************************************
//**********************************************************************************
// METHOD     : key (2)
// DESCRIPTION: Computes Morton-ordered key
//              GMORTON_TYPE is defined as follows:
//              Let X = x7 x6 x5 x4 x3 x2 x1 x0  and
//                  Y = y7 y6 y5 y4 y3 y2 y1 y0 
//              be the X,Y values of a point, where x0, ... x7 represent the bits 
//              0-7 of an 8-bit float, and the same for Y.
//              Then GMORTON_INTERLEAVE means that the computed key will be
//                  KEY = y7x7 y6x6 ... x2x2 y1x1 y0x0
//              while for GMORTON_STACKED the key is s.t.:
//                  KEY = y7y6y5...y2y1y0x7x6...x2x1x0.
//              with an obvious extension to 3D
// ARGUMENTS  : 
//              id   : Array of length 'n' containing the resultant keys. This
//                     must be allocated by caller and be of size 'n'.
//              point: Vector of spatial point locations
//
// RETURNS    : none.
//**********************************************************************************
template<typename TK, typename TF>
void GMorton_KeyGen<TK,TF>::key(GTVector<TK> &id, GTVector<GTPoint<TF>> &point)
{
	GEOFLOW_TRACE();
  GINT        idsz=sizeof(TK);
  GTPoint<TF> p;

  if ( id.size() != point.size() ) {
    std::cout << "GMorton_KeyGen<TK,TF>::key (3): incompatible GTVector dimensions" << std::endl;
    exit(1);
  }
  key(id.data(), point.data(), id.size());

} // end of method key (2)


//**********************************************************************************
//**********************************************************************************
// METHOD     : key (3)
// DESCRIPTION: Computes Morton-ordered key
//              GMORTON_TYPE is defined as follows:
//              Let X = x7 x6 x5 x4 x3 x2 x1 x0  and
//                  Y = y7 y6 y5 y4 y3 y2 y1 y0 
//              be the X,Y values of a point, where x0, ... x7 represent the bits 
//              0-7 of an 8-bit float, and the same for Y.
//              Then GMORTON_INTERLEAVE means that the computed key will be
//                  KEY = y7x7 y6x6 ... x2x2 y1x1 y0x0
//              while for GMORTON_STACKED the key is s.t.:
//                  KEY = y7y6y5...y2y1y0x7x6...x2x1x0.
//              with an obvious extension to 3D
// ARGUMENTS  : 
//              id   : Array of length 'n' containing the resultant keys. This
//                     must be allocated by caller and be of size 'n'.
//              x    : Array of coordinates with each being a vector
//
// RETURNS    : none.
//**********************************************************************************
template<typename TK, typename TF>
void GMorton_KeyGen<TK,TF>::key(GTVector<TK> &id, GTVector<GTVector<TF>> &x)
{
	GEOFLOW_TRACE();
  assert(id.size() == x[0].size() && "GMorton_KeyGen::key(3): incompatible dimensions ");
  GTPoint<TF> p(x.size());

  for ( GSIZET j=0; j<x[0].size(); j++ ) {
    for ( GSIZET i=0; i<p.dim(); i++ ) p[i] = x[i][j];
    key(&id[j], &p, 1);
  }

} // end of method key (3)


//**********************************************************************************
//**********************************************************************************
// METHOD     : key (4)
// DESCRIPTION: Computes Morton-ordered key
//              GMORTON_TYPE is defined as follows:
//              Let X = x7 x6 x5 x4 x3 x2 x1 x0  and
//                  Y = y7 y6 y5 y4 y3 y2 y1 y0 
//              be the X,Y values of a point, where x0, ... x7 represent the bits 
//              0-7 of an 8-bit float, and the same for Y.
//              Then GMORTON_INTERLEAVE means that the computed key will be
//                  KEY = y7x7 y6x6 ... x2x2 y1x1 y0x0
//              while for GMORTON_STACKED the key is s.t.:
//                  KEY = y7y6y5...y2y1y0x7x6...x2x1x0.
//              with an obvious extension to 3D
// ARGUMENTS  : 
//              id   : Array of length 'n' containing the resultant keys. This
//                     must be allocated by caller and be of size 'n'.
//              x    : Array of coordinates with each being a vector
//
// RETURNS    : none.
//**********************************************************************************
template<typename TK, typename TF>
void GMorton_KeyGen<TK,TF>::key(GTVector<TK> &id, GTVector<GTVector<TF>*> &x)
{
	GEOFLOW_TRACE();
  assert(id.size() == x[0]->size() && "GMorton_KeyGen::key(4): incompatible dimensions");
  GTPoint<TF> p(x.size());

  for ( GSIZET j=0; j<x[0]->size(); j++ ) {
    for ( GSIZET i=0; i<p.dim(); i++ ) p[i] = (*x[i])[j];
    key(&id[j], &p, 1);
  }

} // end of method key (4)


//**********************************************************************************
//**********************************************************************************
// METHOD     : key (5)
// DESCRIPTION: Computes Morton-ordered key
//              GMORTON_TYPE is defined as follows:
//              Let X = x7 x6 x5 x4 x3 x2 x1 x0  and
//                  Y = y7 y6 y5 y4 y3 y2 y1 y0 
//              be the X,Y values of a point, where x0, ... x7 represent the bits 
//              0-7 of an 8-bit float, and the same for Y.
//              Then GMORTON_INTERLEAVE means that the computed key will be
//                  KEY = y7x7 y6x6 ... x2x2 y1x1 y0x0
//              while for GMORTON_STACKED the key is s.t.:
//                  KEY = y7y6y5...y2y1y0x7x6...x2x1x0.
//              with an obvious extension to 3D
// ARGUMENTS  : 
//              id   : Array of length 'n' containing the resultant keys. This
//                     must be allocated by caller and be of size 'n'.
//              x    : Array of coordinates with each being a vector
//              ix   : array of indices into x to provide indices, id. There 
//                     will ix.size() indices provided in id
//
// RETURNS    : none.
//**********************************************************************************
template<typename TK, typename TF>
void GMorton_KeyGen<TK,TF>::key(GTVector<TK> &id, GTVector<GTVector<TF>> &x, 
                                GTVector<GINT> &ix)
{
	GEOFLOW_TRACE();
  assert(id.size() == ix.size() && "GMorton_KeyGen::key(5): incompatible dimensions ");
  GTPoint<TF> p(x.size());

  for ( GSIZET j=0; j<ix.size(); j++ ) {
    for ( GSIZET i=0; i<p.dim(); i++ ) p[i] = x[i][ix[j]];
    key(&id[j], &p, 1);
  }

} // end of method key (5)


