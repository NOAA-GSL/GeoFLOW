//==================================================================================
// Module      : gmorton_keygen
// Date        : 7/1/2018 (DLR)
// Description : Encapsulates the access methods and data associated with
//               defining a Morton-type key-generator, as used in GASpAR.
// Copyright   : Copyright 2018. Colorado State University. All rights reserved
// Derived From: GKeyGen
//==================================================================================
#if !defined(MORTON_KEYGEN_HPP)
#define MORTON_KEYGEN_HPP

#include "gtypes.h"
#include <stdlib.h>
#include "gtpoint.hpp"
#include "gkeygen.hpp"
#include "gbitblock.hpp"

enum GMORTON_TYPE {GMORTON_INTERLEAVE=0, GMORTON_STACKED};

// Template args: TK is type for the key; TF is the type for the float
template<typename TK, typename TF> class GMorton_KeyGen: public GKeyGen<TK, TF>
{
public:
                           GMorton_KeyGen();
                          ~GMorton_KeyGen();

         void              setType(GMORTON_TYPE type);                              // Set Morton-order method
//       void              setOrigin(GTPoint<TF> &P0);                              // Set coord origin
//       void              setBox(GTPoint<TF> &P0, GTPoint<TF> &P1);                // Set grid bding box
         void              setIntegralLen(GTPoint<TF> &inP0, GTPoint<TF> &dX);      // Set integral length
         void              setDoLog(GBOOL bDoLog);                                  // Set btakelog_ flag
         void              key(GTVector<TK> &id, GTVector<GTPoint<TF>> &point);     // Use GTPoint
         void              key(TK id[], GTPoint<TF> point[], GINT  n);              // Use GTPoint
         void              key(GTVector<TK> &id, GTVector<GTVector<TF>> &x);        // Use GTVector for points
         void              key(GTVector<TK> &id, GTVector<GTVector<TF>> &x, 
                           GTVector<GINT> &ix);                                     // Use ix elements of GTVector x for points


private:

         // Member data:
         GMORTON_TYPE      itype_;
         GBOOL             btakelog_;
         GBOOL             bintlenset_;
         TF                logFact_;
         TF                delmax_;
         TF                idelmax_;
         TF                idel_;
         TF                ttiny_;
         GTPoint<TF>       P0_;
         GTPoint<TF>       P1_;
         GTPoint<TF>       idX_;
         GTPoint<TF>       dX_;
         GBitBlock         *bb_;
};

#include "gmorton_keygen.ipp"

#endif
