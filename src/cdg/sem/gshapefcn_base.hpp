//==================================================================================
// Module       : gshapefcn_base.hpp
// Date         : 9/19/18 (DLR)
// Description  : Forms abstract base class for iso-parametric shape functions, N_i.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none
//==================================================================================

#if !defined(_GSHAPEFCN_BASE_HPP)
#define _GSHAPEFCN_BASE_HPP
#include "gnbasis.hpp"
#include "gtvector.hpp"

class GShapeFcn_base
{

public:

                          GShapeFcn_base() {};
                          GShapeFcn_base(GTVector<GNBasis<GCTYPE,GFTYPE>*> &b) {};
                          GShapeFcn_base(const GShapeFcn_base &) {};
virtual                  ~GShapeFcn_base() {};

virtual  void             operator=(const GShapeFcn_base &obj)
                          {gbasis_ = obj.gbasis_;}
                   
// Methods:
virtual void              Ni(GTVector<GINT> &I, GTVector<GTVector<GFTYPE>*> &xi, GTVector<GFTYPE> &Ni)=0;                                // shape function,  N_I
virtual void              dNdXi(GTVector<GINT> &I, GINT j, GTVector<GTVector<GFTYPE>*> &xi,
                                                 GTVector<GFTYPE> &out)=0; // dN_I/Xi^j for I, j
virtual void              set_basis(GTVector<GNBasis<GCTYPE,GFTYPE>*> &b)  // set basis functions
                          { gbasis_.resize(b.size()); gbasis_ = b; }

protected:

GTVector<GNBasis<GCTYPE,GFTYPE>*>   gbasis_       ; // expansion basis: one for each direction

};
#endif
