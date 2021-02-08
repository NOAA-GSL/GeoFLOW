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

template<typename T>
class GShapeFcn_base
{

public:

                          GShapeFcn_base(GINT dim=2):dim_(dim) {  
                            assert(dim_>=1 && dim_<=3);
                            this->gbasis_.resize(dim_);
                            this->gbasis_ = NULLPTR; }
                          GShapeFcn_base(GTVector<GNBasis<GCTYPE,GFTYPE>*> &b, GINT dim=2):dim_(dim) {
                            assert(dim_>=1 && dim_<=3); set_basis(b);}
                          GShapeFcn_base(const GShapeFcn_base &) {};
virtual                  ~GShapeFcn_base() {};

virtual  void             operator=(const GShapeFcn_base &obj)
                          {gbasis_ = obj.gbasis_;}
                   
// Methods:
virtual void              Ni(GTVector<GINT> &I, GTVector<GTVector<T>*> &xi, GTVector<T> &Ni)=0;                                // shape function,  N_I
virtual void              dNdXi(GTVector<GINT> &I, GINT j, GTVector<GTVector<T>*> &xi,
                                                 GTVector<T> &out)=0; // dN_I/Xi^j for I, j
virtual void              set_basis(GTVector<GNBasis<GCTYPE,GFTYPE>*> &b)  // set basis functions
                          { gbasis_.resize(b.size()); gbasis_ = b; }

protected:

GINT                                dim_    ; // problem dim (1, 2, 3)     
GTVector<GNBasis<GCTYPE,GFTYPE>*>   gbasis_ ; // expansion basis: one for each direction

};
#endif
