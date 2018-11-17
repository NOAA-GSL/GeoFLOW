//==================================================================================
// Module       : ggrid
// Date         : 8/31/18 (DLR)
// Description  : GeoFLOW grid object. Is primarily a container for
//                a list of elements, as an indirection buffer
//                that breaks down elements by type.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#if !defined(GGRID_HPP)
#define GGRID_HPP 

#include <iostream>
#include "gtvector.hpp"
#include "gtmatrix.hpp"
#include "gelem_base.hpp"


typedef GTVector<GElem_base*> GElemList;


class GGrid 
{
public:
                             GGrid();
                            ~GGrid();

        void                do_typing(); // classify into types
        GElemList           &elems() { return gelems_; }              // get elem list
        GTVector<GSIZET> 
                            &ntype() { return ntype_; }               // no. elems of each type
        GTVector<GTVector<GSIZET>> 
                            &itype() { return itype_; }               // indices for all types
        GTVector<GSIZET>    &itype(GElemType i) { return itype_[i]; } // indices for type i    
        void                 print(GString fname, GBOOL bdof=FALSE);
        GSIZET               ndof();                                  // compute total number elem dof
        GSIZET               nsurfdof();                              // compute total number elem surf dof
        GFTYPE               minsep();                                // find min nodal sep
        GFTYPE               maxsep();                                // find max nodal sep
inline  GTMatrix<GTVector<GFTYPE>>
                           &dXidX(){ return dXidX_; };                // global Rij = dXi^j/dX^i
inline  GTVector<GFTYPE>   *dXidX(GINT i,GINT j){ return &dXidX_(i,j);}

inline  GTVector<GFTYPE>   &Jac(){ return Jac_; }                     // global Jacobian
inline  GTVector<GTVector<GFTYPE>>
                           &faceJac(){ return faceJac_; }             // global face Jacobian
inline  GTVector<GFTYPE>   &faceJac(GINT i){ return faceJac_[i];}     // global face Jacobian
inline  GTVector<GTVector<GTVector<GFTYPE>>>
                           &faceNormal(){ return faceNormal_;}        // global face Jacobian
inline  GTVector<GTVector<GFTYPE>>
                           &faceNormal(GINT i){ return faceNormal_[i];}  // global face Jacobian


friend  std::ostream&        operator<<(std::ostream&, GGrid &);    // Output stream operator
 

protected:
       

private:

GElemList                   gelems_;        // element list
GTVector<GTVector<GSIZET>>  itype_;         // indices in elem list of each type
GTVector<GSIZET>            ntype_;         // no. elems of each type on grid
GTMatrix<GTVector<GFTYPE>>  dXidX_;         // matrix Rij = dXi^j/dX^i, global
GTVector<GFTYPE>            Jac_;           // interior Jacobian
GTVector<GTVector<GFTYPE>>  faceJac_;       // face Jacobian
GTVector<GTVector<GFTYPE>>  
                            faceNormal_;    // normal to face at each node point (2d & 3d), global


};

#endif
