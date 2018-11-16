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
        GSIZET               ndof();                                  // compute total number local dof
        GFTYPE               minsep();                                // find min nodal sep
        GFTYPE               maxsep();                                // find max nodal sep

friend  std::ostream&        operator<<(std::ostream&, GGrid &);    // Output stream operator
 

protected:
       

private:

GElemList                   gelems_;  // element list
GTVector<GTVector<GSIZET>>  itype_;   // indices in elem list of each type
GTVector<GSIZET>            ntype_;   // no. elems of each type on grid

};

#endif
