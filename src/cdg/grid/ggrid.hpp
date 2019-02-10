//==================================================================================
// Module       : ggrid
// Date         : 8/31/18 (DLR)
// Description  : GeoFLOW grid object. Is primarily a container for
//                a list of elements, as an indirection buffer
//                that breaks down elements by type.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#if !defined(_GGRID_HPP)
#define _GGRID_HPP 

#include <iostream>
#include "gtvector.hpp"
#include "gtmatrix.hpp"
#include "gcomm.hpp"
#include "gelem_base.hpp"


typedef GTVector<GElem_base*> GElemList;


class GGrid 
{
public:
                             GGrid(GC_COMM comm);
                            ~GGrid();

        void                do_typing(); // classify into types
        GElemList           &elems() { return gelems_; }              // get elem list
        GSIZET               nelems() { return gelems_.size(); }      // local num elems
        GTVector<GSIZET> 
                            &ntype() { return ntype_; }               // no. elems of each type
        GTVector<GTVector<GSIZET>> 
                            &itype() { return itype_; }               // indices for all types
        GTVector<GSIZET>    &itype(GElemType i) { return itype_[i]; } // indices for type i    
        GElemType            gtype() { return gtype_; }               // get unique elem type on grid       
        GFTYPE               integrate(GTVector<GFTYPE> &u,
                                       GTVector<GFTYPE> &tmp);        // spatial integration of global vector
        void                 print(GString fname, GBOOL bdof=FALSE);
        GSIZET               ndof();                                  // compute total number elem dof
        GSIZET               size() { return ndof(); }
        GSIZET               nsurfdof();                              // compute total number elem surf dof
        void                 init();                                  // inititlize class
        GFTYPE               minlength();                             // find min elem length
        GFTYPE               maxlength();                             // find max elem length
        GTVector<GFTYPE>    &minnodedist()                            // find min node distance
                            {return minnodedist_; }

        GTMatrix<GTVector<GFTYPE>>
                            &dXidX();                                 // global Rij = dXi^j/dX^i
        GTVector<GFTYPE>    &dXidX(GSIZET i,GSIZET j);                // Rij matrix element 
        GTvector<GTVector<GFTYPE>>
                            &xNodes() { return xNodes_; }             // get all nodes coords (Cart)
        GTvector<GFTYPE>    &xNodes(GSIZET i) { return xNodes_[i]; }  // get all nodes coords (Cart)
                            

        GTVector<GFTYPE>    &Jac();                                    // global Jacobian
        GTVector<GFTYPE>
                            &faceJac();                                // global face Jacobian
        GTVector<GTVector<GFTYPE>>
                            &faceNormal();                             // global face normals
        GTVector<GTVector<GSIZET>>
                            &igbdy() { return igbdy_;}                 // global dom bdy indices into u for eacb GBdyType


friend  std::ostream&        operator<<(std::ostream&, GGrid &);      // Output stream operator
 

protected:
       
        void                 def_init();                              // iniitialze deformed elems
        void                 reg_init();                              // initialize regular elems

private:
         
void                        find_min_dist(); 
void                        init_bc_info(); 

GBOOL                       bInitialized_;  // object initialized?
GC_COMM                     comm_;          // communicator
GElemList                   gelems_;        // element list
GTVector<GTVector<GSIZET>>  itype_;         // indices in elem list of each type
GTVector<GSIZET>            ntype_;         // no. elems of each type on grid
GTMatrix<GTVector<GFTYPE>>  dXidX_;         // matrix Rij = dXi^j/dX^i, global
GTVector<GTVector<GFTYPE>>  xNodes_;        // Cart coords of all node points
GTVector<GFTYPE>            Jac_;           // interior Jacobian, global
GTVector<GFTYPE>            faceJac_;       // face Jacobian, global
GTVector<GTVector<GFTYPE>>  faceNormal_;    // normal to face at each node point (2d & 3d), global
GTVector<GFTYPE>            minnodedist_;   // min node length array (for each elem)
GTVector<GTVector<GSIZET>>  igbdy_;         // index into global field indicating a domain bdy
GTVector<GBdyType>          igbdytypes_;    // global domain bdy types for each igbdy index

};

#endif
