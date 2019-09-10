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
#include "gcomm.hpp"
#include "gtvector.hpp"
#include "gtmatrix.hpp"
#include "gnbasis.hpp"
#include "gelem_base.hpp"
#include "gdd_base.hpp"
#include "ggrid.hpp"
#include "gcomm.hpp"
#include "ggfx.hpp"
#include "tbox/property_tree.hpp"


using namespace geoflow::tbox;
using namespace std;

class GMass;

typedef GTVector<GElem_base*> GElemList;


class GGrid 
{
public:
                             GGrid() = delete;
                             GGrid(const geoflow::tbox::PropertyTree &ptree, GTVector<GNBasis<GCTYPE,GFTYPE>*> &b, GC_COMM &comm);

virtual                       ~GGrid();
//static                       GGrid *build(geoflow::tbox::PropertyTree &ptree, GTVector<GNBasis<GCTYPE,GFTYPE>*> &b, GC_COMM comm);

virtual void                 do_elems() = 0;                            // compute grid for irank
virtual void                 do_elems(GTMatrix<GINT> &p,
                               GTVector<GTVector<GFTYPE>> &xnodes) = 0; // compute grid on restart
virtual void                 set_partitioner(GDD_base *d) = 0;         // set and use GDD object

#if 0
virtual void                 set_bdy_callback(
                             std::function<void(GElemList &)> *callback) 
                             {bdycallback_ =  callback; }              // set bdy-set callback
#endif

virtual void                 print(const GString &filename){}          // print grid to file


        void                 grid_init();                             // initialize class
        void                 grid_init(GTMatrix<GINT> &p, 
                               GTVector<GTVector<GFTYPE>> &xnodes);   // initialize class for restart
        void                 do_typing(); // classify into types
        GElemList           &elems() { return gelems_; }              // get elem list
        GSIZET               nelems() { return gelems_.size(); }      // local num elems
        GTVector<GSIZET> 
                            &ntype() { return ntype_; }               // no. elems of each type
        GTVector<GTVector<GSIZET>> 
                            &itype() { return itype_; }               // indices for all types
        GTVector<GSIZET>    &itype(GElemType i) { return itype_[i]; } // indices for type i    
        GElemType            gtype() { return gtype_; }               // get unique elem type on grid       
        void                 deriv(GTVector<GFTYPE> &u, GINT idir, GTVector<GFTYPE> &tmp,
                                   GTVector<GFTYPE> &du );            // derivative of global vector
        GFTYPE               integrate(GTVector<GFTYPE> &u,
                                       GTVector<GFTYPE> &tmp, 
				       GBOOL bglobal = TRUE);         // spatial integration of global vector
        void                 print(GString fname, GBOOL bdof=FALSE);
        GSIZET               ndof();                                  // compute total number elem dof
        GSIZET               size() { return ndof(); }
        GSIZET               nfacedof();                              // compute total number elem face dof
        GSIZET               nbdydof();                               // compute total number elem bdy dof
        GFTYPE               minlength();                             // find min elem length
        GFTYPE               maxlength();                             // find max elem length
        GFTYPE               minnodedist()         
                             {return minnodedist_;}                   // find min node distance
        GFTYPE               volume()         
                             {return volume_;}                        // get grid volume
        GFTYPE               ivolume()         
                             {return ivolume_;}                       // get nverse of grid volume
        GTMatrix<GTVector<GFTYPE>>
                            &dXidX();                                 // global Rij = dXi^j/dX^i
        GTVector<GFTYPE>    &dXidX(GSIZET i,GSIZET j);                // Rij matrix element 
        GTVector<GTVector<GFTYPE>>
                            &xNodes() { return xNodes_; }             // get all nodes coords (Cart)
        GTVector<GFTYPE>    &xNodes(GSIZET i) { return xNodes_[i]; }  // get all nodes coords (Cart)
                            

        GMass               &massop();                                 // global mass op
        GTVector<GFTYPE>    &Jac();                                    // global Jacobian
        GTVector<GFTYPE>
                            &faceJac();                                // global face Jacobian
        GTVector<GTVector<GFTYPE>>
                            &faceNormal();                             // global face normals
        GTVector<GTVector<GSIZET>>
                            &igface() { return igface_;}               // global dom face indices into u for each elem face index
        GTVector<GTVector<GFTYPE>>
                            &bdyNormal();                              // global bdy normals
        GTVector<GTVector<GSIZET>>
                            &igbdy_binned() { return igbdy_binned_;}   // global dom bdy indices binned into GBdyType
        GTVector<GTVector<GSIZET>>
                            &igbdy_byface() { return igbdy_byface_;}   // global dom bdy indices for each face
        GTVector<GTVector<GBdyType>>
                            &igbdyt_byface(){ return igbdyt_byface_;}  // global dom bdy indices for each face
        GTVector<GSIZET>
                            &igbdy() { return igbdy_;}                 // global dom bdy indices into u

        GC_COMM              get_comm() { return comm_; }              // get communicator

virtual void                 config_bdy(const PropertyTree &ptree, 
                             GTVector<GTVector<GSIZET>>   &igbdy, 
                             GTVector<GTVector<GBdyType>> &igbdyt)=0;  // config bdy

        GGFX<GFTYPE>        &get_ggfx() { return *ggfx_; }             // get GGFX op
        void                 set_ggfx(GGFX<GFTYPE> &ggfx) 
                               { ggfx_ = &ggfx; }                      // set GGFX op    

friend  std::ostream&        operator<<(std::ostream&, GGrid &);       // Output stream operator
 

protected:
       
        void                        init_local_face_info();           // get local face info
        void                        init_bc_info();                   // configure bdys
        void                        def_init();                       // iniitialze deformed elems
        void                        reg_init();                       // initialize regular elems
        GFTYPE                      find_min_dist(); 

        GBOOL                       bInitialized_;  // object initialized?
        GBOOL                       is_bdy_time_dep_; // time-dep bdy vals?
        GElemType                   gtype_;         // element types comprising grid
        GINT                        irank_;         // MPI task id
        GINT                        nprocs_;        // number of MPI tasks
        GC_COMM                     comm_;          // communicator
        GFTYPE                      minnodedist_;   // min node length array (for each elem)
	GFTYPE                      volume_;        // grid volume
	GFTYPE                      ivolume_;       // 1 / grid volume
        GElemList                   gelems_;        // element list
        GTVector<GFTYPE>            etmp_;          // elem-level tmp vector
        GTVector<GTVector<GSIZET>>  itype_;         // indices in elem list of each type
        GTVector<GSIZET>            ntype_;         // no. elems of each type on grid
        GTMatrix<GTVector<GFTYPE>>  dXidX_;         // matrix Rij = dXi^j/dX^i, global
        GTVector<GTVector<GFTYPE>>  xNodes_;        // Cart coords of all node points
        GMass                      *mass_;          // mass operator
        GTVector<GFTYPE>            Jac_;           // interior Jacobian, global
        GTVector<GFTYPE>            faceJac_;       // face Jacobian, global
        GTVector<GTVector<GFTYPE>>  faceNormal_;    // normal to eleme faces each face node point (2d & 3d), global
        GTVector<GTVector<GSIZET>>  igface_;        // index into global field indicating elem face node
        GTVector<GTVector<GFTYPE>>  bdyNormal_;     // normal to surface at each bdy node point (2d & 3d), global
        GTVector<GTVector<GSIZET>>  igbdy_binned_;  // index into global field indicating a domain bdy--by type
        GTVector<GSIZET>            igbdy_;         // index into global field indicating a domain bdy
        GTVector<GTVector<GSIZET>>  igbdy_byface_;  // index into global field indicating a domain bdy
        GTVector<GBdyType>          igbdyt_;        // global domain bdy types for each igbdy index
        GTVector<GTVector<GBdyType>>
                                    igbdyt_byface_; // global domain bdy types for each igbdy index
        PropertyTree                ptree_;         // main prop tree
        GGFX<GFTYPE>               *ggfx_;          // connectivity operator

};

#endif
