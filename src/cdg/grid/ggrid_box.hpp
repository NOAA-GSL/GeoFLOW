//==================================================================================
// Module       : ggrid_box
// Date         : 11/11/18 (DLR)
// Description  : Object defining a (global) 2d or 3d box grid. Builds
//                elements that base class then uses to build global
//                computational data structures.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : GGrid.
//==================================================================================
#if !defined(_GGRID_BOX_HPP)
#define _GGRID_BOX_HPP

#include "gtypes.h"
#include <functional>
#include "gtvector.hpp"
#include "gtmatrix.hpp"
#include "gnbasis.hpp"
#include "gelem_base.hpp"
#include "gdd_base.hpp"
#include "ggrid.hpp"
#include "gshapefcn_linear.hpp"
#include "polygon.h"
#include "gtpoint.hpp"
#include "ggrid.hpp"

typedef GTMatrix<GFTYPE> GFTMatrix;

class GGridBox : public GGrid
{

public:
        // Box grid traits:
        struct Traits {
          GTPoint<GFTYPE>     P0;        // global lower point
          GTPoint<GFTYPE>     P1;        // global upper point
          GTVector<GBdyType>  bdyTypes;  // global bdy types
        };

                            GGridBox(const geoflow::tbox::PropertyTree &ptree, GTVector<GNBasis<GCTYPE,GFTYPE>*> &b, GC_COMM &comm);

                           ~GGridBox();

        void                do_elems();                                      // compute elems
        void                do_elems(GTMatrix<GINT> &p,
                              GTVector<GTVector<GFTYPE>> &xnodes);           // compute elems from restart data
        void                do_face_normals();                               // compute normals to elem faces
        void                do_bdy_normals ();                               // compute normals to doimain bdy
        void                set_partitioner(GDD_base *d);                    // set and use GDD object
        void                set_basis(GTVector<GNBasis<GCTYPE,GFTYPE>*> &b); // set element basis
        void                periodize();                                     // periodize coords, if allowed
        void                unperiodize();                                   // un-periodize coords, if allow
         void               config_bdy(const PropertyTree &ptree,
                            GTVector<GTVector<GSIZET>>   &igbdy, 
                            GTVector<GTVector<GBdyType>> &igbdyt);          // config bdy

        void                print(const GString &filename);                  // print grid to file

friend  std::ostream&       operator<<(std::ostream&, GGridBox &);           // Output stream operator
 

private:
         void               init2d();                                       // initialize base icosahedron for 2d grid
         void               init3d();                                       // initialize for 3d grid


         void               do_elems2d();                                   // do 2d grid
         void               do_elems3d();                                   // do 3d grid
         void               do_elems2d(GTMatrix<GINT> &p, 
                              GTVector<GTVector<GFTYPE>> &xnodes);          // do 2d grid restart
         void               do_elems3d(GTMatrix<GINT> &p, 
                              GTVector<GTVector<GFTYPE>> &xnodes);          // do 3d grid restart
         void               find_subdomain();                               // find task's default subdomain
         void               find_bdy_ind2d(GINT, GBOOL, GTVector<GSIZET> &);// find global bdy indices for specified bdy in 2d
         void               find_bdy_ind3d(GINT, GBOOL, GTVector<GSIZET> &);// find global bdy indices for specified bdy in 3d
         GBOOL              is_global_vertex(GTPoint<GFTYPE> &pt);          // pt on global vertex?
         GBOOL              on_global_edge(GINT iface, GTPoint<GFTYPE> &pt);// pt on global edge?

         GINT                ndim_;          // grid dimensionality (2 or 3)
         GFTYPE              eps_;           // float epsilon for comparisons
         GDD_base           *gdd_;           // domain decomposition/partitioning object
         GShapeFcn_linear   *lshapefcn_;     // linear shape func to compute 2d coords
         GTPoint<GFTYPE>     P0_;            // P0 = starting point of box origin
         GTPoint<GFTYPE>     P1_;            // P1 = diagonally-opposing box point
         GTPoint<GFTYPE>     dP_;            // box size
         GTVector<GTPoint<GFTYPE>>
                             gverts_;        // global bdy vertices
         GTVector<GTPoint<GFTYPE>>
                             ftcentroids_;   // centroids of finest elements
         GTVector<GNBasis<GCTYPE,GFTYPE>*> 
                             gbasis_;        // directional bases
         GTVector<GQuad<GFTYPE>> 
                             qmesh_;         // list of vertices for each 2d (quad) element
         GTVector<GHex<GFTYPE>>  
                             hmesh_;         // list of vertices for each 3d (hex) element
         GTVector<GINT>      ne_;            // # elems in each coord direction in 3d
         GTVector<GBOOL>     bPeriodic_;     // is periodic in x, y, or z?
         GTVector<GSIZET>    periodicids_;   // node ids that were made periodic in 1 or more dirs
         GTVector<GUINT>     periodicdirs_;  // integer with bits 1, 2, 3 set for direction of periodiscity
         GTVector<GFTYPE>    Lbox_;          // length of box edges (x, y, [and z])

};

#endif
