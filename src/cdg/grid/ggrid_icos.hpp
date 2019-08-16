//==================================================================================
// Module       : ggrid_icos
// Date         : 8/31/18 (DLR)
// Description  : Object defining a (global) icosahedral grid, that in 2d
//                uses (extrinsic) gnomonic projections to locate element vertices.
//                Vertices always reside on sphere, so centroids will 
//                not (will reside within). In 3d, the base is computed from
//                the same procedure as in 2d, but we use isoparameteric
//                representation on the sphere.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : GGrid.
//==================================================================================
#if !defined(_GGRID_ICOS_HPP)
#define _GGRID_ICOS_HPP

#include "gtypes.h"
#include <functional>
#include "gtvector.hpp"
#include "gtmatrix.hpp"
#include "gnbasis.hpp"
#include "gelem_base.hpp"
#include "gdd_base.hpp"
#include "ggrid.hpp"
#include "gshapefcn_linear.hpp"
#include "gshapefcn_embed.hpp"
#include "polygon.h"

// GICOS_BASE refers to the refined, projected triangular
//   'base' frame which are then partitioned into quad/hex elements
// GICOS_ELEMENTAL refers the computational grid, which is the 
//   final partitioned grid that the elements are constructed with:
// GICOS_CART: print in 3d cartesian coords;
// GICOS_LATLONG: print in r-theta-phi coords in 3d, and 
//    theta-phi in 2d
enum GICOSPTYPE {GICOS_BASE, GICOS_ELEMENTAL}; 
enum GCOORDSYST {GICOS_CART, GICOS_LATLONG}; 

typedef GTMatrix<GFTYPE> GFTMatrix;

class GGridIcos : public GGrid
{
public:
        // ICOS & sphere grid traits:
        struct Traits {
          GINT                ilevel;     // refine level if doing 2D ICOS
          GFTYPE              radiusi;    // inner radius (or just radius if doing ICOS)
          GFTYPE              radiuso;    // outer radius if doing 3D
          GTVector<GBdyType>  bdyTypes  ; // global bdy types (inner outer surf in 3D only)
        };

                            GGridIcos(const geoflow::tbox::PropertyTree &ptree, GTVector<GNBasis<GCTYPE,GFTYPE>*> &b, GC_COMM &comm);
#if 0
#endif
                           ~GGridIcos();

        void                do_elems();                                   // compute grid for irank
        void                do_elems(GTMatrix<GINT> &p,
                              GTVector<GTVector<GFTYPE>> &xnodes);        // compute elems from restart data)
        void                set_partitioner(GDD_base *d);                 // set and use GDD object
        GTVector<GTriangle<GFTYPE>> 
                           &get_tmesh(){ return tmesh_;}                  // get complete triang. mesh
        GTVector    <GHex<GFTYPE>> 
                           &get_hmesh(){ return hmesh_;}                  // get complete hex  mesh
        void                print(const GString &filename, 
                            GCOORDSYST icoord=GICOS_LATLONG);             // print grid to file
        void                config_bdy(const PropertyTree &ptree,
                                       GTVector<GSIZET> &igbdy,
                                       GTVector<GSIZET> &igbdyt);         // config dy cond

friend  std::ostream&       operator<<(std::ostream&, GGridIcos &);       // Output stream operator
 

  private:
         void               init2d();                                       // initialize base icosahedron for 2d grid
         void               init3d();                                       // initialize for 3d grid
         void               project2sphere(GTVector<GTriangle<GFTYPE>> &, 
                                           GFTYPE rad);                     // project Cart mesh to sphere
         void               project2sphere(GTVector<GTPoint<GFTYPE>> &, 
                                           GFTYPE rad);                     // project points to sphere
         void               project2sphere(GTVector<GTVector<GFTYPE>> &, 
                                           GFTYPE rad);                     // project points to sphere
         void               spherical2xyz(GTVector<GTPoint<GFTYPE>*> &);    // (r,lat,long) to (x,y,z)
         void               spherical2xyz(GTVector<GTPoint<GFTYPE>>  &);    // (r,lat,long) to (x,y,z)
         void               spherical2xyz(GTVector<GTVector<GFTYPE>> &);    // (r,lat,long) to (x,y,z)
         void               xyz2spherical(GTVector<GTPoint<GFTYPE>*> &);    // (x,y,z) to (r, lat, long) 
         void               xyz2spherical(GTVector<GTVector<GFTYPE>> &);    // (x,y,z) to (r, lat, long) 
         void               xyz2spherical(GTVector<GTPoint<GFTYPE>>  &);    // (x,y,z) to (r, lat, long) 
         void               cart2gnomonic(GTVector<GTPoint<GFTYPE>> &, GFTYPE, GFTYPE, GFTYPE, 
                                          GTVector<GTPoint<GFTYPE>> &);     // transform to gnomonic space
         void               gnomonic2cart(GTVector<GTVector<GFTYPE>> &, GFTYPE, GFTYPE, GFTYPE, 
                                          GTVector<GTVector<GFTYPE>> &);    // transform from gnomonic space
         void               gnomonic2cart(GTVector<GTPoint<GFTYPE>> &, GFTYPE, GFTYPE, GFTYPE, 
                                          GTVector<GTPoint<GFTYPE>> &);     // transform from gnomonic space
         void               reorderverts2d(GTVector<GTPoint<GFTYPE>> &, GTVector<GSIZET>&,
                                           GTVector<GTPoint<GFTYPE>> &);    // make verts consis with shapefcns

         void               lagrefine();                                    // do 'Lagrange poly'-type refinement of base icos
         void               lagvert(GTPoint<GFTYPE> &a, 
                                    GTPoint<GFTYPE> &b, 
                                    GTPoint<GFTYPE> &c,
                                    GINT I, GTVector<GTPoint<GFTYPE>> &R);   // get list of points, R at index I

         void               interleave(GTVector<GTPoint<GFTYPE>> &R0,           // interleave rows of points for trianlges
                                    GTVector<GTPoint<GFTYPE>> &R1,
                                    GINT I, GTVector<GTPoint<GFTYPE>> &Rz);
         void               order_latlong2d(GTVector<GFPoint> &verts);       // order vertics via lat-long
         void               order_triangles(GTVector<GTriangle<GFTYPE>> &);    // order triangle verts

       
         void               do_elems2d(GINT rank);              // do 2d grid
         void               do_elems3d(GINT rank);              // do 3d grid
         void               do_elems2d(GTMatrix<GINT> &p,
                              GTVector<GTVector<GFTYPE>> &xnodes); // do 2d grid restart
         void               do_elems3d(GTMatrix<GINT> &p,
                              GTVector<GTVector<GFTYPE>> &xnodes); // do 3d grid restart
         void               config_bdy(const PropertyTree &ptree,
                            GTVector<GTVector<GSIZET>>   &igbdy,
                            GTVector<GTVector<GBdyType>> &igbdyt); // configure bdy
         void               find_bdy_ind3d(GFTYPE radius, 
                                           GTVector<GSIZET> &ibdy);// find bdy indices for bdy=radius

         GINT               ilevel_;        // refinement level (>= 0)
         GINT               ndim_;          // grid dimensionality (2 or 3)
         GSIZET             nradelem_;      // # radial elements
         GFTYPE             radiusi_;       // inner radius
         GFTYPE             radiuso_;       // outer radius (=radiusi in 2d)
         GDD_base          *gdd_;           // domain decomposition/partitioning object
         GShapeFcn_linear  *lshapefcn_;     // linear shape func to compute 2d coords
         GTVector<GINT>     iup_;           // triangle pointing 'up' flag

         GTVector<GTriangle<GFTYPE>>    
                            tmesh_;         // array of final mesh triangles
         GTVector<GTPoint<GFTYPE>>
                            ftcentroids_ ;  // centroids of finest triangles/faces/ or hexes
         GTVector<GTriangle<GFTYPE>>     
                             tbase_;        // array of base triangles
         GTVector<GNBasis<GCTYPE,GFTYPE>*> 
                             gbasis_;       // directional bases
         GTVector<GHex<GFTYPE>>  
                             hmesh_;        // list of vertices for each 3d (hex) element
         GTMatrix<GFTYPE>    fv0_;          // vertex list for base icosahedron
         GIMatrix            ifv0_;         // indices into fv0_ for each face of base icosahedron 

};

#endif
