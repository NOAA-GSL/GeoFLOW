//==================================================================================
// Module       : ggrid_icos.hpp
// Date         : 8/31/18 (DLR)
// Description  : Object defining a (global) icosahedral grid, that in 2d
//                uses (extrinsic) gnomonic projections to locate element vertices.
//                Vertices always reside on sphere, so centroids will 
//                not (will reside within). In 3d, the base is computed from
//                the same procedure as in 2d, but we use isoparameteric
//                representation on the sphere.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved.
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
#include "gshapefcn_linear.hpp"
#include "gshapefcn_embed.hpp"
#include "polygon.h"
#include "ggrid.hpp"
#include "gmtk.hpp"
#include "gutils.hpp"


// GICOS_BASE refers to the refined, projected triangular
//   'base' frame which are then partitioned into quad/hex elements
// GICOS_ELEMENTAL refers the computational grid, which is the 
//   final partitioned grid that the elements are constructed with:
// GICOS_CART: print in 3d cartesian coords;
// GICOS_LATLONG: print in r-theta-phi coords in 3d, and 
//    theta-phi in 2d
enum GICOSPTYPE      {GICOS_BASE, GICOS_ELEMENTAL}; 
enum GCOORDSYST      {GICOS_CART, GICOS_LATLONG}; 

typedef GTMatrix<GFTYPE> GFTMatrix;
typedef GFTYPE GTICOS;

using namespace geoflow::pdeint;
using namespace std;

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
        void                do_face_normals();                            // compute normals to elem faces
        void                do_bdy_normals(
                              GTMatrix<GTVector<GFTYPE>>    &dXdXi,
                              GTVector<GTVector<GSIZET>>    &igbdy_face,
                              GTVector<GTVector<GFTYPE>>    &normals,
                              GTVector<GINT>                &idepComp);   // compute normals to domain bdy

        void                set_partitioner(GDD_base<GTICOS> *d);         // set and use GDD object
        GTVector<GTriangle<GTICOS>> 
                           &get_tmesh(){ return tmesh_;}                  // get complete triang. mesh
        GTVector    <GHex<GTICOS>> 
                           &get_hmesh(){ return hmesh_;}                  // get complete hex  mesh
        void                print(const GString &filename, 
                            GCOORDSYST icoord=GICOS_LATLONG);             // print grid to file
        void                config_bdy(const geoflow::tbox::PropertyTree &ptree,
                                       GTVector<GSIZET> &igbdy,
                                       GTVector<GSIZET> &igbdyt);         // config dy cond

friend  std::ostream&       operator<<(std::ostream&, GGridIcos &);       // Output stream operator
 

  private:
         void               init2d();                                       // initialize base icosahedron for 2d grid
         void               init3d();                                       // initialize for 3d grid
         template<typename T>
         void               project2sphere(GTVector<GTriangle<T>> &, 
                                           T rad);                          // project Cart mesh to sphere
         template<typename T>
         void               project2sphere(GTVector<GTPoint<T>> &, 
                                           T rad);                          // project points to sphere
         template<typename T>
         void               project2sphere(GTVector<GTVector<T>> &, 
                                           T  rad);                        // project points to sphere
         template<typename T>
         void               spherical2xyz(GTVector<GTPoint<T>*> &);        // (r,lat,long) to (x,y,z)
         template<typename T>
         void               spherical2xyz(GTVector<GTPoint<T>>  &);        // (r,lat,long) to (x,y,z)
         template<typename T>
         void               spherical2xyz(GTVector<GTVector<T>> &);        // (r,lat,long) to (x,y,z)
         template<typename T>
         void               xyz2spherical(GTVector<GTPoint<T>*> &);        // (x,y,z) to (r, lat, long) 
         template<typename T>
         void               xyz2spherical(GTVector<GTVector<T>> &);         // (x,y,z) to (r, lat, long) 
         template<typename T>
         void               xyz2spherical(GTVector<GTPoint<T>>  &);         // (x,y,z) to (r, lat, long) 
         template<typename T>
         void               cart2gnomonic(GTVector<GTPoint<T>> &, T, T, T,
                                          GTVector<GTPoint<T>> &);          // transform to gnomonic space
         template<typename T>
         void               gnomonic2cart(GTVector<GTVector<T>> &, T, T, T,
                                          GTVector<GTVector<T>> &);         // transform from gnomonic space
         template<typename T>
         void               gnomonic2cart(GTVector<GTPoint<T>> &, T, T, T, 
                                          GTVector<GTPoint<T>> &);          // transform from gnomonic space
         template<typename T>
         void               reorderverts2d(GTVector<GTPoint<T>> &, 
                                           GTVector<GTPoint<T>>&, 
                                           GTVector<GSIZET>&,
                                           GTVector<GTPoint<T>> &);         // make verts consis with shapefcns
         template<typename TF, typename TT>
         void               copycast(GTVector<GTVector<TF>> &from, GTVector<GTVector<TT>> &to);
         template<typename TF, typename TT>
         void               copycast(GTVector<GTVector<TF>*> &from, GTVector<GTVector<TT>*> &to);
         template<typename TF, typename TT>
         void               copycast(GTPoint<TF> &from, GTPoint<TT> &to);

         void               lagrefine();                                    // do 'Lagrange poly'-type refinement of base icos
         template<typename T>
         void               lagvert(GTPoint<T> &a, 
                                    GTPoint<T> &b, 
                                    GTPoint<T> &c,
                                    GINT I, GTVector<GTPoint<T>> &R);      // get list of points, R at index I

         template<typename T>
         void               interleave(GTVector<GTPoint<T>> &R0,           // interleave rows of points for trianlges
                                    GTVector<GTPoint<T>> &R1,
                                    GINT I, GTVector<GTPoint<T>> &Rz);
         template<typename T>
         void               order_latlong2d(GTVector<GTPoint<T>> &verts);  // order vertics via lat-long
         template<typename T>
         void               order_triangles(GTVector<GTriangle<T>> &);     // order triangle verts

       
         void               do_elems2d(GINT rank);              // do 2d grid
         void               do_elems3d(GINT rank);              // do 3d grid
         void               do_elems2d(GTMatrix<GINT> &p,
                              GTVector<GTVector<GFTYPE>> &xnodes); // do 2d grid restart
         void               do_elems3d(GTMatrix<GINT> &p,
                              GTVector<GTVector<GFTYPE>> &xnodes); // do 3d grid restart
         void               config_bdy(const geoflow::tbox::PropertyTree &ptree,
                            GTVector<GTVector<GSIZET>>   &igbdy,
                            GTVector<GTVector<GBdyType>> &igbdyt); // configure bdy
         void               find_bdy_ind3d(GFTYPE radius, 
                                           GTVector<GSIZET> &ibdy);// find bdy indices for bdy=radius
         void               do_face_normals2d();                   // compute normals to elem faces in 2d
         void               do_face_normals3d();                   // compute normals to elem faces in 2d
         void               do_bdy_normals3d(
                              GTMatrix<GTVector<GFTYPE>>    &dXdXi,
                              GTVector<GSIZET>              &igbdy,
                              GINT                           iface,
                              GTVector<GTVector<GFTYPE>>    &normals,
                              GTVector<GINT>              &idepComp);// compute normals to doimain bdy in 3d

         GString            sreftype_;      // subdivision/refinement type
         GINT               ilevel_;        // refinement level (>= 0)
         GINT               nrows_;         // # rows in refine level ilevel
         GINT               ndim_;          // grid dimensionality (2 or 3)
         GSIZET             nradelem_;      // # radial elements
         GFTYPE             radiusi_;       // inner radius
         GFTYPE             radiuso_;       // outer radius (=radiusi in 2d)
         GDD_base<GTICOS>  *gdd_;           // domain decomposition/partitioning object
         GShapeFcn_linear<GTICOS>
                           *lshapefcn_;     // linear shape func to compute 2d coords
         GTVector<GINT>     iup_;           // triangle pointing 'up' flag

         GTVector<GTriangle<GTICOS>>    
                            tmesh_;         // array of final mesh triangles
         GTVector<GTPoint<GTICOS>>
                            ftcentroids_ ;  // centroids of finest triangles/faces/ or hexes
         GTVector<GTriangle<GTICOS>>     
                             tbase_;        // array of base triangles
         GTVector<GNBasis<GCTYPE,GFTYPE>*> 
                             gbasis_;       // directional bases
         GTVector<GHex<GTICOS>>  
                             hmesh_;        // list of vertices for each 3d (hex) element
         GTMatrix<GTICOS>    fv0_;          // vertex list for base icosahedron
         GIMatrix            ifv0_;         // indices into fv0_ for each face of base icosahedron 

};

#include "ggrid_icos.ipp"

#endif
