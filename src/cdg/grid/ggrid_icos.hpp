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

#include "ggrid.hpp"



// GICOS_BASE refers to the refined, projected triangular
//   'base' frame which are then partitioned into quad/hex elements
// GICOS_ELEMENTAL refers the computational grid, which is the 
//   final partitioned grid that the elements are constructed with:
// GICOS_CART: print in 3d cartesian coords;
// GICOS_LATLONG: print in r-theta-phi coords in 3d, and 
//    theta-phi in 2d


using namespace geoflow;
using namespace geoflow::pdeint;  
using namespace geoflow::tbox;
using namespace std;


template<typename TypePack>
class GGridIcos : public GGrid<TypePack>
{

public:
                             enum GICOSPTYPE      {GICOS_BASE, GICOS_ELEMENTAL}; 
                             enum GCOORDSYST      {GICOS_CART, GICOS_LATLONG}; 

                             using Types          = TypePack;
//                           using Grid           = typename Types::Grid;
                             using Mass           = typename Types::Mass;
                             using Ftype          = typename Types::Ftype;
                             using IBdyVol        = typename Types::IBdyVol;
                             using TBdyVol        = typename Types::TBdyVol;
                             using GElemList      = typename Types::GElemList;

                             using Operator       = typename Types::Operator;
                             using Preconditioner = typename Types::Preconditioner;
                             using ConnectivityOp = typename Types::ConnectivityOp;

                             using UpdateBase     = UpdateBdyBase<Types>;
                             using UpdateBasePtr  = std::shared_ptr<UpdateBase>;
                             using BdyUpdateList  = GTVector<GTVector<UpdateBasePtr>>;

                             using VVecFtype      = GTVector<GTVector<Ftype>>;

                             typedef GTMatrix<Ftype> GFTMatrix;
                             typedef Ftype GTICOS;

        // ICOS & sphere grid traits:
        struct Traits {
          GINT                ilevel;     // refine level if doing 2D ICOS
          Ftype               radiusi;    // inner radius (or just radius if doing ICOS)
          Ftype               radiuso;    // outer radius if doing 3D
          GTVector<GBdyType>  bdyTypes  ; // global bdy types (inner outer surf in 3D only)
        };

                            GGridIcos(const geoflow::tbox::PropertyTree &ptree, GTVector<GNBasis<GCTYPE,Ftype>*> &b, GC_COMM &comm);
#if 0
#endif
                           ~GGridIcos();

        void                do_elems();                                   // compute grid for irank
        void                do_elems(GTMatrix<GINT> &p,
                              GTVector<GTVector<Ftype>> &xnodes);        // compute elems from restart data)

        void                set_partitioner(GDD_base<GTICOS> *d);         // set and use GDD object
        GTVector<GTriangle<GTICOS>> 
                           &get_tmesh(){ return tmesh_;}                  // get complete triang. mesh
        GTVector    <GHex<GTICOS>> 
                           &get_hmesh(){ return hmesh_;}                  // get complete hex  mesh
        void                print(const GString &filename, 
                            GCOORDSYST icoord=GICOS_LATLONG);             // print grid to file
        void                config_gbdy(const geoflow::tbox::PropertyTree &ptree,
                            GBOOL                         bterrain,
                            GTVector<GTVector<GSIZET>>   &igbdyf,
                            GTVector<GTVector<GBdyType>> &igbdyft,
                            GTVector<GSIZET>             &igbdy,
                            GTVector<GUINT>              &degbdy);        // configure domain/global bdy
        void                elem_gbdy_data(
                              const GTMatrix<GTVector<Ftype>> &dXdXi,
                              GTVector<GSIZET>                 &igeface,
                              GTVector<Ftype>                 &face_mass,
                              GTVector<GTVector<Ftype>>       &normals);  // compute faces data

  std::size_t       max_duplicates() const;                          // Max duplicate points in grid

//friend  std::ostream&       operator<<(std::ostream&, GGridIcos<Types> &);         // Output stream operator
 

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
                              GTVector<GTVector<Ftype>> &xnodes); // do 2d grid restart
         void               do_elems3d(GTMatrix<GINT> &p,
                              GTVector<GTVector<Ftype>> &xnodes); // do 3d grid restart
         void               config_gbdy(const geoflow::tbox::PropertyTree &ptree,
                            GTVector<GTVector<GSIZET>>   &igbdyf,
                            GTVector<GTVector<GBdyType>> &igbdyt,
                            GTVector<GSIZET>             &igbdy);  // configure bdy
         void               find_gbdy_ind3d(Ftype radius,
                               GTVector<GSIZET> &igbdy,
                               GTVector<GUINT>  &debdy);                    // 
        void                elem_face_data(
                              GTMatrix<GTVector<Ftype>>       &dXdXi,
                              GTVector<GSIZET>                &gieface,
                              GTVector<Ftype>                 &face_mass,
                              GTVector<GTVector<Ftype>>       &normals); // compute elem face data entry point
        void                elem_face_data2d(
                              GTMatrix<GTVector<Ftype>>       &dXdXi,
                              GTVector<GSIZET>                &gieface,
                              GTVector<Ftype>                 &face_mass,
                              GTVector<GTVector<Ftype>>       &normals); // compute 2d elem face data
        void                elem_face_data3d(
                              GTMatrix<GTVector<Ftype>>       &dXdXi,
                              GTVector<GSIZET>                &gieface,
                              GTVector<Ftype>                 &face_mass,
                              GTVector<GTVector<Ftype>>       &normals); // compute 3d elem face data

         void               do_gbdy_normals(
                               const GTMatrix<GTVector<Ftype>> &dXdXi,
                               const GTVector<GSIZET>          &igbdy,
                               const GTVector<GUINT>           &debdy,
                               VVecFtype                       &normals,
                               GTVector<VVecFtype>             &tangents);
         void               do_gbdy_normals3d(
                               const GTMatrix<GTVector<Ftype>> &dXdXi,
                               const GTVector<GSIZET>          &igbdy,
                               const GTVector<GUINT>           &debdy,
                               VVecFtype                       &normals,
                               GTVector<VVecFtype>             &tangents);

         GString            sreftype_;      // subdivision/refinement type
         GINT               ilevel_;        // refinement level (>= 0)
         GINT               nrows_;         // # rows in refine level ilevel
         GINT               ndim_;          // grid dimensionality (2 or 3)
         GSIZET             nesurf_;        // # surface elements
         GSIZET             nnsurf_;        // # surface nodes
         GSIZET             nradelem_;      // # radial elements
         Ftype              radiusi_;       // inner radius
         Ftype              radiuso_;       // outer radius (=radiusi in 2d)
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
         GTVector<GNBasis<GCTYPE,Ftype>*> 
                             gbasis_;       // directional bases
         GTVector<GHex<GTICOS>>  
                             hmesh_;        // list of vertices for each 3d (hex) element
         GTMatrix<GTICOS>    fv0_;          // vertex list for base icosahedron
         GIMatrix            ifv0_;         // indices into fv0_ for each face of base icosahedron 

};


#include "ggrid_icos.ipp"



#endif
