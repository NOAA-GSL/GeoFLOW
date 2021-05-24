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

#include "ggrid.hpp"


using namespace geoflow;
using namespace geoflow::pdeint;
using namespace geoflow::tbox;
using namespace std;


template<typename TypePack>
class GGridBox : public GGrid<TypePack>
{

public:

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



        // Box grid traits:
        struct Traits {
          GTPoint<Ftype>      P0;        // global lower point
          GTPoint<Ftype>      P1;        // global upper point
          GTVector<GBdyType>  bdyTypes;  // global bdy types
        };

                            GGridBox(const geoflow::tbox::PropertyTree &ptree, GTVector<GNBasis<GCTYPE,Ftype>*> &b, GC_COMM &comm);

                           ~GGridBox();

        void                do_elems();                                      // compute elems
        void                do_elems(GTMatrix<GINT> &p,
                              GTVector<GTVector<Ftype>> &xnodes);           // compute elems from restart data
        void                set_partitioner(GDD_base<Ftype> *d);            // set and use GDD object
        void                set_basis(GTVector<GNBasis<GCTYPE,Ftype>*> &b); // set element basis
        void                periodize();                                     // periodize coords, if allowed
        void                unperiodize();                                   // un-periodize coords, if allow
        void                config_gbdy(const geoflow::tbox::PropertyTree &ptree,
                              GBOOL                         bterrain,
                              GTVector<GTVector<GSIZET>>   &igbdyf, 
                              GTVector<GTVector<GBdyType>> &igbdyt,
                              GTVector<GSIZET>             &igbdy,
                              GTVector<GUINT>              &degbdy);         // config bdy
        void                elem_face_data(
                              GTMatrix<GTVector<Ftype>>       &dXdXi,
                              GTVector<GSIZET>                 &gieface,
                              GTVector<Ftype>                 &face_mass,
                              GTVector<GTVector<Ftype>>       &normals);    // compute elem face data
const    GTPoint<Ftype>   &getP0() {return P0_; }                           // get blob bdy point 
const    GTPoint<Ftype>   &getP1() {return P1_; }                           // get blob bdy point 

         std::size_t       max_duplicates() const;                          // Max duplicate points in grid
         void               print(const GString &filename);                  // print grid to file


#if 0
friend  std::ostream&       operator<<(std::ostream&, GGridBox<Types> &);           // Output stream operator
#endif
 

private:
         void               init2d();                                       // initialize base icosahedron for 2d grid
         void               init3d();                                       // initialize for 3d grid


         void               do_elems2d();                                   // do 2d grid
         void               do_elems3d();                                   // do 3d grid
         void               do_elems2d(GTMatrix<GINT> &p, 
                              GTVector<GTVector<Ftype>> &xnodes);          // do 2d grid restart
         void               do_elems3d(GTMatrix<GINT> &p, 
                              GTVector<GTVector<Ftype>> &xnodes);          // do 3d grid restart
         void               find_gbdy_ind(GINT bdyid, GBOOL bunique,
                               GTVector<GSIZET> &ikeep,
                               GTVector<GSIZET> &ibdy,
                               GTVector<GUINT>  &debdy);                    // find glob indices for domain bdys
         void               find_gbdy_ind2d(GINT bdyid, GBOOL bunique,
                               GTVector<GSIZET> &ikeep,
                               GTVector<GSIZET> &ibdy,
                               GTVector<GUINT>  &debdy);                    // find glob indices for 2d domain bdy

         void               find_gbdy_ind3d(GINT bdyid, GBOOL bunique,
                               GTVector<GSIZET> &ikeep,
                               GTVector<GSIZET> &ibdy,
                               GTVector<GUINT>  &debdy);                    // find glob indices for 3d domain bdy

        void                do_gbdy_normals  (
                               const GTMatrix<GTVector<Ftype>> &dXdXi,
                               const GTVector<GSIZET>          &igbdy,
                               const GTVector<GUINT>           &degbdy,
                                     GTVector<GTVector<Ftype>> &normals,
                               GTVector<GINT>                  &idepComp); // compute normals entry point
         void               do_gbdy_normals2d(
                              const GTMatrix<GTVector<Ftype>> &dXdXi,
                              const GTVector<GSIZET>          &igbdy,
                              const GTVector<GUINT>           &degbdy,
                              GTVector<GTVector<Ftype>>       &normals,
                              GTVector<GINT>                  &idepComp);  // compute normals to doimain bdy in 2d
         void               do_gbdy_normals3d(
                              const GTMatrix<GTVector<Ftype>> &dXdXi,
                              const GTVector<GSIZET>          &igbdy,
                              const GTVector<GUINT>           &degbdy,
                              GTVector<GTVector<Ftype>>       &normals,
                              GTVector<GINT>                  &idepComp);       // compute normals to doimain bdy in 3d
        void                find_rank_subdomain();                               // find task's default subdomain

        void                elem_face_data2d(
                              GTMatrix<GTVector<Ftype>>       &dXdXi,
                              GTVector<GSIZET>                &gieface,
                              GTVector<Ftype>                 &face_mass,
                              GTVector<GTVector<Ftype>>       &normals);          // compute 2d elem face data
        void                elem_face_data3d(
                              GTMatrix<GTVector<Ftype>>       &dXdXi,
                              GTVector<GSIZET>                &gieface,
                              GTVector<Ftype>                 &face_mass,
                              GTVector<GTVector<Ftype>>       &normals);          // compute 3d elem face data
         GBOOL              is_global_vertex(GTPoint<Ftype> &pt);          // pt on global vertex?
         GBOOL              on_global_edge(GINT iface, GTPoint<Ftype> &pt);// pt on global edge?



         GINT                ndim_;          // grid dimensionality (2 or 3)
         GDD_base<Ftype>    *gdd_;           // domain decomposition/partitioning object
         GShapeFcn_linear<Ftype> 
                            *lshapefcn_;     // linear shape func to compute 2d coords
         GTPoint<Ftype>      P0_;            // P0 = starting point of box origin
         GTPoint<Ftype>      P1_;            // P1 = diagonally-opposing box point
         GTPoint<Ftype>      dP_;            // box size
         GTVector<GTPoint<Ftype>>
                             gverts_;        // global bdy vertices
         GTVector<GTPoint<Ftype>>
                             ftcentroids_;   // centroids of finest elements
         GTVector<GNBasis<GCTYPE,Ftype>*> 
                             gbasis_;        // directional bases
         GTVector<GQuad<Ftype>> 
                             qmesh_;         // list of vertices for each 2d (quad) element
         GTVector<GHex<Ftype>>  
                             hmesh_;         // list of vertices for each 3d (hex) element
         GTVector<GINT>      ne_;            // # elems in each coord direction in 3d
         GTVector<GBOOL>     bPeriodic_;     // is periodic in x, y, or z?
         GTVector<GSIZET>    periodicids_;   // node ids that were made periodic in 1 or more dirs
         GTVector<GUINT>     periodicdirs_;  // integer with bits 1, 2, 3 set for direction of periodiscity
         GTVector<Ftype>     Lbox_;          // length of box edges (x, y, [and z])


};


#include "ggrid_box.ipp"

#endif
