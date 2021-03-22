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

typedef GTMatrix<GFTYPE> GFTMatrix;

using namespace geoflow::pdeint;
using namespace std;

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
        void                set_partitioner(GDD_base<GFTYPE> *d);            // set and use GDD object
        void                set_basis(GTVector<GNBasis<GCTYPE,GFTYPE>*> &b); // set element basis
        void                periodize();                                     // periodize coords, if allowed
        void                unperiodize();                                   // un-periodize coords, if allow
        void                config_gbdy(const geoflow::tbox::PropertyTree &ptree,
                              GBOOL                         bterrain,
                              GTVector<GTVector<GSIZET>>   &igbdyf, 
                              GTVector<GTVector<GBdyType>> &igbdyt,
                              GTVector<GSIZET>             &igbdy,
                              GTVector<GUINT>              &degbdy);         // config bdy
        void                elem_face_data(
                              const GTMatrix<GTVector<GFTYPE>> &dXdXi,
                              GTVector<GSIZET>                 &gieface,
                              GTVector<GFTYPE>                 &face_mass,
                              GTVector<GTVector<GFTYPE>>       &normals);    // compute elem face data
const    GTPoint<GFTYPE>   &getP0() {return P0_; }                           // get blob bdy point 
const    GTPoint<GFTYPE>   &getP1() {return P1_; }                           // get blob bdy point 

         void               print(const GString &filename);                  // print grid to file


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
                               const GTMatrix<GTVector<GFTYPE>> &dXdXi,
                               const GTVector<GSIZET>           &igbdy,
                               const GTVector<GUINT>            &degbdy,
                                     GTVector<GTVector<GFTYPE>> &normals,
                               GTVector<GINT>                   &idepComp); // compute normals entry point
         void               do_gbdy_normals2d(
                              const GTMatrix<GTVector<GFTYPE>> &dXdXi,
                              const GTVector<GSIZET>           &igbdy,
                              const GTVector<GUINT>            &degbdy,
                              GTVector<GTVector<GFTYPE>>       &normals,
                              GTVector<GINT>                   &idepComp);  // compute normals to doimain bdy in 2d
         void               do_gbdy_normals3d(
                              const GTMatrix<GTVector<GFTYPE>> &dXdXi,
                              const GTVector<GSIZET>           &igbdy,
                              const GTVector<GUINT>            &degbdy,
                              GTVector<GTVector<GFTYPE>>       &normals,
                              GTVector<GINT>                   &idepComp);       // compute normals to doimain bdy in 3d
        void                find_rank_subdomain();                               // find task's default subdomain

        void                elem_face_data2d(
                              const GTMatrix<GTVector<GFTYPE>> &dXdXi,
                              GTVector<GSIZET>                 &gieface,
                              GTVector<GFTYPE>                 &face_mass,
                              GTVector<GTVector<GFTYPE>>       &normals);          // compute 2d elem face data
        void                elem_face_data3d(
                              const GTMatrix<GTVector<GFTYPE>> &dXdXi,
                              GTVector<GSIZET>                 &gieface,
                              GTVector<GFTYPE>                 &face_mass,
                              GTVector<GTVector<GFTYPE>>       &normals);          // compute 3d elem face data
         GBOOL              is_global_vertex(GTPoint<GFTYPE> &pt);          // pt on global vertex?
         GBOOL              on_global_edge(GINT iface, GTPoint<GFTYPE> &pt);// pt on global edge?



         GINT                ndim_;          // grid dimensionality (2 or 3)
         GFTYPE              eps_;           // float epsilon for comparisons
         GDD_base<GFTYPE>    *gdd_;           // domain decomposition/partitioning object
         GShapeFcn_linear<GFTYPE> 
                            *lshapefcn_;     // linear shape func to compute 2d coords
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
