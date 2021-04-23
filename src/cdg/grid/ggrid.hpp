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


#include "gtypes.h"
#include <iostream>
#include <memory>
#include <cmath>
#include <limits>
#include <typeinfo>
#include "gcomm.hpp"
#include "gtvector.hpp"
#include "gtmatrix.hpp"
#include "gnbasis.hpp"
#include "gllbasis.hpp"
#include "gelem_base.hpp"
#include "gdd_base.hpp"
#include "gcblas.hpp"
#include "gutils.hpp"
#include "gcg.hpp"
#include "gmtk.hpp"
#include "ghelmholtz.hpp"
#include "glinop_base.hpp"

#include "polygon.h"
#include "gshapefcn_linear.hpp"
#include "gshapefcn_embed.hpp"
#include "gupdatebdy_factory.hpp"

#include "pdeint/update_bdy_base.hpp"
#include "pdeint/lin_solver_base.hpp"
#include "tbox/error_handler.hpp"
#include "tbox/tracer.hpp"
#include "tbox/property_tree.hpp"


using namespace geoflow;
using namespace geoflow::tbox;
using namespace geoflow::pdeint;
using namespace std;



template<typename TypePack>
class GGrid 
{
public:
                             enum GDerivType {GDV_VARP=0, GDV_CONSTP}; 
                             struct CGTypePack { // define terrain typepack
                                     using Operator         = class GHelmholtz<TypePack>;
                                     using Preconditioner   = GLinOpBase<TypePack>;
                                     using State            = typename TypePack::State;
                                     using StateComp        = typename TypePack::StateComp;
                                     using Grid             = GGrid<TypePack>;
                                     using Ftype            = typename TypePack::Ftype;
                                     using ConnectivityOp   = GGFX<GFTYPE>;
                             };

                             using Types          = TypePack;
                             using EqnBase        = typename Types::EqnBase;
                             using EqnBasePtr     = typename Types::EqnBasePtr;
                             using State          = typename Types::State;
                             using StateComp      = typename Types::StateComp;
//                           using Grid           = typename Types::Grid;
                             using Mass           = typename Types::Mass;
                             using Ftype          = typename Types::Ftype;
                             using Size           = typename Types::Size;
                             using Derivative     = typename Types::Derivative;
                             using Time           = typename Types::Time;
                             using CompDesc       = typename Types::CompDesc;
                             using Jacobian       = typename Types::Jacobian;

                             using IBdyVol        = GTVector<GSIZET>;
                             using TBdyVol        = GTVector<GBdyType>;
                             using GElemList      = GTVector<GElem_base*>;

                             using TerrainOp      = typename CGTypePack::Operator;
                             using TerrainPrecond = typename CGTypePack::Preconditioner;
                             using ConnectivityOp = GGFX<Ftype>;

                             using UpdateBase     = UpdateBdyBase<Types>;
                             using UpdateBasePtr  = std::shared_ptr<UpdateBase>;
                             using BdyUpdateList  = GTVector<GTVector<UpdateBasePtr>>;

                             typedef GTVector<GTVector<GTVector<GSIZET>>>   BinnedBdyIndex;
                             typedef GTVector<GTVector<GTVector<GBdyType>>> BinnedBdyType;

                             static_assert(std::is_same<State,GTVector<GTVector<Ftype>*>>::value,
                               "State is of incorrect type");
                             static_assert(std::is_same<StateComp,GTVector<Ftype>>::value,
                               "StateComp is of incorrect type");
                             static_assert(std::is_same<Derivative,GTVector<GTVector<Ftype>*>>::value,
                               "Derivative is of incorrect type");
                             static_assert(std::is_same<Ftype,Ftype>::value,
                               "Ftype is of incorrect type");


                             GGrid() = delete;
                             GGrid(const geoflow::tbox::PropertyTree &ptree, GTVector<GNBasis<GCTYPE,Ftype>*> &b, GC_COMM &comm);

virtual                     ~GGrid();

virtual void                 do_elems() = 0;                           // compute grid for irank
virtual void                 do_elems(GTMatrix<GINT> &p,
                               GTVector<GTVector<Ftype>> &xnodes) = 0;// compute grid on restart
//virtual void               set_partitioner(GDD_base<GTICOS> *d) = 0; // set and use GDD object

#if 0
virtual void                 set_bdy_callback(
                             std::function<void(GElemList &)> *callback) 
                             {bdycallback_ =  callback; }              // set bdy-set callback
#endif

virtual void                 print(const GString &filename){}          // print grid to file


        void                 grid_init();                             // initialize class
        void                 grid_init(GTMatrix<GINT> &p, 
                               GTVector<GTVector<Ftype>> &xnodes);   // initialize class for restart
        void                 add_terrain(const State &xb, State &tmp);// add terrain 
        void                 do_typing(); // classify into types
        GCBLAS::cuMatBlockDat 
                            &get_cudat() { return cudat_; }           // get CUDA data

        GElemList           &elems() { return gelems_; }              // get elem list
        GSIZET               nelems() { return gelems_.size(); }      // local num elems
        GSIZET               ngelems() { return ngelems_; }           // global num elems
        GTVector<GSIZET> 
                            &ntype() { return ntype_; }               // no. elems of each type
        GTVector<GTVector<GSIZET>> 
                            &itype() { return itype_; }               // indices for all types
        GTVector<GSIZET>    &itype(GElemType i) { return itype_[i]; } // indices for type i    
        GElemType            gtype() { return gtype_; }               // get unique elem type on grid       
        GBOOL                ispconst();                              // is order constant?
        void                 dealias(StateComp &v1, StateComp &v2, 
                                     StateComp &prod);                // dealias for quadratic nonlinearity
        void                 deriv(GTVector<Ftype> &u, GINT idir, GTVector<Ftype> &tmp,
                                   GTVector<Ftype> &du );             // derivative of global vector
        void                 deriv(GTVector<Ftype> &u, GINT idir, GBOOL dotrans, GTVector<Ftype> &tmp,
                                   GTVector<Ftype> &du );            // derivative of global vector
       void                  wderiv(GTVector<Ftype> &q, GINT idir, GBOOL bwghts, 
                                   GTVector<Ftype> &utmp, GTVector<Ftype> &du); // weak derivative
                           


        void                 set_apply_bdy_callback(
                             std::function<void(const Time &t, State &u,
                                         State &ub)> callback)
                                         { bdy_apply_callback_ = callback;
                                           bapplybc_ = TRUE; }        // set bdy-application callback

        Ftype                integrate(GTVector<Ftype> &u,
                                       GTVector<Ftype> &tmp, 
				       GBOOL bglobal = TRUE);         // spatial integration of global vector
        void                 print(GString fname, GBOOL bdof=FALSE);
        GSIZET               ndof();                                  // compute total number elem dof
        GSIZET               size() { return ndof(); }
        GSIZET               nfacedof();                              // compute total number elem face dof
        GSIZET               nbdydof();                               // compute total number elem bdy dof
        Ftype                minlength(GTVector<Ftype> *dx=NULLPTR); // find min elem length
        Ftype                maxlength(GTVector<Ftype> *dx=NULLPTR); // find max elem length
        Ftype                avglength();                             // find avg elem length
        Ftype                minnodedist()         
                             {return minnodedist_;}                   // get min node distance
        GTVector<Ftype>     &dxmin()         
                             {return dxmin_;}                         // get min node distance for each elem
        Ftype                volume()         
                             {return volume_;}                        // get grid volume
        Ftype                ivolume()         
                             {return ivolume_;}                       // get nverse of grid volume
        GTMatrix<GTVector<Ftype>>
                            &dXidX();                                 // global Rij = dXi^j/dX^i
        GTVector<Ftype>     &dXidX(GSIZET i,GSIZET j);                // Rij matrix element 
        GTVector<GTVector<Ftype>>
                            &xNodes() { return xNodes_; }             // get all nodes coords (Cart)
        GTVector<Ftype>     &xNodes(GSIZET i) { return xNodes_[i]; }  // get all nodes coords (Cart)
                            

        Mass                &massop();                                 // global mass op
        Mass                &imassop();                                // global inv.mass op
        GTVector<Ftype>     &Jac();                                    // global Jacobian
        GTVector<Ftype>
                            &faceJac();                                // global face Jacobian
        BdyUpdateList       &bdy_update_list() 
                             { return bdy_update_list_; }              // bdy_update_list pointers
        GTVector<GTVector<Ftype>>
                            &faceNormals()
                             { return faceNormals_; }                  // elem face normals
        GTVector<Ftype>
                            &faceMass()
                             { return faceMass_; }                     // elem face Mass * Jac
        GTVector<GSIZET>
                            &gieface() { return gieface_;}             // elem face indices into glob u for all elem faces
        BinnedBdyIndex
                            &igbdy_binned() { return igbdy_binned_;}   // global dom bdy indices binned into GBdyType
        GTVector<GSIZET>
                            &igbdy() { return igbdy_;}                 // global dom bdy indices into u
        GTVector<GTVector<GSIZET>>  
                            &igbdy_bdyface() { return igbdy_bdyface_;} // bdy ind for each face node
        GTVector<GTVector<GBdyType>>  
                            &igbdyt_bdyface() { return igbdyt_bdyface_;}// bdy types for each face node
        GTVector<GTVector<Ftype>>
                            &bdyNormals() { return bdyNormals_; }      // bdy normals
        GTVector<GINT>      &idepComp  () { return idepComp_; }        // dependent vector components on gbdy 
        GC_COMM              get_comm() { return comm_; }              // get communicator
        void                 set_derivtype(GDerivType gt);             // set deriv. method
        GDerivType           get_derivtype() { return gderivtype_; }   // return deriv. method
        GBOOL                usebdydata() { return do_face_normals_;}  // get usebdydata flag
        void                 set_usebdydata(GBOOL bflag) 
                             {do_face_normals_= bflag;}                // set usebdydata flag
 


        GGFX<Ftype>         &get_ggfx() { return *ggfx_; }             // get GGFX op
        void                 set_ggfx(GGFX<Ftype> &ggfx) 
                               { ggfx_ = &ggfx; }                      // set GGFX op    
        GTVector<Ftype>     &get_mask() { return mask_; }              // get mask

        void                 smooth(GTVector<Ftype> &tmp, 
                                    GTVector<Ftype> &u);              // H1-smoothing operatrion     
        void                 compute_grefderiv(GTVector<Ftype> &u, GTVector<Ftype> &etmp,
                                               GINT idir, GBOOL dotrans, GTVector<Ftype> &du);

        void                 compute_grefderivW(GTVector<Ftype> &u, GTVector<Ftype> &etmp,
                                                GINT idir, GBOOL dotrans, GTVector<Ftype> &du);

        void                 compute_grefderivs(GTVector<Ftype> &u, GTVector<Ftype> &etmp,
                                                GBOOL btrans, GTVector<GTVector<Ftype>*> &du);
        void                 compute_grefderivsW(GTVector<Ftype> &u, GTVector<Ftype> &etmp,
                                                 GBOOL btrans, GTVector<GTVector<Ftype>*> &du);
        void                 compute_grefdiv(GTVector<GTVector<Ftype>*> &u, GTVector<Ftype> &etmp,
                                              GBOOL btrans, GTVector<Ftype> &divu);

         GTVector<GINT>     &testid() { return testid_; }
         GTVector<GINT>     &testty() { return testty_; }


#if 0
friend  std::ostream&        operator<<(std::ostream&, Grid &);       // Output stream operator
#endif
 

protected:
       
        void                 grefderiv_varp  (GTVector<Ftype> &u, GTVector<Ftype> &etmp,
                                              GINT idir, GBOOL dotrans, GTVector<Ftype> &du);
        void                 grefderiv_constp(GTVector<Ftype> &u, GTVector<Ftype> &etmp,
                                              GINT idir, GBOOL dotrans, GTVector<Ftype> &du);
virtual void                 config_gbdy(const PropertyTree &ptree, 
                               GBOOL                         bterrain,
                               GTVector<GTVector<GSIZET>>   &igbdyf, 
                               GTVector<GTVector<GBdyType>> &igbdyft,
                               GTVector<GSIZET>             &igbdy,
                               GTVector<GUINT>              &debdy)=0; // config bdy
virtual void                 elem_face_data(
                               GTMatrix<GTVector<Ftype>>       &dXdXi,
                               GTVector<GSIZET>                 &igeface,
                               GTVector<Ftype>                 &face_mass,
                               GTVector<GTVector<Ftype>>       &normals)=0;// compute elem face data

                                                                       // compute bdy normals entry point

//      void                        init_local_face_info();            // get local face info
        void                        globalize_coords();                // create global coord vecs from elems
        void                        init_bc_info(GBOOL bterrain=FALSE);// configure bdys
        void                        init_qdealias();                   // create quadratic dealias data
        void                        def_geom_init();                   // iniitialze deformed elems
        void                        reg_geom_init();                   // initialize regular elems
        void                        do_face_data();                    // compute normals to elem faces 
        Ftype                       find_min_dist(); 
        void                        find_min_dist(GTVector<Ftype> &dx); 

        GBOOL                       bInitialized_;     // object initialized?
        GBOOL                       bapplybc_;         // bc apply callback set
        GBOOL                       do_face_normals_;  // compute elem face normals for fluxes?
        GBOOL                       bpconst_;          // is p const among elems?
        GBOOL                       usebdydat_;        // flag to set to use bdy data
        GBOOL                       do_gbdy_test_;     // create data required to test gbdys?
        GBOOL                       bInitQDealias_;    // quadratic dealias data initialized?
        GBOOL                       doQDealias_;       // do quadratic dealiasing?
        GINT                        nstreams_;         // no. CUDA streams
        
        GDerivType                  gderivtype_;       // ref. deriv method type

        GElemType                   gtype_;            // element types comprising grid
        GINT                        irank_;            // MPI task id
        GINT                        nprocs_;           // number of MPI tasks
        GSIZET                      ngelems_;          // global number of elements
        GC_COMM                     comm_;             // communicator
        Ftype                       minnodedist_;      // min node length array (for each elem)
        Ftype                       eps_;              // spatial epsilon
	Ftype                       volume_;           // grid volume
	Ftype                       ivolume_;          // 1 / grid volume
        GElemList                   gelems_;           // element list
        GTVector<Ftype>             etmp_;             // elem-level tmp vector
        GTVector<Ftype>             dxmin_;            // elem-based min node dist
        GTVector<GTVector<GSIZET>>  itype_;            // indices in elem list of each type
        GTVector<GSIZET>            ntype_;            // no. elems of each type on grid
        GTMatrix<GTVector<Ftype>>   dXidX_;            // matrix Rij = dXi^j/dX^i, global
        GTMatrix<GTVector<Ftype>>   dXdXi_;            // matrix Bij = dX^j/dXi^i, global, used for constructing normals
        GTVector<GTVector<Ftype>>   xNodes_;           // Cart coords of all node points
        Mass                       *mass_;             // mass operator
        Mass                       *imass_;            // inverse of mass operator
        GTVector<Ftype>             Jac_;              // interior Jacobian, global
        GTVector<Ftype>             faceJac_;          // face Jacobians, global
        GTVector<Ftype>             faceMass_;         // elem face mass * Jacobians
        GTVector<GTVector<Ftype>>   faceNormals_;      // normal to elem faces each face node point (2d & 3d), global
        GTVector<GSIZET>            gieface_;          // index into global field indicating elem face node
        GTVector<GTVector<Ftype>>   bdyNormals_;       // normal to surface at each bdy node point (2d & 3d), global
        GTVector<GINT>              idepComp_;         // dependent component index at each bdy point
        GTVector<GUINT>             degbdy_;           // gbdy node descriptorsgbdy node descriptors
        BinnedBdyIndex              igbdy_binned_;     // index into global field indicating a domain bdy--by type
        GTVector<GTVector<GSIZET>>  igbdy_bdyface_;    // volumbe index for each bdy node on each face
        GTVector<GTVector<GBdyType>> igbdyt_bdyface_;  // bdy type for each igbdt_bdyface_
        GTVector<GSIZET>            igbdy_;            // index into global field indicating a domain bdy
        BdyUpdateList               bdy_update_list_;
                                                       // bdy update class list
        GTVector<Ftype>             mask_;             // bdy mask
        GTVector<GTMatrix<Ftype>>   IQPdealias_;       // 1d quadratic dealias interp mats P->Q
        GTVector<GTMatrix<Ftype>>   IQPdealiasT_;      // transp of 1d quadratic dealias interp mats 
        GTVector<GINT>              qdN_;              // # nodes in each elem for quadratic dealiasing
        GTVector<Ftype>             iWp_;              // inv. of p-based weights (no JAc) for dealiasing
        GTVector<Ftype>             qW_;               // quadratic dealias weights
        std::vector<GINT>           pqdealias_;        // order of quadratic dealias basis in each direction
        GTVector<Ftype>             tptmp_;            // tensor product tmp space
        GTVector<GTVector<Ftype>>   qdtmp_;            // quadratic dealias tmp space
        PropertyTree                ptree_;            // main prop tree
        GGFX<Ftype>                *ggfx_;             // connectivity operator
        typename LinSolverBase<CGTypePack>::Traits
                                    cgtraits_;         // GCG operator traits

        std::function<void(const Time &t, State &u, State &ub)>
                                    bdy_apply_callback_;            
        GCBLAS::cuMatBlockDat       cudat_;            // CUDA data structure

        GTVector<GINT>              testid_;
        GTVector<GINT>              testty_;

};


#include "ggrid.ipp"


#endif
