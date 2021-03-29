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
#include "gstateinfo.hpp"
#include "gcomm.hpp"
#include "gtvector.hpp"
#include "gtmatrix.hpp"
#include "gnbasis.hpp"
#include "gelem_base.hpp"
#include "gdd_base.hpp"
#include "gcomm.hpp"
#include "ggfx.hpp"
#include "glinop.hpp"
#include "ghelmholtz.hpp"
#include "gcblas.hpp"
#include "pdeint/update_bdy_base.hpp"
#include "pdeint/lin_solver_base.hpp"
#include "tbox/property_tree.hpp"


using namespace geoflow::tbox;
using namespace geoflow::pdeint;
using namespace std;

class GMass;

struct BdyTypePack { // define bdy update typepack
        using State            = GTVector<GTVector<GFTYPE>*>;
        using StateComp        = GTVector<GFTYPE>;
        using StateInfo        = GStateInfo;
        using Grid             = GGrid;
        using Value            = GFTYPE;
        using Time             = GFTYPE;
        using Derivative       = GTVector<GTVector<GFTYPE>*>;;
        using CompDesc         = GTVector<GStateCompType>;
        using IBdyVol          = GTVector<GSIZET>;
        using TBdyVol          = GTVector<GBdyType>;
};

struct CGTypePack { // define terrain typepack
        using Types            = CGTypePack;
        using Operator         = class GHelmholtz;
        using Preconditioner   = LinSolverBase<Types>;
        using State            = GTVector<GTVector<GFTYPE>*>;
        using StateComp        = GTVector<GFTYPE>;
        using StateInfo        = GStateInfo;
        using Grid             = GGrid;
        using Value            = GFTYPE;
        using Time             = GFTYPE;
        using ConnectivityOp   = GGFX<GFTYPE>;
};

typedef UpdateBdyBase<BdyTypePack> UpdateBase;
typedef std::shared_ptr<UpdateBase> UpdateBasePtr;
typedef GTVector<GTVector<UpdateBasePtr>> BdyUpdateList;

typedef GLinOp                                 PrecondType;
typedef GTVector<GTVector<GFTYPE>*>            State;
typedef GTVector<GFTYPE>                       StateComp;
typedef GStateInfo                             StateInfo;
typedef GGrid                                  Grid;
typedef GFTYPE                                 Value;
typedef GFTYPE                                 Time;
typedef GTVector<GTVector<GTVector<GSIZET>>>   BinnedBdyIndex;
typedef GTVector<GTVector<GTVector<GBdyType>>> BinnedBdyType;

typedef GTVector<GElem_base*> GElemList;
typedef GFTYPE                Time;
typedef GStateInfo            StateInfo;

class GGrid 
{
public:
                             enum GDerivType {GDV_VARP=0, GDV_CONSTP}; 

                             using CGTypes        = CGTypePack;
                             using Operator       = typename CGTypes::Operator;
                             using Preconditioner = typename CGTypes::Preconditioner;
                             using State          = typename CGTypes::State;
                             using StateComp      = typename CGTypes::StateComp;
                             using Grid           = typename CGTypes::Grid;
                             using Value          = typename CGTypes::Value;
                             using ConnectivityOp = typename CGTypes::ConnectivityOp;

                             GGrid() = delete;
                             GGrid(const geoflow::tbox::PropertyTree &ptree, GTVector<GNBasis<GCTYPE,GFTYPE>*> &b, GC_COMM &comm);

virtual                     ~GGrid();
//static                     GGrid *build(geoflow::tbox::PropertyTree &ptree, GTVector<GNBasis<GCTYPE,GFTYPE>*> &b, GC_COMM comm);

virtual void                 do_elems() = 0;                           // compute grid for irank
virtual void                 do_elems(GTMatrix<GINT> &p,
                             GTVector<GTVector<GFTYPE>> &xnodes) = 0;  // compute grid on restart
//virtual void               set_partitioner(GDD_base<GTICOS> *d) = 0; // set and use GDD object

#if 0
virtual void                 set_bdy_callback(
                             std::function<void(GElemList &)> *callback) 
                             {bdycallback_ =  callback; }              // set bdy-set callback
#endif

virtual void                 print(const GString &filename){}          // print grid to file


        void                 grid_init();                             // initialize class
        void                 grid_init(GTMatrix<GINT> &p, 
                               GTVector<GTVector<GFTYPE>> &xnodes);   // initialize class for restart
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
        void                 deriv(GTVector<GFTYPE> &u, GINT idir, GTVector<GFTYPE> &tmp,
                                   GTVector<GFTYPE> &du );            // derivative of global vector
        void                 deriv(GTVector<GFTYPE> &u, GINT idir, GBOOL dotrans, GTVector<GFTYPE> &tmp,
                                   GTVector<GFTYPE> &du );            // derivative of global vector
       void                  wderiv(GTVector<GFTYPE> &q, GINT idir, GBOOL bwghts, 
                                   GTVector<GFTYPE> &utmp, GTVector<GFTYPE> &du); // weak derivative
                           


        void                 set_apply_bdy_callback(
                             std::function<void(const Time &t, State &u,
                                         State &ub)> callback)
                                         { bdy_apply_callback_ = callback;
                                           bapplybc_ = TRUE; }        // set bdy-application callback

        GFTYPE               integrate(GTVector<GFTYPE> &u,
                                       GTVector<GFTYPE> &tmp, 
				       GBOOL bglobal = TRUE);         // spatial integration of global vector
        void                 print(GString fname, GBOOL bdof=FALSE);
        GSIZET               ndof();                                  // compute total number elem dof
        GSIZET               size() { return ndof(); }
        GSIZET               nfacedof();                              // compute total number elem face dof
        GSIZET               nbdydof();                               // compute total number elem bdy dof
        GFTYPE               minlength(GTVector<GFTYPE> *dx=NULLPTR); // find min elem length
        GFTYPE               maxlength(GTVector<GFTYPE> *dx=NULLPTR); // find max elem length
        GFTYPE               avglength();                             // find avg elem length
        GFTYPE               minnodedist()         
                             {return minnodedist_;}                   // get min node distance
        GTVector<GFTYPE>    &dxmin()         
                             {return dxmin_;}                         // get min node distance for each elem
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
        GMass               &imassop();                                // global inv.mass op
        GTVector<GFTYPE>    &Jac();                                    // global Jacobian
        GTVector<GFTYPE>
                            &faceJac();                                // global face Jacobian
        BdyUpdateList       &bdy_update_list() 
                             { return bdy_update_list_; }              // bdy_update_list pointers
        GTVector<GTVector<GFTYPE>>
                            &faceNormals()
                             { return faceNormals_; }                  // elem face normals
        GTVector<GFTYPE>
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
        GTVector<GTVector<GFTYPE>>
                            &bdyNormals() { return bdyNormals_; }      // bdy normals
        GTVector<GINT>      &idepComp  () { return idepComp_; }        // dependent vector components on gbdy 
        GC_COMM              get_comm() { return comm_; }              // get communicator
        void                 set_derivtype(GDerivType gt);             // set deriv. method
        GDerivType           get_derivtype() { return gderivtype_; }   // return deriv. method
        GBOOL                usebdydata() { return do_face_normals_;}  // get usebdydata flag
        void                 set_usebdydata(GBOOL bflag) 
                             {do_face_normals_= bflag;}                // set usebdydata flag
 


        GGFX<GFTYPE>        &get_ggfx() { return *ggfx_; }             // get GGFX op
        void                 set_ggfx(GGFX<GFTYPE> &ggfx) 
                               { ggfx_ = &ggfx; }                      // set GGFX op    
        GTVector<GFTYPE>    &get_mask() { return mask_; }              // get mask

        void                 smooth(GTVector<GFTYPE> &tmp, 
                                    GTVector<GFTYPE> &u);              // H1-smoothing operatrion     
        void                 compute_grefderiv(GTVector<GFTYPE> &u, GTVector<GFTYPE> &etmp,
                                               GINT idir, GBOOL dotrans, GTVector<GFTYPE> &du);

        void                 compute_grefderivW(GTVector<GFTYPE> &u, GTVector<GFTYPE> &etmp,
                                                GINT idir, GBOOL dotrans, GTVector<GFTYPE> &du);

        void                 compute_grefderivs(GTVector<GFTYPE> &u, GTVector<GFTYPE> &etmp,
                                                GBOOL btrans, GTVector<GTVector<GFTYPE>*> &du);
        void                 compute_grefderivsW(GTVector<GFTYPE> &u, GTVector<GFTYPE> &etmp,
                                                 GBOOL btrans, GTVector<GTVector<GFTYPE>*> &du);
        void                 compute_grefdiv(GTVector<GTVector<GFTYPE>*> &u, GTVector<GFTYPE> &etmp,
                                              GBOOL btrans, GTVector<GFTYPE> &divu);

         GTVector<GINT>     &testid() { return testid_; }
         GTVector<GINT>     &testty() { return testty_; }


friend  std::ostream&        operator<<(std::ostream&, GGrid &);       // Output stream operator
 

protected:
       
        void                 grefderiv_varp  (GTVector<GFTYPE> &u, GTVector<GFTYPE> &etmp,
                                              GINT idir, GBOOL dotrans, GTVector<GFTYPE> &du);
        void                 grefderiv_constp(GTVector<GFTYPE> &u, GTVector<GFTYPE> &etmp,
                                              GINT idir, GBOOL dotrans, GTVector<GFTYPE> &du);
        GBOOL                ispconst();
virtual void                 config_gbdy(const PropertyTree &ptree, 
                               GBOOL                         bterrain,
                               GTVector<GTVector<GSIZET>>   &igbdyf, 
                               GTVector<GTVector<GBdyType>> &igbdyft,
                               GTVector<GSIZET>             &igbdy,
                               GTVector<GUINT>              &debdy)=0; // config bdy
virtual void                 elem_face_data(
                               GTMatrix<GTVector<GFTYPE>>       &dXdXi,
                               GTVector<GSIZET>                 &igeface,
                               GTVector<GFTYPE>                 &face_mass,
                               GTVector<GTVector<GFTYPE>>       &normals)=0;// compute elem face data

                                                                       // compute bdy normals entry point

//      void                        init_local_face_info();            // get local face info
        void                        globalize_coords();                // create global coord vecs from elems
        void                        init_bc_info(GBOOL bterrain=FALSE);// configure bdys
        void                        def_geom_init();                   // iniitialze deformed elems
        void                        reg_geom_init();                   // initialize regular elems
        void                        do_face_data();                    // compute normals to elem faces 
        GFTYPE                      find_min_dist(); 
        void                        find_min_dist(GTVector<GFTYPE> &dx); 

        GBOOL                       bInitialized_;     // object initialized?
        GBOOL                       bapplybc_;         // bc apply callback set
        GBOOL                       do_face_normals_;  // compute elem face normals for fluxes?
        GBOOL                       bpconst_;          // is p const among elems?
        GBOOL                       usebdydat_;        // flag to set to use bdy data
        GBOOL                       do_gbdy_test_;     // create data required to test gbdys?
        GINT                        nstreams_;         // no. CUDA streams
        GDerivType                  gderivtype_;       // ref. deriv method type

        GElemType                   gtype_;            // element types comprising grid
        GINT                        irank_;            // MPI task id
        GINT                        nprocs_;           // number of MPI tasks
        GSIZET                      ngelems_;          // global number of elements
        GC_COMM                     comm_;             // communicator
        GFTYPE                      minnodedist_;      // min node length array (for each elem)
        GFTYPE                      eps_;              // spatial epsilon
	GFTYPE                      volume_;           // grid volume
	GFTYPE                      ivolume_;          // 1 / grid volume
        GElemList                   gelems_;           // element list
        GTVector<GFTYPE>            etmp_;             // elem-level tmp vector
        GTVector<GFTYPE>            dxmin_;            // elem-based min node dist
        GTVector<GTVector<GSIZET>>  itype_;            // indices in elem list of each type
        GTVector<GSIZET>            ntype_;            // no. elems of each type on grid
        GTMatrix<GTVector<GFTYPE>>  dXidX_;            // matrix Rij = dXi^j/dX^i, global
        GTMatrix<GTVector<GFTYPE>>  dXdXi_;            // matrix Bij = dX^j/dXi^i, global, used for constructing normals
        GTVector<GTVector<GFTYPE>>  xNodes_;           // Cart coords of all node points
        GMass                      *mass_;             // mass operator
        GMass                      *imass_;            // inverse of mass operator
        GTVector<GFTYPE>            Jac_;              // interior Jacobian, global
        GTVector<GFTYPE>            faceJac_;          // face Jacobians, global
        GTVector<GFTYPE>            faceMass_;         // elem face mass * Jacobians
        GTVector<GTVector<GFTYPE>>  faceNormals_;      // normal to elem faces each face node point (2d & 3d), global
        GTVector<GSIZET>            gieface_;          // index into global field indicating elem face node
        GTVector<GTVector<GFTYPE>>  bdyNormals_;       // normal to surface at each bdy node point (2d & 3d), global
        GTVector<GINT>              idepComp_;         // dependent component index at each bdy point
        GTVector<GUINT>             degbdy_;           // gbdy node descriptorsgbdy node descriptors
        BinnedBdyIndex              igbdy_binned_;     // index into global field indicating a domain bdy--by type
        GTVector<GTVector<GSIZET>>  igbdy_bdyface_;    // volumbe index for each bdy node on each face
        GTVector<GTVector<GBdyType>> igbdyt_bdyface_;  // bdy type for each igbdt_bdyface_
        GTVector<GSIZET>            igbdy_;            // index into global field indicating a domain bdy
        BdyUpdateList               bdy_update_list_;
                                                       // bdy update class list
        GTVector<GFTYPE>            mask_;             // bdy mask
        PropertyTree                ptree_;            // main prop tree
        GGFX<GFTYPE>               *ggfx_;             // connectivity operator
        LinSolverBase<CGTypes>::Traits
                                    cgtraits_;         // GCG operator traits

        std::function<void(const Time &t, State &u, State &ub)>
                                    bdy_apply_callback_;            
        GCBLAS::cuMatBlockDat       cudat_;            // CUDA data structure

        GTVector<GINT>              testid_;
        GTVector<GINT>              testty_;
};

#endif
