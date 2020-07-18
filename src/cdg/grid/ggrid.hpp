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
#include "ggrid.hpp"
#include "gcomm.hpp"
#include "ggfx.hpp"
#include "glinop.hpp"
#include "ghelmholtz.hpp"
#include "gcg.hpp"
#include "gupdatebdy_factory.hpp"
#include "tbox/property_tree.hpp"


using namespace geoflow::tbox;
using namespace std;


typename PrecondType         = GLinOp;
typename State               = GTVector<GTVector<GFTYPE>*>;
typename StateComp           = GTVector<GFTYPE>;
typename StateInfo           = GStateInfo;
typename Grid                = GGrid;
typename Value               = GFTYPE;
typename Time                = GFTYPE;
typename BinnedBdyIndex      = GTVector<GTVector<GTVector<GSIZET>>>;
typename BinnedBdyType       = GTVector<GTVector<GTVector<GBdyType>>>;


struct BdyTypePack { // define bdy update typepack
        using Types            = BdyTypePack;
        using State            = GTVector<GTVector<GFTYPE>*>;
        using StateComp        = GTVector<GFTYPE>;
        using StateInfo        = GStateInfo;
        using Grid             = GGrid;
        using Value            = GFTYPE;
        using Time             = GFTYPE;
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

typename UpdateBase          = UpdateBdyBase<BdyTypePack>;
typename UpdateBasePtr       = std::shared_ptr<UpdateBase>;
typename BdyUpdateList       = GTVector<GTVector<UpdateBasePtr>>;

class GMass;

typedef GTVector<GElem_base*> GElemList;
typedef GFTYPE                Time;
typedef GStateInfo            StateInfo;

class GGrid 
{
public:
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

virtual                       ~GGrid();
//static                       GGrid *build(geoflow::tbox::PropertyTree &ptree, GTVector<GNBasis<GCTYPE,GFTYPE>*> &b, GC_COMM comm);

virtual void                 do_elems() = 0;                            // compute grid for irank
virtual void                 do_elems(GTMatrix<GINT> &p,
                               GTVector<GTVector<GFTYPE>> &xnodes) = 0; // compute grid on restart
//virtual void                 set_partitioner(GDD_base<GTICOS> *d) = 0;   // set and use GDD object

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

        void                 set_bdy_update_callback(
                             std::function<void(const geoflow::tbox::PropertyTree &ptree,
                                            GString &supdate, Grid &grid,
                                            StateInfo &stinfo, const Time &t, 
                                            State &utmp, State &u, State &ub)> callback)
                                           { bdy_update_callback_ =  callback;
                                             bupdatebc_ = TRUE; }       // set bdy-update callback
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
        GFTYPE               minlength();                             // find min elem length
        GFTYPE               maxlength();                             // find max elem length
        GFTYPE               avglength();                             // find avg elem length
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
        GMass               &imassop();                                // global inv.mass op
        GTVector<GFTYPE>    &Jac();                                    // global Jacobian
        GTVector<GFTYPE>
                            &faceJac();                                // global face Jacobian
        BdyUpdateList       &bdy_update_list() 
                             { return bdy_update_list_; }              // bdy_update_list pointers
        GTVector<GTVector<GFTYPE>>
                            &faceNormal();                             // global face normals
        GTVector<GSIZET>
                            &gieface() { return gieface_;}             // elem face indices into glob u for all elem faces
        BinnedBdyIndex
                            &igbdy_binned() { return igbdy_binned_;}   // global dom bdy indices binned into GBdyType
#if 0
        BinnedBdyIndex
                            &ilbdy_binned() { return ilbdy_binned_;}   // indices into bdy arrays binned into GBdyType
#endif
        GTVector<GSIZET>
                            &igbdy() { return igbdy_;}                 // global dom bdy indices into u
        GTVector<GTVector<GFTYPE>>
                            &bdyNormals() { return bdyNormals_; }      // bdy normals
        GTVector<GINT>      &idepComp  () { return idepComp_; }        // dependent vector components on bdy 
        GC_COMM              get_comm() { return comm_; }              // get communicator

virtual void                 config_bdy(const PropertyTree &ptree, 
                             GTVector<GTVector<GSIZET>>   &igbdy, 
                             GTVector<GTVector<GBdyType>> &igbdyt)=0;  // config bdy

        GGFX<GFTYPE>        &get_ggfx() { return *ggfx_; }             // get GGFX op
        void                 set_ggfx(GGFX<GFTYPE> &ggfx) 
                               { ggfx_ = &ggfx; }                      // set GGFX op    
        GTVector<GFTYPE>    &get_mask() { return mask_; }              // get mask

friend  std::ostream&        operator<<(std::ostream&, GGrid &);       // Output stream operator
 

protected:
       
virtual void                        do_face_normals()=0;              // compute normals to elem faces 
virtual void                        do_bdy_normals(GTMatrix<GTVector<GFTYPE>>
                                                                       &dXdXi,
                                                   GTVector<GTVector<GSIZET>>  
                                                                       &igbdy,
                                                   GTVector<GTVector<GFTYPE>> 
                                                                       &normals,
                                                   GTVector<GINT>      &idepComp)=0;
                                                                      // compute bdy normals entry point

        void                        init_local_face_info();           // get local face info
        void                        globalize_coords();               // create global coord vecs from elems
        void                        init_bc_info();                   // configure bdys
        void                        def_geom_init();                  // iniitialze deformed elems
        void                        reg_geom_init();                  // initialize regular elems
        void                        do_normals();                     // compute normals to elem faces and domain bdy
        GFTYPE                      find_min_dist(); 

        GBOOL                       bInitialized_;  // object initialized?
        GBOOL                       bapplybc_;      // bc apply callback set
        GBOOL                       bupdatebc_;     // bc update callback set
        GBOOL                       do_face_normals_; // compute elem face normals for fluxes?
        GElemType                   gtype_;         // element types comprising grid
        GINT                        irank_;         // MPI task id
        GINT                        nprocs_;        // number of MPI tasks
        GSIZET                      ngelems_;       // global number of elements
        GC_COMM                     comm_;          // communicator
        GFTYPE                      minnodedist_;   // min node length array (for each elem)
	GFTYPE                      volume_;        // grid volume
	GFTYPE                      ivolume_;       // 1 / grid volume
        GElemList                   gelems_;        // element list
        GTVector<GFTYPE>            etmp_;          // elem-level tmp vector
        GTVector<GTVector<GSIZET>>  itype_;         // indices in elem list of each type
        GTVector<GSIZET>            ntype_;         // no. elems of each type on grid
        GTMatrix<GTVector<GFTYPE>>  dXidX_;         // matrix Rij = dXi^j/dX^i, global
        GTMatrix<GTVector<GFTYPE>>  dXdXi_;         // matrix Bij = dX^j/dXi^i, global, used for constructing normals
        GTVector<GTVector<GFTYPE>>  xNodes_;        // Cart coords of all node points
        GMass                      *mass_;          // mass operator
        GMass                      *imass_;         // inverse of mass operator
        GTVector<GFTYPE>            Jac_;           // interior Jacobian, global
        GTVector<GFTYPE>            faceJac_;       // face Jacobians, global
        GTVector<GTVector<GFTYPE>>  faceNormals_;   // normal to eleme faces each face node point (2d & 3d), global
        GTVector<GSIZET>            gieface_;       // index into global field indicating elem face node
        GTVector<GTVector<GFTYPE>>  bdyNormals_;    // normal to surface at each bdy node point (2d & 3d), global
        GTVector<GINT>              idepComp_;      // dependent component index at each bdy point
        GTVector<GTVector<GSIZET>>  igbdy_binned_;  // index into global field indicating a domain bdy--by type
//      GTVector<GTVector<GSIZET>>  ilbdy_binned_;  // index into bdy arrays--by type
//      GTVector<GTVector<GSIZET>>  igbdy_bdyface_;// volumbe index for each bdy node on each face
        GTVector<GTVector<GBdyType>> igbdyt_bdyface_;// bdy type for each igbdt_bdyface_
        GTVector<GSIZET>            igbdy_;         // index into global field indicating a domain bdy
        BdyUpdateList               bdy_update_list_;
                                                    // bdy update class list
        GTVector<GFTYPE>            mask_;          // bdy mask
        PropertyTree                ptree_;         // main prop tree
        GGFX<GFTYPE>               *ggfx_;          // connectivity operator
        LinSolverBase<CGTypes>::Traits
                                    cgtraits_;      // GCG operator traits

        std::function<void(const Time &t, State &u, State &ub)>
                                    bdy_apply_callback_;            
                                                    // bdy apply callback
        std::function<void(const geoflow::tbox::PropertyTree& ptree,
                      GString &supdate, Grid &grid,
                      StateInfo &stinfo, const Time &t, 
                      State &utmp,  State &u, State &ub)> 
                                    bdy_update_callback_;            
                                                    // bdy update callback

};

#endif
