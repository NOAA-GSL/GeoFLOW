//==================================================================================
// Module       : hackathon_driver.cpp
// Date         : 10/21/20 (DLR)
// Description  : GeoFLOW mini-app for tensor-product performance 
//                improvement
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : 
//==================================================================================


#include "gexec.h"
#include "gtypes.h"
#include "gcomm.hpp"
#include "ggfx.hpp"
#include "gllbasis.hpp"
#include "gmass.hpp"
#include "gmorton_keygen.hpp"
#include "gmtk.hpp"
#include "ggrid_factory.hpp"
#include "pdeint/observer_base.hpp"
#include "pdeint/observer_factory.hpp"
#include "tbox/pio.hpp"
#include "tbox/property_tree.hpp"
#include "tbox/mpixx.hpp"
#include "tbox/global_manager.hpp"
#include "tbox/input_manager.hpp"
#include "tbox/tracer.hpp"

#include <cassert>
//#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <random>
#include <unistd.h>


using namespace geoflow::tbox;
using namespace std;



template< // Complete typepack
typename StateType     = GTVector<GTVector<GFTYPE>*>,
typename StateCompType = GTVector<GFTYPE>,
typename StateInfoType = GStateInfo,
typename GridType      = GGrid,
typename MassOpType    = GMass,
typename ValueType     = GFTYPE,
typename DerivType     = StateType,
typename TimeType      = ValueType,
typename CompType      = GTVector<GStateCompType>,
typename JacoType      = StateType,
typename SizeType      = GSIZET
>
struct TypePack {
        using State      = StateType;
        using StateComp  = StateCompType;
        using StateInfo  = StateInfoType;
        using Grid       = GridType;
        using Mass       = MassOpType;
        using Value      = ValueType;
        using Derivative = DerivType;
        using Time       = TimeType;
        using CompDesc   = CompType;
        using Jacobian   = JacoType;
        using Size       = SizeType;
};
using MyTypes       = TypePack<>;           // Define grid types used
using Grid          = GGrid;
using IOBaseType    = IOBase<MyTypes>;          // IO Base type
using IOBasePtr     = std::shared_ptr<IOBaseType>;// IO Base ptr
using ObsTraitsType = ObserverBase<MyTypes>::Traits;

Grid      *grid_ = NULLPTR;
IOBasePtr  pIO_  = NULLPTR; // ptr to IOBase operator


#define NMETH 2


int main(int argc, char **argv)
{
	GEOFLOW_TRACE_INITIALIZE();

    GINT    errcode=0 ;
    GINT    nc=GDIM; // no. coords
    GFTYPE  eps=100.0*std::numeric_limits<GFTYPE>::epsilon();
    GFTYPE  dnorm, told=0.0, tnew=0.0;
    GString sgrid;// name of JSON grid object to use
    GTVector<GINT>
            pvec;
    std::vector<GINT> 
            pstd(GDIM);  
    GC_COMM comm = GC_COMM_WORLD;

    // Initialize comm:
    GComm::InitComm(&argc, &argv);
    GINT myrank  = GComm::WorldRank();
    GINT nprocs  = GComm::WorldSize();
    pio::initialize(myrank);
    GEOFLOW_TRACE();

    // Read main prop tree; may ovewrite with
    // certain command line args:
    PropertyTree ptree;     // main prop tree
    PropertyTree gridptree; // grid prop tree
    PropertyTree polyptree; // poly_test prop tree

    /////////////////////////////////////////////////////////////////
    /////////////////////// Initialize system ///////////////////////
    /////////////////////////////////////////////////////////////////
    ptree.load_file("hack_input.jsn");       // main param file structure

    // Create other prop trees for various objects:
    sgrid       = ptree.getValue<GString>("grid_type");
    pstd        = ptree.getArray<GINT>("exp_order");
    gridptree   = ptree.getPropertyTree(sgrid);
    polyptree   = ptree.getPropertyTree("poly_test");
    nc          = GDIM;
    assert(sgrid == "grid_box"); // use only Cartesian grids
    assert(pstd.size() >= GDIM); 

    pvec.resize(pstd.size()); pvec = pstd; pvec.range(0,GDIM-1);

    // Create basis:
    GTVector<GNBasis<GCTYPE,GFTYPE>*> gbasis(GDIM);
    for ( auto k=0; k<GDIM; k++ ) {
      gbasis [k] = new GLLBasis<GCTYPE,GFTYPE>(pstd[k]);
    }
    
    // Create grid:
    ObsTraitsType binobstraits;
    grid_ = GGridFactory<MyTypes>::build(ptree, gbasis, pIO_, binobstraits, comm);

    /////////////////////////////////////////////////////////////////
    ////////////////////// Allocate arrays //////////////////////////
    /////////////////////////////////////////////////////////////////

    // Create state and tmp space:
    State     utmp (3);
    State     da   (nc);
    State     du   (nc);
    StateComp diff, dunew, duold, u;
    
    for ( auto j=0; j<utmp .size(); j++ ) utmp [j] = new StateComp(grid_->size());
    for ( auto j=0; j<da   .size(); j++ ) da   [j] = new StateComp(grid_->size());
    u    .resize(grid_->size());
    diff .resize(grid_->size());
    dunew.resize(grid_->size());
    duold.resize(grid_->size());
    du = NULLPTR;

    /////////////////////////////////////////////////////////////////
    ////////////////////// Compute solutions/////////////////////////
    /////////////////////////////////////////////////////////////////

    // Initialize u: set p, q, r exponents;
    //   u = x^p y^q z^r:
    std::vector<GFTYPE> 
            vpoly = polyptree.getArray<GFTYPE>("poly");
    
    GINT    ncyc  = polyptree.getValue<GINT>("ncycles",100);
    GINT    idir  = polyptree.getValue<GINT>("idir",2);
    GFTYPE  p, q, r, x, y, z;
    GTVector<GTVector<GFTYPE>> 
           *xnodes = &grid_->xNodes();   

    assert(vpoly.size() >= GDIM); 
    assert(idir >= 1 && idir <= GDIM);
    p = vpoly[0]; q = vpoly[1]; r = GDIM == 3 ? vpoly[2] : 0.0;

    // Set scalar field, and analytic derivative:
    dnorm = 0.0;
    for ( auto j=0; j<(*xnodes)[0].size(); j++ ) {
      x = (*xnodes)[0][j];
      y = (*xnodes)[1][j];
      z = GDIM == 3 ? (*xnodes)[2][j] : 1.0;
        u     [j] = pow(x,p)*pow(y,q)*pow(z,r);
      (*da[0])[j] = p==0 ? 0.0 : p*pow(x,p-1)*pow(y,q)*pow(z,r);
      (*da[1])[j] = q==0 ? 0.0 : q*pow(x,p)*pow(y,q-1)*pow(z,r);
      if ( GDIM == 3 ) (*da[2])[j] = r==0 ? 0.0 
                        : r*pow(x,p)*pow(y,q)*pow(z,r-1);
    } // end, loop over grid points

    for ( auto j=0; j<nc; j++ ) dnorm = MAX(dnorm, da[j]->amax());

    // Compute numerical derivs of u in each direction, using
    // different methods:
    grid_->set_derivtype(GGrid::GDV_VARP); // variable order
    GEOFLOW_TRACE_START("old_deriv");
    for ( auto n=0; n<ncyc; n++ ) {
       grid_->deriv(u, idir, *utmp[0], duold);
    }
    GEOFLOW_TRACE_STOP();

    grid_->set_derivtype(GGrid::GDV_CONSTP); // const order
    GEOFLOW_TRACE_START("new_deriv");
    for ( auto n=0; n<ncyc; n++ ) {
       grid_->deriv(u, idir, *utmp[0], dunew);
    }
    GEOFLOW_TRACE_STOP();

//cout << "da_y  =" << *da   [idir-1] << endl;
//cout << "dnew_y=" <<  dunew << endl;

    // Find inf-norm and L2-norm errors for each method::
    GTMatrix<GFTYPE> errs(NMETH,2); // for each method, Linf and L2 errs
    StateComp        lnorm(2), gnorm(2);
    std::string      smethod[NMETH] = {"old", "new"};

    /////////////////////////////////////////////////////////////////
    //////////////////////// Compute Errors /////////////////////////
    /////////////////////////////////////////////////////////////////
    du[0] = &duold; du[1] = &dunew;
    for ( auto n=0; n<NMETH; n++ ) { // over old and new methods
      diff     = (*da[idir-1]) - (*du[n]);
     *utmp [0] = diff;                   // for inf-norm
     *utmp [1] = diff; utmp[1]->rpow(2); // for L2 norm

      lnorm[0] = utmp[0]->infnorm (); 
      gnorm[1] = sqrt(grid_->integrate(*utmp[1],*utmp[2]))/dnorm;

      // Accumulate to find global inf-norm:
      GComm::Allreduce(lnorm.data(), gnorm.data(), 1, T2GCDatatype<GFTYPE>(), GC_OP_MAX, comm);
      errs(n,0) = gnorm[0]; errs(n,1) = gnorm[1];
      if ( myrank == 0 ) {
        if ( errs(n,1) > eps ) {
          std::cout << "main: ---------------------------derivative FAILED: " << errs(n,1) << " : direction=" << idir << " method: " << smethod[n]  << std::endl;
          errcode += 1;
        } else {
          std::cout << "main: ---------------------------derivative OK: " << errs(n,1) << " : direction=" << idir << " method: " << smethod[n] << std::endl;
          errcode += 0;
        }
      }

    } // end, new-old method loop


    // Print convergence data to file:
    std::ifstream itst;
    std::ofstream ios;
    itst.open("deriv_err.txt");
    ios.open("deriv_err.txt",std::ios_base::app);

    // Write header, if required:
    if ( itst.peek() == std::ofstream::traits_type::eof() ) {
    ios << "#  idir   p      num_elems  ncyc  inf_err_old  L2_err_old  t_old   inf_err_new  L2_err_new   t_new" << std::endl;
    }
    itst.close();

    ios << idir << "  " << pvec[0] ;
        for ( auto j=1; j<GDIM; j++ ) ios << " " <<  pvec[j]; 
    ios << "  " << grid_->ngelems() 
        << "  " << ncyc
        << "  " << errs(0,0) << "  " << errs(0,1) << "  " << told
        << "  " << errs(1,0) << "  " << errs(1,1) << "  " << tnew
        << std::endl;
    ios.close();
 
    pio::finalize();
    GEOFLOW_TRACE_STOP();
    GComm::TermComm();
    GEOFLOW_TRACE_FINALIZE();

    return( errcode );

} // end, main

