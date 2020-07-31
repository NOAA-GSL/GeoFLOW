//==================================================================================
// Module       : gtest_ylm.cpp
// Date         : 10/24/18 (DLR)
// Description  : GeoFLOW test of spherical harmonics implementation
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

#include "gexec.h"
#include "gtypes.h"
#include <cstdio>
#include <unistd.h>
#include <iostream>
#include "pdeint/observer_base.hpp"
#include "pdeint/io_base.hpp"
#include "pdeint/observer_factory.hpp"
#include "pdeint/null_observer.hpp"
#include "tbox/property_tree.hpp"
#include "tbox/mpixx.hpp"
#include "tbox/global_manager.hpp"
#include "tbox/input_manager.hpp"
#include "tbox/error_handler.hpp"
#include "gcomm.hpp"
#include "gmass.hpp"
#include "ggfx.hpp"
#include "gllbasis.hpp"
#include "gmorton_keygen.hpp"
#include "ggrid_box.hpp"
#include "ggrid_factory.hpp"
#include "gio_observer.hpp"
#include "gio.hpp"
#if defined(_G_USE_GPTL)
  #include "gptl.h"
#endif


using namespace geoflow::pdeint;
using namespace geoflow::tbox;
using namespace std;

struct stStateInfo {
  GINT        sttype  = 1;       // state type index (grid=0 or state=1)
  GINT        gtype   = 0;       // check src/cdg/include/gtypes.h
  GSIZET      index   = 0;       // time index
  GSIZET      nelems  = 0;       // num elems
  GSIZET      cycle   = 0;       // continuous time cycle
  GFTYPE      time    = 0.0;     // state time
  std::vector<GString>
              svars;             // names of state members
  GTVector<GStateCompType>
              icomptype;         // encoding of state component types    
  GTMatrix<GINT>
              porder;            // if ivers=0, is 1 X GDIM; else nelems X GDIM;
  GString     idir;              // input directory
  GString     odir;              // output directory
};

template< // default template arg types
typename StateType     = GTVector<GTVector<GFTYPE>*>,
typename StateInfoType = stStateInfo,
typename GridType      = GGrid,
typename ValueType     = GFTYPE,
typename DerivType     = StateType,
typename TimeType      = ValueType,
typename CompType      = GTVector<GStateCompType>,
typename JacoType      = StateType,
typename SizeType      = GSIZET 
>
struct EquationTypes {
        using State      = StateType;
        using StateInfo  = StateInfoType;
        using Grid       = GridType;
        using Value      = ValueType;
        using Derivative = DerivType;
        using Time       = TimeType;
        using CompDesc   = CompType;
        using Jacobian   = JacoType;
        using Size       = SizeType;
};
using MyTypes       = EquationTypes<>;           // Define types used
using EqnBase       = EquationBase<MyTypes>;     // Equation Base type
using EqnBasePtr    = std::shared_ptr<EqnBase>;  // Equation Base ptr
using IOBaseType    = IOBase<MyTypes>;           // IO Base type
using IOBasePtr     = std::shared_ptr<IOBaseType>;// IO Base ptr
using Grid          = GGrid;


GGrid       *grid_;
GC_COMM      comm_=GC_COMM_WORLD ;      // communicator

void init_ggfx(PropertyTree &ptree, Grid &grid, GGFX<GFTYPE> &ggfx);


int main(int argc, char **argv)
{
    GString   serr ="main: ";
    GINT      errcode, gerrcode;
    GINT      itest, lmax;
    GSIZET    ncheck;
    GFTYPE    eps=1e4*std::numeric_limits<GFTYPE>::epsilon();
    GFTYPE    cexcl, radius, rexcl;
    GFTYPE    err;
    GFTYPE    zero;
    GTVector<GFTYPE> vylm(2);
    IOBasePtr pIO;         // ptr to IOBase operator
    const char * const stest[] = { "orthonormality", "1-derivative"};

    typename ObserverBase<MyTypes>::Traits
                   binobstraits;

    EH_MESSAGE("main: Starting...");

    if ( argc > 1 ) {
      cout << "No arguments accepted" << endl;
      exit(1);
    }
    errcode = 0;

    // Initialize comm:
    GComm::InitComm(&argc, &argv);
    mpixx::environment  env(argc,argv); // init GeoFLOW comm
    mpixx::communicator world;
    GlobalManager::initialize(argc,argv);
    GlobalManager::startup();



    GINT myrank  = GComm::WorldRank();
    GINT nprocs  = GComm::WorldSize();

#if defined(_G_USE_GPTL)
    // Set GTPL options:
    GPTLsetoption (GPTLcpu, 1);

    // Initialize GPTL:
    GPTLinitialize();
#endif
    EH_MESSAGE("main: Read prop tree...");

    // Get minimal property tree:
    PropertyTree gtree, ptree; 
    GString sgrid;
    std::vector<GINT> pstd(GDIM);

    ptree.load_file("ylm_input.jsn");
    sgrid       = ptree.getValue<GString>("grid_type");
    pstd        = ptree.getArray<GINT>("exp_order");
    gtree       = ptree.getPropertyTree(sgrid);
    lmax        = ptree.getValue<GINT>("lmax",4);
    zero        = ptree.getValue<GFTYPE>("zero",1.0e-10);
    eps         = ptree.getValue<GFTYPE>("eps",1.0e-10);
    cexcl       = ptree.getValue<GFTYPE>("excl_angle",1.0); //exclude angle (deg)

    assert(sgrid == "grid_icos" && "Must use ICOS grid for now");

    radius = gtree.getValue<GFTYPE>("radius",1.0);
    rexcl  = cexcl * PI/180.0;

    // Create basis:
    GTVector<GNBasis<GCTYPE,GFTYPE>*> gbasis(GDIM);
    for ( auto k=0; k<GDIM; k++ ) {
      gbasis [k] = new GLLBasis<GCTYPE,GFTYPE>(pstd[k]);
    }

    EH_MESSAGE("main: Generate grid...");

    // Generate grid:
#if defined(_G_USE_GPTL)
    GPTLstart("gen_grid");
#endif

    ObserverFactory<MyTypes>::get_traits(ptree, "gio_observer", binobstraits);
    grid_ = GGridFactory<MyTypes>::build(ptree, gbasis, pIO, binobstraits, comm_);

#if defined(_G_USE_GPTL)
    GPTLstop("gen_grid");
#endif

    EH_MESSAGE("main: Initialize gather-scatter...");

    // Initialize gather/scatter operator:
    GGFX<GFTYPE> ggfx;
    init_ggfx(ptree, *grid_, ggfx);

    GLONG            na = pow(lmax*lmax,2);
    GTVector<GLONG> al (na);
    GTVector<GLONG> am (na);
    GTVector<GLONG> alp(na);
    GTVector<GLONG> amp(na);
    GTVector<GFTYPE> delta(na);
    GTVector<GFTYPE> rylm (grid_->ndof());
    GTVector<GFTYPE> rylmp(grid_->ndof());
    GTVector<GTVector<GFTYPE>*> ytmp(4);
    GTVector<GTVector<GFTYPE>> *xnodes = &grid_->xNodes();
   
    for ( auto j=0; j<ytmp.size(); j++ ) {
      ytmp[j] = new GTVector<GFTYPE>(grid_->ndof());
    }

    EH_MESSAGE("main: Check orthonormality...");

    GFTYPE integral;

    itest = 0;
    // Check orthonormality of real Ylm:
    GLONG l, m, lp, mp, n;
    for ( l=0, n=0; l<lmax; l++ ) {
      for ( m=-l; m<=l; m++ ) {
        for ( lp=0; lp<lmax; lp++ ) {
          for ( mp=-lp; mp<=lp; mp++, n++ ) {
            GMTK::rYlm_cart <GFTYPE>(l , m , *xnodes, *ytmp[0], rylm)  ; // rYlm basis
            GMTK::rYlm_cart <GFTYPE>(lp, mp, *xnodes, *ytmp[0], rylmp) ; // rYlm' basis
            rylm *= rylmp;
            integral = grid_->integrate(rylm, *ytmp[0])/(radius*radius);
//          cout << "main: l=" << l << " m=" << m << " lp=" << lp << " mp=" << mp << " integral=" << integral << endl;
            al   [n] = l;
            am   [n] = m;
            alp  [n] = lp;
            amp  [n] = mp;
            delta[n] = integral;
          }
        }
      }
    }


    GLONG  nbad=0;
    GTVector<GLONG> ibad(delta.size());
    for ( auto j=0; j< delta.size(); j++ ) {
      if ( al[j] == alp[j] && am[j] == amp[j] ) {
        if ( !FUZZYEQ(1.0,delta[j],eps) ) ibad[nbad++] = j;
      } 
      else {
        if ( fabs(delta[j]) > zero ) ibad[nbad++] = j;
      }
    }

   errcode = nbad == 0 ? 0 : 1;
   if ( nbad > 0 ) goto prerr;

   EH_MESSAGE("main: Check 1-derivative...");

   // Check first derivative: Integrate over dOMega:
   //   I_theta =Integral dYlm/dtheta sin theta dtheta dphi
   //           -Integral cot theta Ylm sin theta dtheta dphi
   itest = 1;


   // Reset lp, mp:
   alp = 999;
   amp = 999;

   GFTYPE x, y, z, colat, lon, r;
   GFTYPE ci, cotc;
   for ( l=0, n=0; l<=lmax; l++ ) {
     for ( m=-l; m<=l; m++, n++ ) {

        // Compute ulm(t=0):
        GMTK::rYlm_cart <GFTYPE>(l, m, *xnodes, *ytmp[0], rylm)  ; // rYlm basis
        GMTK::drYlm_cart<GFTYPE>(l, m, *xnodes,  1, cexcl, ytmp, rylmp) ; // drYlm/dth 
       *ytmp[0] = 0.0;
        for ( auto j=0; j<(*xnodes)[0].size(); j++ ) {
          x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
          r = sqrt(x*x + y*y + z*z); colat = acos(z/r); lon = atan2(y,x);
          colat = acos(z/r);
          if ( colat < rexcl || (PI-colat) < rexcl ) continue;
          cotc  = cos(colat)/sin(colat);
          if ( !FUZZYEQ(0.0,colat,eps) ) { vylm[0] = sin(colat)*rylm[j]; }
          if ( !FUZZYEQ(PI ,colat,eps) ) { vylm[1] = sin(colat)*rylm[j]; }
          
         (*ytmp[0])[j] =  -cotc*rylm[j]; // integrand
        }
        integral = grid_->integrate(*ytmp[0], *ytmp[1])/(radius*radius)
                 + 2.0*PI*(vylm[1] - vylm[0]);
        ci       = grid_->integrate(rylmp  , *ytmp[1])/(radius*radius); // Integral of deriv

        al   [n] = l;
        am   [n] = m;
        delta[n] = fabs(ci-integral) / integral; // rel error
        delta[n] = fabs(ci-integral); 
  cout << "main: l=" << l << " m=" << m << " integ=" << integral << " cinteg=" << ci << endl;
     }
   }
   ncheck = n;

   cout << "main: number of derivatives to check: " << ncheck << endl;

   nbad   = 0;
   for ( auto j=0; j< ncheck; j++ ) {
     if ( delta[j] > eps ) ibad[nbad++] = j;
   }
   
   errcode = nbad == 0 ? 0 : 2;
   if ( nbad > 0 ) goto prerr;


prerr:

    EH_MESSAGE("main: Write to file...");

    std::ifstream itst;
    std::ofstream ios;
    itst.open("ylm_err.txt");
    ios.open("ylm_err.txt",std::ios_base::app);

    if ( itst.peek() == std::ofstream::traits_type::eof() ) {
      ios << "test" << "  ";
      ios << "#elems" << "  ";
      for ( auto j=0; j<GDIM; j++ ) ios << "p" << j+1 << "  ";
      ios <<  "ndof    l     lp     m      mp      delta" << std::endl;
    }
    itst.close();

    // Get global no elements and dof:
    GLONG j;
    GTVector<GSIZET> lsz(2), gsz(2);
    lsz[0] = grid_->nelems();
    lsz[1] = grid_->ndof();
    GComm::Allreduce(lsz.data(), gsz.data(), 2, T2GCDatatype<GSIZET>() , GC_OP_SUM, comm_);

    for ( auto jj=0; jj<nbad; jj++ ) {
      j = ibad[jj];
      ios <<   stest[itest]  << "  ";
      ios <<   gsz[0] << "  ";
      for ( auto k=0; k<GDIM; k++ ) ios << pstd[k] << "  ";
      ios << gsz[1]     << "  "
          << al [j]     << "  " 
          << alp[j]     << "  " 
          << am [j]     << "  " 
          << amp[j]     << "  " 
          << delta[j]   << "  "
          << std::endl;
    }
    ios.close();
    

    // Accumulate errors:
    GComm::Allreduce(&errcode, &gerrcode, 1, T2GCDatatype<GINT>() , GC_OP_MAX, comm_);

 
    if ( gerrcode != 0 ) {
      cout << serr << " Error: test: " << stest[itest] << " code=" << errcode << endl;
      cout << serr << " Error: nbad: " << nbad << endl;
    }
    else {
      cout << serr << " Success!" << endl;
    }

#if defined(_G_USE_GPTL)
    GPTLpr_file("timing.txt");
    GPTLfinalize();
#endif

    GComm::TermComm();

    return(gerrcode);

} // end, main


//==================================================================================
// Module       : gtools.cpp
// Date         : 5/5/19 (DLR)
// Description  : GeoFLOW initialization tools
// Copyright    : Copyright 2019. Colorado State University. All rights reserved.
// Derived From : 
//==================================================================================
//#include "gtools.h"

//**********************************************************************************
//**********************************************************************************
// METHOD: init_ggfx
// DESC  : Initialize gather/scatter operator
// ARGS  : ptree   : main property tree
//         grid    : Grid object
//         ggfx    : gather/scatter op, GGFX
//**********************************************************************************
void init_ggfx(PropertyTree &ptree, Grid &grid, GGFX<GFTYPE> &ggfx)
{
  GString                        serr = "init_ggfx: ";
  GFTYPE                         delta;
  GFTYPE                         rad;
  GMorton_KeyGen<GNODEID,GFTYPE> gmorton;
  GTPoint<GFTYPE>                dX, porigin, P0;
  GTVector<GNODEID>              glob_indices;
  GTVector<GTVector<GFTYPE>>    *xnodes;
  GString                        sgrid;
  std::vector<GFTYPE>            pstd;
  PropertyTree                   gtree;

  sgrid = ptree.getValue<GString>("grid_type");
  gtree = ptree.getPropertyTree(sgrid);

  if ( sgrid == "grid_box" ) {
    pstd = gtree.getArray<GFTYPE>("xyz0");
    P0   = pstd;
  }
  if ( sgrid == "grid_icos" ) {
    rad   = gtree.getValue<GFTYPE>("radius");
    P0.x1 = -rad ;
    P0.x2 = -rad ;
    P0.x3 = -rad ;
  }
  if ( sgrid == "grid_sphere" ) {
    rad   = gtree.getValue<GFTYPE>("radiuso");
    P0.x1 = -rad ;
    P0.x2 = -rad ;
    P0.x3 = -rad ;
  }

  // First, periodize coords if required to, 
  // before labeling nodes:
  if ( typeid(grid) == typeid(GGridBox) ) { 
    static_cast<GGridBox*>(&grid)->periodize();
  }

  delta  = grid.minnodedist();
  dX     = 0.1*delta;
  xnodes = &grid.xNodes();
  glob_indices.resize(grid.ndof());

  // Integralize *all* internal nodes
  // using Morton indices:
  gmorton.setIntegralLen(P0,dX);

  gmorton.key(glob_indices, *xnodes);

  // Initialize gather/scatter operator:
  GBOOL bret;
  bret = ggfx.init(glob_indices);
  assert(bret && "Initialization of GGFX operator failed");

  // Unperiodize nodes now that connectivity map is
  // generated, so that coordinates mean what they should:
  if ( typeid(grid) == typeid(GGridBox) ) { // periodize coords
    static_cast<GGridBox*>(&grid)->unperiodize();
  }

} // end method init_ggfx

