//==================================================================================
// Module       : gtest_derivs.cpp
// Date         : 2/24/19 (DLR)
// Description  : GeoFLOW test of derivatives
// Copyright    : Copyright 2019. Colorado State University. All rights reserved.
// Derived From : 
//==================================================================================

#include "gexec.h"
#include "gtypes.h"
#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <gptl.h>
#include <memory>
#include <cstdlib>
#include <cassert>
#include <random>
#include "gcomm.hpp"
#include "ggfx.hpp"
#include "gllbasis.hpp"
#include "gmass.hpp"
#include "gmorton_keygen.hpp"
#include "ggrid_factory.hpp"
#include "gmtk.hpp"
#include "tbox/property_tree.hpp"
#include "tbox/mpixx.hpp"
#include "tbox/global_manager.hpp"
#include "tbox/input_manager.hpp"

using namespace geoflow::tbox;
using namespace std;

#define  _DO_REFDERIV
#undef   _DO_REFDERIVW
#if defined(_DO_REFDERIVW) && defined(_DO_REFDERIV)
  #error "Cannot define both _DO_REFDERIVW AND _DO_REFDERIV"
#endif

GGrid *grid_ = NULLPTR;
void init_ggfx(GGrid &grid, GGFX<GFTYPE> &ggfx);

int main(int argc, char **argv)
{

    GString serr ="main: ";
    GBOOL  bret;
    GINT   errcode, iopt;
    GINT   ilevel=0;// 2d ICOS refinement level
    GINT   np=1;    // elem 'order'
    GINT   nc=GDIM; // no. coords
    std::vector<GINT> ne(3); // # elements in each direction in 3d
    std::vector<GINT> pstd(GDIM);  
    GString sgrid;// name of JSON grid object to use
    GC_COMM comm = GC_COMM_WORLD;

    // Initialize comm:
    GComm::InitComm(&argc, &argv);
    GINT myrank  = GComm::WorldRank();
    GINT nprocs  = GComm::WorldSize();

    // Read main prop tree; may ovewrite with
    // certain command line args:
    PropertyTree ptree;     // main prop tree
    PropertyTree gridptree; // grid prop tree
    PropertyTree polyptree; // poly_test prop tree

    ptree.load_file("input.jsn");       // main param file structure
    // Create other prop trees for various objects:
    sgrid       = ptree.getValue<GString>("grid_type");
    pstd        = ptree.getArray<GINT>("exp_order");
    gridptree   = ptree.getPropertyTree(sgrid);
    polyptree   = ptree.getPropertyTree("poly_test");

    ne          = gridptree.getArray<GINT>("num_elems");  // may be modified by command line

    nc = sgrid == "grid_icos" ? 3: GDIM;
    

    // If using GridBox, then P0, dP are used as expected, 
    // they aren't used for now:
    GTPoint<GFTYPE> P0, dP;
    std::vector<GFTYPE> vstd;
    if ( sgrid != "grid_icos" ) {
      vstd = gridptree.getArray<GFTYPE>("xyz0");
      P0   = vstd;
      vstd = gridptree.getArray<GFTYPE>("delxyz");
      dP   = vstd;
      }

#if 1

    // Parse command line. ':' after char
    // option indicates that it takes an argument.
    // Note: -i reserved for InputManager:
    while ((iopt = getopt(argc, argv, "i:j:k:l:m:p:h")) != -1) {
      switch (iopt) {
      case 'i': // handled by InputManager
          break;
      case 'j': // get # elements in r/x
          ne[0] = atoi(optarg);
          gridptree.setArray<GINT>("num_elems",ne);
          break;
      case 'k': // get # elements in lat/y
          ne[1] = atoi(optarg);
          gridptree.setArray<GINT>("num_elems",ne);
          break;
      case 'm': // get # elements in long/z
          ne[2] = atoi(optarg);
          gridptree.setArray<GINT>("num_elems",ne);
          break;
      case 'l': // # 2d refinement level
          ilevel = atoi(optarg);
          gridptree.setValue<GINT>("ilevel",ilevel);
          break;
      case 'p': // get nodal exp order
          np = atoi(optarg);
          pstd.assign(GDIM,np);
          ptree.setValue<GINT>("exp_order",np);
          break;
      case 'h': // help
          std::cout << "usage: " << std::endl <<
          argv[0] << " [-h] [-j #Elems in x/r] [-k #Elems in y/lat]  [-m #Elems in z/long] [-l refine level] -p expansion order] " << std::endl;
          std::cout << "Note: Icos grid only requires refine level to specify number of elements. " << std::endl;
          exit(1); 
          break;
      case ':': // missing option argument
          std::cout << argv[0] << ": option " << optopt << " requires an argument" << std::endl;
          exit(1); 
          break;
      case '?':
      default: // invalid option
          std::cout << argv[0] << ": option " << optopt << " invalid" << std::endl;
          exit(1);
          break;
      }
    }

#endif

    // Set GTPL options:
    GPTLsetoption (GPTLcpu, 1);

    // Initialize GPTL:
    GPTLinitialize();


    // Create basis:
    GTVector<GNBasis<GCTYPE,GFTYPE>*> gbasis(GDIM);
    for ( GSIZET k=0; k<GDIM; k++ ) {
      gbasis [k] = new GLLBasis<GCTYPE,GFTYPE>(pstd[k]);
    }
    
    GPTLstart("gen_grid");
    // Create grid:
    grid_ = GGridFactory::build(ptree, gbasis, comm);
    GPTLstop("gen_grid");



    GPTLstart("do_gather_op");
    // Initialize gather/scatter operator:
    GGFX<GFTYPE> ggfx;
    init_ggfx(*grid_, ggfx);
    GPTLstop("do_gather_op");


    // Create state and tmp space:
    GTVector<GTVector<GFTYPE>*> utmp(4);
    GTVector<GTVector<GFTYPE>*> u   (1);
    GTVector<GTVector<GFTYPE>*> du (nc);
    GTVector<GTVector<GFTYPE>*> da (nc);
    
    for ( GSIZET j=0; j<utmp.size(); j++ ) utmp[j] = new GTVector<GFTYPE>(grid_->size());
    for ( GSIZET j=0; j<u   .size(); j++ ) u   [j] = new GTVector<GFTYPE>(grid_->size());
    for ( GSIZET j=0; j<du  .size(); j++ ) du  [j] = new GTVector<GFTYPE>(grid_->size());
    for ( GSIZET j=0; j<da  .size(); j++ ) da  [j] = new GTVector<GFTYPE>(grid_->size());

    // Initialize u: set p, q, r exponents
    // (Can set up to read from input file):
    GFTYPE p=polyptree.getValue("xpoly",1);
    GFTYPE q=polyptree.getValue("ypoly",0);
    GFTYPE r=polyptree.getValue("zpoly",0);
    GFTYPE x, y, z=1.0;
    GTVector<GFTYPE> etmp1;
    GTVector<GTVector<GFTYPE>> *xnodes = &grid_->xNodes();   
    GTVector<GFTYPE>           *jac    = &grid_->Jac();   
    GTMatrix<GTVector<GFTYPE>> *dXidX  = &grid_->dXidX();   
    GMass                       mass(*grid_);
    GSIZET nxy = (*xnodes)[0].size();

    assert(p>=0 && q>=0 && r>=0 && "Polynomial order must be >= 0");
    for ( GSIZET j=0; j<nxy; j++ ) {
      x = (*xnodes)[0][j];
      y = (*xnodes)[1][j];
      if ( xnodes->size() > 2 ) z = (*xnodes)[2][j];
      (*u [0])[j] = pow(x,p)*pow(y,q)*pow(z,r);
      (*da[0])[j] = p==0 ? 0.0 : p*pow(x,p-1)*pow(y,q)*pow(z,r);
      (*da[1])[j] = q==0 ? 0.0 : q*pow(x,p)*pow(y,q-1)*pow(z,r);
      if ( xnodes->size() > 2 ) (*da[2])[j] = r==0 ? 0.0 : r*pow(x,p)*pow(y,q)*pow(z,r-1);
    }
#if defined(_DO_REFDERIVW)
    // Compute nc derivs on u, with weightsr; compare solutions
    // later with integrated analytic solution: 
    GMTK::compute_grefderivsW(*grid_, *u[0], etmp1, FALSE, du);
    for ( GSIZET j=0; j<nc; j++ ) {  
      du[j]->pointProd((*dXidX)(j,0)); // do chain rule for box grid
      du[j]->pointProd(*jac);
    }
#elif defined(_DO_REFDERIV)
    // Compute nc derivs on u, without weights: 
    for ( GSIZET j=0; j<nc; j++ ) {  // do chain rule
      if ( grid_->gtype() == GE_REGULAR ) {
        GMTK::compute_grefderiv (*grid_, *u[0], etmp1, j+1, FALSE, *du[j]);
        du[j]->pointProd((*dXidX)(j,0));
      }
      else {
        *du[j] = 0.0;
        for ( GSIZET k=0; k<nc; k++ ) {
#if 0
          GMTK::compute_grefderiv (*grid_, *u[0], etmp1, j+1, FALSE, *utmp[0]);
          utmp[0]->pointProd((*dXidX)(j,k));
          *du[j] += *utmp[0];
#else
          grid_->deriv(*u[0], j+1, *utmp[0], *du[j]);
#endif
        }
      }
    }
#else
    for ( GSIZET j=0; j<nc; j++ ) {
      grid_->deriv(*u[0], j+1, *utmp[0], *du[j]);
    }
#endif


    GSIZET lnelems=grid_->nelems();
    GSIZET gnelems;
    GFTYPE ftmp, nnorm, eps=10.0*std::numeric_limits<GFTYPE>::epsilon();
    GFTYPE x0, x1, y0, y1, z0, z1;
    std::ifstream itst;
    std::ofstream ios;
    GTVector<GFTYPE> da_int(nc), du_int(3), lnorm(3), gnorm(3), maxerror(3);

    if ( grid_->gtype() == GE_REGULAR ) {
      x0 = P0.x1; x1 = P0.x1+dP.x1;
      y0 = P0.x2; y1 = P0.x2+dP.x2;
      z0 = P0.x3; z1 = P0.x3+dP.x3;
    // So we don't try to take 0^0 = pow(0,0)
    z0 = z0 == 0 ? 1.0 : z0;
    assert(x0 != 0.0 
        && y0 != 0.0
        && z0 != 0.0
        && "Don't allow zero domain endpoint!");
    }

    maxerror = 0.0;

#if defined(_DO_REFDERIVW)
    assert(grid_->gtype()!=GE_2DEMBEDDED && "Do not set _DO_REFDERIVW with ICOS grid");

    // Compute integral of analytic solutions, compare
    if ( grid_->gtype() == GE_REGULAR ) {
      if ( nc == 2 ) {
        da_int[0] = (pow(x1,p  )-pow(x0,p  )) * (pow(y1,q+1)-pow(y0,q+1))/(q+1);
        da_int[1] = (pow(x1,p+1)-pow(x0,p+1)) * (pow(y1,q  )-pow(y0,q  ))/(p+1);
      } else if ( nc == 3 ) {
        da_int[0] = (pow(x1,p  )-pow(x0,p  )) * (pow(y1,q+1)-pow(y0,q+1)) * (pow(z1,r+1)-pow(z0,r+1))/((q+1)*(r+1));
        da_int[1] =  (pow(x1,p+1)-pow(x0,p+1)) * (pow(y1,q  )-pow(y0,q  )) * (pow(z1,r+1)-pow(z0,r+1))/((p+1)*(r+1));
        da_int[2] =  (pow(x1,p+1)-pow(x0,p+1)) * (pow(y1,q+1)-pow(y0,q+1)) * (pow(z1,r  )-pow(z0,r  ))/((p+1)*(q+1));
      }
    } 
    for ( GSIZET j=0; j<da.size(); j++ ) { // integral errors:
      cout << "main: error da[" << j << "]=" << *da[j] << endl; 
      cout << "main: error du[" << j << "]=" << *du[j] << endl; 
      ftmp = du[j]->sum();
      GComm::Allreduce(&ftmp  , du_int.data()+j , 1, T2GCDatatype<GFTYPE>() , GC_OP_SUM, comm);
      maxerror[j] = fabs(da_int[j] - du_int[j]) / (da_int[j]+1.0e-15);
      cout << "main: da_int[" << j << "]=" << da_int[j] 
           <<      " du_int[" << j << "]=" << du_int[j] << endl;
      cout << "main: error da_int-du_int[" << j << "]=" << maxerror[j] << endl;
    }
    // Print convergence data to file:
    itst.open("deriv_err.txt");
    ios.open("deriv_err.txt",std::ios_base::app);

    // Write header, if required:
    if ( itst.peek() == std::ofstream::traits_type::eof() ) {
    ios << "# p  num_elems   rel_err_x   rel_err_y  rel_err_z " << std::endl;
    }
    itst.close();

    ios << np  << "  "  << "  " << gnelems << "  "
        << "  " << maxerror[0] << "  " << maxerror[1] << "  " << maxerror[2]
        << std::endl;
    ios.close();

#else
    // Compute collocated  analytic solution, do comparisons:
    maxerror = 0.0;
    for ( GSIZET j=0; j<nc; j++ ) { //local errors
      for ( GSIZET i=0; i<da[j]->size(); i++ ) (*utmp[0])[i] = pow((*da[j])[i],2);
      nnorm = grid_->integrate(*utmp[0], *utmp[1]);
      nnorm = nnorm > std::numeric_limits<GFTYPE>::epsilon() ? nnorm : 1.0;
cout << "main: nnorm=" << nnorm << endl;
cout << "main: da[" << j << "]=" << *da[j] << endl;
cout << "main: du[" << j << "]=" << *du[j] << endl;
      *utmp[0] = *du[j] - *da[j];
      for ( GSIZET i=0; i<da[j]->size(); i++ ) (*utmp[1])[i] = fabs((*utmp[0])[i]);
      for ( GSIZET i=0; i<da[j]->size(); i++ ) (*utmp[2])[i] = pow((*utmp[0])[i],2);
      lnorm[0]  = utmp[0]->infnorm (); // inf-norm
      gnorm[1]  = grid_->integrate(*utmp[1],*utmp[0])/sqrt(nnorm);
      gnorm[2]  = sqrt(grid_->integrate(*utmp[2],*utmp[0])/nnorm);
       
      // Accumulate to find global errors for this field:
      GComm::Allreduce(lnorm.data()  , gnorm.data()  , 1, T2GCDatatype<GFTYPE>() , GC_OP_MAX, comm);
cout << "main: gnorm[" << j << "]=" << gnorm << endl;
      // now find max errors of each type for each field:
      for ( GSIZET i=0; i<maxerror.size(); i++ ) maxerror[i] = MAX(maxerror[i],gnorm[i]);
    }

    cout << "main: maxerror = " << maxerror << endl;
   
    GComm::Allreduce(&lnelems, &gnelems, 1, T2GCDatatype<GSIZET>() , GC_OP_SUM, comm);
    if ( maxerror[2] > 100.0*std::numeric_limits<GFTYPE>::epsilon() ) {
      std::cout << "main: -------------------------------------derivative FAILED" << std::endl;
      errcode = 1;
    } else {
      std::cout << "main: -------------------------------------derivative OK" << std::endl;
      errcode = 0;
    }

    // Print convergence data to file:
    itst.open("deriv_err.txt");
    ios.open("deriv_err.txt",std::ios_base::app);

    // Write header, if required:
    if ( itst.peek() == std::ofstream::traits_type::eof() ) {
    ios << "# p  num_elems   inf_err   L1_err   L2_err" << std::endl;
    }
    itst.close();

    ios << np  << "  "  << "  " << gnelems << "  "
        << "  " << maxerror[0] << "  " << maxerror[1] << "  " << maxerror[2]
        << std::endl;
    ios.close();
#endif
 
    GPTLpr_file("timing.txt");
    GPTLfinalize();

    GComm::TermComm();
    if ( grid_ != NULLPTR ) delete grid_;

    return( errcode );

} // end, main


//**********************************************************************************
//**********************************************************************************
// METHOD: init_ggfx
// DESC  : Initialize gather scatter operator
// ARGS  : grid    : GGrid object
//         ggfx    : gather/scatter op, GGFX
//**********************************************************************************
void init_ggfx(GGrid &grid, GGFX<GFTYPE> &ggfx)
{
  GFTYPE                         delta;
  GMorton_KeyGen<GNODEID,GFTYPE> gmorton;
  GTPoint<GFTYPE>                dX, porigin, P0;
  GTVector<GNODEID>              glob_indices;
  GTVector<GTVector<GFTYPE>>    *xnodes;

  delta  = grid.minnodedist();
  dX     = 0.5*delta;
  xnodes = &grid.xNodes();
  glob_indices.resize(grid.ndof());

  // Integralize *all* internal nodes
  // using Morton indices:
  gmorton.setIntegralLen(P0,dX);
  gmorton.key(glob_indices, *xnodes);

  // Initialize gather/scatter operator:
  GBOOL bret = ggfx.init(glob_indices);
  assert(bret && "Initialization of GGFX operator failed");

} // end method init_ggfx








