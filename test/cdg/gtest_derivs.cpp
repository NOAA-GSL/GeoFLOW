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
#include "gtools.h"
#include "gio.h"

using namespace geoflow::tbox;
using namespace std;

#define   _DO_REFDERIV
#undef   _DO_REFDERIVW
#if defined(_DO_REFDERIVW) && defined(_DO_REFDERIV)
  #error "Cannot define both _DO_REFDERIVW AND _DO_REFDERIV"
#endif

GGrid *grid_ = NULLPTR;

void shape_deriv(GGrid &grid, GTVector<GFTYPE> &u, GINT jder, GTVector<GFTYPE> &utmp, GTVector<GFTYPE> &du);
void print(GGrid &grid, GTVector<GFTYPE> &u, const GString fpref, GC_COMM &comm);

int main(int argc, char **argv)
{

    GString serr ="main: ";
    GBOOL  bret;
    GINT   errcode, iopt;
    GINT   ilevel=0;// 2d ICOS refinement level
    GINT   np=1;    // elem 'order'
    GINT   nc=GDIM; // no. coords
    GFTYPE den, phi, psi, rad, theta;
    GFTYPE dphidx, dphidy;
    GFTYPE dthdx , dthdy , dthdz ;
    GFTYPE eps=10.0*std::numeric_limits<GFTYPE>::epsilon();
    std::vector<GINT> ne(3); // # elements in each direction in 3d
    std::vector<GINT> pstd(GDIM);  
    GString sgrid;// name of JSON grid object to use
    GString stmp ;
    std::stringstream ss;
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
    init_ggfx(ptree, *grid_, ggfx);
    GPTLstop("do_gather_op");

#if 0
    GTMatrix<GTVector<GFTYPE>> *dXidX = &grid_->dXidX();
    // Smooth dXidX elements:
    for ( GSIZET j=0; j<dXidX->size(2); j++ ) {
      for ( GSIZET i=0; i<dXidX->size(1); i++ ) {
        ggfx.doOp((*dXidX)(i,j), GGFX_OP_SMOOTH);
      }
    }
#endif

    // Create state and tmp space:
    GTVector<GTVector<GFTYPE>*> utmp(4);
    GTVector<GTVector<GFTYPE>*> u   (1);
    GTVector<GTVector<GFTYPE>*> du (nc);
    GTVector<GTVector<GFTYPE>*> da (nc);
    GTVector<GTVector<GFTYPE>*> diff(nc);
    
    for ( GSIZET j=0; j<utmp.size(); j++ ) utmp[j] = new GTVector<GFTYPE>(grid_->size());
    for ( GSIZET j=0; j<u   .size(); j++ ) u   [j] = new GTVector<GFTYPE>(grid_->size());
    for ( GSIZET j=0; j<du  .size(); j++ ) du  [j] = new GTVector<GFTYPE>(grid_->size());
    for ( GSIZET j=0; j<da  .size(); j++ ) da  [j] = new GTVector<GFTYPE>(grid_->size());
    for ( GSIZET j=0; j<diff.size(); j++ ) diff[j] = new GTVector<GFTYPE>(grid_->size());

    for ( GSIZET j=0; j<da  .size(); j++ ) *da  [j] = 0.0;

    // Initialize u: set p, q, r exponents
    // (Can set up to read from input file):
    GFTYPE p = polyptree.getValue<GFTYPE>("xpoly",1);
    GFTYPE q = polyptree.getValue<GFTYPE>("ypoly",1);
    GFTYPE r = polyptree.getValue<GFTYPE>("zpoly",1);
    GFTYPE sig=polyptree.getValue<GFTYPE>("sigma",0.05);
    GFTYPE dthetaex=polyptree.getValue<GFTYPE>("excise_angle",5.0); // excise angle (deg)
    GFTYPE x, y, z=1.0;
    GTVector<GFTYPE> etmp1, mask;
    GTVector<GTVector<GFTYPE>> *xnodes = &grid_->xNodes();   
    GTVector<GFTYPE>           *jac    = &grid_->Jac();   
    GMass                       mass(*grid_);
    GSIZET nxy = (*xnodes)[0].size();

    dthetaex *= (PI/180.0);

    if ( sgrid == "grid_icos" )
      assert(p>=0 && q>=0 && r>=0 && "Polynomial order p>=0, and q>=0");
    else 
      assert(p>=0 && q>=0 && r>=0 && "Polynomial order must be >= 0");

    mask.resize(nxy);
    mask = 1.0;
    for ( GSIZET j=0; j<nxy; j++ ) {
      x = (*xnodes)[0][j];
      y = (*xnodes)[1][j];
      if ( xnodes->size() > 2 ) z = (*xnodes)[2][j];
      if ( sgrid != "grid_icos" ) {
        (*u [0])[j] = pow(x,p)*pow(y,q)*pow(z,r);
        (*da[0])[j] = p==0 ? 0.0 : p*pow(x,p-1)*pow(y,q)*pow(z,r);
        (*da[1])[j] = q==0 ? 0.0 : q*pow(x,p)*pow(y,q-1)*pow(z,r);
        if ( xnodes->size() > 2 ) (*da[2])[j] = r==0 ? 0.0 : r*pow(x,p)*pow(y,q)*pow(z,r-1);

      }
      else {
        rad         = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
        theta       = asin(z/rad);
        phi         = atan2(y,x);
        (*u [0])[j] = pow(cos(theta),p)*cos(q*phi);
        mask[j]     = 0.5*PI - abs(theta) < dthetaex ? 0 : 1;
        dthdx       =-sin(theta)*cos(phi)/rad;
        dthdy       =-sin(theta)*sin(phi)/rad;
        dthdz       = cos(theta)/rad;
        dphidx      =-sin(phi)/(rad*cos(theta));
        dphidy      = cos(phi)/(rad*cos(theta));
        (*da[0])[j] = -p*pow(cos(theta),p-1)*cos(q*phi)*sin(theta)*dthdx  
                    -  q*pow(cos(theta),p  )*sin(q*phi)           *dphidx;
        (*da[1])[j] = -p*pow(cos(theta),p-1)*cos(q*phi)*sin(theta)*dthdy 
                    -  q*pow(cos(theta),p  )*sin(q*phi)           *dphidy;
        (*da[2])[j] = -p*pow(cos(theta),p-1)*cos(q*phi)*sin(theta)*dthdz;
      }
    } // end, loop over grid points

    for ( GSIZET j=0; j<da.size(); j++ ) {
      da[j]->pointProd(mask);
    }

    print(*grid_, *u[0], "u", comm);
    for ( GSIZET j=0; j<da.size(); j++ ) {
      ss << "da" << j;
      print(*grid_, *da[j], ss.str(), comm); 
      ss.str("");
    }

#if defined(_DO_REFDERIVW)
    assert(grid_->gtype() != GE_2DEMBEDDED && 
           "_DO_REFDERIVW not allowed with this grid");
    // Compute nc derivs on u, with weightsr; compare solutions
    // later with integrated analytic solution: 
    GMTK::compute_grefderivsW(*grid_, *u[0], etmp1, FALSE, du);
    for ( GSIZET j=0; j<nc; j++ ) {  
      du[j]->pointProd((*dXidX)(j,0)); // do chain rule for box grid
      du[j]->pointProd(*jac);
    }
#elif defined(_DO_REFDERIV)
    // Compute nc derivs on u, without weights: 
    for ( GSIZET j=0; j<du.size(); j++ ) {  // do chain rule
        grid_->deriv(*u[0], j+1, *utmp[0], *du[j]);
    } // end, j-loop
#else
    for ( GSIZET j=0; j<du.size(); j++ ) {
      shape_deriv(*grid_, *u[0], j+1, *utmp[0], *du[j]);
    }
#endif

    for ( GSIZET j=0; j<du.size(); j++ ) {
      du[j]->pointProd(mask);
//    ggfx.doOp(*du[j], GGFX_OP_SMOOTH);
    }


// GMTK::project2sphere<GFTYPE>(*grid_, du, du);

    GSIZET lnelems=grid_->nelems();
    GSIZET gnelems;
    GFTYPE ftmp, nnorm;
    GFTYPE dxmin, lmin;
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
    dxmin = grid_->minnodedist();
    lmin  = grid_->minlength();

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
    ios << "# p  num_elems  dx_min  lmin   rel_err_x   rel_err_y  rel_err_z " << std::endl;
    }
    itst.close();

    ios << np  << "  "  << "  " << gnelems << "  " << dxmin << "  " << lmin
        << "  " << maxerror[0] << "  " << maxerror[1] << "  " << maxerror[2]
        << std::endl;
    ios.close();

#else

cout << "main: u=" << *u[0] << endl;

    // Compute collocated  analytic solution, do comparisons:
    maxerror = 0.0;
    for ( GSIZET j=0; j<du.size(); j++ ) { //local errors
     *diff[j]  = (*da[j]) - (*du[j]);
     *utmp[0]  = *da[j]; utmp[0]->pow(2);
#if 0
      nnorm    = grid_->integrate(*utmp[0], *utmp[1]);
      nnorm    = nnorm > std::numeric_limits<GFTYPE>::epsilon() ? nnorm : 1.0;
#else
      nnorm    = 1.0;
#endif

cout << "main: nnorm=" << nnorm << endl;
     *utmp[0]  = *diff[j];
     *utmp[1]  = *diff[j]; utmp[1]->abs();
     *utmp[2]  = *diff[j]; utmp[2]->pow(2);
      lnorm[0] = utmp[0]->infnorm (); // inf-norm

      gnorm[1] = grid_->integrate(*utmp[1],*utmp[0])/sqrt(nnorm);
      gnorm[2] = sqrt(grid_->integrate(*utmp[2],*utmp[0])/nnorm);

       
      // Accumulate to find global errors for this field:
      GComm::Allreduce(lnorm.data()  , gnorm.data()  , 1, T2GCDatatype<GFTYPE>() , GC_OP_MAX, comm);
cout << "main: gnorm[" << j << "]=" << gnorm << endl;
      // now find max errors of each type for each field:
      for ( GSIZET i=0; i<maxerror.size(); i++ ) maxerror[i] = MAX(maxerror[i],gnorm[i]);

      ss << "diff" << j;
      print(*grid_, *diff[j], ss.str(), comm); 
      ss.str("");
      ss << "du" << j;
      print(*grid_, *du[j], ss.str(), comm); 
      ss.str("");
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
    ios << "# p  num_elems   dx_min   lmin  inf_err   L1_err   L2_err" << std::endl;
    }
    itst.close();

    ios << np  << "  "  << "  " << gnelems << "  " << dxmin << "  " << lmin
        << "  " << maxerror[0] << "  " << maxerror[1] << "  " << maxerror[2]
        << std::endl;
    ios.close();
#endif
 
    GPTLpr_file("timing.txt");
    GPTLfinalize();

    GComm::TermComm();

    return( errcode );

} // end, main

//**********************************************************************************
//**********************************************************************************
// METHOD: shape_deriv
// DESC  : Compute derivative in j direction using shape function interface
// ARGS  : grid    : GGrid object
//         u       : field
//         jder    : which direction to take deriv in
//         utmp    : tmp space
//         da      : return derivative
//**********************************************************************************
void shape_deriv(GGrid &grid, GTVector<GFTYPE> &u, GINT jder, GTVector<GFTYPE> &utmp, GTVector<GFTYPE> &du)
{
  GString   serr = "shape_deriv: ";
//GINT      nxy = grid.gtype()==GE_2DEMBEDDED ? GDIM+1 : GDIM;
  GSIZET    ibeg, iend, nref;
  GElemList *gelems = &grid.elems();
  GShapeFcn_embed *gshapefcn;
  GTMatrix<GTVector<GFTYPE>> *dXidX = &grid.dXidX();
  GTVector<GNBasis<GCTYPE,GFTYPE>*> gbasis(GDIM);
  GTVector<GINT>     N(GDIM), I(GDIM);
  GTVector<GFTYPE>   dNi(grid.ndof());
  GTVector<GTVector<GFTYPE>>   *xNodes = &grid.xNodes();
  GTVector<GTVector<GFTYPE>*>  xi_ev(GDIM);

  gshapefcn = new GShapeFcn_embed();

  nref = grid.gtype()==GE_2DEMBEDDED ? GDIM+1 : GDIM;

  du = 0.0;
  for ( GSIZET e=0; e<grid.elems().size(); e++ ) {

    ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
    u.range(ibeg, iend); // restrict global vecs to local range
    du.range(ibeg, iend);
    utmp.range(ibeg, iend);
    for ( GSIZET j=0; j<dXidX->size(2); j++ ) {
      for ( GSIZET i=0; i<dXidX->size(1); i++ ) { 
        (*dXidX)(i,j).range(ibeg,iend);
      }
    }
    for ( GSIZET k=0; k<GDIM ; k++ ) {
      gbasis[k] = (*gelems)[e]->gbasis(k);
      N[k]= (*gelems)[e]->size(k);
      xi_ev[k] = gbasis[k]->getXiNodes();
    }
    gshapefcn->set_basis(gbasis);
    // Compute derivative in reference space, D_k then
    // compute du/dx_j = Sum_k dXi_k/dx^j D_k u,
    for ( GINT r=0; r<nref; r++ ) {

      if ( grid.gtype() == GE_REGULAR  && r != (jder-1) ) continue;

      utmp  = 0.0;
#if defined(_G_IS2D)
      for ( GINT j=0, n=0; j<N[1]; j++ ) { 
        for (GINT i=0; i<N[0]; i++, n++ ) { 
          I[0] = i; I[1] = j;
          gshapefcn->dNdXi(I, r+1, xi_ev, dNi);
          dNi  *= u[n];  
          utmp += dNi;
        } // i-loop
      } // j-loop
#elif defined(_G_IS3D)
      for ( GINT k=0, n=0; k<N[2]; k++ ) { 
        for ( GINT j=0; j<N[1]; j++ ) { 
          for (GINT i=0; i<N[0]; i++, n++ ) { 
            I[0] = i; I[1] = j; I[2] = k;
            gshapefcn->dNdXi(I, r+1, xi_ev, dNi);
            dNi  *= u[n];  
            utmp += dNi;
          } // i-loop
        } // k-loop
      } // k-loop
#endif
//cout << serr << " ref_deriv[" << r << "]=" << utmp << endl;
      if ( grid.gtype() == GE_REGULAR ) {
//cout << serr << " dXidX[" << jder << "]=" << (*dXidX)(jder-1, 0) << endl;
        utmp.pointProd((*dXidX)(jder-1, 0));
        du += utmp;
      }
      else {
        utmp.pointProd((*dXidX)(r,jder-1));
//cout << serr << " du_partial[" << r << "]=" << utmp << endl;
        du += utmp;
      }

    } // end, reference deriv loop

  } // end, elem loop
  u.range_reset();
  du.range_reset();
  utmp.range_reset();
  for ( GSIZET j=0; j<dXidX->size(2); j++ ) {
    for ( GSIZET i=0; i<dXidX->size(1); i++ ) { 
      (*dXidX)(i,j).range_reset();
    }
  }


  delete gshapefcn;

} // end, shape_deriv

void print(GGrid &grid, GTVector<GFTYPE> &u, const GString fpref, GC_COMM &comm)
{
  GINT                         nxy;
  std::vector<GString>         scoord = {"x","y","z"};
  GTVector<GINT>               istate(1);
  GTVector<GString>            snames(1);
  GTVector<GString>            gnames(GDIM+1);
  GTVector<GTVector<GFTYPE>*>  ustate(1);
  GElemList                   *gelems = &grid.elems();
  GIOTraits                    iot;
  std::stringstream            ss;

  ustate[0] = &u;
  snames[0] = fpref;

  iot.nelems = grid_->nelems();
  iot.gtype  = grid_->gtype();
  iot.porder.resize(1,GDIM);

  istate = 0;

  for ( GSIZET j=0; j<GDIM; j++ ) iot.porder(0,j) = (*gelems)[0]->size(j)-1;

  // Write state:
  gio_write_state(iot, grid, ustate, istate, snames, comm);

  // Write grid:
  nxy = grid.gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;
  for ( GSIZET j=0; j<nxy; j++ ) {
    ss << scoord[j].c_str() << "grid";
    gnames[j] = ss.str();
    ss.str("");
  }
  gio_write_grid(iot, grid, gnames, comm);

  GTMatrix<GTVector<GFTYPE>> *dxidx = &grid.dXidX();
  for ( GSIZET j=0; j<dxidx->size(2); j++ ) {
    for ( GSIZET i=0; i<dxidx->size(1); i++ ) {
      ss << "dXidX" << i+1 << j+1;
      snames[0] = ss.str();
      ustate[0] = &(*dxidx)(i,j);
      gio_write_state(iot, grid, ustate, istate, snames, comm);
      ss.str("");
    }
  }

  GTVector<GFTYPE> *jac= &grid.Jac();
  ss << "Jac";
  snames[0] = ss.str();
  ustate[0] = jac;
  gio_write_state(iot, grid, ustate, istate, snames, comm);
  ss.str("");

} // end, print method
