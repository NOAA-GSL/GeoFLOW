//==================================================================================
// Module       : gtest_helm.cpp
// Date         : 2/24/19 (DLR)
// Description  : GeoFLOW test of Helmholtz operator
// Copyright    : Copyright 2019. Colorado State University. All rights reserved.
// Derived From : 
//==================================================================================

#include "gexec.h"
#include "gtypes.h"
#include <cstdio>
#include <unistd.h>
#include <iostream>
#if defined(GEOFLOW_USE_GPTL)
  #include "gptl.h"
#endif
#include <memory>
#include <cstdlib>
#include <cassert>
#include <random>
#include "gcomm.hpp"
#include "ggfx.hpp"
#include "gllbasis.hpp"
#include "gmorton_keygen.hpp"
#include "ggrid_factory.hpp"
#include "gmtk.hpp"
#include "ghelmholtz.hpp"
#include "tbox/property_tree.hpp"
#include "tbox/mpixx.hpp"
#include "tbox/global_manager.hpp"
#include "tbox/input_manager.hpp"
#include "gtools.h"

using namespace geoflow::tbox;
using namespace std;


GGrid *grid_ = NULLPTR;

int main(int argc, char **argv)
{

    GString serr ="main: ";
    GBOOL  bret;
    GINT   iopt;
    GINT   ilevel=0;// 2d ICOS refinement level
    GINT   np=1;    // elem 'order'
    GINT   nstate=GDIM;  // number 'state' arrays
    GSIZET maxSteps;
    GFTYPE radiusi=1, radiuso=2;
    std::vector<GINT> ne(3); // # elements in each direction in 3d
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
    PropertyTree polyptree; // polynomial prop tree

    ptree.load_file("input.jsn");       // main param file structure
    // Create other prop trees for various objects:
    sgrid       = ptree.getValue<GString>("grid_type");
    np          = ptree.getValue<GINT>("exp_order");
    
    assert(sgrid=="grid_box" && "Box grid must be specified");
    gridptree   = ptree.getPropertyTree(sgrid);
    polyptree   = ptree.getPropertyTree("poly_test");

    GTPoint<GFTYPE> P0(3), dP(3);
    ne          = gridptree.getArray<GINT>("num_elems");  // may be modified by command line
    std::vector<GFTYPE> xyz0  = gridptree.getArray<GFTYPE>("xyz0");
    std::vector<GFTYPE> dxyz0 = gridptree.getArray<GFTYPE>("delxyz");
    P0 = xyz0;
    dP = dxyz0;

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

#if defined(GEOFLOW_USE_GPTL)
    // Set GTPL options:
    GPTLsetoption (GPTLcpu, 1);

    // Initialize GPTL:
    GPTLinitialize();
#endif


    // Create basis:
    GTVector<GNBasis<GCTYPE,GFTYPE>*> gbasis(GDIM);
    for ( GSIZET k=0; k<GDIM; k++ ) {
      gbasis [k] = new GLLBasis<GCTYPE,GFTYPE>(np);
    }
    

#if defined(GEOFLOW_USE_GPTL)
    GPTLstart("gen_grid");
#endif
    // Create grid:
    grid_ = GGridFactory::build(gridptree, gbasis, comm);

#if defined(GEOFLOW_USE_GPTL)
    GPTLstop("gen_grid");
#endif


#if defined(GEOFLOW_USE_GPTL)
    GPTLstart("do_gather_op");
#endif

    // Initialize gather/scatter operator:
    GGFX ggfx;
    init_ggfx(ptree, *grid_, ggfx);

#if defined(GEOFLOW_USE_GPTL)
    GPTLstop("do_gather_op");
#endif


    // Create state and tmp space:
    GTVector<GTVector<GFTYPE>*> utmp(5); // tmp space
    GTVector<GTVector<GFTYPE>*> u   (1); // field
    GTVector<GTVector<GFTYPE>*> du  (1); // computed solutn 
    GTVector<GTVector<GFTYPE>*> da  (1); // analytic Laplacian
    
    for ( GSIZET j=0; j<utmp.size(); j++ ) utmp[j] = new GTVector<GFTYPE>(grid_->size());
    for ( GSIZET j=0; j<u   .size(); j++ ) u   [j] = new GTVector<GFTYPE>(grid_->size());
    for ( GSIZET j=0; j<du  .size(); j++ ) du  [j] = new GTVector<GFTYPE>(grid_->size());
    for ( GSIZET j=0; j<da  .size(); j++ ) da  [j] = new GTVector<GFTYPE>(grid_->size());
    // Instantiate Helmholtz operator:
    GINT       nt = utmp.size();
    GHelmholtz helm(*grid_);
//  GMass      mass(*grid_);

    // Initialize u = x^p y^q z^r:
    GFTYPE p = polyptree.getValue<GFTYPE>("xpoly"); // exponent on x
    GFTYPE q = polyptree.getValue<GFTYPE>("ypoly"); // exponent on y
    GFTYPE r = polyptree.getValue<GFTYPE>("zpoly"); // exponent on z
    GFTYPE x, y, z=1.0;
    GTVector<GFTYPE> etmp1;
    GTVector<GTVector<GFTYPE>> *xnodes = &grid_->xNodes();   
    GSIZET nxy = (*xnodes)[0].size();

    assert( p>=0 && q>=0 && r >=0 && "Exponents must be >=0");

    // Initialize field, u = x^p y^q z^r;
    // Also, initialize collocated analytic Lap u:
    if ( GDIM < 3 ) r = 0;
    for ( GSIZET j=0; j<nxy; j++ ) {
      x = (*xnodes)[0][j] == 0.0 ? 1.0 : (*xnodes)[0][j];
      y = (*xnodes)[1][j] == 0.0 ? 1.0 : (*xnodes)[1][j];
      if ( GDIM == 3 ) z = (*xnodes)[2][j] == 0.0 ? 1.1 : (*xnodes)[2][j];
      (*u [0])[j] = pow(x,p)*pow(y,q)*pow(z,r);
      (*da[0])[j] = 0.0;
      if ( p> 0 ) (*da[0])[j] += p*(p-1)*pow(x,p-2)*pow(y,q)*pow(z,r);
      if ( q> 0 ) (*da[0])[j] += q*(q-1)*pow(x,p)*pow(y,q-2)*pow(z,r);
      if ( r> 0 ) (*da[0])[j] += r*(r-1)*pow(x,p)*pow(y,q)*pow(z,r-2);
    }

    helm.opVec_prod(*u [0], utmp, *du[0]);

    cout << "main: u =" << *u[0] << endl;
    cout << "main: Hu=" << *du[0] << endl;

    // Compute analytic integral of Laplacian acting on u:
    GFTYPE da_int, err;
#if 1
    GFTYPE x0=P0.x1, x1=P0.x1+dP.x1;
    GFTYPE y0=P0.x2, y1=P0.x2+dP.x2;
    GFTYPE z0=P0.x3, z1=P0.x3+dP.x3;
    da_int = 0.0;
    if ( GDIM == 2 ) {
      if ( p>0 ) da_int += p/(q+1)*(pow(x1,p-1)-pow(x0,p-1)) * (pow(y1,q+1)-pow(y0,q+1));
      if ( q>0 ) da_int += q/(p+1)*(pow(x1,p+1)-pow(x0,p+1)) * (pow(y1,q-1)-pow(y0,q-1));
    } else if ( GDIM == 3 ){
      if ( p>0 ) da_int += q/((p+1)*(r+1))*(pow(x1,p+1)-pow(x0,p+1)) * (pow(y1,q-1)-pow(y0,q-1)) * (pow(z1,r+1)-pow(z0,r+1));
      if ( q>0 ) da_int += q/((p+1)*(r+1))*(pow(x1,p+1)-pow(x0,p+1)) * (pow(y1,q-1)-pow(y0,q-1)) * (pow(z1,r+1)-pow(z0,r+1));
      if ( r>0 ) da_int += r/((p+1)*(q+1))*(pow(x1,p+1)-pow(x0,p+1)) * (pow(y1,q+1)-pow(y0,q+1)) * (pow(z1,r-1)-pow(z0,r-1));
    }
#else
    // Semi-analytic integral of collocated Laplacian u:
    da_int = grid_->integrate(*da[0],*utmp[0]);
#endif

    // Compute numerical integral of (weak) Laplacian 
    // (Already contains mass and Jacobian):
    GFTYPE eps=10.0*std::numeric_limits<GFTYPE>::epsilon();
    GTVector<GFTYPE> lnorm(1), gnorm(1);
    lnorm[0] = du[0]->sum();
    GComm::Allreduce(lnorm.data()  , gnorm.data()  , 1, T2GCDatatype<GFTYPE>() , GC_OP_SUM, comm);
    cout << "main: analyt_int=" << da_int << " comp_int=" << gnorm[0] << endl;

    err = (da_int-gnorm[0])/(da_int+eps);
    cout << "main: .................error=" << err << endl;
    

    // Get add'l info for error file:
    GSIZET lnelems=grid_->nelems();
    GSIZET gnelems;
    GComm::Allreduce(&lnelems, &gnelems, 1, T2GCDatatype<GSIZET>() , GC_OP_SUM, comm);

    // Print convergence data to file:
    std::ifstream itst;
    std::ofstream ios;
    itst.open("helm_err.txt");
    ios.open("helm_err.txt",std::ios_base::app);

    // Write header, if required:
    if ( itst.peek() == std::ofstream::traits_type::eof() ) {
    ios << "# p       num_elems        Rel_err " << std::endl;
    }
    itst.close();

    ios << np  << "     "  << gnelems << "     " << err << std::endl;
    ios.close();

 
#if defined(GEOFLOW_USE_GPTL)
    GPTLpr_file("timing.txt");
    GPTLfinalize();
#endif

    GComm::TermComm();
    if ( grid_ != NULLPTR ) delete grid_;

    return(0);

} // end, main

