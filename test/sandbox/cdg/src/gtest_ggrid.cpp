//==================================================================================
// Module       : gtest_ggrid.cpp
// Date         : 7/24/18 (DLR)
// Description  : GeoFLOW tests of GGrid classes
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

#include "gexec.h"
#include "gtypes.h"
#include <cstdio>
#include <unistd.h>
#include <iostream>
#include "gptl.h"
#include <random>
#include "gcomm.hpp"
#include "gllbasis.hpp"
#include "ggrid_icos.hpp"


int main(int argc, char **argv)
{
    // Initialize comm:
    GComm::InitComm(&argc, &argv);

    GString serr ="main: ";
    GBOOL  bret;
    GINT   iopt;
    GINT   ilevel=0;// 2d ICOS refinement level
    GINT   errcode=0, gerrcode;
    GINT   nlat=10, nlong=20;
    GINT   np=1;    // elem 'order'
    GFTYPE radiusi=1, radiuso=2;
    GFTYPE xnoise = 1.0e-2; // noise level in points around centroids
    GTVector<GINT> ne[3]; // # elements in each direction in 3d
    GC_COMM comm = GC_COMM_WORLD;

    for ( GSIZET j=0; j<3; j++ ) ne[j] = 10;


#if 1

    // Parse command line. ':' after char
    // option indicates that it takes an argument:
    while ((iopt = getopt(argc, argv, "i:j:k;l:n:p:q:r:h")) != -1) {
      switch (iopt) {
      case 'i': // get # elements in r
          ne[0] = atoi(optarg);
          break;
      case 'j': // get # elements in lat
          ne[1] = atoi(optarg);
          break;
      case 'k': // get # elements in long
          ne[2] = atoi(optarg);
          break;
      case 'l': // # 2d refinement level
          ilevel = atoi(optarg);
          break;
      case 'n': // # 2d refinement level
          xnoise = atof(optarg);
          break;
      case 'p': // get nodal exp order
          np = atoi(optarg);
          break;
      case 'q': // inner radius for 2d/3d
          radiusi = atoi(optarg);
          break;
      case 'r': // outer radius for 3d
          radiuso = atoi(optarg);
          break;
      case 'h': // help
          std::cout << "usage: " << std::endl <<
          argv[0] << " [-h] [-i #Elems in r] [-j #Elems in lat]  [-k #Elems in long] [-l refine level] [-n noise_level] [-p expansion order] [-q rad_inner] [-r rad_outer]  " << std::endl;
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

    GINT myrank  = GComm::WorldRank();
    GINT nprocs  = GComm::WorldSize();

    // Set GTPL options:
    GPTLsetoption (GPTLcpu, 1);

    // Initialize GPTL:
    GPTLinitialize();


    // Create basis:
    GTVector<GNBasis<GCTYPE,GFTYPE>*> gbasis(GDIM);
    for ( GSIZET k=0; k<GDIM; k++ ) gbasis[k] = new GLLBasis<GCTYPE,GFTYPE>(np);
    
    GGrid grid(comm);

    GPTLstart("gen_base_grid");

    GGridIcos::Traits traits;
    #if defined(_G_IS2D)
    traits.radiusi = radiusi;
    traits.ilevel= ilevel;
    GGridIcos gen_icos(traits, gbasis, nprocs);
    #elif defined(_G_IS3D)
    traits.radiusi = radiusi;
    traits.radiuso = radius0;
    GGridIcos gen_icos(traits, ne, gbasis, nprocs);
    #endif

    GPTLstop("gen_base_grid");

    // Test some of the coord transformation methods:
    GFTYPE xlat, xlatc, xlong, xlongc;
    GFTYPE eps = std::numeric_limits<GFTYPE>::epsilon();
    GTPoint<GFTYPE> pdiff;
    GTVector<GTPoint<GFTYPE>> cart(1), cref(1), sph(1), tcart(1), tsph(1);
    GTVector<GTPoint<GFTYPE>> gnom(1);

    // Gnomonic points are only 2 dimensional:
    for ( GSIZET j=0; j<gnom.size(); j++ ) gnom[j].resize(2);

    eps *= 1000.0;
    GPP(comm,"eps=" << eps);

    // Check a range of reference locations (xlatc, xlongc) to
    // check gnomonic transform:
    GSIZET nbad=0;
    GFTYPE dlat  = PI/(nlat-1);
    GFTYPE dlong = 2.0*PI/(nlong-1);
    GFTYPE error,  emin=std::numeric_limits<GFTYPE>::max(), emax=0.0; // max error
    std::default_random_engine rgen; 
    std::uniform_real_distribution<GFTYPE> dist(-xnoise,xnoise);
    for ( GSIZET j=0; j<nlat; j++ ) {
      xlatc = PI/2.0 - j*dlat; // go to next reference latitude
      for ( GSIZET i=0; i<nlong; i++ ) {
        xlongc = i*dlong; // go to next references longitude
        cref[0].x1 = radiusi; cref[0].x2 = xlatc; cref[0].x3 = xlongc; // reference point
        cart[0].x1 = radiusi;
        cart[0].x2 = MIN( MAX(cref[0].x2 + dist(rgen),-PI/2.0), PI/2.0);
        cart[0].x3 = MIN( MAX(cref[0].x3 + dist(rgen),    0.0), 2.0*PI);
        cref[0].x2 *= 180.0/PI; cref[0].x3 *= 180.0/PI;
        gen_icos.spherical2xyz(cart);
        sph[0] = cart[0];
        gen_icos.xyz2spherical(sph);
        sph[0].x2 *= 180.0/PI; sph[0].x3 *= 180.0/PI;
        gen_icos.cart2gnomonic(cart, radiusi, xlatc, xlongc, gnom);
        gen_icos.gnomonic2cart(gnom, radiusi, xlatc, xlongc, tcart);
        tsph[0] = tcart[0];
        gen_icos.xyz2spherical(tsph);
        tsph[0].x2 *= 180.0/PI; tsph[0].x3 *= 180.0/PI;
        pdiff = tcart[0] - cart[0];
        error = pdiff.norm();
        emax = !std::isnan(error)? MAX(emax, error) : emax;
        emin = !std::isnan(error)? MIN(emax, error) : emin;
        if (std::isnan(error) || error > eps || std::isnan(error) ) { 
          errcode = 1;
          std::cout << "...............ilong=" << i << " jlat=" << j << ":  sph =" << sph [0] << std::endl;
          std::cout << "cref =" << cref[0] << std::endl;
          std::cout << "cart =" << cart[0] << std::endl;
          std::cout << "gnom =" << gnom[0] << std::endl;
          std::cout << "tcart=" << tcart[0] << std::endl;
          std::cout << "tsph =" << tsph [0] << std::endl;
          std::cout << "diff norm=" << error   << std::endl;
          nbad++;
if ( std::isnan(error) ) exit(1);
        }     
      }
    }

    GPP(comm,"noise level=" << xnoise);
    GPP(comm,"ntotal = " << nlat*nlong <<  " nbad =" << nbad);
    GPP(comm,"max error = " << emax);
    GPP(comm,"min error = " << emin);
    if ( errcode > 0 ) {
      GPP(comm,"main: ----------------------------gnomonic transform  FAILED");
      exit(errcode);
    } else {
      GPP(comm,"main: ----------------------------gnomonic transform OK");
    }

    char    sbuff[1024];
    GString supp, fname;

    GPTLstart("print_base_grid");
    // Print base grid:
    #if defined(_G_IS2D)
    sprintf(sbuff, "level%d_p%d", ilevel, np);
    supp = sbuff;
    fname = "tgrid_" + supp + ".txt";
    gen_icos.print(fname,GICOS_CART);
    #endif
    GPTLstop("print_base_grid");


    GPTLstart("gen_elem_grid");
    // Generate grid:
    gen_icos.do_grid(grid, myrank);
    GPP(comm,"main: ----------------------------gelem.size=" << grid.elems().size());
    GPTLstop("gen_elem_grid");

   GPP(comm,"main: ----------------------------gelem.minlen=" << grid.minlength());
   GPP(comm,"main: ----------------------------gelem.maxlen=" << grid.maxlength());

    // Print grid:
    sprintf(sbuff, "level%d_p%d", ilevel, np);
    supp = sbuff;
    fname = "grid_wire_" + supp + ".txt";
    GPTLstart("print_wire_grid");
    if ( myrank == 0 ) {
      grid.print(fname); // print vertices only
    }
    GPTLstop("print_wire_grid");

    GPTLstart("print_elem_grid");
    fname = "grid_" + supp + ".txt";
    if ( myrank == 0 ) {
      grid.print(fname, TRUE); // print internal dof too
    }
    GPTLstop("print_elem_grid");

    // Accumulate errors:
    GComm::Allreduce(&errcode, &gerrcode, 1, T2GCDatatype<GINT>() , GC_OP_MAX, comm);

 
    GPTLpr_file("timing.txt");
    GPTLfinalize();


term:
    if ( gerrcode != 0 ) {
      GPP(comm,serr << " Error: code=" << errcode);
    }
    else {
      GPP(comm,serr << "     Success!");
    }

    GComm::TermComm();
    return(gerrcode);
} // end, main
