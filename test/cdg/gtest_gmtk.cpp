//==================================================================================
// Module       : gtest_gmtk.cpp
// Date         : 7/24/18 (DLR)
// Description  : GeoFLOW test of GMTK
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <vector>
#include "gexec.h"
#include "gcomm.hpp"
#include "gtvector.hpp"
#include "gtmatrix.hpp"
#include "gtpoint.hpp"
#include "gmtk.hpp"
#include "gllbasis.hpp"

#if defined(_G_USE_GPTL)
    #include <gptl.h>
#endif

int main(int argc, char **argv)
{
    // Initialize comm:
    GComm::InitComm(&argc, &argv);

    GString serr ="main: ";
    GBOOL  bret;
    GINT   nrpt=1;
    GINT   iopt;
    GINT   errcode=0, gerrcode;
    GSIZET ne=4;    
    GSIZET np=1;    // elem 'order'
    GC_COMM comm = GC_COMM_WORLD;


#if 1

    // Parse command line. ':' after char
    // option indicates that it takes an argument:
    while ((iopt = getopt(argc, argv, "i:n:p:h")) != -1) {
      switch (iopt) {
      case 'i': // get iteration count
          nrpt = atoi(optarg);
          break;
      case 'p': // get array/vector size 
          np = atoi(optarg);
          break;
      case 'n': // # elems per task
          ne = atoi(optarg);
          break;
      case 'h': // help
          std::cout << "usage: " << std::endl <<
          argv[0] << " [-h] [-n #Elems per task] [-p expansion order] [-i #Repeat] " << std::endl;
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

    GTVector<GSIZET> N(3);
    N[0] = np+1; N[1] = np+2; N[2] = np+3;
    // Check push_back, reserve:
    GTVector<GDOUBLE>  u2d (N[0]*N[1]);
    GTVector<GDOUBLE>  u2da(N[0]*N[1]);
    GTVector<GDOUBLE>  E2  (N[0]*N[1]);
    GTVector<GDOUBLE>  y2  (N[0]*N[1]);
    GTVector<GDOUBLE>  u3d (N[0]*N[1]*N[2]);
    GTVector<GDOUBLE>  u3da(N[0]*N[1]*N[2]);
    GTVector<GDOUBLE>  E3  (N[0]*N[1]*N[2]);
    GTVector<GDOUBLE>  y3  (N[0]*N[1]*N[2]);
    GTMatrix<GDOUBLE>  D1  (N[0],N[0]);
    GTMatrix<GDOUBLE>  I1  (N[0],N[0]);
    GTMatrix<GDOUBLE>  D2  (N[1],N[1]);
    GTMatrix<GDOUBLE>  I2  (N[1],N[1]);
    GTMatrix<GDOUBLE>  D2T (N[1],N[1]);
    GTMatrix<GDOUBLE>  D3  (N[2],N[2]);
    GTMatrix<GDOUBLE>  I3  (N[2],N[2]);
    GTMatrix<GDOUBLE>  D3T (N[2],N[2]);

    I1.createIdentity();
    I2.createIdentity();
    I3.createIdentity();

    // Set 1d matrices, solution vector:
    for ( GSIZET j=0; j<N[0]; j++ ) {
      for ( GSIZET i=0; i<N[0]; i++ ) {
        D1(i,j) = 2.0*static_cast<GDOUBLE>(i)+N[0]*static_cast<GDOUBLE>(j);
      }
    }
    for ( GSIZET j=0; j<N[1]; j++ ) {
      for ( GSIZET i=0; i<N[1]; i++ ) {
        D2(i,j) = 2.0*static_cast<GDOUBLE>(i)+N[0]*static_cast<GDOUBLE>(j);
      }
    }
    for ( GSIZET j=0; j<N[2]; j++ ) {
      for ( GSIZET i=0; i<N[2]; i++ ) {
        D3(i,j) = 2.0*static_cast<GDOUBLE>(i)+N[0]*static_cast<GDOUBLE>(j);
      }
    }
    for ( GSIZET j=0; j<N[1]; j++ ) {
      for ( GSIZET i=0; i<N[0]; i++ ) {
        u2d[i+j*N[0]] = 2.0*static_cast<GDOUBLE>(i)+N[0]*static_cast<GDOUBLE>(j);
      }
    }
    for ( GSIZET k=0; k<N[2]; k++ ) {
      for ( GSIZET j=0; j<N[1]; j++ ) {
        for ( GSIZET i=0; i<N[0]; i++ ) {
          u3d[i+j*N[0]+k*N[0]*N[1]] = 2.0*static_cast<GDOUBLE>(i)+N[0]*static_cast<GDOUBLE>(j);
        }
      }
    }
    D2.transpose(D2T);
    D3.transpose(D3T);

std::cout << "D1=" << D1 << std::endl;
std::cout << "D2T=" << D2T << std::endl;


    // Create solution vector for 2d:
    // First, create full matrix:
    GSIZET l, m, n;
    GTMatrix<GDOUBLE> BIG2(N[0]*N[1],N[0]*N[1]);
    GTVector<GDOUBLE> tmp(N[0]*N[1]);
    for ( GSIZET j=0; j<N[1]; j++ ) {
      for ( GSIZET i=0; i<N[1]; i++ ) {

        for ( GSIZET m=0; m<N[0]; m++ ) {
          for ( GSIZET l=0; l<N[0]; l++ ) {
            BIG2(l+i*N[0],m+j*N[0]) = D2(i,j)*D1(l,m);
          }  
        }  

      }
    }
    u2da = BIG2 * u2d;

    GMTK::D2_X_D1(D1, D2T, u2d, tmp, y2); 
    E2 = y2;
    E2 -= u2da;
  
#if 0
std::cout << "BIG2 = " << BIG2 << std::endl;
std::cout << "y2=" << y2 << std::endl; 
std::cout << "y2a  = " << u2da << std::endl;
#endif
 
   if ( E2.Eucnorm() > 0 ) {
      std::cout << "main: -------------------------------------D2_X_D1 FAILED" << std::endl;
      errcode = 1;
   } else {
      std::cout << "main: -------------------------------------D2_X_D1 OK" << std::endl;
   }

    for ( GSIZET j=0; j<N[1]; j++ ) {
      for ( GSIZET i=0; i<N[1]; i++ ) {

        for ( GSIZET m=0; m<N[0]; m++ ) {
          for ( GSIZET l=0; l<N[0]; l++ ) {
            BIG2(l+i*N[0],m+j*N[0]) = I2(i,j)*D1(l,m);
          }  
        }  

      }
    }
    u2da = BIG2 * u2d;

    GMTK::I2_X_D1(D1, u2d, N[0], N[1], y2); 
    E2 = y2;
    E2 -= u2da;
  
#if 1
std::cout << "BIG2 = " << BIG2 << std::endl;
std::cout << "y2=" << y2 << std::endl; 
std::cout << "y2a  = " << u2da << std::endl;
#endif
 
   if ( E2.Eucnorm() > 0 ) {
      std::cout << "main: -------------------------------------I2_X_D1 FAILED" << std::endl;
      errcode = 1;
   } else {
      std::cout << "main: -------------------------------------I2_X_D1 OK" << std::endl;
   }

    for ( GSIZET j=0; j<N[1]; j++ ) {
      for ( GSIZET i=0; i<N[1]; i++ ) {

        for ( GSIZET m=0; m<N[0]; m++ ) {
          for ( GSIZET l=0; l<N[0]; l++ ) {
            BIG2(l+i*N[0],m+j*N[0]) = D2(i,j)*I1(l,m);
          }  
        }  

      }
    }
    u2da = BIG2 * u2d;

    GMTK::D2_X_I1(D2T, u2d, N[0], N[1], y2); 
    E2 = y2;
    E2 -= u2da;
  
#if 0
std::cout << "BIG2 = " << BIG2 << std::endl;
std::cout << "y2=" << y2 << std::endl; 
std::cout << "y2a  = " << u2da << std::endl;
#endif
 
   if ( E2.Eucnorm() > 0 ) {
      std::cout << "main: -------------------------------------D2_X_I1 FAILED" << std::endl;
      errcode = 1;
   } else {
      std::cout << "main: -------------------------------------D2_X_I1 OK" << std::endl;
   }

    // Create solution vector for 3d:
    // First, create full matrix:
    GTMatrix<GDOUBLE> BIG3(N[0]*N[1]*N[2],N[0]*N[1]*N[2]);
    for ( GSIZET j=0; j<N[2]; j++ ) {
      for ( GSIZET i=0; i<N[2]; i++ ) {

        for ( GSIZET m=0; m<N[0]*N[1]; m++ ) {
          for ( GSIZET l=0; l<N[0]*N[1]; l++ ) {
            BIG3(l+i*N[0]*N[1],m+j*N[0]*N[1]) = D3(i,j)*BIG2(l,m);
          }
        }

      }
    }
    u3da = BIG3 * u3d;

    GMTK::D3_X_D2_X_D1(D1, D2T, D3T, u3d, tmp, y3);
    E3 = y3;
    E3 -= u3da;

#if 0
std::cout << "BIG3 = " << BIG3 << std::endl;
std::cout << "y3   = " << y3   << std::endl;
std::cout << "y3a  = " << u3da << std::endl;
#endif
   for ( GSIZET j=0; j<y3.size(); j++ ) {
     if ( y3[j] != u3da[j] ) std::cout << j << " :  del u=" << y3[j]-u3da[j] << std::endl;
   }

   if ( E3.Eucnorm() > 0 ) {
      std::cout << "main: --------------------------------D3_X_D2_X_D1 FAILED" << std::endl;
      std::cout << "                                      Eucnorm=" << E3.Eucnorm() << std::endl;
      errcode = 2;
   } else {
      std::cout << "main: --------------------------------D3_X_D2_X_D1 OK" << std::endl;
   }


#if 1
   GLLBasis<GCTYPE,GFTYPE> gbasis(N[0]-1);
   GLLBasis<GCTYPE,GFTYPE> gobasis(N[0]+1);
   GTVector<GFTYPE> u, v, t;
   GTMatrix<GFTYPE> M1, M2, M2T;
   M1.resize(gobasis.getOrder()+1,gbasis.getOrder()+1);
   M2.resize(gobasis.getOrder()+1,gbasis.getOrder()+1);
   M2T.resize(M2.size(2),M2.size(1));

   gbasis.evalBasis(*(gobasis.getXiNodes()),M1);
   gbasis.evalBasis(*(gobasis.getXiNodes()),M2);
   M2.transpose(M2T);
   u.resize((gbasis.getOrder()+1)*(gbasis.getOrder()+1));
   v.resize((gobasis.getOrder()+1)*(gobasis.getOrder()+1));
   t.resize((gobasis.getOrder()+1)*(gobasis.getOrder()+1));
   u = 3.0;
   GMTK::D2_X_D1(M1, M2T, u, t, v); 
   std::cout << "main: M1=" << M1 << std::endl;
   std::cout << "main: M2=" << M2 << std::endl;
   std::cout << "main: u =" << u  << std::endl;
   std::cout << "main: --------------------------------interp_chk: " << v <<  std::endl;
   std::cout << "main: --------------------------------M1*M2: " << M1*M2 <<  std::endl;
#endif
  


  gerrcode = errcode;
  
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
