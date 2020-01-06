//==================================================================================
// Module       : gtest_blas.cpp
// Date         : 6/30/18 (DLR)
// Description  : GeoFLOW BLAS test
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

#include "gexec.h"
#include "gtypes.h"
#include <cstdio>
#include <unistd.h>
#include <iostream>
#include "gtvector.hpp"
#include "gtmatrix.hpp"

#if defined(_G_USE_GPTL)
    #include "gptl.h"
#endif

int main(int argc, char **argv)
{
    GString serr ="main: ";
    GBOOL   bret;
    GINT    errcode=0;

#if defined(_G_USE_GPTL)
    // Set GTPL options:
    GPTLsetoption (GPTLcpu, 1);

    // Initialize GPTL:
    GPTLinitialize();
#endif

    // Set matrices:
    GFLOAT a[16] = { 1 , 5 , 9 , 13,
                     2 , 6 , 10, 14, 
                     3 , 7 , 11, 15, 
                     4 , 8 , 12, 16 };
    GFLOAT b[16] = { 17, 21, 25, 29,
                     18, 22, 26, 30,
                     19, 23, 27, 31,
                     20, 24, 28, 32 };

    GTMatrix<GFLOAT> A(4,4);
    GTMatrix<GFLOAT> A34(3,4);
    GTMatrix<GFLOAT> B(4,4);
    GTMatrix<GFLOAT> C(4,4);
    GTMatrix<GFLOAT> D(4,4);
    GTMatrix<GFLOAT> D34(3,4);
    GTMatrix<GFLOAT> D43(4,3);
    GTMatrix<GFLOAT> E(4,4);
    GTMatrix<GFLOAT> E34(3,4);
    GTMatrix<GFLOAT> E43(4,3);
    GTMatrix<GFLOAT> F(4,4);
    GTMatrix<GFLOAT> F34(3,4);
    GTMatrix<GDOUBLE> Adbl(4,4);

    GTMatrix<GFLOAT> A3(3,3);
    GTMatrix<GFLOAT> B3(3,3);
    GTMatrix<GFLOAT> C3(3,3);
    GTMatrix<GFLOAT> D3(3,3);
    GTMatrix<GFLOAT> E3(3,3);

    GTVector<GFLOAT> d(4);
    GTVector<GFLOAT> e(4);
    GTVector<GFLOAT> u(4);
    GTVector<GFLOAT> v(4);

    A3(0,0) = 1 ; A3(0,1) = 0 ; A3(0,2) = 2 ; 
    A3(1,0) =-1 ; A3(1,1) = 5 ; A3(1,2) = 0 ;
    A3(2,0) = 0 ; A3(2,1) = 3 ; A3(2,2) = -9;

    A(0,0) = 1 ; A(0,1) = 2 ; A(0,2) = 3 ; A(0,3) = 4 ; 
    A(1,0) = 5 ; A(1,1) = 6 ; A(1,2) = 7 ; A(1,3) = 8 ; 
    A(2,0) = 9 ; A(2,1) = 10; A(2,2) = 11; A(2,3) = 12; 
    A(3,0) = 13; A(3,1) = 14; A(3,2) = 15; A(3,3) = 16; 

    A34(0,0) = 1 ; A34(0,1) = 2 ; A34(0,2) = 3 ; A34(0,3) = 4 ; 
    A34(1,0) = 5 ; A34(1,1) = 6 ; A34(1,2) = 7 ; A34(1,3) = 8 ; 
    A34(2,0) = 9 ; A34(2,1) = 10; A34(2,2) = 11; A34(2,3) = 12; 

    B(0,0) = 17; B(0,1) = 18; B(0,2) = 19; B(0,3) = 20; 
    B(1,0) = 21; B(1,1) = 22; B(1,2) = 23; B(1,3) = 24; 
    B(2,0) = 25; B(2,1) = 26; B(2,2) = 27; B(2,3) = 28; 
    B(3,0) = 29; B(3,1) = 30; B(3,2) = 31; B(3,3) = 32; 

    v  [0] = 3 ; v  [1] = 5 ; v  [2] = 7 ; v  [3] = 9 ;

#if defined(_G_USE_GPTL)
    GPTLstart("blas_stuff");
#endif

    // A + B  solution:
    D(0,0) = 18; D(0,1) = 20; D(0,2) = 22; D(0,3) = 24; 
    D(1,0) = 26; D(1,1) = 28; D(1,2) = 30; D(1,3) = 32; 
    D(2,0) = 34; D(2,1) = 36; D(2,2) = 38; D(2,3) = 40; 
    D(3,0) = 42; D(3,1) = 44; D(3,2) = 46; D(3,3) = 48; 

    C = A + B;
    E = C - D;
    if ( E.data().Eucnorm() > 0 ) {
      std::cout << "main: -------------------------------------A+B FAILED" << std::endl;
      errcode = 1;
    } else {
      std::cout << "main: -------------------------------------A+B OK" << std::endl;
    }

    // A - B solution:
    D.set(-16.0);

    C = A - B;
    E = C - D;
    if ( E.data().Eucnorm() > 0 ) {
      std::cout << "main: -------------------------------------A-B FAILED" << std::endl;
      errcode = 2;
    } else {
      std::cout << "main: -------------------------------------A-B OK" << std::endl;
    }
    
    // A += 3 solution:
    C = A;
    D(0,0) = 4 ; D(0,1) = 5 ; D(0,2) = 6 ; D(0,3) = 7 ; 
    D(1,0) = 8 ; D(1,1) = 9 ; D(1,2) = 10; D(1,3) = 11; 
    D(2,0) = 12; D(2,1) = 13; D(2,2) = 14; D(2,3) = 15; 
    D(3,0) = 16; D(3,1) = 17; D(3,2) = 18; D(3,3) = 19; 
    C += 3.0;
    E = C - D;
    if ( E.data().Eucnorm() > 0 ) {
      std::cout << "main: -------------------------------------A+=3 FAILED" << std::endl;
      errcode = 3;
    } else {
      std::cout << "main: -------------------------------------A+=3 OK" << std::endl;
    }

    // A *= 3 solution:
    C = A;
    D(0,0) = 3 ; D(0,1) = 6 ; D(0,2) = 9 ; D(0,3) = 12; 
    D(1,0) = 15; D(1,1) = 18; D(1,2) = 21; D(1,3) = 24; 
    D(2,0) = 27; D(2,1) = 30; D(2,2) = 33; D(2,3) = 36; 
    D(3,0) = 39; D(3,1) = 42; D(3,2) = 45; D(3,3) = 48; 
    C *= 3.0;
    E = C - D;
    if ( E.data().Eucnorm() > 0 ) {
      std::cout << "main: -------------------------------------A*=3 FAILED" << std::endl;
      errcode = 5;
    } else {
      std::cout << "main: -------------------------------------A*=3 OK" << std::endl;
    }

    // D = 4A - 5B solution:
    D(0,0) =-81; D(0,1) =-82; D(0,2) =-83; D(0,3) =-84; 
    D(1,0) =-85; D(1,1) =-86; D(1,2) =-87; D(1,3) =-88; 
    D(2,0) =-89; D(2,1) =-90; D(2,2) =-91; D(2,3) =-92; 
    D(3,0) =-93; D(3,1) =-94; D(3,2) =-95; D(3,3) =-96; 
    C = A*4 - B*5;
    E = C - D;
    if ( E.data().Eucnorm() > 0 ) {
      std::cout << "main: -------------------------------------A-B FAILED" << std::endl;
      errcode = 6;
    } else {
      std::cout << "main: -------------------------------------A-B OK" << std::endl;
    }

    // A += B solution:
    C = A;
    D(0,0) = 18; D(0,1) = 20; D(0,2) = 22; D(0,3) = 24; 
    D(1,0) = 26; D(1,1) = 28; D(1,2) = 30; D(1,3) = 32; 
    D(2,0) = 34; D(2,1) = 36; D(2,2) = 38; D(2,3) = 40; 
    D(3,0) = 42; D(3,1) = 44; D(3,2) = 46; D(3,3) = 48; 
    C += B;
    E = C - D;
    if ( E.data().Eucnorm() > 0 ) {
      std::cout << "main: -------------------------------------A+=B FAILED" << std::endl;
      errcode = 7;
    } else {
      std::cout << "main: -------------------------------------A+=B OK" << std::endl;
    }

    // C = A B solution:
    D(0,0) = 250 ; D(0,1) = 260 ; D(0,2) = 270 ; D(0,3) = 280 ; 
    D(1,0) = 618 ; D(1,1) = 644 ; D(1,2) = 670 ; D(1,3) = 696 ; 
    D(2,0) = 986 ; D(2,1) = 1028; D(2,2) = 1070; D(2,3) = 1112; 
    D(3,0) = 1354; D(3,1) = 1412; D(3,2) = 1470; D(3,3) = 1528; 
    C = A * B;
std::cout << "A  =" << A << std::endl;
std::cout << "B  =" << B << std::endl;
std::cout << "A*B=" << C << std::endl;
    E = C - D;
    if ( E.data().Eucnorm() > 0 ) {
      std::cout << "main: -------------------------------------A*B FAILED" << std::endl;
      errcode = 7;
    } else {
      std::cout << "main: -------------------------------------A*B OK" << std::endl;
    }

    // u = A v solution:
    d [0] = 70; d[1] = 166; d[2] = 262; d[3] = 358;
    u = A * v;
    e = u - d;
    if ( e.Eucnorm() > 0 ) {
      std::cout << "main: -------------------------------------A*v FAILED" << std::endl;
      errcode = 7;
    } else {
      std::cout << "main: -------------------------------------A*v OK" << std::endl;
    }

    // A-Transose solution:
    D(0,0) = 1 ; D(0,1) = 5 ; D(0,2) = 9 ; D(0,3) = 13; 
    D(1,0) = 2 ; D(1,1) = 6 ; D(1,2) = 10; D(1,3) = 14; 
    D(2,0) = 3 ; D(2,1) = 7 ; D(2,2) = 11; D(2,3) = 15; 
    D(3,0) = 4 ; D(3,1) = 8 ; D(3,2) = 12; D(3,3) = 16; 
    A.transpose(C);
std::cout << "Transpose(A)=" << C << std::endl;

    E = C - D;
    if ( E.data().Eucnorm() > 0 ) {
      std::cout << "main: -------------------------------------A Transpose FAILED" << std::endl;
      errcode = 7;
    } else {
      std::cout << "main: -------------------------------------A Transpose OK" << std::endl;
    }

    // A-Transose in place solution:
    A34(0,0) = 1 ; A34(0,1) = 2 ; A34(0,2) = 3 ; A34(0,3) = 4 ; 
    A34(1,0) = 5 ; A34(1,1) = 6 ; A34(1,2) = 7 ; A34(1,3) = 8 ; 
    A34(2,0) = 9 ; A34(2,1) = 10; A34(2,2) = 11; A34(2,3) = 12; 

#if 0
    F34 = A34;
    D43(0,0) = 1 ; D43(0,1) = 5 ; D43(0,2) = 9 ; 
    D43(1,0) = 2 ; D43(1,1) = 6 ; D43(1,2) = 10; 
    D43(2,0) = 3 ; D43(2,1) = 7 ; D43(2,2) = 11; 
    D43(3,0) = 4 ; D43(3,1) = 8 ; D43(3,2) = 12; 
    F34.transpose(); // transpose in place

    E43 = F34 - D43;

    std::cout << "F34=" << F34 << std::endl;
    if ( E43.data().Eucnorm() > 0 ) {
      std::cout << "main: -------------------------------------A Transpose-in-place FAILED" << std::endl;
      errcode = 8;
    } else {
      std::cout << "main: -------------------------------------A Transpose-in-place OK" << std::endl;
    }
#endif

#if 1
    // Inverse solution (not used now; not accurate enough):
//  D3(0,0) = 0.8824 ; D3(0,1) =-0.1176 ; D3(0,2) = 0.1961 ; 
//  D3(1,0) = 0.1765 ; D3(1,1) = 0.1765 ; D3(1,2) = 0.0392 ; 
//  D3(2,0) = 0.0588 ; D3(2,1) = 0.0588 ; D3(2,2) =-0.0980; 
    A3.inverse(B3);
    // check:
    D3 = A3 * B3;
    std::cout << "A3 * A3^-1=" << D3 << std::endl;
    C3.createIdentity();
    E3 = C3 - D3;
    if ( E3.data().Eucnorm() > 1.0e-7 ) {
      std::cout << "main: -------------------------------------A^-1 FAILED" << std::endl;
      errcode = 9;
    } else {
      std::cout << "main: -------------------------------------A^-1 OK" << std::endl;
    }
#endif

#if defined(_G_USE_GPTL)
    GPTLstop("blas_stuff");

    GPTLpr_file("timing.txt");
    GPTLfinalize();
#endif


term:
    if ( errcode != 0 ) {
      std::cout << serr << " Error: code=" << errcode <<  std::endl;
    }
    else {
      std::cout << serr << "     Success!" << std::endl;
    }

    return(errcode);
} // end, main
