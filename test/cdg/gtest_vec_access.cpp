//==================================================================================
// Module       : gtest_vec_access.cpp
// Date         : 1/1/18 (DLR)
// Description  : GeoFLOW test of vector access times, including STL
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

#include <cstdio>
#include <unistd.h>
#include <iostream>
#include <vector>
#include "gptl.h"
#include "gtvector.hpp"

int main(int argc, char **argv)
{
    GINT   nrpt=1000;
    GINT   iopt;
    GSIZET n=100000;
    GSIZET i,j;
    std::vector<GDOUBLE>  vec1,  vec2,  vec3;
    GTVector<GDOUBLE>    gvec1, gvec2, gvec3;

    // Parse command line. ':' after char
    // option indicates that it takes an argument:
    while ((iopt = getopt(argc, argv, "i:n:h")) != -1) {
      switch (iopt) {
      case 'i': // get iteration count
          nrpt = atoi(optarg);
          break;
      case 'n': // get array/vector size 
          n = atoi(optarg);
          break;
      case 'h': // help
          std::cout << "usage: " << std::endl <<
          argv[0] << " [-h] [-n #Size] [-i #Repeat] " << std::endl;
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


    // Set GTPL options:
    GPTLsetoption (GPTLcpu, 1);
//  GPTLsetoption (GPTLwall, 0);
//  GPTLsetoption (GPTLnanotime,1);
    
    // Initialize GPTL:
    GPTLinitialize();

    // Initialize vectors:
    gvec1.resize(n);
    gvec2.resize(n);
    gvec3.resize(n);
    for ( i=0; i<n; i++ )
    {
      vec1.push_back((GDOUBLE)(2*i));
      vec2.push_back((GDOUBLE)( i));
      gvec1[i] = vec1[i];
      gvec2[i] = vec2[i];
    }
    
    vec3 .resize(n);

    //////////////////////////////////////////////////////////////
    // Check vectorization/optimization using []-access methods:
    //////////////////////////////////////////////////////////////
    #pragma acc data copy(vec3[0:n-1]),copyin(vec1[0:n-1],vec2[0:n-1])
    GPTLstart("std:[]access_mthd");
    for ( j=0; j<nrpt; j++ ) 
    {
      #pragma acc kernels
      for ( i=0; i<n; i++ )
      {
        vec3[i] = vec1[i] + vec2[i];
      }
    }
    GPTLstop("std:[]access_mthd");

    std::cout << "After check after using std access methods:" << std::endl;
    std::vector<GDOUBLE>::iterator vi3=vec3.begin();
    std::vector<GDOUBLE>::iterator vi1=vec1.begin();
    std::vector<GDOUBLE>::iterator vi2=vec2.begin();
    while ( vi3 < vec3.begin()+10 )
    {
      std::cout << *vi3 << " ";
      vi3++; 
    }  
    std::cout << std::endl;

    vec3.assign(n,0);

    #pragma acc data copy(gvec3[0:n-1]),copyin(gvec1[0:n-1],gvec2[0:n-1])
    GPTLstart("gtvec:[]access_mthd");
    for ( j=0; j<nrpt; j++ ) 
    {
      #pragma acc kernels
//    #pragma forceinline recursive
      for ( i=0; i<n; i++ )
      {
        gvec3[i] = gvec1[i] + gvec2[i];
      }
    }
    GPTLstop("gtvec:[]access_mthd");

    std::cout << "After check after using gtvector[] data access:" << std::endl;


    //////////////////////////////////////////////////////////////
    // Check vectorization/optimization using direct data access:
    //////////////////////////////////////////////////////////////
    GDOUBLE *d1 = vec1.data();
    GDOUBLE *d2 = vec2.data();
    GDOUBLE *d3 = vec3.data();

    #pragma acc data copy(vec3[0:n-1]),copyin(vec1[0:n-1],vec2[0:n-1])
    GPTLstart("std:direct_access");
    for ( j=0; j<nrpt; j++ ) 
    {
      #pragma acc kernels
      for ( i=0; i<n; i++ )
      {
        d3[i] = d1[i] + d2[i];
      }
    }
    GPTLstop("std:direct_access");

    std::cout << "After check after using direct data access:" << std::endl;
    vi3=vec3.begin();
    vi1=vec1.begin();
    vi2=vec2.begin();
    while ( vi3 < vec3.begin()+10 )
    {
      std::cout << *vi3 << " ";
      vi3++; 
    }  
    std::cout << std::endl;

    vec3.assign(n,0);

    d1 = gvec1.data();
    d2 = gvec2.data();
    d3 = gvec3.data();

    #pragma acc data copy(gvec3[0:n-1]),copyin(gvec1[0:n-1],gvec2[0:n-1])
    GPTLstart("gtvec:direct_access");
    for ( j=0; j<nrpt; j++ )
    {
      #pragma acc kernels
      for ( i=0; i<n; i++ )
      {
        d3[i] = d1[i] + d2[i];
      }
    }
    GPTLstop("gtvec:direct_access");


    //////////////////////////////////////////////////////////////
    // Check vectorization/optimization using iterators.
    // Admittedly, this is somewhat artifical, but worth a check:
    //////////////////////////////////////////////////////////////

    #pragma acc data copy(vec3[0:n-1]),copyin(vec1[0:n-1],vec2[0:n-1])

    GPTLstart("std:iterators");
    for ( j=0; j<nrpt; j++ ) 
    {
      vi1=vec1.begin();
      vi2=vec2.begin();
      vi3=vec3.begin();
//    NOTE: must delcare vi3 private below to prevent 
      #pragma acc kernels loop private(vi3)
      for ( i=0 ; i<n; i++,vi1++,vi2++,vi3++ )
      {
        *vi3 = (*vi1) + (*vi2);
      }
    }
    GPTLstop("std:iterators");

    std::cout << "After check with iterators:" << std::endl;
    vi3=vec3.begin();
    while ( vi3 < vec3.begin()+10 )
    {
      std::cout << *vi3 << " ";
      vi3++; 
    }  
    std::cout << std::cout;

    GPTLpr_file("gptl.txt");
    GPTLfinalize();

    return(0);
} // end, main
