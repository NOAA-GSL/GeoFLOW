//==================================================================================
// Module       : gtest_allreduce.cpp
// Date         : 12/20/19 (DLR)
// Description  : GeoFLOW MPI 'allreduce' test. 
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#include <sys/stat.h>
#include <cstdio>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include "gptl.h"
#include <mpi.h>

using namespace std;

int main(int argc, char **argv)
{
    bool  bret;
    int   nrpt=1;
    int   nnodes=0;
    int   impierr, iopt, irank, ntasks;
    int   bscale=0, i, iret, j, m, nv=24;
    size_t n=100000;
    size_t lsz[2];
    size_t gsz[2];
    double         ngb=1.0, ti, tmin, tmax, tsum;

    tmax=numeric_limits<double>::min(), 
    tmin=numeric_limits<double>::max();


    // Parse command line. ':' after char
    // option indicates that it takes an argument:
    while ((iopt = getopt(argc, argv, "i:n:d:sh")) != -1) {
      switch (iopt) {
      case 'i': // get iteration count
          nrpt = atoi(optarg);
          break;
      case 'n': // get no. GB per task
          ngb = (double)atof(optarg);
          break;
      case 's': // scale input data size by no. MPP tasks
          bscale = 1;
          break;
      case 'h': // help
          cout << "usage: " << endl <<
          argv[0] << " [-h] [-s] [-d #Nodes] [-n #Size] [-i #Repeat] " << endl;
          exit(1); 
          break;
      case ':': // missing option argument
          cout << argv[0] << ": option " << optopt << " requires an argument" << endl;
          exit(1); 
          break;
      case '?':
      default: // invalid option
          cout << argv[0] << ": option " << optopt << " invalid" << endl;
          exit(1);
          break;
      }
    }

    impierr = MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &irank);
    MPI_Comm_size(comm, &ntasks);


    if ( bscale == 1 )
    {
      ngb = ((double)n ) / ((double)ntasks);
    }

    if ( ngb <= 0 ) 
    {
      cout << irank << ": main: illegal data size : ngb=" << n << endl;
      cout << irank << ": main: number of tasks is: " << ntasks << endl;
      exit(1);
    }

    cout << irank << ": main: Cycling over ntasks= " << ntasks << endl;

    // Set GTPL options:
//  GPTLsetoption (GPTLcpu, 1);
//  GPTLsetoption (GPTLwall, 0);
//  GPTLsetoption (GPTLnanotime,1);
    
    // Initialize GPTL:
    GPTLinitialize();

    // Convert # GB per task to vec size, n:
    n = ngb/nv/((double)sizeof(double)) * 1024e6;

    cout << "n = " << n << endl;

    // Initialize vectors:
    vector<double> gsend[nv];

    // 

    for ( auto j=0; j<nv; j++ ) gsend[j].resize(n);

    lsz[0] = n;
    lsz[1] = n;


    tsum = 0.0;
    for ( auto j=0; j<nrpt; j++ ) {
      GPTLstart("ti");
      MPI_Allreduce(lsz, gsz, 2, MPI_UNSIGNED_LONG_LONG, MPI_MAX, comm);
      GPTLstop("ti");
      GPTLget_wallclock("ti",-1,&ti);
      tmin   = min<double>(tmin,ti);
      tmax   = max<double>(tmax,ti);
      tsum  += ti;
    }
    

//  GPTLget_count    ("total",-1,&count[0]);
//  GPTLget_wallclock("total",-1,&tmean]);

    ofstream os;
    bool bnew;
    string fn = "times.dat";

    if ( irank == 0 ) {
      bnew = access(fn.c_str(), F_OK) < 0;

      cout << "main: bnew= " << bnew << endl;

      os.open(fn.c_str(),ios::out|ios::app);
      if ( bnew )
      {
        os << "#Tasks" << "  " << "Size (GB)"  << "  " << "#Iter" << "  " << "Tmean" << "  " << "Tmin" << "  " << "Tmax"  << endl;
      }
    
      os << ntasks << "  " << ngb << "  " << nrpt << "  " << tsum/nrpt << "  " << tmin << "  " << tmax << endl;
      os.close();
    }

//  GPTLpr_file("gptl.txt");
    GPTLfinalize();
    MPI_Finalize();

    cout << irank << ": main: done." << endl;

    return(0);
} // end, main
