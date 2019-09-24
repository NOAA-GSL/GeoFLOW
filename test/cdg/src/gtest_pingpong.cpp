//==================================================================================
// Module       : gtest_pingpong.cpp
// Date         : 3/1/18 (DLR)
// Description  : GeoFLOW MPI 'ping-pong' test. Really, it's an all-to-all
//                test with fixed or scalable (to the number of MPI tasks)
//                data.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#include <sys/stat.h>
#include <cstdio>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include "gptl.h"
#include <mpi.h>
#include "gcomm.hpp"
#include "gtvector.hpp"

using namespace std;

int main(int argc, char **argv)
{
    GBOOL  bret;
    GINT   nrpt=1000;
    GINT   nnodes=0;
    GINT   iopt, irank, ntasks;
    GINT   bscale=0, i, iret, j, m;
    GSIZET n=100000;

    // Parse command line. ':' after char
    // option indicates that it takes an argument:
    while ((iopt = getopt(argc, argv, "i:n:d:sh")) != -1) {
      switch (iopt) {
      case 'i': // get iteration count
          nrpt = atoi(optarg);
          break;
      case 'n': // get array/vector size 
          n = atoi(optarg);
          break;
      case 'd': // get array/vector size 
          nnodes = atoi(optarg);
          break;
      case 's': // scale input data size by no. MPP tasks
          bscale = 1;
          break;
      case 'h': // help
          std::cout << "usage: " << endl <<
          argv[0] << " [-h] [-s] [-d #Nodes] [-n #Size] [-i #Repeat] " << endl;
          exit(1); 
          break;
      case ':': // missing option argument
          std::cout << argv[0] << ": option " << optopt << " requires an argument" << endl;
          exit(1); 
          break;
      case '?':
      default: // invalid option
          std::cout << argv[0] << ": option " << optopt << " invalid" << endl;
          exit(1);
          break;
      }
    }

    GComm::InitComm(&argc, &argv);
    irank = GComm::WorldRank(); 
    ntasks = GComm::WorldSize();

    if ( bscale == 1 )
    {
      n = ((GDOUBLE)n ) / ((GDOUBLE)ntasks);
    }

    if ( n <= 0 ) 
    {
      std::cout << irank << ": main: illegal data size : n=" << n << endl;
      std::cout << irank << ": main: number of tasks is: " << ntasks << endl;
      exit(1);
    }

    std::cout << irank << ": main: Cycling over ntasks= " << ntasks << endl;

    // Set GTPL options:
//  GPTLsetoption (GPTLcpu, 1);
//  GPTLsetoption (GPTLwall, 0);
//  GPTLsetoption (GPTLnanotime,1);
    
    // Initialize GPTL:
    GPTLinitialize();

    // Initialize vectors:
    GTVector<GDOUBLE> gsend, rbuff[ntasks-1];

    gsend.resize(n);
    for ( i=0; i<n; i++ )
    {
      gsend[i] = irank*n+i;
    }
    for ( i=0; i<ntasks-1; i++ )
    {
      rbuff[i].resize(n);
    }


    MPI_Request *rhandle = new MPI_Request [ntasks-1];
    MPI_Request *shandle = new MPI_Request [ntasks-1];
    
    GPTLstart("total");
    
    for ( j=0; j<nrpt; j++ ) 
    {
      // Post receives:
      GPTLstart("post_recvs");
      for ( i=0, m=0 ; i<ntasks; i++ )
      {
        if ( i == irank ) continue;  // don't recv from self
        bret = GComm::IRecv(rbuff[m].data(), n, GC_GDOUBLE, i, &rhandle[m]);
        if ( !bret )
        {
          std::cout << irank<< ": main: IRecv error: post from task: " << i << endl;
          exit(1);
        }
        m++;
      }
      GPTLstop("post_recvs");

      // Do sends:
      GPTLstart("sends");
      for ( i=0, m=0 ; i<ntasks; i++ )
      {
        if ( i == irank ) continue; // don't send to self
        bret = GComm::ISend(gsend.data(), n, GC_GDOUBLE, i, &shandle[m]);
        if ( !bret )
        {
          std::cout << irank << ": main: BSend error: to task: " << i << endl;
          exit(1);
        }
        m++;
      }
      GPTLstop("sends");

      // Wait:
      GPTLstart("waitall");
      bret = GComm::BWaitAll(shandle, ntasks-1);
      if ( !bret )
      {
        std::cout << irank << ": main: BWaitAll send error" << endl;
        exit(1);
      }

      bret = GComm::BWaitAll(rhandle, ntasks-1);
      if ( !bret )
      {
        std::cout << irank << ": main: BWaitAll recv error" << endl;
        exit(1);
      }
      GPTLstop("waitall");
    }
    GPTLstop("total");

    std::vector   <GINT> count(4);
    std::vector<GDOUBLE> time (4);
    GPTLget_count    ("total",-1,&count[0]);
    GPTLget_wallclock("total",-1,&time [0]);
    GPTLget_count    ("post_recvs",-1,&count[1]);
    GPTLget_wallclock("post_recvs",-1,&time [1]);
    GPTLget_count    ("sends",-1,&count[2]);
    GPTLget_wallclock("sends",-1,&time [2]);
    GPTLget_count    ("waitall",-1,&count[3]);
    GPTLget_wallclock("waitall",-1,&time [3]);

    ofstream os;
    GBOOL bnew;
    GString fn = "times.dat";

    if ( irank == 0 ) {
      bnew = access(fn.c_str(), F_OK) < 0;

      cout << "main: bnew= " << bnew << endl;

      os.open(fn.c_str(),ios::out|ios::app);
      if ( bnew )
      {
        os << "#Nodes" << "  " << "#Tasks" << "  " << "Size"  << "  " << "#Iter" << "  " << "TTotal" << "  " << "TRecv"  << "  " << "TSend"  << "  " << "TWait" << endl;
      }
    
      os << nnodes << "  " << ntasks << "  " << n << "  " << nrpt << "  " << time[0] << "  " << time[1] << "  " << time[2] << "  " << time[3] << endl;
      os.close();
    }

//  GPTLpr_file("gptl.txt");
    GPTLfinalize();
    GComm::TermComm();

    std::cout << irank << ": main: done." << endl;

    return(0);
} // end, main
