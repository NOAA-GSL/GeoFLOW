//==================================================================================
// Module       : gtest_gen_type.cpp
// Date         : 2/8/18 (BF)
// Description  : GeoFLOW test of contiguity of data in 'elements'
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#include <cstdint>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <random>
#include <type_traits>
#include <vector>
#include <list>
#include <iterator>
#include "gptl.h"
#if defined(USE_PAPI)
#include <papi.h>
#endif
#include "gtypes.h"
#include "gtvector.hpp"
#include "cff_blas.h"

#define  RAT(a,b) ( ((GDOUBLE)a) / ((GDOUBLE)b) )

#define DO_PAPI
#define DO_MATMAT         // Do Mat-Mat mult as 'work'

using namespace std;


#define PORDER          MIN(8,G_MAX_ORDER)
                          // Vector Length within Element (not changeable)
#define CACHE_SIZE      8 // Cache-blocking factor default
#define NWORD          12 // Bytes extra in each Element (offset)


/**
 * Element with fixed size memory
 */
class ElementA {

	static const GSIZET size_ = PORDER*PORDER;

public:

	ElementA(GSIZET icsz=CACHE_SIZE):
        icsz_ (icsz) {
		for(GSIZET i = 0; i < size_; ++i){
			p_[i] = i;
			q_[i] = 2*i;
			r_[i] = 3*i;
			s_[i] = 0.0;
		}
	}

	ElementA(ElementA const& EA){
               icsz_ = EA.icsz_;
		for(GSIZET i = 0; i < size_; ++i){
			p_[i] = EA.p_[i];
			q_[i] = EA.q_[i];
			r_[i] = EA.r_[i];
			s_[i] = EA.s_[i];
		}
	}

	~ElementA(){
	}

	GSIZET get_size() const{
		return size_;
	}

	void do_work(){
                #if !defined(DO_MATMAT)
		for(GSIZET i = 0; i < size_; ++i){
			s_[i] = p_[i] * q_[i] * r_[i];
		}
                #else
                GSIZET nr=PORDER;
                GSIZET nc=PORDER;
                dmxm(s_, p_, &nr, &nc, r_, &nr, &nc, &icsz_); 
                #endif
	}

private:
	GCHAR   dummy[NWORD]; // other vars
        GSIZET   icsz_;
	GDOUBLE p_[size_];
	GDOUBLE q_[size_];
	GDOUBLE r_[size_];
	GDOUBLE s_[size_];
};


/**
 * Element with dynamic allocated memory
 */
class ElementB {

//     static const GSIZET size_ = PORDER*PORDER;

public:

	ElementB(GSIZET psz, GSIZET icsz=CACHE_SIZE):
        psz_  (psz),
        size_ (psz*psz),
        icsz_ (icsz) {
		p_ = new GDOUBLE[size_];
		q_ = new GDOUBLE[size_];
		r_ = new GDOUBLE[size_];
		s_ = new GDOUBLE[size_];
		for(GSIZET i = 0; i < size_; ++i){
			p_[i] = i;
			q_[i] = 2*i;
			r_[i] = 3*i;
			s_[i] = 0.0;
		}
	}

	ElementB(ElementB const& EB){
                psz_  = EB.psz_;
                size_ = EB.size_;

		p_ = new GDOUBLE[size_];
		q_ = new GDOUBLE[size_];
		r_ = new GDOUBLE[size_];
		s_ = new GDOUBLE[size_];
                icsz_ = EB.icsz_;
		for(GSIZET i = 0; i < size_; ++i){
			p_[i] = EB.p_[i];
			q_[i] = EB.r_[i];
			r_[i] = EB.q_[i];
			s_[i] = EB.s_[i];
		}
	}

	~ElementB(){
		delete[] p_;
		delete[] q_;
		delete[] r_;
		delete[] s_;
	}

	GSIZET get_size() const{
			return size_;
	}

	void do_work(){
                #if !defined(DO_MATMAT)
		for(GSIZET i = 0; i < size_; ++i){
			s_[i] = p_[i] * q_[i] * r_[i];
		}
                #else
                GSIZET nr=psz_;
                GSIZET nc=psz_;
                dmxm(s_, p_, &nr, &nc, r_, &nr, &nc, &icsz_);
                #endif

	}

private:
	GCHAR    dummy[NWORD];
        GSIZET   psz_;
        GSIZET   size_;
        GSIZET   icsz_;
	GDOUBLE* p_;
	GDOUBLE* q_;
	GDOUBLE* r_;
	GDOUBLE* s_;
};

/**
 * Element with pointers to memory
 */
class ElementC {

	static const GSIZET size_ = PORDER*PORDER;

public:

	ElementC(GSIZET icsz=CACHE_SIZE) :
		p_(NULLPTR),
		q_(NULLPTR),
		r_(NULLPTR),
		s_(NULLPTR),
                icsz_ (icsz) {
	}

	ElementC(GDOUBLE* p, GDOUBLE* q, GDOUBLE*r, GDOUBLE* s, GSIZET icsz=CACHE_SIZE) :
		p_(p), q_(q), r_(r), s_(s), icsz_(icsz){
		for(GSIZET i = 0; i < size_; ++i){
			p_[i] = i;
			q_[i] = 2*i;
			r_[i] = 3*i;
			s_[i] = 0.0;
		}
	}

	ElementC(ElementC const& EC) :
		p_(EC.p_), q_(EC.q_), r_(EC.r_), s_(EC.s_), icsz_(EC.icsz_) {
	}

	~ElementC(){
		p_ = NULLPTR;
		q_ = NULLPTR;
		r_ = NULLPTR;
		s_ = NULLPTR;
	}

	GSIZET get_size() const{
		return size_;
	}

	void do_work(){
                #if !defined(DO_MATMAT)
		for(GSIZET i = 0; i < size_; ++i){
			s_[i] = p_[i] * q_[i] * r_[i];
		}
                #else
                GSIZET nr=PORDER;
                GSIZET nc=PORDER;
                dmxm(s_, p_, &nr, &nc, r_, &nr, &nc, &icsz_);
                #endif

	}

private:
	GCHAR    dummy[NWORD]; // other vars
        GSIZET   icsz_;
	GDOUBLE* p_;
	GDOUBLE* q_;
	GDOUBLE* r_;
	GDOUBLE* s_;
};

GINT main(int argc, char *argv[]){
        GString swork;
        GINT    iopt;
        GSIZET  nelem  = 1000000; //MIN(100000,G_MAX_ELEMS); 
        GSIZET  nrpt   = 100;
        GSIZET  icsz   = CACHE_SIZE;
        GSIZET  psize  = PORDER; // 'exp' order to dyn allocation
 
        #if defined(DO_MATMAT)
          swork = "_matmat";
        #else
          swork = "_3prod";
        #endif

        // Parse command line. ':' after char
        // option indicates that it takes an argument:
        while ((iopt = getopt(argc, argv, "i:n:c:p:h")) != -1) {
          switch (iopt) {
          case 'i': // get nrpt/iteration count
              nrpt = atoi(optarg);
              break;
          case 'c': // get cache-blocking factor
              icsz = atoi(optarg);
              break;
          case 'n': // num elems
              nelem = atoi(optarg);
              break;
          case 'p': // 'expansion order'
              psize = atoi(optarg);
              break;
          case 'h': // help
              cout << "usage: " << std::endl <<
              argv[0] << " [-h] [-n #Elems] [-p #PolyOrder] [-i #Repeat] " << std::endl;
              exit(1);
              break;
          case ':': // missing option argument
              cout << argv[0] << ": option " << optopt << " requires an argument" << std::endl;
              exit(1);
              break;
          case '?':
          default: // invalid option
              cout << argv[0] << ": option " << optopt << " invalid" << std::endl;
              exit(1);
              break;
          }
        }

        // Initialize GPTL:
        GINT iret;

        #if defined(USE_PAPI) && defined(DO_PAPI)
        iret=GPTLsetoption (PAPI_L1_DCM  , 1);
        iret=GPTLsetoption (PAPI_L2_DCM  , 1);
        iret=GPTLsetoption (PAPI_L2_DCA  , 1);
        iret=GPTLsetoption (PAPI_TOT_CYC , 1);
        iret=GPTLsetoption (PAPI_TOT_INS , 1);
        iret=GPTLsetoption (PAPI_LST_INS , 1);
        iret=GPTLsetoption (PAPI_REF_CYC , 1);

//      iret=GPTLsetoption (PAPI_L1_LDM  , 1);
//      iret=GPTLsetoption (PAPI_L2_LDM  , 1);
//      iret=GPTLsetoption (PAPI_L3_DCM  , 1);
//      iret=GPTLsetoption (PAPI_L3_DCA  , 1);
//      iret=GPTLsetoption (PAPI_PRF_DM  , 1);
  
        #endif

        if ( (iret=GPTLinitialize ())  != 0) {
          std::cout << "main: GPTLinitialize failure" << std::endl;
          exit(1);
        }

	// Fixed Elements
	std::vector<ElementA> VectorA;
	for(GSIZET i = 0; i < nelem; ++i){
		VectorA.push_back( ElementA(icsz) );
	}


	// Dynamic Elements
	std::vector<ElementB> VectorB;
	for(GSIZET i = 0; i < nelem; ++i){
		VectorB.push_back( ElementB(psize,icsz) );
	}

	// Pointer Elements (Do not time !!!)
	std::vector<ElementC> VectorC;
        std::vector<ElementC> VectorD;
	std::vector<GDOUBLE> pqrs (4*PORDER*PORDER*nelem);
        std::vector<GDOUBLE> pqrs2(4*PORDER*PORDER*nelem);
        std::list<ElementC>  Elist;

	// Memory Layout
	// [p1, p2, ..., pN, q1, q2, ..., qN, ... s1, s2, ... sN]
	for(GSIZET i = 0; i < nelem; ++i){
		VectorC.push_back( ElementC(
				pqrs.data()+(PORDER*PORDER*nelem*0 + i*PORDER*PORDER),
				pqrs.data()+(PORDER*PORDER*nelem*1 + i*PORDER*PORDER),
				pqrs.data()+(PORDER*PORDER*nelem*2 + i*PORDER*PORDER),
				pqrs.data()+(PORDER*PORDER*nelem*3 + i*PORDER*PORDER,icsz)) );
	}


	// Memory Layout
	// [p1, q1, r1, s1, p2, q2, r2, s2, ... pN, qN, rN, sN]
	for(GSIZET i = 0; i < nelem; ++i){
		VectorD.push_back( ElementC(
                                pqrs2.data()+(4*i+0)*PORDER*PORDER,
                                pqrs2.data()+(4*i+1)*PORDER*PORDER,
                                pqrs2.data()+(4*i+2)*PORDER*PORDER,
                                pqrs2.data()+(4*i+3)*PORDER*PORDER,icsz));
	}

	// Memory Layout
	// [p1, p2, ..., pN, q1, q2, ..., qN, ... s1, s2, ... sN]
	for(GSIZET i = 0; i < nelem; ++i){
		Elist.push_back( ElementC(
				pqrs.data()+(PORDER*PORDER*nelem*0 + i*PORDER*PORDER),
				pqrs.data()+(PORDER*PORDER*nelem*1 + i*PORDER*PORDER),
				pqrs.data()+(PORDER*PORDER*nelem*2 + i*PORDER*PORDER),
				pqrs.data()+(PORDER*PORDER*nelem*3 + i*PORDER*PORDER),icsz) );
	}
	// -------------------------------------------------------------
        cout << "--- Testing ---" << std::endl;



	std::mt19937 rng(12345);
	std::uniform_int_distribution<GINT> gen(0,4);


        GINT  last_path = -1;
	for(GSIZET s = 0; s < nrpt; ++s){
		GINT random_path = gen(rng);
                if( last_path == random_path ){
                  random_path -= 1;
                  if( random_path == -1 ){
                    random_path = 5;
                  }
                }
                last_path = random_path;
                  
                cout << "main: starting iteration " << s << " random_path " << random_path << "..." << std::endl;
		if( random_path == 0 ){

                        GPTLstart("ElemA");
			for(GSIZET i = 0; i < nelem; ++i){
				VectorA[i].do_work();
			}

                        GPTLstop("ElemA");

		}
		else if( random_path == 1 ){

                        GPTLstart("ElemB");
			for(GSIZET i = 0; i < nelem; ++i){
				VectorB[i].do_work();
			}
                        GPTLstop("ElemB");

		}
		else if( random_path == 2 ) {

                        GPTLstart("ElemC");
			for(GSIZET i = 0; i < nelem; ++i){
				VectorC[i].do_work();
			}
                        GPTLstop("ElemC");

		}
                else if( random_path == 3 ) {
                        GPTLstart("ElemD");
			for(GSIZET i = 0; i < nelem; ++i){
				VectorD[i].do_work();
			}
                        GPTLstop("ElemD");
                }
#if !defined(_G_USE_OPENACC)
                else if( random_path == 4 ) {
                        GPTLstart("ElemE");
                        std::list<ElementC>::iterator ie=Elist.begin();
	                for(; ie!=Elist.end(); ++ie){ 
				ie->do_work();
			}
                        GPTLstop("ElemE");
                }
#endif

                cout << "main: ...iteration " << s << " random_path " << random_path << " done." << std::endl;
	}

        // Print timer data to appended file for easier access:
        std::vector   <GINT> count(6);
        std::vector<GDOUBLE> time (6);
        GPTLget_count    ("ElemA",-1,&count[0]);
        GPTLget_wallclock("ElemA",-1,&time [0]);
        GPTLget_count    ("ElemB",-1,&count[1]);
        GPTLget_wallclock("ElemB",-1,&time [1]);
        GPTLget_count    ("ElemC",-1,&count[2]);
        GPTLget_wallclock("ElemC",-1,&time [2]);
        GPTLget_count    ("ElemD",-1,&count[3]);
        GPTLget_wallclock("ElemD",-1,&time [3]);
        #if !defined(_G_USE_OPENACC)
        GPTLget_count    ("ElemE",-1,&count[4]);
        GPTLget_wallclock("ElemE",-1,&time [4]);
        #endif

        // Get PAPI data if available:
        GTVector<GDOUBLE> l1m(6), l2m(6), l2a(6), l3m(6), 
                          ins(6), lst(6), cyc(6), prf(6), 
                          lm1(6), lm2(6), ref(6);
        l1m = 0.0; l2m = 0.0; l2a = 0.0; l3m = 0.0; ins = 0.0; 
        lst = 0.0; cyc = 0.0; prf = 0.0; lm1 = 0.0; lm2 = 0.0; 
        ref = 0.0;
#if defined(USE_PAPI) && defined(DO_PAPI)
std::cout << "main: USE_PAPI..." << std::endl;
        GPTLget_eventvalue("ElemA"  ,"PAPI_L1_DCM" ,-1,&l1m[0]);
        GPTLget_eventvalue("ElemA"  ,"PAPI_L2_DCM" ,-1,&l2m[0]);
        GPTLget_eventvalue("ElemA"  ,"PAPI_L2_DCA" ,-1,&l2a[0]);
        GPTLget_eventvalue("ElemA"  ,"PAPI_LST_INS",-1,&lst[0]);
        GPTLget_eventvalue("ElemA"  ,"PAPI_TOT_CYC",-1,&cyc[0]);
        GPTLget_eventvalue("ElemA"  ,"PAPI_TOT_INS",-1,&ins[0]);
        GPTLget_eventvalue("ElemA"  ,"PAPI_REF_CYC",-1,&ref[0]);

        GPTLget_eventvalue("ElemB"  ,"PAPI_L1_DCM" ,-1,&l1m[1]);
        GPTLget_eventvalue("ElemB"  ,"PAPI_L2_DCM" ,-1,&l2m[1]);
        GPTLget_eventvalue("ElemB"  ,"PAPI_L2_DCA" ,-1,&l2a[1]);
        GPTLget_eventvalue("ElemB"  ,"PAPI_LST_INS",-1,&lst[1]);
        GPTLget_eventvalue("ElemB"  ,"PAPI_TOT_CYC",-1,&cyc[1]);
        GPTLget_eventvalue("ElemB"  ,"PAPI_TOT_INS",-1,&ins[1]);
        GPTLget_eventvalue("ElemB"  ,"PAPI_REF_CYC",-1,&ref[1]);

        GPTLget_eventvalue("ElemC"  ,"PAPI_L1_DCM" ,-1,&l1m[2]);
        GPTLget_eventvalue("ElemC"  ,"PAPI_L2_DCM" ,-1,&l2m[2]);
        GPTLget_eventvalue("ElemC"  ,"PAPI_L2_DCA" ,-1,&l2a[2]);
        GPTLget_eventvalue("ElemC"  ,"PAPI_LST_INS",-1,&lst[2]);
        GPTLget_eventvalue("ElemC"  ,"PAPI_TOT_CYC",-1,&cyc[2]);
        GPTLget_eventvalue("ElemC"  ,"PAPI_TOT_INS",-1,&ins[2]);
        GPTLget_eventvalue("ElemC"  ,"PAPI_REF_CYC",-1,&ref[2]);

        GPTLget_eventvalue("ElemD"  ,"PAPI_L1_DCM" ,-1,&l1m[3]);
        GPTLget_eventvalue("ElemD"  ,"PAPI_L2_DCM" ,-1,&l2m[3]);
        GPTLget_eventvalue("ElemD"  ,"PAPI_L2_DCA" ,-1,&l2a[3]);
        GPTLget_eventvalue("ElemD"  ,"PAPI_LST_INS",-1,&lst[3]);
        GPTLget_eventvalue("ElemD"  ,"PAPI_TOT_CYC",-1,&cyc[3]);
        GPTLget_eventvalue("ElemD"  ,"PAPI_TOT_INS",-1,&ins[3]);
        GPTLget_eventvalue("ElemD"  ,"PAPI_REF_CYC",-1,&ref[3]);
  #if !defined(_G_USE_OPENACC)
        GPTLget_eventvalue("ElemE"  ,"PAPI_L1_DCM" ,-1,&l1m[4]);
        GPTLget_eventvalue("ElemE"  ,"PAPI_L2_DCM" ,-1,&l2m[4]);
        GPTLget_eventvalue("ElemE"  ,"PAPI_L2_DCA" ,-1,&l2a[4]);
        GPTLget_eventvalue("ElemE"  ,"PAPI_LST_INS",-1,&lst[4]);
        GPTLget_eventvalue("ElemE"  ,"PAPI_TOT_CYC",-1,&cyc[4]);
        GPTLget_eventvalue("ElemE"  ,"PAPI_TOT_INS",-1,&ins[4]);
        GPTLget_eventvalue("ElemE"  ,"PAPI_REF_CYC",-1,&ref[4]);
  #endif
#endif

        ofstream os[6];
        std::ostringstream fn;
        
        fn << "elemA_fixed" << swork << ".txt";
        os[0].open(fn.str().c_str(),ios::out|ios::app);
        os[0] << nelem << " " << PORDER << " " << icsz << " " << nrpt << " " << time[0]/count[0] << " " << l1m[0] << " " << l2m[0] << " " << l2a[0] <<  " " << cyc[0] << " " <<  ins[0] << " " << ref[0]  << " " << lst[0] << std::endl;
        os[0].close();

        fn.str("");
        fn << "elemB_dyn" << swork << ".txt";
        os[1].open(fn.str().c_str(),ios::out|ios::app);
        os[1] << nelem << " " << PORDER << " " << icsz << " " << nrpt << " " << time[1]/count[1] << " " << l1m[1] << " " << l2m[1] << " " << l2a[1] << " " << cyc[1] << " " <<  ins[1] << " " << ref[1]  << " " << lst[1] << std::endl;
        os[1].close();

        fn.str("");
        fn << "elemC_p1p2p3" << swork << ".txt";
        os[2].open(fn.str().c_str(),ios::out|ios::app);
        os[2] << nelem << " " << PORDER << " " << icsz << " " << nrpt << " " << time[2]/count[2] << " " << l1m[2] << " " << l2m[2] << " " << l2a[2] << " " << cyc[2] << " " <<  ins[2] << " " << ref[2]  << " " << lst[2] << std::endl;
        os[2].close();

        fn.str("");
        fn << "elemD_p1q2r3" << swork << ".txt";
        os[3].open(fn.str().c_str(),ios::out|ios::app);
        os[3] << nelem << " " << PORDER << " " << icsz << " " << nrpt << " " << time[3]/count[3] << " " << l1m[3] << " " << l2m[3] << " " << l2a[3] << " " << cyc[3] << " " << ins[3] << " " << ref[3]  << " " << lst[3] << std::endl;
        os[3].close();

        #if !defined(_G_USE_OPENACC)
        fn.str("");
        fn << "elemE_p1p2p3" << swork << ".txt";
        os[4].open(fn.str().c_str(),ios::out|ios::app);
        os[4] << nelem << " " << PORDER << " " << icsz << " " << nrpt << " " << time[4]/count[4] << " " << l1m[4] << " " << l2m[4] << " " << l2a[4] << " " << cyc[4] << " " <<  ins[4] << " " << ref[4] << " " << lst[4] << std::endl;
        os[4].close();
        #endif

        cout << "main: test done. " << std::endl;
        GPTLpr_file("gptl.txt");

        GPTLfinalize();


};
