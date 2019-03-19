//==================================================================================
// Module       : gout_simple_observer.cpp
// Date         : 3/18/19 (DLR)
// Description  : Observer object for carrying out simple POSIX-based  
//                binary output.
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#include "gout_simple.hpp"
#include "ggio.h"

using namespace std;

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Instantiate with Traits
// ARGS   : traits: Traits sturcture
//**********************************************************************************
template<typename T>
GOutSimpleObserver<T>::GOutSimpleObserver(Traits &traits)
: Observer_base(traits),
bgrid_printed_        (FALSE),
cycle_                    (0),
cycle_last_               (0),
time_last_              (0.0)
{ 
} // end of constructor (1) method


//**********************************************************************************
//**********************************************************************************
// METHOD     : observe
// DESCRIPTION: Prints state to files specified by traits. Format is:
//                  var1.CCCCCC.TTTT.out,
//              where CCCCCC represents a cycle number, and TTTT represents
//              the mpi task doing the writing.
//              NOTE: an internal cycle counter is maintained, as this 
/                     observer, like all others,  should be called at 
//                    each time step.
//
// ARGUMENTS  : t    : time, t^n, for state, uin=u^n
//              u    : state
//               
// RETURNS    : none.
//**********************************************************************************
template<typename EquationType>
void GOutSimpleObserver<T>::observe(const Time t, const State &u)
{
  init(u);
   
  if ( (traits_.itype == OBS_CYCLE && (cycle-cycle_last_)%traits_.cycle_interval == 0)
  ||   (traits_.itype == OBS_TIME  && t-time_last_ >= traits_time_interval) ) {
    gio(*grid_, u, traits.istate, cycle_, state_names_, bprgrid_);
    cycle_last_ = cycle_;
    time_last_  = t;
  }
  cycle_++;
  
} // end of method observe


//**********************************************************************************
//**********************************************************************************
// METHOD     : init
// DESCRIPTION: 
// ARGUMENTS  : u  : state variable
// RETURNS    : none.
//**********************************************************************************
template<typename EquationType>
void GOutSimpleObserver<T>::init(State &u)
{
   if ( state_names_.size() > 0 ) return;

   char    stmp[1024];
 
   if ( state_names_.size() == 0 ) {
     if ( traits_.state_names.size() == 0 ) {
       for ( auto j=0; j<u.size(); j++ ) {
         sprintf(stmp, '%s%d', "u", j+1);
         state_names_.push_back(stmp); 
       } 
     } 
     else {
       for ( auto j=0; j<u.size(); j++ ) {
         state_names_.push_back(traits_.state_names[j].data()); 
       } 
     }
   }

   if ( state_index_.size() == 0 ) {
     if ( traits_.state_index.size() == 0 ) {
       for ( auto j=0; j<u.size(); j++ ) {
         state_names_.push_back(0); 
       } 
     } 
     else {
       for ( auto j=0; j<u.size(); j++ ) {
         state_names_.push_back(traits_.state_index[j]); 
       } 
     }
   }

} // end of method init


