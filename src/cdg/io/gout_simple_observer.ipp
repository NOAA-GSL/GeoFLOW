//==================================================================================
// Module       : gout_simple_observer.ipp
// Date         : 3/18/19 (DLR)
// Description  : Observer object for carrying out simple POSIX-based  
//                binary output.
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Instantiate with Traits
// ARGS   : traits: Traits sturcture
//**********************************************************************************
template<typename EquationType>
GOutSimpleObserver<EquationType>::GOutSimpleObserver(typename ObserverBase<EquationType>::Traits &traits, Grid &grid):
bprgrid_        (TRUE),
cycle_          (0),
ocycle_         (1),
cycle_last_     (0),
time_last_      (0.0)
{ 
  this->traits_ = traits;
  this->grid_   = &grid;
} // end of constructor (1) method


//**********************************************************************************
//**********************************************************************************
// METHOD     : observe_impl
// DESCRIPTION: Prints state to files specified by traits. Format is:
//                  var1.CCCCCC.TTTT.out,
//              where CCCCCC represents a cycle number, and TTTT represents
//              the mpi task doing the writing.
//              NOTE: an internal cycle counter is maintained, as this 
//                    observer, like all others,  should be called at 
//                    each time step.
//
// ARGUMENTS  : t    : time, t^n, for state, uin=u^n
//              u    : state
//               
// RETURNS    : none.
//**********************************************************************************
template<typename EquationType>
void GOutSimpleObserver<EquationType>::observe_impl(const Time &t, const State &u)
{
  init(t,u);

  mpixx::communicator comm;
   
  if ( (this->traits_.itype == ObserverBase<EquationType>::OBS_CYCLE 
        && (cycle_-cycle_last_) == this->traits_.cycle_interval)
    || (this->traits_.itype == ObserverBase<EquationType>::OBS_TIME  
        &&  t-time_last_ >= this->traits_.time_interval) ) {
    gio(*(this->grid_), u, state_index_, ocycle_, t, state_names_, sdir_, comm, bprgrid_);
    bprgrid_ = FALSE;
    cycle_last_ = cycle_;
    time_last_  = t;
    ocycle_++;
  }
  cycle_++;
  
} // end of method observe_impl


//**********************************************************************************
//**********************************************************************************
// METHOD     : init
// DESCRIPTION: Fill member index and name data based on traits
// ARGUMENTS  : t  : state time
//              u  : state variable
// RETURNS    : none.
//**********************************************************************************
template<typename EquationType>
void GOutSimpleObserver<EquationType>::init(const Time t, const State &u)
{

   char    stmp[1024];

   sdir_ = this->traits_.dir;
 
   if ( icycle_ == 0 ) {
     time_last_ = t; 
   }
 
   // Set state names member data, if not already set:
   if ( state_names_.size()  <= 0 ) {
     if ( this->traits_.state_names.size() == 0 ) {
       for ( auto j=0; j<u.size(); j++ ) {
         sprintf(stmp, "%s%d", "u", j+1);
         state_names_.push_back(stmp); 
       } 
     } 
     else {
       for ( auto j=0; j<u.size(); j++ ) {
         state_names_.push_back(this->traits_.state_names[j].data()); 
       } 
     }
   }

   // Set state index member data, if not already set:
   if ( state_index_.size()  <= 0 ) {
     if ( this->traits_.state_index.size() == 0 ) {
       for ( auto j=0; j<state_names_.size(); j++ ) {
         state_index_.push_back(j); 
       } 
     } 
     else {
       for ( auto j=0; j<this->traits_.state_index.size(); j++ ) {
         state_index_.push_back(this->traits_.state_index[j]); 
       } 
     }
  }

} // end of method init

;
