//==================================================================================
// Module       : gglobaldiag_basic.ipp
// Date         : 3/28/19 (DLR)
// Description  : Observer object for carrying out L2 & extrema diagnostics for
//                kinetic quantities
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : ObserverBase.
//==================================================================================

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Instantiate with Traits
// ARGS   : traits: Traits sturcture
//**********************************************************************************
template<typename EquationType>
GGlobalDiag_basic<EquationType>::GGlobalDiag_basic(const EqnBasePtr &equation, Grid &grid, typename ObserverBase<EquationType>::Traits &traits):
ObserverBase<EquationType>(equation, grid, traits),
bInit_          (FALSE),
cycle_          (0),
ocycle_         (1),
cycle_last_     (0),
time_last_      (0.0),
grid_           (&grid)
{ 
  traits_ = traits;
  utmp_   = static_cast<GTVector<GTVector<GFTYPE>*>*>(utmp_);
  myrank_ = GComm::WorldRank(grid.get_comm());
} // end of constructor (1) method


//**********************************************************************************
//**********************************************************************************
// METHOD     : observe_impl
// DESCRIPTION: Compute energy, enstrophy, helicity, and energy injection, 
//              and output to one file. Compute max of energy, enstrophy,
//              and output to another file.
//
// ARGUMENTS  : t    : time, t^n, for state, uin=u^n
//              u    : state
//              uf   : forcing
//               
// RETURNS    : none.
//**********************************************************************************
template<typename EquationType>
void GGlobalDiag_basic<EquationType>::observe_impl(const Time &t, const State &u, const State &uf)
{
  init(t,u);

  mpixx::communicator comm;

  if ( (traits_.itype == ObserverBase<EquationType>::OBS_CYCLE 
        && (cycle_-cycle_last_+1) >= traits_.cycle_interval)
    || (traits_.itype == ObserverBase<EquationType>::OBS_TIME  
        &&  t-time_last_ >= traits_.time_interval) ) {

    do_kinetic_L2 (t, u, uf, "gbalance.txt");
    do_kinetic_max(t, u, uf, "gmax.txt");
    cycle_last_ = cycle_+1;
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
void GGlobalDiag_basic<EquationType>::init(const Time t, const State &u)
{
   assert(utmp_ != NULLPTR && this->utmp_->size() > 1
       && "tmp space not set, or is insufficient");

   if ( bInit_ ) return;

   sidir_ = traits_.idir;
   sodir_ = traits_.odir;
 
   time_last_  = this->traits_.start_time ;
   ocycle_     = this->traits_.start_ocycle;


   // Find State's kinetic components:
   assert(this->eqn_ptr_ != NULL && "Equation implementation must be set");

   GSIZET   *iwhere=NULLPTR;
   GSIZET    nwhere=0;
   CompDesc *icomptype = &(this->eqn_ptr_->comptype());
   icomptype->contains(GSC_KINETIC, iwhere, nwhere);
   for ( GSIZET j=0; j<nwhere; j++ ) ikinetic_.push_back(iwhere[j]);

   if ( iwhere != NULLPTR ) delete [] iwhere;

   ku_.resize(ikinetic_.size());
   

   bInit_ = TRUE;
 
} // end of method init


//**********************************************************************************
//**********************************************************************************
// METHOD     : do_kinetic_L2
// DESCRIPTION: Compute integrated diagnostic quantities, and output to file
// ARGUMENTS  : t  : state time
//              uu : state variable
//              uf : forcing
//              fname: file name
// RETURNS    : none.
//**********************************************************************************
template<typename EquationType>
void GGlobalDiag_basic<EquationType>::do_kinetic_L2(const Time t, const State &u, const State &uf, const GString fname)
{
  assert(utmp_ != NULLPTR && utmp_->size() > 3
      && "tmp space not set, or is insufficient");

  
  GBOOL   isreduced= FALSE;
  GBOOL   ismax    = FALSE;
  GINT    ndim = grid_->gtype() == GE_2DEMBEDDED ? 3 : GDIM;
  GFTYPE absu, absw, ener, enst, hel, fv, rhel;
  GTVector<GFTYPE> lmax(5), gmax(5);

  // Make things a little easier:
  GTVector<GTVector<GFTYPE>*> utmp(5);
  for ( GINT j=0; j<5; j++ ) utmp[j] = (*utmp_)[j];

  // Find kinetic components to operate on:
  for ( GINT j=0; j<ikinetic_.size(); j++ ) ku_[j] = u[ikinetic_[j]];


  // Energy = <u^2>/2
  lmax[0] = GMTK::energy(*grid_, ku_, utmp, isreduced, ismax);
 

  // Enstrophy = <omega^2>/2
  lmax[1] = 0.0;
  if ( ku_.size() == 1 ) {
    for ( GINT j=0; j<ndim; j++ ) {
      GMTK::grad<GFTYPE>(*grid_, *ku_[0], j+1, utmp, *utmp[2]);
      utmp[2]->pow(2);
      lmax[1] += grid_->integrate(*utmp[2],*utmp[0], isreduced); 
    }
    lmax[1] *= 0.5*grid_->ivolume();
  }
  else {
    lmax[1] = GMTK::enstrophy(*grid_, ku_, utmp, isreduced, ismax);
  }

  // Energy injection = <f.u>
  lmax[2] = GMTK::energyinj(*grid_, ku_, uf, utmp, isreduced, ismax);

  // Helicity = <u.omega>
  lmax[3] = GMTK::helicity(*grid_, ku_, utmp, isreduced, ismax);

  // Relative helicity = <u.omega/(|u|*|omega|)>
  lmax[4] = GMTK::relhelicity(*grid_, ku_, utmp, isreduced, ismax);

  // Gather final max's:
  GComm::Allreduce(lmax.data(), gmax.data(), 5, T2GCDatatype<GFTYPE>(), GC_OP_SUM, grid_->get_comm());
  ener = gmax[0]; enst = gmax[1]; fv = gmax[2]; hel = gmax[3]; rhel = gmax[4];

  // Print data to file:
  std::ifstream itst;
  std::ofstream ios;
  GString       fullfile = sodir_ + "/" + fname;

  if ( myrank_ == 0 ) {
    itst.open(fullfile);
    ios.open(fullfile,std::ios_base::app);

    // Write header, if required:
    if ( itst.peek() == std::ofstream::traits_type::eof() ) {
    ios << "#time    KE     Enst     f.v    hel     rhel " << std::endl;
    }
    itst.close();

    ios << t  
        << "    " << ener  << "    "  << enst 
        << "    " << fv    << "    "  << hel
        << "    " << rhel  
        << std::endl;
    ios.close();
  }
 
} // end of method do_kinetic_L2


//**********************************************************************************
//**********************************************************************************
// METHOD     : do_kinetic_max
// DESCRIPTION: Compute max quantities, and output to file
// ARGUMENTS  : t    : state time
//              u    : state variable
//              uf   : forcing
//              fname: file name
// RETURNS    : none.
//**********************************************************************************
template<typename EquationType>
void GGlobalDiag_basic<EquationType>::do_kinetic_max(const Time t, const State &u, const State &uf, const GString fname)
{
  assert(utmp_ != NULLPTR && utmp_->size() > 5
      && "tmp space not set, or is insufficient");

  GBOOL   isreduced= FALSE;
  GBOOL   ismax    = TRUE;
  GINT   ndim = grid_->gtype() == GE_2DEMBEDDED ? 3 : GDIM;
  GFTYPE absu, absw, ener, enst, hel, fv, rhel;
  GTVector<GFTYPE> lmax(5), gmax(5);

  // Make things a little easier:
  GTVector<GTVector<GFTYPE>*> utmp(6);
  for ( GINT j=0; j<6; j++ ) utmp[j] = (*utmp_)[j];

  // Find kinetic components to operate on:
  for ( GINT j=0; j<ikinetic_.size(); j++ ) ku_[j] = u[ikinetic_[j]];

  // Energy = <u^2>/2
  lmax[0] = GMTK::energy(*grid_, ku_, utmp, isreduced, ismax);
 
  // Enstrophy = <omega^2>/2
  lmax[1] = 0.0;
  if ( ku_.size() == 1 ) {
    for ( GINT j=0; j<ndim; j++ ) {
      GMTK::grad<GFTYPE>(*grid_, *ku_[0], j+1, utmp, *utmp[2]);
      utmp[2]->pow(2);
      lmax[1] += grid_->integrate(*utmp[2],*utmp[0], isreduced); 
    }
    lmax[1] *= 0.5*grid_->ivolume();
  }
  else {
    lmax[1] = GMTK::enstrophy(*grid_, ku_, utmp, isreduced, ismax);
  }

  // Energy injection = <f.u>
  lmax[2] = GMTK::energyinj(*grid_, ku_, uf, utmp, isreduced, ismax);

  // Helicity = <u.omega>
  lmax[3] = GMTK::helicity(*grid_, ku_, utmp, isreduced, ismax);

  // Relative helicity = <u.omega/(|u|*|omega|)>
  lmax[4] = GMTK::relhelicity(*grid_, ku_, utmp, isreduced, ismax);


  // Gather final max's:
  GComm::Allreduce(lmax.data(), gmax.data(), 5, T2GCDatatype<GFTYPE>(), GC_OP_MAX, grid_->get_comm());
  ener = gmax[0]; enst = gmax[1]; fv = gmax[2]; hel = gmax[3]; rhel = gmax[4];
  

  // Print data to file:
  std::ifstream itst;
  std::ofstream ios;
  GString       fullfile = sodir_ + "/" + fname;

  if ( myrank_ == 0 ) {
    itst.open(fullfile);
    ios.open(fullfile,std::ios_base::app);

    // Write header, if required:
    if ( itst.peek() == std::ofstream::traits_type::eof() ) {
    ios << "#time    KE     Enst     f.v    hel     rhel " << std::endl;
    }
    itst.close();

    ios << t  
        << "    " << ener  << "    "  << enst 
        << "    " << fv    << "    "  << hel
        << "    " << rhel  
        << std::endl;
    ios.close();
  }
 
} // end of method do_kinetic_max

;
