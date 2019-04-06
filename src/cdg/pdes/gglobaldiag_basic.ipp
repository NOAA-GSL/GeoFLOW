//==================================================================================
// Module       : gglobaldiag_basic.ipp
// Date         : 3/28/19 (DLR)
// Description  : Observer object for carrying out L2 & extrema diagnostics for
//                Burgers equation.
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
GGlobalDiag_basic<EquationType>::GGlobalDiag_basic(typename ObserverBase<EquationType>::Traits &traits, Grid &grid):
bInit_          (FALSE),
cycle_          (0),
ocycle_         (1),
cycle_last_     (1),
time_last_      (0.0),
ivol_           (1.0)
{ 
  traits_ = traits;
  grid_   = &grid;
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

    do_L2 (t, u, uf, "gbalance.txt");
    do_max(t, u, uf, "gmax.txt");
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
void GGlobalDiag_basic<EquationType>::init(const Time t, const State &u)
{
   assert(utmp_ != NULLPTR && this->utmp_->size() > 1
       && "tmp space not set, or is insufficient");

   if ( bInit_ ) return;

   sidir_ = traits_.idir;
   sodir_ = traits_.odir;
 
   time_last_  = this->traits_.start_time ;
   ocycle_     = this->traits_.start_ocycle;

   *(*utmp_)[0] = 1.0;
   GFTYPE vol = grid_->integrate(*(*utmp_)[0],*(*utmp_)[1]);
   assert(vol > 0.0 && "Invalid volume integral");
   ivol_ = 1.0/vol;

   bInit_ = TRUE;
 
} // end of method init


//**********************************************************************************
//**********************************************************************************
// METHOD     : do_L2
// DESCRIPTION: Compute integrated diagnostic quantities, and output to file
// ARGUMENTS  : t  : state time
//              uu : state variable
//              uf : forcing
//              fname: file name
// RETURNS    : none.
//**********************************************************************************
template<typename EquationType>
void GGlobalDiag_basic<EquationType>::do_L2(const Time t, const State &u, const State &uf, const GString fname)
{
  assert(utmp_ != NULLPTR && utmp_->size() > 3
      && "tmp space not set, or is insufficient");

  GFTYPE absu, absw, ener, enst, hel, fv, rhel;

  // Make things a little easier:
  GTVector<GTVector<GFTYPE>*> utmp(4);
  for ( GSIZET j=0; j<4; j++ ) utmp[j] = (*utmp_)[j];

  // Energy = <u^2>/2:
  ener = 0.0;
  for ( GSIZET j=0; j<u.size(); j++ ) {
   *utmp[0] = *u[j];
    utmp[0]->pow(2);
    ener += grid_->integrate(*utmp[0],*utmp[1]); 
  }
  ener *= 0.5*ivol_;
 
  // Enstrophy = <omega^2>/2
  enst = 0.0;
  if ( GDIM == 2 && u.size() == 2 ) {
    GMTK::curl<GFTYPE>(*grid_, u, 3, utmp, *utmp[2]);
    utmp[2]->pow(2);
    enst += grid_->integrate(*utmp[2],*utmp[0]); 
  }
  else {
    for ( GINT j=0; j<GDIM; j++ ) {
      GMTK::curl<GFTYPE>(*grid_, u, j+1, utmp, *utmp[2]);
      utmp[2]->pow(2);
      enst += grid_->integrate(*utmp[2],*utmp[0]); 
    }
  }
  enst *= 0.5*ivol_;

  // Energy injection = <f.u>
  fv = 0.0;
  if ( uf.size() > 0 && uf[0] != NULLPTR 
    && uf[0]->size() > 0 ) {
    *utmp[1] = 0.0;
    for ( GINT j=0; j<GDIM; j++ ) {
      *utmp[1] = *uf[j];
      utmp[1]->pointProd(*u[j]);
      fv += grid_->integrate(*utmp[1],*utmp[0]); 
    }
    fv *= ivol_;
  }

  // Helicity = <u.omega>
  hel = 0.0;
  if ( GDIM > 2 || u.size() > 2 ) {
    for ( GINT j=0; j<GDIM; j++ ) {
      GMTK::curl<GFTYPE>(*grid_, u, j+1, utmp, *utmp[2]);
      utmp[2]->pointProd(*u[j]);
      hel += grid_->integrate(*utmp[2],*utmp[0]); 
    }
  }
  hel *= ivol_;

  // Relative helicity = <u.omega/(|u|*|omega|)>
  rhel = 0.0;
  if ( GDIM > 2 || u.size() > 2 ) {
    // Compute |u|:
    *utmp[3] = 0.0;
    for ( GINT j=0; j<GDIM; j++ ) {
     *utmp[1] = *u[j];
      utmp[1]->pow(2);
     *utmp[3] += *utmp[1];
    }
    utmp[3]->pow(0.5);
    
    // Compute |curl u| = |omega|:
    *utmp[4] = 0.0;
    for ( GINT j=0; j<GDIM; j++ ) {
      GMTK::curl<GFTYPE>(*grid_, u, j+1, utmp, *utmp[2]);
      utmp[2]->pow(2);
     *utmp[4] += *utmp[2];
    }
    utmp[4]->pow(0.5);

    // Create 1/|u| |omega| :
    GFTYPE tiny = std::numeric_limits<GFTYPE>::epsilon();
    for ( GSIZET k=0; k<u[0]->size(); k++ )  
      (*utmp[3])[k] = 1.0/( (*utmp[3])[k] * (*utmp[4])[k] + tiny );

    // Compute <u.omega / |u| |omega| >:
    for ( GINT j=0; j<GDIM; j++ ) {
      GMTK::curl<GFTYPE>(*grid_, u, j+1, utmp, *utmp[2]);
      utmp[2]->pointProd(*u[j]);
      utmp[2]->pointProd(*utmp[3]);
      rhel += grid_->integrate(*utmp[2],*utmp[0]); 
    }
  }
  rhel *= ivol_;


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
 
} // end of method do_L2


//**********************************************************************************
//**********************************************************************************
// METHOD     : do_max
// DESCRIPTION: Compute max quantities, and output to file
// ARGUMENTS  : t    : state time
//              u    : state variable
//              uf   : forcing
//              fname: file name
// RETURNS    : none.
//**********************************************************************************
template<typename EquationType>
void GGlobalDiag_basic<EquationType>::do_max(const Time t, const State &u, const State &uf, const GString fname)
{
  assert(utmp_ != NULLPTR && utmp_->size() > 5
      && "tmp space not set, or is insufficient");

  GFTYPE absu, absw, ener, enst, hel, fv, rhel;

  // Make things a little easier:
  GTVector<GTVector<GFTYPE>*> utmp(6);
  for ( GSIZET j=0; j<6; j++ ) utmp[j] = (*utmp_)[j];


  // Energy = u^2/2:
  *utmp[1] = 0.0;
  for ( GSIZET j=0; j<u.size(); j++ ) {
   *utmp[0] = *u[j];
    utmp[0]->pow(2);
   *utmp[1] += *utmp[0];
  }
  ener = 0.5*utmp[0]->max();
 
  // Enstrophy = omega^2/2
  *utmp[3] = 0.0;
  if ( GDIM == 2 && u.size() == 2 ) {
    GMTK::curl<GFTYPE>(*grid_, u, 3, utmp, *utmp[2]);
    utmp[2]->pow(2);
   *utmp[3] += *utmp[2];
  }
  else {
    for ( GINT j=0; j<GDIM; j++ ) {
      GMTK::curl<GFTYPE>(*grid_, u, j+1, utmp, *utmp[2]);
      utmp[2]->pow(2);
     *utmp[3] += *utmp[2];
    }
  }
  enst = 0.5*utmp[3]->max();

  // Energy injection = f.u
  fv = 0.0;
  if ( uf.size() > 0 && uf[0] != NULLPTR 
    && uf[0]->size() > 0 ) {
    *utmp[3] = 0.0;
    for ( GINT j=0; j<GDIM; j++ ) {
      *utmp[1] = *uf[j];
      utmp[1]->pointProd(*u[j]);
    }
    fv = utmp[3]->amax();
  }

  // Helicity = u.omega
  *utmp[3] = 0.0;
  if ( GDIM > 2 || u.size() > 2 ) {
    for ( GINT j=0; j<GDIM; j++ ) {
      GMTK::curl<GFTYPE>(*grid_, u, j+1, utmp, *utmp[2]);
      utmp[2]->pointProd(*u[j]);
     *utmp[3] += *utmp[2];
    }
  }
  hel = utmp[3]->amax();

  // Relative helicity = u.omega/(|u|*|omega|)
  *utmp[5] = 0.0;
  if ( GDIM > 2 || u.size() > 2 ) {
    // Compute |u|:
    *utmp[3] = 0.0;
    for ( GINT j=0; j<GDIM; j++ ) {
     *utmp[1] = *u[j];
      utmp[1]->pow(2);
     *utmp[3] += *utmp[1];
    }
    utmp[3]->pow(0.5);
    
    // Compute |curl u| = |omega|:
    *utmp[4] = 0.0;
    for ( GINT j=0; j<GDIM; j++ ) {
      GMTK::curl<GFTYPE>(*grid_, u, j+1, utmp, *utmp[2]);
      utmp[2]->pow(2);
     *utmp[4] += *utmp[2];
    }
    utmp[4]->pow(0.5);

    // Create 1/|u| |omega| :
    GFTYPE tiny = std::numeric_limits<GFTYPE>::epsilon();
    for ( GSIZET k=0; k<u[0]->size(); k++ )  
      (*utmp[3])[k] = 1.0/( (*utmp[3])[k] * (*utmp[4])[k] + tiny );

    // Compute u.omega / |u| |omega|: 
    for ( GINT j=0; j<GDIM; j++ ) {
      GMTK::curl<GFTYPE>(*grid_, u, j+1, utmp, *utmp[2]);
      utmp[2]->pointProd(*u[j]);
      utmp[2]->pointProd(*utmp[3]);
     *utmp[5] += *utmp[2];
    }
  }
  rhel = utmp[5]->amax();

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
 
} // end of method do_max

;
