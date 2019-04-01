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
bprgrid_        (TRUE),
cycle_          (0),
ocycle_         (1),
cycle_last_     (0),
time_last_      (0.0),
bInit_          (FALSE),
ivol_           (1.0)
{ 
  this->traits_ = traits;
  this->grid_   = &grid;
  myrank_       = GComm::WorldRank(grid.get_comm());
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
   
  if ( (this->traits_.itype == ObserverBase<EquationType>::OBS_CYCLE 
        && (cycle_-cycle_last_) == this->traits_.cycle_interval)
    || (this->traits_.itype == ObserverBase<EquationType>::OBS_TIME  
        &&  t-time_last_ >= this->traits_.time_interval) ) {

    do_global(t, u, uf, "gbalance.txt");
    do_max   (t, u, uf, "gmax.txt");
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
void GGlobalDiag_basic<EquationType>::init(const Time t, const State &u)
{
   assert(utmp_ != NULLPTR && utmp_->size() > 1
       && "tmp space not set, or is insufficient");

   sdir_ = this->traits_.dir;
 
   if ( cycle_ == 0 ) {
     time_last_ = t; 
   }

   *(*utmp_)[0] = 1.0;
   ivol_ = 1.0/grid_->integrate(*(*utmp_)[0],*(*utmp_)[1]); 

   bInit_ = TRUE;
 
} // end of method init


//**********************************************************************************
//**********************************************************************************
// METHOD     : do_global
// DESCRIPTION: Compute global quantities, and output to file
// ARGUMENTS  : t  : state time
//              u  : state variable
//              uf : forcing
//              fname: file name
// RETURNS    : none.
//**********************************************************************************
template<typename EquationType>
void GGlobalDiag_basic<EquationType>::do_global(const Time t, const State &u, const State &uf, const GString fname)
{
  assert(utmp_ != NULLPTR && utmp_->size() > 3
      && "tmp space not set, or is insufficient");

  GFTYPE absu, absw, ener, enst, hel, fv, rhel;


  // Energy = <u^2>/2:
  ener = 0.0;
  for ( GSIZET j=0; j<u.size(); j++ ) {
   *(*utmp_)[0] = *u[j];
    (*utmp_)[0]->pow(2);
    ener += grid_->integrate(*(*utmp_)[0],*(*utmp_)[1]); 
  }
  ener *= 0.5*ivol_;
 
  // Enstrophy = <omega^2>/2
  enst = 0.0;
  if ( GDIM == 2 && u.size() == 2 ) {
   *(*utmp_)[0] = *u[j];
    GMTK::curl(*grid_, u, 3, *utmp_, *(*utmp_)[2]);
    (*utmp_)[2]->pow(2);
    enst += grid_->integrate(*(*utmp_)[2],*(*utmp_)[0]); 
  }
  else {
    for ( GSIZET j=0; j<GDIM; j++ ) {
      GMTK::curl(*grid_, u, j+1, *utmp_, *(*utmp_)[2]);
      (*utmp_)[2]->pow(2);
      enst += grid_->integrate(*(*utmp_)[2],*(*utmp_)[0]); 
    }
  }
  enst *= 0.5*ivol_;

  // Energy injection = <f.u>
  fv = 0.0;
  *(*utmp_)[1] = 0.0;
  for ( GSIZET j=0; j<GDIM; j++ ) {
    *(*utmp_)[1] = *uf[j];
    (*utmp_)[1]->pointProd(*u[j]);
    fv += grid_->integrate(*(*utmp_)[1],*(*utmp_)[0]); 
  }
  fv *= ivol_;

  // Helicity = <u.omega>
  hel = 0.0;
  if ( GDIM > 2 || u.size() > 2 ) {
    for ( GSIZET j=0; j<GDIM; j++ ) {
      GMTK::curl(*grid_, u, j+1, *utmp_, *(*utmp_)[2]);
      (*utmp_)[2]->pointProd(*u[j]);
      hel += grid_->integrate(*(*utmp_)[2],*(*utmp_)[0]); 
    }
  }
  hel *= ivol_;

  // Relative helicity = <u.omega/(|u|*|omega|)>
  rhel = 0.0;
  if ( GDIM > 2 || u.size() > 2 ) {
    // Compute |u|:
    *(*utmp_)[3] = 0.0;
    for ( GSIZET j=0; j<GDIM; j++ ) {
     *(*utmp_)[1] = *u[j];
      (*utmp_)[1]->pow(2);
     *(*utmp_)[3] += *(*utmp_)[1]
    }
    (*utmp_)[3]->pow(0.5);
    
    // Compute |curl u| = |omega|:
    *(*utmp_)[4] = 0.0;
    for ( GSIZET j=0; j<GDIM; j++ ) {
      GMTK::curl(*grid_, u, j+1, *utmp_, *(*utmp_)[2]);
      (*utmp_)[2]->pow(2);
     *(*utmp_)[4] += *(*utmp_)[2];
    }
    (*utmp_)[4]->pow(0.5);

    // Create 1/|u| |omega| :
    GFTYPE tiny = std::numeric_limits<GFTYPE>::epsilon();
    for ( GSIZET j=0; j<u[0]->size(); j++ )  
      (*(*utmp_)[3])[k] = 1.0/( (*(*utmp_)[3])[k] * (*(*utmp_)[4])[k] + tiny );

    // Compute <u.omega / |u| |omega| >:
    for ( GSIZET j=0; j<GDIM; j++ ) {
      GMTK::curl(*grid_, u, j+1, *utmp_, *(*utmp_)[2]);
      (*utmp_)[2]->pointProd(u[j]);
      (*utmp_)[2]->pointProd(*(*utmp_)[3]);
      rhel += grid_->integrate(*(*utmp_)[2],*(*utmp_)[0]); 
    }
  }
  rhel *= ivol_;


  // Print data to file:
  std::ifstream itst;
  std::ofstream ios;
  GString       fullfile = sdir_ + fname;

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
        << "    " << rhel  << 
        << std::endl;
    ios.close();
  }
 
} // end of method do_global


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
  assert(utmp_ != NULLPTR && utmp_->size() > 3
      && "tmp space not set, or is insufficient");

  GFTYPE absu, absw, ener, enst, hel, fv, rhel;


  // Energy = u^2/2:
  *(*utmp_)[1] = 0.0;
  for ( GSIZET j=0; j<u.size(); j++ ) {
   *(*utmp_)[0] = *u[j];
    (*utmp_)[0]->pow(2);
   *(*utmp_)[1] += *(*utmp_)[0];
  }
  ener = 0.5*(*utmp_)[0]->max();
 
  // Enstrophy = omega^2/2
  *(*utmp_)[3] = 0.0;
  if ( GDIM == 2 && u.size() == 2 ) {
   *(*utmp_)[0] = *u[j];
    GMTK::curl(*grid_, u, 3, *utmp_, *(*utmp_)[2]);
    (*utmp_)[2]->pow(2);
   *(*utmp_)[3] += *(*utmp_)[2];
  }
  else {
    for ( GSIZET j=0; j<GDIM; j++ ) {
      GMTK::curl(*grid_, u, j+1, *utmp_, *(*utmp_)[2]);
      (*utmp_)[2]->pow(2);
     *(*utmp_)[3] += *(*utmp_)[2];
    }
  }
  enst = 0.5*(*utmp_)[3]->max();

  // Energy injection = f.u
  *(*utmp_)[3] = 0.0;
  for ( GSIZET j=0; j<GDIM; j++ ) {
    *(*utmp_)[1] = *uf[j];
    (*utmp_)[1]->pointProd(*u[j]);
  }
  (*utmp_)[3]->abs();
  fv = (*utmp_)[3]->max();

  // Helicity = u.omega
  *(*utmp_)[3] = 0.0;
  if ( GDIM > 2 || u.size() > 2 ) {
    for ( GSIZET j=0; j<GDIM; j++ ) {
      GMTK::curl(*grid_, u, j+1, *utmp_, *(*utmp_)[2]);
      (*utmp_)[2]->pointProd(*u[j]);
     *(*utmp_)[3] += *(*utmp_)[2];
    }
  }
  (*utmp_)[3]->abs();
  hel = (*utmp_)[3]->max();

  // Relative helicity = u.omega/(|u|*|omega|)
  *(*utmp_)[5] = 0.0;
  if ( GDIM > 2 || u.size() > 2 ) {
    // Compute |u|:
    *(*utmp_)[3] = 0.0;
    for ( GSIZET j=0; j<GDIM; j++ ) {
     *(*utmp_)[1] = *u[j];
      (*utmp_)[1]->pow(2);
     *(*utmp_)[3] += *(*utmp_)[1]
    }
    (*utmp_)[3]->pow(0.5);
    
    // Compute |curl u| = |omega|:
    *(*utmp_)[4] = 0.0;
    for ( GSIZET j=0; j<GDIM; j++ ) {
      GMTK::curl(*grid_, u, j+1, *utmp_, *(*utmp_)[2]);
      (*utmp_)[2]->pow(2);
     *(*utmp_)[4] += *(*utmp_)[2];
    }
    (*utmp_)[4]->pow(0.5);

    // Create 1/|u| |omega| :
    GFTYPE tiny = std::numeric_limits<GFTYPE>::epsilon();
    for ( GSIZET j=0; j<u[0]->size(); j++ )  
      (*(*utmp_)[3])[k] = 1.0/( (*(*utmp_)[3])[k] * (*(*utmp_)[4])[k] + tiny );

    // Compute u.omega / |u| |omega|: 
    for ( GSIZET j=0; j<GDIM; j++ ) {
      GMTK::curl(*grid_, u, j+1, *utmp_, *(*utmp_)[2]);
      (*utmp_)[2]->pointProd(u[j]);
      (*utmp_)[2]->pointProd(*(*utmp_)[3]);
     *(*utmp_)[5] += *(*utmp_)[2];
    }
  }
  (*utmp_)[5]->abs();
  rhel = (*utmp_)[5]->max();


  // Print data to file:
  std::ifstream itst;
  std::ofstream ios;
  GString       fullfile = sdir_ + fname;

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
        << "    " << rhel  << 
        << std::endl;
    ios.close();
  }
 
} // end of method do_max

;
