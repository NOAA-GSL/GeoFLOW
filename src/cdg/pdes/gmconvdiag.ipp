//==================================================================================
// Module       : gmconvdiag.ipp
// Date         : 3/18/21 (DLR)
// Description  : Observer object for carrying out L2 & extrema diagnostics for
//                GMConv solver
// Copyright    : Copyright 2021. Colorado State University. All rights reserved.
// Derived From : ObserverBase.
//==================================================================================

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Instantiate with Traits
// ARGS   : traits: Traits sturcture
//**********************************************************************************
template<typename EquationType>
GMConvDiag<EquationType>::GMConvDiag(EqnBasePtr &equation, Grid &grid, typename ObserverBase<EquationType>::Traits &traits):
ObserverBase<EquationType>(equation, grid, traits),
bInit_          (FALSE),
cycle_          (0),
ocycle_         (0),
cycle_last_     (0),
time_last_      (0.0),
grid_           (&grid)
{ 
  traits_ = traits;
  utmp_   = static_cast<GTVector<GTVector<GFTYPE>*>*>(utmp_);
  myrank_ = GComm::WorldRank(grid.get_comm());

  solver_ = dynamic_cast<GMConv<EquationType>*>(equation.get());
  assert(solver_); // must be for GMConv

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
void GMConvDiag<EquationType>::observe_impl(const Time &t, const State &u, const State &uf)
{
  StateInfo info;

  init_impl(info);

  mpixx::communicator comm;

  if ( ((traits_.itype == ObserverBase<EquationType>::OBS_CYCLE)
        && ((cycle_-cycle_last_+1) >= traits_.cycle_interval))
    || ((traits_.itype == ObserverBase<EquationType>::OBS_TIME)
        &&  (t-time_last_ >= traits_.time_interval))
    || (cycle_ == 0) ) {

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
// METHOD     : init_impl
// DESCRIPTION: Fill member index and name data based on traits
// ARGUMENTS  : info: StateInfo structure
// RETURNS    : none.
//**********************************************************************************
template<typename EquationType>
void GMConvDiag<EquationType>::init_impl(StateInfo &info)
{
   assert(utmp_ != NULLPTR && utmp_->size() > 1
       && "tmp space not set, or is insufficient");

   if ( bInit_ ) return;

   sidir_ = traits_.idir;
   sodir_ = traits_.odir;
 
   time_last_  = this->traits_.start_time ;
   ocycle_     = this->traits_.start_ocycle;

   bInit_ = TRUE;
 
} // end of method init_impl


//**********************************************************************************
//**********************************************************************************
// METHOD     : do_L2
// DESCRIPTION: Compute integrated diagnostic quantities, and output to file
// ARGUMENTS  : t  : state time
//              u  : state variable
//              uf : forcing
//              fname: file name
// RETURNS    : none.
//**********************************************************************************
template<typename EquationType>
void GMConvDiag<EquationType>::do_L2(const Time t, const State &u, const State &uf, const GString fname)
{
  assert(utmp_ != NULLPTR && utmp_->size() > 3
      && "tmp space not set, or is insufficient");

  
  GBOOL   isreduced= FALSE;
  GBOOL   ismax    = FALSE;
  GINT    nd, ndim = grid_->gtype() == GE_2DEMBEDDED ? 3 : GDIM;
  GFTYPE absu, absw, eint, ke, mass;
  GTVector<GFTYPE> *d, *e;
  GTVector<GFTYPE> lmax(3), gmax(3);
  typename GMConv<EquationType>::Traits trsolver;

  trsolver = solver_->get_traits();

  // Make things a little easier:
  GTVector<GTVector<GFTYPE>*> utmp(3);
  for ( auto j=0; j<utmp.size(); j++ ) utmp[j] = (*utmp_)[j];

  // Find internal energy density, <e>:
  e = u[solver_->ENERGY];
  lmax[1] = grid_->integrate(*e, *utmp[0], FALSE);

  // Find local integrated mass:
  d = u[solver_->DENSITY];
 *utmp[1] = *d;
  if ( trsolver.usebase ) *utmp[1] += *u[solver_->BASESTATE];
  lmax[0] = grid_->integrate(*utmp[1], *utmp[0], FALSE);

  // Find kinetic energy density,  <0.5 rho v^2>
  *utmp[0] = *utmp[1]; utmp[0]->rpow(-1); // inverse total density
  solver_->compute_v(u, 1, *utmp[0], *utmp[2]); utmp[2]->rpow(2);
  for ( auto j=1; j<ndim; j++ ) {
    solver_->compute_v(u, j+1, *utmp[0], *utmp[1]);
    utmp[1]->rpow(2);
   *utmp[2] +=  *utmp[1];
  }
 *utmp[1] = *d;
  if ( trsolver.usebase ) *utmp[1] += *u[solver_->BASESTATE];
  utmp[2]->apointProd(0.5, *utmp[1]); // d v^2
  lmax[2] = grid_->integrate(*utmp[2], *utmp[0], FALSE);


  // Gather final sums:
  GComm::Allreduce(lmax.data(), gmax.data(), 3, T2GCDatatype<GFTYPE>(), GC_OP_SUM, grid_->get_comm());
  mass = gmax[0]; eint = gmax[1]/grid_->volume(); ke = gmax[2]/grid_->volume(); 

  // Print data to file:
  GBOOL         doheader=FALSE;
  GSIZET        le, ne;
  GFTYPE        dxmin, elmin, elmax, elavg;
  std::ofstream ios;
  GString       fullfile = sodir_ + "/" + fname;

  if ( geoflow::file_empty(fullfile) ) {
    le    = grid_->nelems();
    dxmin = grid_->minnodedist();
    elmin = grid_->minlength();
    elmax = grid_->maxlength();
    elavg = grid_->avglength();
    GComm::Allreduce(&le, &ne, 1, T2GCDatatype<GSIZET>() , GC_OP_SUM, grid_->get_comm());
    doheader = TRUE;
  }

  if ( myrank_ == 0 ) {
    ios.open(fullfile,std::ios_base::app);
    if ( doheader ) {
      ios << "#nelems=" << ne << " dxmin=" << dxmin << " elmin=" << elmin << " elmax=" << elmax << " elavg=" << elavg << std::endl;
      ios << "#time      Mass       <e>         <d KE>     " << std::endl;
    }

    ios << t  << scientific << setprecision(15) 
        << "    " << mass  << "    "  << eint
        << "    " << ke
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
void GMConvDiag<EquationType>::do_max(const Time t, const State &u, const State &uf, const GString fname)
{
  assert(utmp_ != NULLPTR && utmp_->size() > 3
      && "tmp space not set, or is insufficient");

  GBOOL   isreduced= FALSE;
  GBOOL   ismax    = FALSE;
  GINT    nd, ndim = grid_->gtype() == GE_2DEMBEDDED ? 3 : GDIM;
  GFTYPE absu, absw, eint, ke, mass;
  GTVector<GFTYPE> *d, *e;
  GTVector<GFTYPE> lmax(3), gmax(3);
  typename GMConv<EquationType>::Traits trsolver;

  trsolver = solver_->get_traits();

  // Make things a little easier:
  GTVector<GTVector<GFTYPE>*> utmp(3);
  for ( auto j=0; j<utmp.size(); j++ ) utmp[j] = (*utmp_)[j];

  // Find internal energy density:
  e = u[solver_->ENERGY];
  lmax[1] = e->amax();

  // Find local integrated mass:
  d = u[solver_->DENSITY];
 *utmp[1] = *d;
  if ( trsolver.usebase ) *utmp[1] += *u[solver_->BASESTATE];
  lmax[0] = utmp[1]->amax();

  // Find kinetic energy density:  0.5 rho v^2
  *utmp[0] = *utmp[1]; utmp[0]->rpow(-1); // inverse total density
  solver_->compute_v(u, 1, *utmp[0], *utmp[2]); utmp[2]->rpow(2);
  for ( auto j=1; j<ndim; j++ ) {
    solver_->compute_v(u, j+1, *utmp[0], *utmp[1]);
    utmp[1]->rpow(2);
   *utmp[2] +=  *utmp[1];
  }
 *utmp[1] = *d;
  if ( trsolver.usebase ) *utmp[1] += *u[solver_->BASESTATE];
  utmp[2]->apointProd(0.5, *utmp[1]);;
  lmax[2] = utmp[2]->amax();


  // Gather final extrema:
  GComm::Allreduce(lmax.data(), gmax.data(), 3, T2GCDatatype<GFTYPE>(), GC_OP_MAX, grid_->get_comm());
  mass = gmax[0]; eint = gmax[1]; ke = gmax[2]; 


  // Print data to file:
  GBOOL         doheader=FALSE;
  GSIZET        le, ne;
  GFTYPE        dxmin, elmin, elmax, elavg;
  std::ofstream ios;
  GString       fullfile = sodir_ + "/" + fname;

  if ( geoflow::file_empty(fullfile) ) {
    le    = grid_->nelems();
    dxmin = grid_->minnodedist();
    elmin = grid_->minlength();
    elmax = grid_->maxlength();
    elavg = grid_->avglength();
    GComm::Allreduce(&le, &ne, 1, T2GCDatatype<GSIZET>() , GC_OP_SUM, grid_->get_comm());
    doheader = TRUE;
  }
  if ( myrank_ == 0 ) {
    ios.open(fullfile,std::ios_base::app);
    if ( doheader ) {
      ios << "#nelems=" << ne << " dxmin=" << dxmin << " elmin=" << elmin << " elmax=" << elmax << " elavg=" << elavg << std::endl;
      ios << "#time      den_tot        e           d KE      " << std::endl;
    }

    ios << t  << scientific << setprecision(15) 
        << "    " << mass  << "    "  << eint
        << "    " << ke
        << std::endl;
    ios.close();
  }
 
} // end of method do_max

