//==================================================================================
// Module       : gmconv.ipp
// Date         : 6/11/20 (DLR)
// Description  : Object defining a moist convection solver:
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : EquationBase.
//==================================================================================

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method  (1)
// DESC   : Instantiate with grid + state + tmp. 
//          grid      : grid object
//          traits    : GMConv:Traits struct
// RETURNS: none
//**********************************************************************************
template<typename TypePack>
GMConv<TypePack>::GMConv(Grid &grid, GMConv<TypePack>::Traits &traits) :
EquationBase<TypePack>(),
bInit_                   (FALSE),
bforced_                 (FALSE),
bsteptop_                (FALSE),
bvterm_                  (FALSE),
istage_                      (0),
nevolve_                     (0),
nhydro_                      (0),
nmoist_                      (0),
icycle_                      (0),
nc_                       (GDIM),
gmass_                 (NULLPTR),
gimass_                (NULLPTR),
/*
gflux_                 (NULLPTR),
*/
ghelm_                 (NULLPTR),
gstressen_             (NULLPTR),
gadvect_               (NULLPTR),
gexrk_                 (NULLPTR),
gpdv_                  (NULLPTR),
gdiv_                  (NULLPTR),
grid_                    (&grid),
ggfx_         (&grid.get_ggfx()),
comm_(grid.get_ggfx().getComm()),
steptop_callback_      (NULLPTR)
{

  GGridIcos *icos = dynamic_cast<GGridIcos*>(grid_);

  traits_.iforced.resize(traits.iforced.size());
  traits_ = traits;

  nc_ = grid_->gtype() == GE_2DEMBEDDED ? 3 : GDIM;

  if ( GDIM == 2 && icos ) {
    assert( !traits_.dograv && !traits_.dofallout  
         && "Embedded 2D spherical grid not allowed with preferred directions");
  }

  if ( traits_.dofallout ) {
    assert( !traits_.dodry
         && "Must do moist convection with fallout");
  }

  if ( traits_.usebase ) {
    assert( traits_.dograv
         && "Must dograv with base state");
  }

} // end of constructor method (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor method 
// DESC   : 
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<typename TypePack>
GMConv<TypePack>::~GMConv()
{
//if ( gflux_     != NULLPTR ) delete gflux_;
  if ( ghelm_     != NULLPTR ) delete ghelm_;
  if ( gstressen_ != NULLPTR ) delete gstressen_;
  if ( gadvect_   != NULLPTR ) delete gadvect_;
  if ( gpdv_      != NULLPTR ) delete gpdv_;
  if ( gdiv_      != NULLPTR ) delete gdiv_;
  if ( gexrk_     != NULLPTR ) delete gexrk_;

} // end, destructor


//**********************************************************************************
//**********************************************************************************
// METHOD : dt_impl
// DESC   : Compute time step, assuming a Courant number of 1:
//            dt = min_grid(dx/u)
// ARGS   : t : time
//          u : (full) state
//          dt: timestep, returned
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::dt_impl(const Time &t, State &u, Time &dt)
{
  GString    serr = "GMConv<TypePack>::dt_impl: ";
  Ftype      dtmin, dt1;
  StateComp *csq, *p;
  StateComp *rhoT, *tmp1, *tmp2;

  assert(utmp_.size() >= 6 );

  // This is an estimate. We assume the timestep is
  // is governed by fast sonic waves with speed
  //  |v| + c,
  // where 
  //   c^2 = p/rho; 
  // Then, dt is computed element-by-element from
  //   dt = dx_min/(|v| + c)_max
  // where min and max are computed over the element.
  // Here, approximate |v| + c as sqrt(v^2 + c^2)
   
  dtmin = std::numeric_limits<Ftype>::max();

  // Assign pointers:
  p    = utmp_[utmp_.size()-1];
  rhoT = utmp_[utmp_.size()-2];
  tmp1 = utmp_[utmp_.size()-3];
  tmp2 = utmp_[utmp_.size()-4];

 *rhoT = *u[DENSITY]; 
  if ( traits_.usebase ) *rhoT += *u[BASESTATE];
 *tmp1 = *rhoT; tmp1->rpow(-1.0);
  compute_v(u, *tmp1, v_); 

  compute_cv(u, *tmp1, *tmp2);                     // Cv
  geoflow::compute_temp<Ftype>(*u[ENERGY], *rhoT, *tmp2, *tmp1);  // temperature
  compute_qd(u, *tmp2);                            // dry mass ratio
  geoflow::compute_p<Ftype>(*tmp1, *rhoT, *tmp2, RD, *p); // partial pressure for dry air
  if ( !traits_.dodry ) {
    geoflow::compute_p<Ftype>(*tmp1, *rhoT, *u[VAPOR], RV, *tmp2); // partial pressure for vapor
   *p += *tmp2;
  }

  csq = utmp_[utmp_.size()-4];
  for ( auto j=0; j<p->size(); j++ ) { // sound speed, csq
    (*csq)[j] = (*p)[j] / (*rhoT)[j] ;
  }

   // Compute v^2 + c^2:
  *tmp1 = *v_[0]; tmp1->rpow(2);
   for ( auto k=1; k<v_.size(); k++ ) {    // each advecting v 
     for ( auto j=0; j<v_[0]->size(); j++ ) { // v^2
       (*tmp1)[j] += (*v_[k])[j] * (*v_[k])[j];
     }
   }
   *tmp1 += *csq; // v^2 + c^2

  
   // Compute max(v^2 + c^2) for each element:
   GMTK::maxbyelem<Ftype>(*grid_, *tmp1, maxbyelem_);
   
   // Note: maxbyelem_ is an array with the max of v^2 + c^2 
   //       on each element

   // Find estimate of smallest dt on this task:
   for ( auto e=1; e<dxmin_.size(); e++ ) { // check each element
     dt1 = dxmin_[e] * dxmin_[e] / maxbyelem_[e]; // this dt^2
     dtmin = MIN(dtmin, sqrt(dt1)); 
   }

   // Find minimum dt over all tasks:
   GComm::Allreduce(&dtmin, &dt1, 1, T2GCDatatype<Ftype>() , GC_OP_MIN, comm_);

   // Limit any timestep-to-timestep increae to 10%:
   dt = MIN(dt1*traits_.courant, 1.01*dt);

} // end of method dt_impl


//**********************************************************************************
//**********************************************************************************
// METHOD : dudt_impl
// DESC   : Compute RHS for explicit schemes
// ARGS   : t   : time
//          u   : state. 
//          uf  : forcing tendency state vector
//          ub  : bdy state vector
//          dt  : time step
//          dudt: accelerations, returned
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::dudt_impl(const Time &t, const State &u, const State &uf,  const State &ub, const Time &dt, Derivative &dudt)
{
  if ( traits_.dodry ) {
    dudt_dry(t, u, uf, ub, dt, dudt);
  }
  else {
    dudt_wet(t, u, uf, ub, dt, dudt);
  }
}


//**********************************************************************************
//**********************************************************************************
// METHOD : dudt_dry
// DESC   : Compute RHS for explicit schemes
// ARGS   : t   : time
//          u   : state. 
//          uf  : forcing tendency state vector
//          ub  : bdy state vector
//          dt  : time step
//          dudt: accelerations, returned
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::dudt_dry(const Time &t, const State &u, const State &uf,  const State &ub, const Time &dt, Derivative &dudt)
{
  assert(!traits_.bconserved &&
         "conservation not yet supported"); 

  GString    serr = "GMConv<TypePack>::dudt_dry: ";
  GINT       nice, nliq ;
  StateComp *irhoT, *Ltot;
  StateComp *dd, *e, *p, *rhoT, *T; // energy, den, pressure, temperature
  StateComp *Mass=grid_->massop().data();
  StateComp *tmp1, *tmp2;
  State      g(nc_); 

  // NOTE:
  // Make sure that, in init(), Helmholtz op is using only
  // weak Laplacian (q * mass term isn't being used), or there will 
  // be problems. This is required for explicit schemes, for
  // which this method is called.

  assert( !traits_.bconserved ); // don't allow conservative form yet

  // Set tmp pool for RHS computations:
  assert(urhstmp_.size() >= szrhstmp());
  Ltot  = urhstmp_[urhstmp_.size()-1];
  rhoT  = urhstmp_[urhstmp_.size()-2];
  irhoT = urhstmp_[urhstmp_.size()-3];
  tmp1  = urhstmp_[urhstmp_.size()-4];
  tmp2  = urhstmp_[urhstmp_.size()-5];

#if 0
  // Do filtering, if any:
  FilterList *fl = &this->get_filter_list();
  for ( auto j=0; j<fl->size(); j++ ) {
    if ( (*fl)[j] != NULLPTR ) {
      (*fl)[j]->apply(t, *u[j], urhstmp_);
    }
  }
#endif


  // Get total density and inverse: *rhoT  = *u[DENSITY]; 
 *rhoT = *u[DENSITY];
  if ( traits_.usebase ) *rhoT +=  *ubase_[0];   
 *irhoT = *rhoT; irhoT->rpow(-1.0); // 1/rhoT

  // Compute velocity for timestep:
  compute_v(s_, *irhoT, v_); // stored in v_
  
  // Compute all terms as though they are on the LHS, then
  // change the sign and divide by Mass at the end....


  // *************************************************************
  // Total density RHS:
  // *************************************************************

  gdiv_->apply(*rhoT, v_, urhstmp_, *dudt[DENSITY]); 

  if ( uf[DENSITY] != NULLPTR ) {
cout << "dudt_dry: uf[MOMENTUM]" << endl;
     *dudt[DENSITY] -= *uf[DENSITY];//  += sdot(s_rhoT)
  }
  
  p     = urhstmp_[urhstmp_.size()-3]; // holds pressure
  T     = urhstmp_[urhstmp_.size()-6]; // holds temperature
  e     = u[ENERGY];                   // internal energy density

  // *************************************************************
  // Energy equation RHS:
  // *************************************************************
  compute_cv(u, *tmp1, *tmp2);                     // Cv
  geoflow::compute_temp<Ftype>(*e, *rhoT, *tmp2, *T);     // temperature
  compute_qd  (u, *tmp1);                          // dry mass ratio
  geoflow::compute_p<Ftype>(*T, *rhoT, *tmp1, RD, *p);    // partial pressure for dry air

  GMTK::saxpby<Ftype>(*tmp1, *e, 1.0, *p, 1.0);    // h = p+e, enthalpy density

  gdiv_->apply(*tmp1, v_, urhstmp_, *dudt[ENERGY]); 

  gstressen_->apply(*rhoT, v_, urhstmp_, *tmp1);   // [mu u_i s^{ij}],j
 *dudt[ENERGY] -= *tmp1;                           // -= [mu u^i s^{ij}],j

  gadvect_->apply(*p, v_, urhstmp_, *tmp1);         // v.Grad p 
 *dudt[ENERGY] -= *tmp1;                            // -= v . Grad p

  GMTK::dot<Ftype>(fv_, v_, *tmp2, *tmp1);
 *dudt[ENERGY] += *tmp1;                            // += f_kinetic . v

  if ( uf[ENERGY] != NULLPTR ) {                    
cout << "dudt_dry: uf[ENERGY]" << endl;
    *tmp1 = *uf[ENERGY]; *tmp1 *= *Mass;
    GMTK::saxpby<Ftype>(*dudt[ENERGY], 1.0, *tmp1, -1.0); 
                                                    // -= q_heat
  }

  // *************************************************************
  // Momentum equations RHS:
  // *************************************************************
  dd   = urhstmp_[urhstmp_.size()-6]; // holds density fluctuation
 *dd   = (*rhoT); 
  if ( traits_.usebase ) {
   *dd -= *ubase_[0];                 // density fluctuation
   *p  -= *ubase_[1];                 // pressure fluctuation
  }
  for ( auto j=0; j<v_.size(); j++ ) { // for each component

    gdiv_->apply(*s_[j], v_, urhstmp_, *dudt[j]); 


#if 0
    grid_->wderiv(*p, j+1, TRUE, *tmp2, *tmp1);       // Grad p'
   *dudt[j] -= *tmp1;                                 // -= Grad p'
#else
    grid_->deriv(*p, j+1, *tmp2, *tmp1);              // Grad p'
   *tmp1 *= *Mass;                                    // M Grad p' 
   *dudt[j] += *tmp1;                                 // += Grad p'
#endif

    gstressen_->apply(*rhoT, v_, j+1, urhstmp_, 
                                         *tmp1);      // [mu s^{ij}],j
   *dudt[j] -= *tmp1;                                 // -= [mu s^{ij}],j

    if ( traits_.docoriolis ) {
      GMTK::cross_prod_s(traits_.omega, s_, j+1, *tmp1);
     *tmp1 *= *Mass;             
      GMTK::saxpby<Ftype>(*dudt[j], 1.0, *tmp1, 2.0);  // += 2 Omega X (rhoT v) M J
    }

    if ( traits_.dograv || traits_.usebase ) {
     *tmp1 = -GG; 
      compute_vpref(*tmp1, j+1, *tmp2);               // compute grav component
      tmp2->pointProd(*dd, *tmp2);
     *tmp2 *= *Mass;             
     *dudt[j] -= *tmp2;                               // -= rho' vec{g} M J
    }

    if ( traits_.bforced && uf[j] != NULLPTR ) {                    
cout << "dudt_dry: uf[MOMENTUM]" << endl;
      *tmp1 = *uf[j]; *tmp1 *= *Mass;
      GMTK::saxpby<Ftype>(*dudt[j], 1.0, *tmp1, -1.0); 
                                                      // -= f_v
    }

  } // end, momentum loop

  // Multiply RHS by -M^-1 to (1) place all terms on the 'RHS',
  // and (2) to account for factor of M on du/dt term:
  for ( auto j=0; j<nevolve_; j++ ) {
    dudt[j]->apointProd(-1.0, *gimass_->data());// dudt -> -M^-1 dudt
  }

  istage_++;
  
} // end of method dudt_dry


//**********************************************************************************
//**********************************************************************************
// METHOD : dudt_wet
// DESC   : Compute RHS for explicit schemes
// ARGS   : t   : time
//          u   : state. 
//          uf  : forcing tendency state vector
//          ub  : bdy state vector
//          dt  : time step
//          dudt: accelerations, returned
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::dudt_wet(const Time &t, const State &u, const State &uf,  const State &ub, const Time &dt, Derivative &dudt)
{
  assert(!traits_.bconserved &&
         "conservation not yet supported"); 

  GString    serr = "GMConv<TypePack>::dudt_wet: ";
  GINT       nice, nliq ;
  StateComp *irhoT, *Ltot;
  StateComp *dd, *e, *p, *rhoT, *T; // energy, den, pressure, temperature
  StateComp *Mass=grid_->massop().data();
  StateComp *tmp1, *tmp2;
  State      g(nc_); 

  // NOTE:
  // Make sure that, in init(), Helmholtz op is using only
  // weak Laplacian (q * mass term isn't being used), or there will 
  // be problems. This is required for explicit schemes, for
  // which this method is called.

  assert( !traits_.bconserved ); // don't allow conservative form yet

  // Set tmp pool for RHS computations:
  assert(urhstmp_.size() >= szrhstmp());
  Ltot  = urhstmp_[urhstmp_.size()-1];
  rhoT  = urhstmp_[urhstmp_.size()-2];
  irhoT = urhstmp_[urhstmp_.size()-3];
  tmp1  = urhstmp_[urhstmp_.size()-4];
  tmp2  = urhstmp_[urhstmp_.size()-5];

#if 0
  // Do filtering, if any:
  FilterList *fl = &this->get_filter_list();
  for ( auto j=0; j<fl->size(); j++ ) {
    if ( (*fl)[j] != NULLPTR ) {
      (*fl)[j]->apply(t, *u[j], urhstmp_);
    }
  }
#endif


  // Get total density and inverse: *rhoT  = *u[DENSITY]; 
 *rhoT = *u[DENSITY];
  if ( traits_.usebase ) *rhoT +=  *ubase_[0];   
 *irhoT = *rhoT; irhoT->rpow(-1.0);

  // Compute velocity for timestep:
  compute_v(s_, *irhoT, v_); // stored in v_
  
  // Compute all operators as though they are on the LHS, then
  // change the sign and add Mass at the end....


  // *************************************************************
  // Total density RHS:
  // *************************************************************

  gdiv_->apply(*rhoT, v_, urhstmp_, *dudt[DENSITY]); 

  if ( traits_.dofallout ) {
    compute_falloutsrc(*rhoT, qi_, tvi_, -1, urhstmp_, *Ltot);
    GMTK::saxpby<Ftype>(*dudt[DENSITY], 1.0, *Ltot, 1.0);   // += Ltot
  }
  if ( uf[DENSITY] != NULLPTR ) *dudt[DENSITY] -= *uf[DENSITY];//  += sdot(s_rhoT)
  
  // *************************************************************
  // Mass fraction equations (vapor + all hyrodmeteors) RHS:
  // We solve
  //  dq_i/dt + u.Grad q_i + -div(q_i rhoT W_i)/rhoT 
  //                       - q_i/rhoT Ltot + dot(s_i)/rhoT = 0
  // where Ltot is total fallout source over liq + ice sectors:
  // *************************************************************
  for ( auto j=0; j<nmoist_; j++ ) {
    gadvect_->apply(*qi_[j], v_, urhstmp_, *dudt[VAPOR+j]); // apply advection
    compute_vpref(*tvi_[j], W_);
   *tmp1 = (*qi_[j]) * (*rhoT);                 // q_i rhoT
    gdiv_->apply(*tmp1, W_, urhstmp_, *tmp2); 
   *tmp2 *= *irhoT;                             // Div(q_i rhoT W)/rhoT
   *dudt[VAPOR+j] = *tmp2;                      // += Div(q_i rhoT W)/rhoT
    *tmp1  = (*Ltot) * (*irhoT); *tmp1 *= (*qi_[j]);// q_i/rhoT Ltot
    *dudt[VAPOR+j] -= *tmp1;                    // += -q_i/rhoT Ltot
    if ( uf[VAPOR+j] != NULLPTR ) {             // add in sdot(s_i)/rhoT
      *tmp1 = *uf[VAPOR+1]; *tmp1 *= *irhoT;    // dot(s)/rhoT 
      *tmp1 *= *Mass;
      GMTK::saxpby<Ftype>(*dudt[VAPOR+j], 1.0, *tmp1, -1.0); 
                                               // -= dot(s)/rhoT
    }
  }
 
  p     = urhstmp_[urhstmp_.size()-3]; // holds pressure
  T     = urhstmp_[urhstmp_.size()-6]; // holds temperature
  e     = u[ENERGY];                   // internal energy density

  // *************************************************************
  // Energy equation RHS:
  // *************************************************************
  compute_cv(u, *tmp1, *tmp2);                     // Cv
  geoflow::compute_temp<Ftype>(*e, *rhoT, *tmp2, *T);     // temperature
  compute_qd  (u, *tmp1);                          // dry mass ratio
  geoflow::compute_p<Ftype>(*T, *rhoT, *tmp1, RD, *p);    // partial pressure for dry air
  if ( !traits_.dodry ) {
    geoflow::compute_p<Ftype>(*T, *rhoT, *u[VAPOR], RV, *tmp1); // partial pressure for vapor
   *p += *tmp1;
  }

  GMTK::saxpby<Ftype>(*tmp1, *e, 1.0, *p, 1.0);     // h = p+e, enthalpy density

  gdiv_->apply(*tmp1, v_, urhstmp_, *dudt[ENERGY]); 

  gstressen_->apply(*rhoT, v_, urhstmp_, *tmp1);   // [mu u_i s^{ij}],j
 *dudt[ENERGY] -= *tmp1;                           // -= [mu u^i s^{ij}],j

  if ( traits_.dofallout || !traits_.dodry ) {
    GMTK::paxy(*tmp1, *rhoT, CVL, *T);             // tmp1 = C_liq rhoT T
    compute_falloutsrc(*tmp1, qliq_, tvliq_, -1.0, urhstmp_, *Ltot);
                                                   // liquid fallout src
   *dudt[ENERGY] += *Ltot;                         // += L_liq
    GMTK::paxy(*tmp1, *rhoT, CVI, *T);             // tmp1 = C_ice rhoT T
    compute_falloutsrc(*tmp1, qice_, tvice_, -1.0, urhstmp_, *Ltot); 
                                                   // ice fallout src
   *dudt[ENERGY] += *Ltot;                         // += L_ice
  }

  gadvect_->apply(*p, v_, urhstmp_, *tmp1);         // v.Grad p 
 *dudt[ENERGY] -= *tmp1;                            // -= v . Grad p

  if ( traits_.dograv && traits_.dofallout ) {
    compute_pe(*rhoT, qi_, tvi_, urhstmp_, *tmp1);
   *dudt[ENERGY] += *tmp1;                          // += Sum_i rhoT q_i g.W_i
  }

  GMTK::dot<Ftype>(fv_, v_, *tmp2, *tmp1);
 *dudt[ENERGY] += *tmp1;                            // += f_kinetic . v

  if ( uf[ENERGY] != NULLPTR ) {                    
    *tmp1 = *uf[ENERGY]; *tmp1 *= *Mass;
    GMTK::saxpby<Ftype>(*dudt[ENERGY], 1.0, *tmp1, -1.0); 
                                                    // -= q_heat
  }

  // *************************************************************
  // Momentum equations RHS:
  // *************************************************************
  dd   = urhstmp_[urhstmp_.size()-6]; // holds density fluctuation
 *dd   = (*rhoT); 
  if ( traits_.usebase ) {
   *dd -= *ubase_[0];                 // density fluctuation
   *p  -= *ubase_[1];                 // pressure fluctuation
  }
  for ( auto j=0; j<v_.size(); j++ ) { // for each component

    gdiv_->apply(*s_[j], v_, urhstmp_, *dudt[j]); 

    if ( traits_.dofallout || !traits_.dodry ) {
      compute_falloutsrc(*u[j], qliq_, tvi_,-1.0, urhstmp_, *Ltot);
                                                      // hydrometeor fallout src
     *dudt[j] += *Ltot;                               // += L_tot
    }

#if 0
    grid_->wderiv(*p, j+1, TRUE, *tmp2, *tmp1);       // Grad p'
   *dudt[j] -= *tmp1;                                 // -= Grad p'
#else
    grid_->deriv(*p, j+1, *tmp2, *tmp1);              // Grad p'
   *tmp1 *= *Mass;                                    // M Grad p' 
   *dudt[j] += *tmp1;                                 // += Grad p'
#endif

    gstressen_->apply(*rhoT, v_, j+1, urhstmp_, 
                                         *tmp1);      // [mu s^{ij}],j
   *dudt[j] -= *tmp1;                                 // -= [mu s^{ij}],j

    if ( traits_.docoriolis ) {
      GMTK::cross_prod_s(traits_.omega, s_, j+1, *tmp1);
     *tmp1 *= *Mass;             
      GMTK::saxpby<Ftype>(*dudt[j], 1.0, *tmp1, 2.0);  // += 2 Omega X (rhoT v) M J
    }

    if ( traits_.dograv || traits_.usebase ) {
     *tmp1 = -GG; 
      compute_vpref(*tmp1, j+1, *tmp2);               // compute grav component
      tmp2->pointProd(*dd, *tmp2);
     *tmp2 *= *Mass;             
     *dudt[j] -= *tmp2;                               // -= rho' vec{g} M J
    }

    if ( traits_.bforced && uf[j] != NULLPTR ) {                    
      *tmp1 = *uf[j]; *tmp1 *= *Mass;
      GMTK::saxpby<Ftype>(*dudt[j], 1.0, *tmp1, -1.0); 
                                                      // -= f_v
    }

  } // end, momentum loop

  // Multiply RHS by -M^-1 to (1) place all terms on the 'RHS',
  // and (2) to account for factor of M on du/dt term:
  for ( auto j=0; j<nevolve_; j++ ) {
    dudt[j]->apointProd(-1.0, *gimass_->data());// dudt -> -M^-1 dudt
  }

  istage_++;
  
} // end of method dudt_wet


//**********************************************************************************
//**********************************************************************************
// METHOD : step_impl (1)
// DESC   : Step implementation method entry point
// ARGS   : t   : time
//          uin : input state, modified on output with update
//          uf  : force-tendency vector
//          ub  : bdy vector
//          dt  : time step
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::step_impl(const Time &t, State &uin, State &uf, State &ub, const Time &dt)
{

  GBOOL bret;

  assert(bInit_);

  // If there's a top-of-the-timestep callback, 
  // call it here:
  if ( bsteptop_ ) {
    steptop_callback_(t, uin, dt);
  }

  // Set bdy conditions:
  apply_bc_impl(t, uin, ub);

  // Set evolved state vars from input state.
  // These are not deep copies:
  for ( auto j=0; j<traits_.nsolve; j++ ) uevolve_ [j] = uin[j];
  

  switch ( traits_.isteptype ) {
    case GSTEPPER_EXRK:
      // Assign qi, tvi, qice, qliq, tvice, tvliq:
      assign_helpers(uin, uf);
      istage_ = 0;
      for ( auto j=0; j<uold_.size(); j++ ) *uold_[j] = *uevolve_[j];
      step_exrk(t, uold_, uf, ub, dt, uevolve_);
      break;
    case GSTEPPER_BDFAB:
    case GSTEPPER_BDFEXT:
      step_multistep(t, uevolve_, uf, ub, dt);
      break;
  }

  // Do filtering, if any:
  FilterList *fl = &this->get_filter_list();
  for ( auto j=0; j<fl->size(); j++ ) {
    if ( (*fl)[j] != NULLPTR ) {
      (*fl)[j]->apply(t, *uevolve_[j], urhstmp_);
    }
  }

  apply_bc_impl(t, uin, ub);

  // Check solution for NaN and Inf:
  bret = TRUE;
  for ( auto j=0; j<uevolve_.size() && bret; j++ ) {
//  bret = bret && uevolve_ [j]->isfinite();
  }
  assert(bret && "Solution not finite");

  icycle_++;

} // end of method step_impl (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : step_impl (2)
// DESC   : Step implementation method entry point
// ARGS   : t   : time
//          uin : input state, modified on output with update
//          ub  : bdy vector
//          dt  : time step
//          uout: output state
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::step_impl(const Time &t, const State &uin, State &uf,  State &ub, const Time &dt, State &uout)
{
  assert(FALSE && "step_impl(2) not available");

} // end of method step_impl (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : step_multistep
// DESC   : Carries out multistep update. The time derivative and 
//          advection terms are handlex using a multistep expansion. The
//          advection term is 'extrapolated' to the new time level
//          using known state data so as to obviate the need for a fully 
//          implicit treatment. The dissipation term is handled implicitly.
// ARGS   : t   : time
//          u   : state
//          uf  : force tendency vector
//          ub  : bdy vector
//          dt  : time step
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::step_multistep(const Time &t, State &uin, State &uf, State &ub, const Time &dt)
{
  assert(FALSE && "Multistep methods not yet available");
  
} // end of method step_multistep


//**********************************************************************************
//**********************************************************************************
// METHOD : step_exrk
// DESC   : Take a step using Explicit RK method
// ARGS   : t   : time
//          uin : input state; must not be modified
//          uf  : force tendency vector
//          ub  : bdy vector
//          dt  : time step
//          uout: output/updated state
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::step_exrk(const Time &t, State &uin, State &uf, State &ub, const Time &dt, State &uout)
{
  assert(gexrk_ != NULLPTR && "GExRK operator not instantiated");

  // GExRK stepper steps entire state over one dt:
  gexrk_->step(t, uin, uf, ub, dt, urktmp_, uout);

//GMTK::constrain2sphere(*grid_, uout);

} // end of method step_exrk


//**********************************************************************************
//**********************************************************************************
// METHOD : init_impl
// DESC   : Initialize equation object
// ARGS   : tmp : Array of tmp vector pointers, pointing to vectors
//                of same size as State. Must be MAX(2*DIM+2,iorder+1)
//                vectors
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::init_impl(State &u, State &tmp)
{
  GString serr = "GMConv<TypePack>::init: ";

  GBOOL      bmultilevel = FALSE;
  GSIZET     n;
  GINT       nexcl, nrhstmp;
  GGridIcos *icos = dynamic_cast<GGridIcos*>(grid_);
  CompDesc  *icomptype = &this->stateinfo().icomptype;
 
  assert(tmp.size() >= this->tmp_size());

  v_.resize(nc_); v_ = NULLPTR;
  W_.resize(nc_); W_ = NULLPTR;

  
  // Set space for state velocity, if needed:
  nexcl = 0;
  if ( traits_.usemomden ) {
    for ( auto j=0; j<v_.size(); j++ ) {
      v_[j] = tmp[j];
      nexcl++;
    }
  }

  // Set space for terminal velocity, if needed:
  if ( icos && (traits_.nlsector || traits_.nisector)) {
    for ( auto j=0; j<W_.size(); j++ ) {
      W_[j] = tmp[j+v_.size()];
      nexcl++;
    }
  }

  // Set up tmp pool:
  utmp_.resize(tmp.size()-nexcl); 
  for ( auto j=0; j<utmp_.size(); j++ ) {
    utmp_[j] = tmp[j+nexcl];
  }

  // Check if specified stepper type is valid:
  GBOOL  bfound;
  GSIZET itype;
  valid_types_.push_back("GSTEPPER_EXRK");
/*
  valid_types_.push_back("GSTEPPER_BDFAB");
  valid_types_.push_back("GSTEPPER_BDFEXT");
*/
  bfound = valid_types_.contains(traits_.ssteptype, itype);
  assert( bfound && "Invalid stepping method specified");
  traits_.isteptype = static_cast<GStepperType>(itype);

  if ( traits_.dodry ) traits_.dofallout = FALSE;

  // Set std::vector for traits.iforced:
  stdiforced_.resize(traits_.iforced.size());
  for ( auto j=0; j<stdiforced_.size(); j++ ) {
    stdiforced_[j] = traits_.iforced[j];
  }

  // Set dissipation from traits. Note that
  // this is set to be constant, based on configuration,
  // even though the solver can accommodate spatially
  // variable dissipation:
  nu_  .resize(1);
  nu_  = traits_.nu;
  eta_ .resize(1);
  eta_ = traits_.eta;

  // Find no. state and solve members, and component types:
  for( auto j=0; j<v_.size(); j++ ) icomptype->push_back(GSC_KINETIC); 
  icomptype->push_back(GSC_DENSITYT); 
  icomptype->push_back(GSC_ENERGY); 
  if ( !traits_.dodry ) {
    n = traits_.nlsector + traits_.nisector;
    for( auto j=0; j<traits_.nlsector; j++ ) icomptype->push_back(GSC_MASSFRAC); 
    for( auto j=0; j<traits_.nisector; j++ ) icomptype->push_back(GSC_MASSFRAC); 
  }
  if ( traits_.usebase ) { // base state components
    for( auto j=0; j<2; j++ ) icomptype->push_back(GSC_PRESCRIBED); 
  }
  if ( !traits_.dodry && traits_.dofallout ) { // include terminal velocities:
    for( auto j=0; j<n; j++ ) icomptype->push_back(GSC_PRESCRIBED); 
  }


  // Specify starting indices for each sector
  // (negative if invalid):
  MOMENTUM   = 0;
  ENERGY     = nc_;
  DENSITY    = ENERGY+1;
  VAPOR      = traits_.dodry ? -1 : DENSITY+1;
  LIQMASS    = traits_.dodry ? -1 : VAPOR+1;
  ICEMASS    = traits_.dodry ? -1 : LIQMASS+traits_.nlsector;
  PRESCRIBED = traits_.dodry ? nc_+2 : ICEMASS+traits_.nisector;
  BASESTATE  = PRESCRIBED;
  LIQTERMV   = traits_.dodry ? -1 : BASESTATE+traits_.nbase;
  ICETERMV   = traits_.dodry ? -1 : LIQTERMV+traits_.nlsector;

  

  // Find multistep/multistage time stepping coefficients:
  GMultilevel_coeffs_base<Ftype> *tcoeff_obj=NULLPTR; // time deriv coeffs
  GMultilevel_coeffs_base<Ftype> *acoeff_obj=NULLPTR; // adv op. coeffs


  std::function<void(const Time &t,                    // RHS callback function
                     const State  &uin,
                     const State  &uf,
                     const State  &ub,
                     const Time   &dt,
                     State &dudt)> rhs
                  = [this](const Time &t,           
                     const State  &uin, 
                     const State  &uf, 
                     const State  &ub, 
                     const Time   &dt,
                     State &dudt){dudt_impl(t, uin, uf, ub, dt, dudt);}; 

  std::function<void(const Time &t,                    // Bdy apply callback function
                     State  &uin, 
                     State &ub)> applybc 
                  = [this](const Time &t,              
                     State  &uin, 
                     State &ub){apply_bc_impl(t, uin, ub);}; 

  // Configure time stepping:
  GExRKStepper<GFTYPE>::Traits rktraits;
  switch ( traits_.isteptype ) {
    case GSTEPPER_EXRK:
      rktraits.bSSP   = traits_.bSSP;
      rktraits.norder = traits_.itorder;
      rktraits.nstage = traits_.nstage;
      gexrk_ = new GExRKStepper<Ftype>(rktraits, *grid_);
      gexrk_->setRHSfunction(rhs);
      gexrk_->set_apply_bdy_callback(applybc);
      gexrk_->set_ggfx(ggfx_);
      // Set 'helper' tmp arrays from main pool, utmp_, so that
      // we're sure there's no overlap:
      uold_   .resize(traits_.nsolve); // RK-solution at time level n
      uevolve_.resize(traits_.nsolve); // current RK solution
      ubase_  .resize(traits_.nbase);  // points to base-state components
      urktmp_ .resize(traits_.nsolve*(traits_.itorder+1)+1); // RK stepping work space
      urhstmp_.resize(szrhstmp());     // work space for RHS
      nrhstmp = utmp_.size()-uold_.size()-urktmp_.size();

      assert(nrhstmp >= szrhstmp() && "Invalid rhstmp array size");
      // Make sure there is no overlap between tmp arrays:
      n = 0;
      for ( GSIZET j=0; j<traits_.nsolve ; j++, n++ ) uold_   [j] = utmp_[n];
      for ( GSIZET j=0; j<urktmp_ .size(); j++, n++ ) urktmp_ [j] = utmp_[n];
      for ( GSIZET j=0; j<urhstmp_.size(); j++, n++ ) urhstmp_[j] = utmp_[n];
      break;
/*
    case GSTEPPER_BDFAB:
      dthist_.resize(MAX(traits_.itorder,traits_.inorder));
      tcoeff_obj = new G_BDF<Ftype>(traits_.itorder, dthist_);
      acoeff_obj = new G_AB<Ftype> (traits_.inorder, dthist_);
      tcoeffs_.resize(tcoeff_obj->getCoeffs().size());
      acoeffs_.resize(acoeff_obj->getCoeffs().size());
      tcoeffs_ = tcoeff_obj->getCoeffs(); 
      acoeffs_ = acoeff_obj->getCoeffs();
      uold_   .resize(traits_.nsolve); // solution at time level n
      for ( GSIZET j=0; j<traits_.nsolve; j++ ) uold_[j] = utmp_[j];
      bmultilevel = TRUE;
      break;
    case GSTEPPER_BDFEXT:
      dthist_.resize(MAX(traits_.itorder,traits_.inorder));
      tcoeff_obj = new G_BDF<Ftype>(traits_.itorder, dthist_);
      acoeff_obj = new G_EXT<Ftype>(traits_.inorder, dthist_);
      tcoeffs_.resize(tcoeff_obj->getCoeffs().size());
      acoeffs_.resize(acoeff_obj->getCoeffs().size());
      tcoeffs_ = tcoeff_obj->getCoeffs(); 
      acoeffs_ = acoeff_obj->getCoeffs();
      urhstmp_.resize(utmp_.size()-urktmp_.size());
      for ( GSIZET j=0; j<utmp_.size(); j++ ) urhstmp_[j] = utmp_[j];
      bmultilevel = TRUE;
      break;
*/
    default:
      assert(FALSE && "Invalid stepper type");
  }
  if ( tcoeff_obj != NULLPTR ) delete tcoeff_obj;
  if ( acoeff_obj != NULLPTR ) delete acoeff_obj;
  
  // Instantiate spatial discretization operators:
  gmass_      = &grid_->massop();
//ghelm_      = new GHelmholtz(*grid_);
  
  typename GStressEnOp<TypePack>::Traits trstress;
  trstress.type       = GStressEnOp<TypePack>::GSTRESS_REDUCED;
  trstress.full_colloc= TRUE;
  trstress.Stokes_hyp = traits_.Stokeshyp;
  trstress.indep_diss = traits_.bindepdiss;
  trstress.nu    .resize(nu_  .size());  trstress.nu    = nu_;
  trstress.zeta  .resize(nu_  .size());  trstress.zeta  = traits_.zeta; 
  trstress.eta   .resize(eta_ .size());  trstress.eta   = eta_;
  trstress.lambda.resize(eta_ .size());  trstress.lambda= traits_.lambda;
  gstressen_  = new GStressEnOp<TypePack>(trstress, *grid_);

  if ( traits_.isteptype ==  GSTEPPER_EXRK ) {
    gimass_ = &grid_->imassop();
  }

  // If doing semi-implicit time stepping; handle viscous term 
  // (linear) inplicitly, which implies using full Helmholtz operator:
  if ( traits_.isteptype == GSTEPPER_BDFAB || traits_.isteptype == GSTEPPER_BDFEXT ) {
    assert(FALSE && "Implicit time stepping not yet supported");
  }

  typename GDivOp<TypePack>::Traits trgdiv;
  trgdiv.docollocation = traits_.divopcolloc;
  grid_->set_usebdydata(traits_.usebdydata); // Set usebdydata flag in grid
  if ( traits_.bconserved ) {
    assert(FALSE && "Conservation not yet supported");
//  gpdv_  = new GpdV<TypePack>(*grid_);
    gdiv_  = new GDivOp<TypePack>(trgdiv, *grid_);
    assert( (gmass_     != NULLPTR
//        && ghelm_     != NULLPTR
          && gstressen_ != NULLPTR
//        && gpdv_      != NULLPTR
          && gdiv_      != NULLPTR ) && "1 or more operators undefined");
  }
  else {
    gadvect_ = new GAdvect<TypePack>(*grid_);
//  gpdv_    = new GpdV<TypePack>(*grid_);
    gdiv_    = new GDivOp<TypePack>(trgdiv, *grid_);
    assert( (gmass_     != NULLPTR
//        && ghelm_     != NULLPTR
          && gstressen_ != NULLPTR
//        && gpdv_      != NULLPTR
          && gdiv_      != NULLPTR
          && gadvect_   != NULLPTR) && "1 or more operators undefined");
  }

  // If doing a multi-step method, instantiate (deep) space for 
  // required time levels for state:
  if ( bmultilevel ) {
    ukeep_ .resize(traits_.itorder);
    for ( auto i=0; i<traits_.itorder-1; i++ ) { // for each time level
      ukeep_[i].resize(traits_.nsolve);
      for ( auto j=0; j<traits_.nsolve; j++ ) ukeep_[i][j] = new GTVector<Ftype>(grid_->ndof());
    }
  }

  // Find minimum element edge/face lengths for 
  // timestep computation:
  maxbyelem_.resize(grid_->nelems());
  grid_->minlength(&dxmin_);

  // Set size of mass frac and 
  // misc. helper arrays:
  nhydro_ = traits_.dodry ? 0 : traits_.nlsector + traits_.nisector + 1;
  nmoist_ = traits_.dodry ? 0 : nhydro_ + 1;
  nevolve_ = nc_ + 2 + nmoist_;
  this->stateinfo().nevolve = traits_.nsolve;
  this->stateinfo().npresc  = traits_.nstate - traits_.nsolve;
  qi_   .resize(nmoist_);
  tvi_  .resize(nmoist_);
  qliq_ .resize(traits_.nlsector);
  tvliq_.resize(traits_.nlsector);
  qice_ .resize(traits_.nisector);
  tvice_.resize(traits_.nisector);
  fv_   .resize(nc_); 
  s_   . resize(nc_); 

  qi_   = NULLPTR;
  tvi_  = NULLPTR;
  qliq_ = NULLPTR;
  tvliq_= NULLPTR;
  qice_ = NULLPTR;
  tvice_= NULLPTR;
  fv_   = NULLPTR;
  s_    = NULLPTR;


  compute_base(u);

  bInit_ = TRUE;

} // end of method init_impl


//**********************************************************************************
//**********************************************************************************
// METHOD : cycle_keep
// DESC   : Cycle the mult-level states making sure the most
//          recent is at index 0, the next most recent, at index 1, etc...
// ARGS   : u     : State variable providing most recent state
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::cycle_keep(const State &u)
{

  // Make sure following index map contains the 
  // correct time level information:
  //   ukeep[0] <--> time level n (most recent)
  //   ukeep[1] <--> time level n-1
  //   ukeep[2] <--> time level n-2 ...
  ukeep_ .resize(traits_.itorder);
  for ( auto i=traits_.itorder-1; i>=1; i-- ) ukeep_[i] = ukeep_[i+1];
  ukeep_[0] = u;

} // end of method cycle_keep


#if 0
//**********************************************************************************
//**********************************************************************************
// METHOD : set_nu
// DESC   : Set viscosity quantity. This may be a field, or effectively a
//          scalar (where only element 0 contains valid data). If not set, 
//          Helmholtz op creates a default of nu = 1.
//          nu : viscosity vector (may be of length 1). Will be managed
//               by caller; only a pointer is used by internal operators.
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::set_nu(GTVector<Ftype> &nu)
{
  assert(ghelm_ != NULLPTR && "Init must be called first");
  nu_ = nu; // Not sure this class actually needs this. May be removed later
  ghelm_->set_Lap_scalar(nu_);

} // end of method set_nu

#endif


//**********************************************************************************
//**********************************************************************************
// METHOD : apply_bc_impl
// DESC   : Apply global domain boundary conditions, ub
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::apply_bc_impl(const Time &t, State &u, State &ub)
{
  Time ttime = t;

  BdyUpdateList *updatelist = &grid_->bdy_update_list();;


  // Update bdy values if required to:
  for ( auto k=0; k<updatelist->size(); k++ ) { // foreach grid bdy
    for ( auto j=0; j<(*updatelist)[j].size(); j++ ) { // each update method
      (*updatelist)[k][j]->update(*grid_, this->stateinfo(), ttime, utmp_, u, ub);
    }
  }
} // end of method apply_bc_impl


//**********************************************************************************
//**********************************************************************************
// METHOD : compute_cv
// DESC   : Compute total specific heat at const vol, Cv
//             Cv = Cvd qd + Cvv qv + Sum_i(Cl_i ql_i) + Sum_j(Ci_j qi_j).
//          where ql are the liquid mass fractions, and qi are the ice
//          mass fractions. 
// ARGS   : u    : state
//          utmp : tmp vector
//          cv   : cv field
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::compute_cv(const State &u, StateComp &utmp, StateComp &cv)
{
   GString    serr = "GMConv<TypePack>::compute_cv: ";
   GINT       ibeg;

   // Compute Cv:
   if ( traits_.dodry ) { // if dry dynamics only
     cv = CVD;  
     return;
   }

   utmp     = 1.0;              // running total: subtract qi's  
   utmp    -= (*u[VAPOR]);      // -q_v
   cv       = (*u[VAPOR])*CVV;  // Cv = Cvv * q_vapor
   ibeg     = LIQMASS;
   for ( auto k=ibeg; k<ibeg+traits_.nlsector+1; k++ ) { // liquids
      utmp    -= *u[k];         // subtract in ql_k
      GMTK::saxpby<Ftype>(cv, 0.0, *u[k], CVL); // add in Cvl * ql_k
   }
   ibeg = LIQMASS + traits_.nlsector;;
   for ( auto k=ibeg; k<ibeg+traits_.nisector; k++ ) { // ices
     utmp    -= *u[k];         // subtract in qi_k
     cv       += (*u[k]) * CVI; // add in Cvi * qi_k
     GMTK::saxpby<Ftype>(cv, 0.0, *u[k], CVI); // add in Cvi * qi_k
   }
   // After subtracting q_i, final result is qd=q_dry, so:
   GMTK::saxpby<Ftype>(cv, 0.0, utmp, CVD); // Final Cv += Cvd * qd

} // end of method compute_cv


//**********************************************************************************
//**********************************************************************************
// METHOD : compute_qd
// DESC   : Compute dry mass fraction from other mass fractions:
//             Cv = 1 - Sum_i q_i
// ARGS   : u    : state
//          qd   : dry mass fraction field
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::compute_qd(const State &u, StateComp &qd)
{
   GString    serr = "GMConv<TypePack>::compute_qd: ";
   GINT       ibeg;

   // Compute qd:
   if ( traits_.dodry ) { // if dry dynamics only
     qd = 1.0;  
     return;
   }

   qd  = 1.0;  
   qd -= (*u[VAPOR]);      // running total 1- Sum_k q_k
   ibeg = LIQMASS;
   for ( auto k=ibeg; k<ibeg+traits_.nlsector+1; k++ ) { // liquids
     qd -= *u[k];         // subtract in ql_k
   }
   ibeg = LIQMASS + traits_.nlsector;;
   for ( auto k=ibeg; k<ibeg+traits_.nisector; k++ ) { // ices
     qd  -= *u[k];         // subtract in qi_k
   }

} // end of method compute_qd


//**********************************************************************************
//**********************************************************************************
// METHOD : compute_v
// DESC   : Compute velocity from momentum density in state vector.
//             v_i = s_i/ d
//          where v_i is the member data array, d  is (total) density 
// ARGS   : u    : state
//          id   : 1/denstiy
//          v    : velocity state; components may change on exit
// RETURNS: none. Member data, v_, is set here
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::compute_v(const State &u, StateComp &id, State &v)
{
   GString    serr = "GMConv<TypePack>::compute_v: ";


   if ( !traits_.usemomden ) {
     // State uses velocity form already so 
     // do pointer assignment:
     for ( auto j=0; j<v.size(); j++ ) { // v_i = u_i
       v[j] = u[j];
     }
     return;
   }

   // Find velocity from momentum density:
   
   for ( auto j=0; j<v.size(); j++ ) {
     *v[j]  = *u[j];  // deep copy
     *v[j] *= id[j];  // divide by density
   }

} // end of method compute_v


//**********************************************************************************
//**********************************************************************************
// METHOD : compute_vpref (1)
// DESC   : Transform preferred direction component into vector.
//          Components of W may be NULL.
// ARGS   : tvi   : terminal velocity or gravity magnitude in preferred
//                  direction
//          W     : Cartesian vector field 
// RETURNS: none. 
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::compute_vpref(StateComp &tvi, State &W)
{
   GString    serr = "GMConv<TypePack>::compute_vpref(1): ";
   Ftype      r, x, y, z;
   Ftype      lat, lon;
   GTVector<GTVector<Ftype>> 
             *xnodes = &grid_->xNodes();

   GGridIcos *icos = dynamic_cast<GGridIcos*>(grid_);
   GGridBox  *box  = dynamic_cast <GGridBox*>(grid_);

   assert(W.size() == nc_);

   // Specify components for Cartesian grid:
   if ( box != NULLPTR ) {
     // In Cartesian coords, select the 'z' direction
     // as preferred 'fallout' direction. In 2d, this
     // will be the 2-coord; in 3d, the 3-coord:
    *W[nc_-1] = tvi; 
     return;
   }

   // Specify components for spherical grid. Here,
   // the preferred direction is along radial direction; 
   // we assume that the prescribed term velocity
   // is the radial component only, which must be
   // transformed into Cartesian components:

   for ( auto i=0; i<(*xnodes)[0].size(); i++ ) {
     x = (*xnodes)[0][i]; y = (*xnodes)[1][i]; z = (*xnodes)[2][i];
     r     = sqrt(x*x + y*y + z*z);
     lat   = asin(z/r);
     lon   = atan2(y,x);

     (*W[0])[i] = cos(lat)*cos(lon);
     (*W[1])[i] = cos(lat)*sin(lon);
     (*W[2])[i] = sin(lat);
   }

   for ( auto j=0; j<W_.size(); j++ ) {
     *W[j] *= tvi;
   }


} // end of method compute_vpref (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : compute_vpref (2)
// DESC   : Transform preferred direction component specified in coord
//          direction idir
// ARGS   : tvi   : terminal velocity or gravity magnitude in preferred
//                  direction
//          idir  : vetor component direction, in 1, 2,...GDIM
//          W     : idir-th Cartesian vector field component
//                  sought.
// RETURNS: none. 
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::compute_vpref(StateComp &tvi, GINT idir, StateComp &W)
{
   GString    serr = "GMConv<TypePack>::compute_vpref (2): ";
   Ftype      r, x, y, z;
   Ftype      lat, lon;
   GTVector<GTVector<Ftype>> 
             *xnodes = &grid_->xNodes();

   GGridIcos *icos = dynamic_cast<GGridIcos*>(grid_);
   GGridBox  *box  = dynamic_cast <GGridBox*>(grid_);

   assert(idir > 0 && idir <= nc_);

   // Specify components for Cartesian grid:
   if ( box != NULLPTR ) {
     // In Cartesian coords, select the 'z' direction
     // as preferred 'fallout' direction. In 2d, this
     // will be the 2-coord; in 3d, the 3-coord:
     if ( idir != nc_ ) {
       W = 0.0;
     }
     else {
       W = tvi; 
     }
     return;
   }

   // Specify components for spherical grid. Here,
   // the preferred direction is along radial direction; 
   // we assume that the prescribed term velocity
   // is the radial component only, which must be
   // transformed into Cartesian components:

   if ( idir == 1 ) {
     for ( auto i=0; i<(*xnodes)[0].size(); i++ ) {
       x = (*xnodes)[0][i]; y = (*xnodes)[1][i]; z = (*xnodes)[2][i];
       r     = sqrt(x*x + y*y + z*z);
       lat   = asin(z/r);
       lon   = atan2(y,x);
       W[i]  = tvi[i]*cos(lat)*cos(lon);
     }
   }
   else if ( idir == 2 ) {
     for ( auto i=0; i<(*xnodes)[0].size(); i++ ) {
       x = (*xnodes)[0][i]; y = (*xnodes)[1][i]; z = (*xnodes)[2][i];
       r     = sqrt(x*x + y*y + z*z);
       lat   = asin(z/r);
       lon   = atan2(y,x);
       W[i]  = tvi[i]*cos(lat)*sin(lon);
     }
   }
   else {
     for ( auto i=0; i<(*xnodes)[0].size(); i++ ) {
       x = (*xnodes)[0][i]; y = (*xnodes)[1][i]; z = (*xnodes)[2][i];
       r     = sqrt(x*x + y*y + z*z);
       lat   = asin(z/r);
       lon   = atan2(y,x);
       W[i]  = tvi[i]*sin(lat);
     }
   }

} // end of method compute_vpref (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : compute_falloutsrc
// DESC   : Compute total effect due to fallout of density quantity, g. 
//            r = Sum_i Div [q_i g vec{W_i}] 
//          where vec{W_i} is the fallout terminal velocity for precipitating
//          species i (from either liquid or ice sectors of state). 
// ARGS   : 
//          g    : quantiy to fall out ('flux-out'). May be: total density,
//                 momentum density component, internal energy density
//          qi   : hydrometeor mass fractions
//          tvi  : terminal velocity vector for each qi hydrometeor, from which
//                 vec{W} is computed
//          jexcl: which hydrometeor index to exclude from sum. If jexcl<0,
//                 exclude none. This index must be the 0-starting index of
//                 the hydrometeor in the qi array.
//          utmp : tmp vectors; at least 5 required
//          r    : fallout src field
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::compute_falloutsrc(StateComp &g, State &qi, State &tvi, GINT jexcl, State &utmp, StateComp &r)
{
   GString     serr = "GMConv<TypePack>::compute_falloutsrc: ";
   GINT        nhydro; // no. hydrometeors
   StateComp  *div, *qg; 

   assert(tvi.size() == qi.size());
   assert(utmp.size() >= 5);
   assert(jexcl < 0 || (jexcl >=0 && jexcl < qi.size()));

   // Compute:
   //    r = -Sum_i Div (rho_t q_i vec{W}_i )

   r = 0.0;
   if ( !traits_.dofallout || traits_.dodry ) return;

   nhydro = qi.size(); // traits_.nlsector + traits_.nisector;

   qg    = utmp[3];    // temp
   div   = utmp[4];    // temp
   for ( auto j=0; j<nhydro; j++ ) {
     if ( j == jexcl ) continue;  

     *qg = g; (*qg) *= (*qi[j]); // compute g q_i

     // Convert terminal velocities to required 
     // (Cartesian) components. vterm may point
     // to tvi, or to v_ member data, or may be NULL:
     compute_vpref(*tvi[j], W_);

     // Compute i_th contribution to source term:
     gdiv_->apply(*qg, W_, utmp, *div); 
     r += *div;
   }

} // end of method compute_falloutsrc


//**********************************************************************************
//**********************************************************************************
// METHOD : compute_pe
// DESC   : Compute potential energy of hydrometeors:
//       
//              r = Sum_i rhoT q_i vec{g} . vec{W}_i,
//          where i is over all hydrometeor sectors, q_i & W_i
//          are mass fraction and terminal velocity of ith hydrometeor. 
// ARGS   : 
//          rhoT : total density
//          qi   : hydrometeor mass fractions
//          tvi  : terminal velocity vector for each qi hydrometeor, from which
//                 vec{W} is computed
//          utmp : tmp vectors; nc_+2 required
//          r    : fallout potential energy
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::compute_pe(StateComp &rhoT, State &qi, State &tvi, State &utmp, StateComp &r)
{
   GString     serr = "GMConv<TypePack>::compute_pe: ";
   GINT        nhydro; // no. hydrometeors
   StateComp   *tmp1, *tmp2;
   State       g(nc_); 

   assert(tvi.size() == qi.size());
   assert(utmp.size() >= nc_+2);

   // Compute:
   //    r = Sum_i rhoT q_i vec{g} . vec{W}_i,

   r = 0.0;
   if ( !traits_.dograv || !traits_.dofallout || traits_.dodry ) return;

   nhydro = qi.size(); // traits_.nlsector + traits_.nisector;

   tmp1 = utmp[nc_];
   tmp2 = utmp[nc_+1];

   // Set tmp space for gravity vector:
   for ( auto j=0; j<nc_; j++ ) g[j] = utmp[j];
  
   // Compute Cartersian components of gravity vector field:
  *tmp1 = -GG;
   compute_vpref(*tmp1, g);

   r = 0.0;
   for ( auto j=0; j<nhydro; j++ ) {

     compute_vpref(*tvi[j], W_);       // compute W
     GMTK::dot<Ftype>(g, W_, *tmp1, *tmp2);  // g.W
     tmp2->pointProd(*qi[j]);          // qi * (g.W)
     r += *tmp2;

   }
   r.pointProd(rhoT);                  // rhoT Sum_i q_i vec{g}.vec{W}_i

} // end of method compute_pe


//**********************************************************************************
//**********************************************************************************
// METHOD : compute_base
// DESC   : Compute base state. Is grid-dependent.
// ARGS   : u : state vector
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::compute_base(State &u)
{
   GString    serr = "GMConv<TypePack>::compute_base: ";
   Ftype      x, y, z;
   Ftype      r, T;
   StateComp  dtmp, ptmp;
   GTVector<GTVector<Ftype>> 
             *xnodes = &grid_->xNodes();

   if ( !traits_.usebase ) return;

   GGridIcos *icos = dynamic_cast<GGridIcos*>(grid_);
   GGridBox  *box  = dynamic_cast <GGridBox*>(grid_);

   dtmp.resize(grid_->ndof());
   ptmp.resize(grid_->ndof());
   // Specify components for Cartesian grid:
   if ( box != NULLPTR ) {
     // In Cartesian coords, select the 'z' direction
     // as preferred 'fallout' direction. In 2d, this
     // will be the 2-coord; in 3d, the 3-coord:
     for ( auto j=0; j<(*xnodes)[0].size(); j++ ) {
       x = (*xnodes)[0][j]; y = (*xnodes)[1][j];
       if ( nc_ == 3 ) z = (*xnodes)[2][j];
       r         = nc_ == 3 ? z : y;
       T         = traits_.Ts_base - GG*r/CPD;
       ptmp[j]   = traits_.P0_base*pow(T/traits_.Ts_base,CPD/RD);
       dtmp[j]   = ptmp[j] / ( RD * T );
     }
     *u  [BASESTATE] = dtmp;
     *u[BASESTATE+1] = ptmp;
     return;
   }

   // Specify state for spherical grid. Here,
   // the preferred direction is along radial direction:
   for ( auto j=0; j<(*xnodes)[0].size(); j++ ) {
     x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
     r        = sqrt(x*x + y*y + z*z);
     T        = traits_.Ts_base - GG*r/CPD;
     ptmp[j]  = traits_.P0_base*pow(T/traits_.Ts_base,CPD/RD);
     dtmp[j]  = ptmp[j] / ( RD * T );
   }
   *u  [BASESTATE] = dtmp;
   *u[BASESTATE+1] = ptmp;

} // end of method compute_base


//**********************************************************************************
//**********************************************************************************
// METHOD : assign_helpers
// DESC   : Assigns helper arrays to (full) state components 
//       
// ARGS   : u : incoming _full_ state vector
//          uf: incoming _full_ forcing state vector
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::assign_helpers(const State &u, const State &uf)
{
   GString     serr = "GMConv<TypePack>::assign_helpers: ";

   // Point helper arrays to corresp state components:
   qi_  = NULLPTR;
   tvi_ = NULLPTR;
   for ( auto j=0; j<nhydro_; j++ ) { //
     qi_ [j] = u[nc_+3+j]; // set mass frac vector
     tvi_[j] = u[nc_+3+nhydro_+j]; // set term speed for each qi
   }
   for ( auto j=0; j<nc_; j++ ) fv_[j] = uf[j]; // kinetic forcing vector
   for ( auto j=0; j<traits_.nbase; j++ ) ubase_[j] = u[BASESTATE+j]; // base state

   GINT nliq = traits_.nlsector;
   GINT nice = traits_.nisector;
   for ( auto j=0; j<s_.size(); j++ ) s_    [j] = u[j];
   for ( auto j=0; j<nliq; j++ ) qliq_ [j] = u[nc_+3+j];
   for ( auto j=0; j<nice; j++ ) qice_ [j] = u[nc_+3+nliq+j];
   for ( auto j=0; j<nliq; j++ ) tvliq_[j] = u[nc_+3+  nliq+nice+j];
   for ( auto j=0; j<nice; j++ ) tvice_[j] = u[nc_+3+2*nliq+nice+j];

} // end of method assign_helpers


//**********************************************************************************
//**********************************************************************************
// METHOD : szrhstmp
// DESC   : Get size of tmp array for RHS computations
//       
// ARGS   : none
// RETURNS: GINT size
//**********************************************************************************
template<typename TypePack>
GINT GMConv<TypePack>::szrhstmp()
{
  GINT       maxop, sum = 0;
  GGridBox  *box = dynamic_cast <GGridBox*>(grid_);
   
  maxop = MAX(4,nc_+2); // stressen max
  if ( traits_.dofallout && !traits_.dodry ) maxop = MAX(5,maxop);

  sum += maxop;
  sum += 6;        // size for misc tmp space in dudt_impl

  return sum;

} //end szrhstmp


//**********************************************************************************
//**********************************************************************************
// METHOD : tmp_size_impl
// DESC   : Find required tmp size on GLL grid
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
GINT GMConv<TypePack>::tmp_size_impl()
{
  GINT sum = 0;
 
  sum += nc_;                                 // for v_ 
  sum += traits_.nlsector || traits_.nisector ? nc_ : 0;  // for W_
  sum += traits_.nsolve;                     // old state storage
  sum += traits_.nsolve
       * (traits_.itorder+1)+1;              // RKK storage
  sum += szrhstmp();                         // RHS tmp size
 
  return sum;
  
} // end of method tmp_size_impl


//**********************************************************************************
//**********************************************************************************
// METHOD : compute_derived_impl
// DESC   : Compute quantities derived from state
// ARGS   : u     : input state (full)
//          sop   : operation/quantity to compute
//          utmp  : tmp space
//          uout  : resultant quantity
//          iuout : indices in uout used to contain result
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::compute_derived_impl(const State &u, GString sop,
                                            State &utmp, State &uout, 
                                            std::vector<GINT> &iuout)
{ 
  Ftype  fact1, fact2, tb;
  State  tu(1);

  if      ( "temp"     == sop ) { // temperature
    assert(uout .size() >= 1   && "Incorrect no. output components");
    assert(utmp .size() >= 2   && "Incorrect no. tmp components");
    compute_cv(u, *utmp[0], *utmp[1]);                      // Cv
   *utmp[0] = *u[DENSITY];
    if ( traits_.usebase ) *utmp[0] += *u[BASESTATE];       // total density
    geoflow::compute_temp<Ftype>(*u[ENERGY], *utmp[0], *utmp[1], *uout[0]);  // temperature
    iuout.resize(1); iuout[0] = 0;
  }
  else if ( "press"    == sop ) { // pressure (total)
    assert(uout .size() >= 1   && "Incorrect no. output components");
    assert(utmp .size() >= 4   && "Incorrect no. tmp components");
    compute_cv(u, *utmp[0], *utmp[1]);                     // Cv
   *utmp[2] = *u[DENSITY];
    if ( traits_.usebase ) *utmp[2] += *u[BASESTATE];       // total density
    geoflow::compute_temp<Ftype>(*u[ENERGY], *utmp[2], *utmp[1], *utmp[3]);  // temperature
    compute_qd(u, *utmp[0]);
    geoflow::compute_p<Ftype>(*utmp[3], *utmp[2], *utmp[0], RD, *uout[0]);
    if ( !traits_.dodry ) {
      geoflow::compute_p<Ftype>(*utmp[3], *utmp[2], *u[VAPOR], RV, *utmp[0]);
     *uout[0] += *utmp[0];
    }
    iuout.resize(1); iuout[0] = 0;
  }
  else if ( "dptemp"    == sop ) { // potential temp fluctuation
    assert(uout .size() >= 1   && "Incorrect no. output components");
    assert(utmp .size() >= 4   && "Incorrect no. tmp components");
    tu[0] = utmp[3];
    this->compute_derived_impl(u, "press", utmp, uout, iuout);
    this->compute_derived_impl(u, "temp" , utmp, tu  , iuout);
    fact1 = -RD / CPD;
    fact2 = 1.0/traits_.P0_base;
    // Compute theta from T = theta * (P/P0)^(R/C_p)
    // Note: Do we need to change definition for moist dynamics?
    for ( auto j=0; j<uout[0]->size(); j++ ) {
      (*uout[0])[j] = (*tu[0])[j] * pow(fact2*(*uout[0])[j], fact1);
    }
    if ( traits_.usebase ) { // subtract base potl temp
      // theta_base = P0/(RD *rho_base) (P_base/P0)^(1-R/Cp)
      for ( auto j=0; j<uout[0]->size(); j++ ) {
        tb = (traits_.P0_base/(RD*(*u[BASESTATE])[j])) 
           * pow((*u[BASESTATE+1])[j]/traits_.P0_base,1.0+fact1);
        (*uout[0])[j] -= tb;
      }
    }

  }
  else if ( "ptemp"    == sop ) { // potential temp
    assert(uout .size() >= 1   && "Incorrect no. output components");
    assert(utmp .size() >= 4   && "Incorrect no. tmp components");
    tu[0] = utmp[3];
    this->compute_derived_impl(u, "press", utmp, uout, iuout);
    this->compute_derived_impl(u, "temp" , utmp, tu  , iuout);
    fact1 = -RD / CPD;
    fact2 = 1.0/traits_.P0_base;
    // Compute theta from T = theta * (P/P0)^(R/C_p)
    // Note: Do we need to change definition for moist dynamics?
    for ( auto j=0; j<uout[0]->size(); j++ ) {
      (*uout[0])[j] = (*tu[0])[j] * pow(fact2*(*uout[0])[j], fact1);
    }
  }
  else if ( "den"      == sop ) { // density (total)
    assert(uout .size() >= 1   && "Incorrect no. output components");
    if ( traits_.usebase ) {
      GMTK::saxpby(*uout[0], *u[DENSITY], 1.0, *u[BASESTATE], 1.0); 
    }
    else {
      *uout[0] = *u[DENSITY];
    }
    iuout.resize(1); iuout[0] = 0;
  }
  else if ( "vel"      == sop ) { // x-velocity
    assert(uout .size() >= nc_  && "Incorrect no. output components");
    assert(utmp .size() >= 1    && "Incorrect no. tmp components");
    iuout.resize(nc_); 
    if ( !traits_.usemomden ) {
     for ( auto j=0; j<nc_; j++ ) {
       *uout[j] = *u[j];
       iuout[j] = j;
     }
    }
    else {
     *utmp[0] = *u[DENSITY];
      if ( traits_.usebase ) *utmp[0] += *u[BASESTATE];
      utmp[0]->rpow(-1.0); // 1/d
      for ( auto j=0; j<nc_; j++ ) {
       *uout[j] = *u[j]; *uout[j] *= *utmp[0];
       iuout[j] = j;
     }
    }
  }
  else if ( "ke"      == sop ) { // kinetic energy
    assert(uout .size() >= 1    && "Incorrect no. output components");
    assert(utmp .size() >= 2    && "Incorrect no. tmp components");
   *utmp[1] = *u[DENSITY];
    utmp[1]->rpow(-1);
    compute_v(u, *utmp[1], v_);
    GMTK::domathop(*grid_, v_, "vmag", utmp, uout, iuout);
    uout[0]->pointProd(*u[DENSITY]);
   *uout[0] *= 0.5;
  }
  else {
    assert(FALSE);
  }

} // end of method compute_derived_impl

