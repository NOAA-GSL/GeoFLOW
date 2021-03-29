//==================================================================================
// Module       : gstress.hpp
// Date         : 09/05/20 (DLR)
// Description  : 
//               
//                      
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none
//==================================================================================


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Default constructor
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<typename TypePack>
GStressEnOp<TypePack>::GStressEnOp(Traits &traits, Grid &grid)
:
traits_               (traits),
grid_                  (&grid),
massop_       (&grid.massop()),
nu_                  (NULLPTR),
zeta_                (NULLPTR),
eta_               (NULLPTR),
lambda_              (NULLPTR)
{
  assert(grid_->ntype().multiplicity(0) == GE_MAX-1 
        && "Only a single element type allowed on grid");
  
  if  ( traits_.Stokes_hyp ) {
    assert(traits_.nu.size() > 0);
    traits_.zeta  .resize(traits_.nu.size());
    traits_.zeta  = traits_.nu  ; traits_.zeta   *= -2.0/GDIM;
    nu_     = &traits_.nu;
    zeta_   = &traits_.zeta;
    eta_    = &traits_.nu;
    lambda_ = &traits_.zeta;
  }
  else if ( traits_.indep_diss ) {// mom & en visocisities spec independently
    assert(traits_.nu .size() > 0 && traits_.zeta  .size() > 0 );
    assert(traits_.eta.size() > 0 && traits_.lambda.size() > 0 );
    nu_     = &traits_.nu;
    zeta_   = &traits_.zeta;
    eta_    = &traits_.eta;
    lambda_ = &traits_.lambda;
  } 
  else { // energy coeffs same as for mom
    assert(traits_.nu.size() > 0 && traits_.zeta.size() > 0 );
    // From Eyink 2018 PRX 8:011022:
    traits_.zeta -= traits_.nu * (1.0/static_cast<Ftype>(GDIM));    
    nu_     = &traits_.nu;
    zeta_   = &traits_.zeta;
    eta_    = &traits_.nu;
    lambda_ = &traits_.zeta;
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
GStressEnOp<TypePack>::~GStressEnOp()
{
} // end, destructor


//**********************************************************************************
//**********************************************************************************
// METHOD : apply (1)
// DESC   : Compute application of this operator to input velocity vector:
//            so_i = [ mu (u_j,i + u_i,j) ] + zeta Div u delta_ij],j
//          or
//            so_i = [ mu (u_j,i ],j
//          where
//            mu = d nu,
//          and traits select the operator used.
//           
// ARGS   : d   : density
//          u   : input vector field
//          idir: which momentum component we're computing for
//          utmp: array of tmp arrays. At least 3 required.
//          so  : output (result) vector component, idir
//             
// RETURNS:  none
//**********************************************************************************
template<typename TypePack>
void GStressEnOp<TypePack>::apply(StateComp &d, State &u, GINT idir, State &utmp, StateComp &so)
{

  if      ( GSTRESS_FULL    == traits_.type ) {
    if ( traits_.full_colloc )
      mom_update_full_coll(d, u, idir, utmp, so);
    else
      mom_update_full_cons(d, u, idir, utmp, so);
  }
  else if ( GSTRESS_REDUCED == traits_.type ) {
    mom_update_reduced(d, u, idir, utmp, so); // only a colloc method
  }
  else {
    assert(FALSE);
  }

} // enbd, apply (1) (momementum)  


//**********************************************************************************
//**********************************************************************************
// METHOD : apply (2)
// DESC   : Compute application of this operator to input energy:
//            eo = [ eta  u_i (Del_i u_j + Del_j u_i)],j 
//               + [ lambda u_i (Div u delta_ij) ],j 
//            eta  = d eta; lambda = d zeta'
//          or
//            eo = [ eta  u_i ( Del_j u_i)],j 
//            eta  = d eta; 
//          and traits select the operator used.
//
// ARGS   : d   : density
//          u   : input vector field
//          utmp: array of tmp arrays. At least 4 required.
//          eo  : output (result) vector component, idir
//             
// RETURNS:  none
//**********************************************************************************
template<typename TypePack>
void GStressEnOp<TypePack>::apply(StateComp &d, State &u, State &utmp, StateComp &eo) 
{

  if      ( GSTRESS_FULL    == traits_.type ) {
    if ( traits_.full_colloc )
      energy_update_full_coll(d, u, utmp, eo);
    else
      energy_update_full_coll(d, u, utmp, eo);
  }
  else if ( GSTRESS_REDUCED == traits_.type ) {
    energy_update_reduced(d, u, utmp, eo); // only a colloc version
  }
  else {
    assert(FALSE);
  }

} // enbd, apply (2) (energy)  


//**********************************************************************************
//**********************************************************************************
// METHOD : mom_update_full_coll
// DESC   : Compute application of this operator to input momentum vector:
//            so_i = [ mu (u_j,i + u_i,j) ] + zeta Div u delta_ij],j
//            mu = d nu; zeta = d zeta',
//          where nu and zeta' are the kinemtic quantities. This is version
//          is based on a 'collocation' discretization of the divergence.
//          
// ARGS   : d   : density
//          u   : input vector field
//          idir: which momentum component we're computing for
//          utmp: array of tmp arrays. At least 3 required.
//          so  : output (result) vector component, idir
//             
// RETURNS:  none
//**********************************************************************************
template<typename TypePack>
void GStressEnOp<TypePack>::mom_update_full_coll(StateComp &d, State &u, GINT idir, State &utmp, StateComp &so) 
{

  assert( utmp.size() >= 4
       && "Insufficient temp space specified");

  GBOOL      usebdy = grid_->usebdydata();
  GINT       nxy = grid_->gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;

  assert( idir > 0 && idir <= nxy );

  // so = D^{j} [mu d [D_i u_j + Dj u_i) + Dk zeta d u_k delta_ij ]:
  // Below, i = idir:

  // Do D^{j} [mu (D_i u_j) ] terms:
  so = 0.0;
  for ( auto j=0; j<nxy; j++ ) { 
    grid_->deriv(*u[j], idir, *utmp[0], *utmp[1]);
    // Point-multiply by mu before taking 'divergence':
    utmp[1]->pointProd(d);
    utmp[1]->pointProd(*nu_);
    grid_->deriv(*utmp[1]  , j+1, *utmp[0], *utmp[2]);
    so += *utmp[2];
  }

  // Do D^{j} [mu d (D_j u_i) ] terms:
  for ( auto j=0; j<nxy; j++ ) { 
    grid_->deriv(*u[idir-1], j+1, *utmp[0], *utmp[1]);
    // Point-multiply by mu before taking 'divergence':
    utmp[1]->pointProd(d);
    utmp[1]->pointProd(*nu_);
    grid_->deriv(*utmp[1]  , j+1, *utmp[0], *utmp[2]);
    so += *utmp[2];
  }

  // Compute dilitation term:
  //   D^{j} (zeta d (Div u) delta_ij):
  grid_->deriv(*u[0]  , 1, *utmp[0], *utmp[1]); // store Div in utmp[1]]
  for ( auto j=1; j<nxy; j++ ) { 
    grid_->deriv(*u[j]  , j+1, *utmp[0], *utmp[2]); 
    *utmp[1] += *utmp[2];
  }

  utmp[1]->pointProd(d);
  utmp[1]->pointProd(*zeta_);  // zeta Div u

  grid_->deriv(*utmp[1], idir, *utmp[0], *utmp[2]);
  so += *utmp[2];
 
} // end of method mom_update_full_coll


//**********************************************************************************
//**********************************************************************************
// METHOD : energy_update_full_coll
// DESC   : Compute application of this operator to input energy:
//            eo = [ eta  u_i (Del_i u_j + Del_j u_i)],j 
//               + [ lambda u_i (Div u delta_ij) ],j 
//            eta  = d eta; lambda = d lambda'
//          where eta and lambda' are the kinematic quantities This is version
//          is based on a 'collocation' discretization of the divergence.
// ARGS   : d   : density
//          u   : input vector field
//          utmp: array of tmp arrays. At least 4 required.
//          eo  : output (result) vector component, idir
//             
// RETURNS:  none
//**********************************************************************************
template<typename TypePack>
void GStressEnOp<TypePack>::energy_update_full_coll(StateComp &d, State &u, State &utmp, StateComp &eo) 
{

  assert( utmp.size() >= 4
       && "Insufficient temp space specified");

  GBOOL      usebdy = grid_->usebdydata();
  GINT       nxy = grid_->gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;
  GTVector<GSIZET>          *ieface  = &grid_->gieface() ;
  GTVector<GTVector<Ftype>> *normals = &grid_->faceNormals();
  StateComp                 *bmass   = &grid_->faceMass();

  // eo = D^{j} [ eta d  u^i [D_i u_j + Dj u_i) 
  //    + lambda d u^i (Dk u^k) delta_ij ]

  // D^{j} [ eta d u^i (D_i u_j) ] terms:
  eo = 0.0;
  for ( auto j=0; j<nxy; j++ ) { 
   *utmp[1] = 0.0;
    for ( auto i=0; i<nxy; i++ ) {
       grid_->deriv(*u[j], i+1, *utmp[0], *utmp[2]);
       utmp[2]->pointProd(*u[i]);
      *utmp[1] += *utmp[2];
    }
    // Point-multiply by eta before taking 'divergence':
    utmp[1]->pointProd(d);
    utmp[1]->pointProd(*eta_);
    grid_->deriv(*utmp[1], j+1, *utmp[0], *utmp[2]);
    eo += *utmp[2];
  }

  // = D^{j} [ eta d u^i (D_j u_i) ] terms:
  for ( auto j=0; j<nxy; j++ ) { 
   *utmp[1] = 0.0;
    for ( auto i=0; i<nxy; i++ ) {
       grid_->deriv(*u[i], j+1, *utmp[0], *utmp[2]);
       utmp[2]->pointProd(*u[i]);
       *utmp[1] += *utmp[2];
    }
    // Point-multiply by eta before taking 'divergence':
    utmp[1]->pointProd(d);
    utmp[1]->pointProd(*eta_);
    grid_->deriv(*utmp[1], j+1, *utmp[0], *utmp[2]);
    eo += *utmp[2];
  }

  // Compute dilitation term:
  //   = D^{j} (lambda d (Div u) delta_ij):
  //   ... First, compute Div u:
  // (NOTE: we'll use MTK to compute Div u eventually):

  // eo = [ eta  d u_i (Del_i u_j + Del_j u_i)],j 
  //    + [ lambda d u_i (Div u delta_ij) ],j 

  grid_->deriv(*u[0]  , 1, *utmp[0], *utmp[1]); // store Div in utmp[1]]
  for ( auto j=1; j<nxy; j++ ) { 
    grid_->deriv(*u[j], j+1, *utmp[0], *utmp[2]); 
    *utmp[1] += *utmp[2];
  }

  utmp[1]->pointProd(d);
  utmp[1]->pointProd(*lambda_);

  // Now compute
  //  = D^{j} [lambda d u^i (Div u) delta_ij]:
  for ( auto j=0; j<nxy; j++ ) { 
    u[j]->pointProd(*utmp[1],*utmp[2]); 
    grid_->deriv(*utmp[2], j+1, *utmp[0], *utmp[3]); 
    eo += *utmp[3];
  }

} // end of method energy_update_full_coll


//**********************************************************************************
//**********************************************************************************
// METHOD : mom_update_full_cons 
// DESC   : Compute application of this operator to input momentum vector:
//            so_i = [ mu (u_j,i + u_i,j) ] + zeta Div u delta_ij],j
//            mu = d nu; zeta = d zeta',
//          where nu and zeta' are the kinemtic quantities. This is version
//          is based on a 'conservative' discretization of the divergence.
//          
// ARGS   : d   : density
//          u   : input vector field
//          idir: which momentum component we're computing for
//          utmp: array of tmp arrays. At least 3 required.
//          so  : output (result) vector component, idir
//             
// RETURNS:  none
//**********************************************************************************
template<typename TypePack>
void GStressEnOp<TypePack>::mom_update_full_cons(StateComp &d, State &u, GINT idir, State &utmp, StateComp &so) 
{

  assert( utmp.size() >= 4
       && "Insufficient temp space specified");

  GBOOL      usebdy = grid_->usebdydata();
  GINT       nxy = grid_->gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;
  GINT                       isgn;
  GSIZET                     k;
  Ftype                      fsgn;
  GTVector<GSIZET>          *ieface  = &grid_->gieface() ;
  GTVector<GTVector<Ftype>> *normals = &grid_->faceNormals();
  StateComp                 *bmass   = &grid_->faceMass();

  assert( idir > 0 && idir <= nxy );

#if defined(DO_COMPRESS_MODES_ONLY)
  tfact_.resizem(u[0]->size());
#endif

  // so = -D^{T,j} [mu d [D_i u_j + Dj u_i) + Dk zeta d u_k delta_ij ]:
  //    + bdy surface terms:
  // Below, i = idir:

  // Do -D^{T,j} [mu (D_i u_j) ] terms:
  so = 0.0;
  for ( auto j=0; j<nxy; j++ ) { 
    grid_->deriv(*u[j], idir, *utmp[0], *utmp[1]);
    // Point-multiply by mu before taking 'divergence':
    utmp[1]->pointProd(d);
    utmp[1]->pointProd(*nu_);
    grid_->wderiv(*utmp[1]  , j+1, TRUE, *utmp[0], *utmp[2]);
    so -= *utmp[2];

    // Compute bdy terms for this component, j:
    for ( auto b=0; usebdy && b<ieface->size(); b++ ) {
      k = (*ieface)[b];
      so[k] += (*utmp[1])[k] * (*normals)[j][b] * (*bmass)[b];
    }
  }

  // Do -D^{T,j} [mu d (D_j u_i) ] terms:
  for ( auto j=0; j<nxy; j++ ) { 
    grid_->deriv(*u[idir-1], j+1, *utmp[0], *utmp[1]);
    // Point-multiply by mu before taking 'divergence':
    utmp[1]->pointProd(d);
    utmp[1]->pointProd(*nu_);
    grid_->wderiv(*utmp[1]  , j+1, TRUE, *utmp[0], *utmp[2]);
    so -= *utmp[2];

    // Compute surface terms for this component, j:
    for ( auto b=0; usebdy && b<ieface->size(); b++ ) {
      k = (*ieface)[b];
      so[k] += (*utmp[1])[k] * (*normals)[j][b] * (*bmass)[b];
    }
  }

  // Compute dilitation term:
  //   -D^{T,j} (zeta d (Div u) delta_ij):
  grid_->deriv(*u[0]  , 1, *utmp[0], *utmp[1]); // store Div in utmp[1]]
  for ( auto j=1; j<nxy; j++ ) { 
    grid_->deriv(*u[j]  , j+1, *utmp[0], *utmp[2]); 
    *utmp[1] += *utmp[2];
  }

  utmp[1]->pointProd(d);
  utmp[1]->pointProd(*zeta_);  // zeta Div u

#if defined(DO_COMPRESS_MODES_ONLY)
  for ( auto i=0; i<utmp[1]->size(); i++ ) {
    isgn      = sgn<Ftype>((*utmp[1])[i]);
    fsgn      = static_cast<Ftype>(isgn);
    tfact_[i] = isgn == 0 ? 0.0 : 0.5*(1.0-fsgn);
  }
#endif

  grid_->wderiv(*utmp[1], idir, TRUE, *utmp[0], *utmp[2]);
#if defined(DO_COMPRESS_MODES_ONLY)
  utmp[2]->pointProd(tfact_);
#endif
  so -= *utmp[2];
 
  // Compute surface terms for
  //  Integral zeta d (Div u) delta_ij.n^j dV:
  // Use kernel above, for i=idir:
  for ( auto b=0; usebdy && b<ieface->size(); b++ ) {
    k = (*ieface)[b];
  #if defined(DO_COMPRESS_MODES_ONLY)
    so[k] += (*utmp[1])[k]*tfact_[k] * (*normals)[idir-1][b] * (*bmass)[b];
  #else
    so[k] += (*utmp[1])[k] * (*normals)[idir-1][b] * (*bmass)[b];
  #endif
  }

} // end of method mom_update_full_cons


//**********************************************************************************
//**********************************************************************************
// METHOD : energy_update_full_cons
// DESC   : Compute application of this operator to input energy:
//            eo = [ eta  u_i (Del_i u_j + Del_j u_i)],j 
//               + [ lambda u_i (Div u delta_ij) ],j 
//            eta  = d eta; lambda = d lambda'
//          where eta and lambda' are the kinematic quantities This is version
//          is based on a 'conservative' discretization of the divergence.
// ARGS   : d   : density
//          u   : input vector field
//          utmp: array of tmp arrays. At least 4 required.
//          eo  : output (result) vector component, idir
//             
// RETURNS:  none
//**********************************************************************************
template<typename TypePack>
void GStressEnOp<TypePack>::energy_update_full_cons(StateComp &d, State &u, State &utmp, StateComp &eo) 
{

  assert( utmp.size() >= 4
       && "Insufficient temp space specified");

  GBOOL      usebdy = grid_->usebdydata();
  GINT       nxy = grid_->gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;
  Ftype                      isgn;
  GSIZET                     k;
  Ftype                      fsgn;
  GTVector<GSIZET>          *ieface  = &grid_->gieface() ;
  GTVector<GTVector<Ftype>> *normals = &grid_->faceNormals();
  StateComp                 *bmass   = &grid_->faceMass();

#if defined(DO_COMPRESS_MODES_ONLY)
  tfact_.resizem(u[0]->size());
#endif

  // eo -= D^{T,j} [ eta d  u^i [D_i u_j + Dj u_i) 
  //    + lambda d u^i (Dk u^k) delta_ij ]
  //    + surface terms:

  // - D^{T,j} [ eta d u^i (D_i u_j) ] terms:
  eo = 0.0;
  for ( auto j=0; j<nxy; j++ ) { 
   *utmp[1] = 0.0;
    for ( auto i=0; i<nxy; i++ ) {
       grid_->deriv(*u[j], i+1, *utmp[0], *utmp[2]);
       utmp[2]->pointProd(*u[i]);
      *utmp[1] += *utmp[2];
    }
    // Point-multiply by eta before taking 'divergence':
    utmp[1]->pointProd(d);
    utmp[1]->pointProd(*eta_);
    grid_->wderiv(*utmp[1], j+1, TRUE, *utmp[0], *utmp[2]);
    eo -= *utmp[2];

    // Do the surface terms for jth component of normal:
    for ( auto b=0; usebdy && b<ieface->size(); b++ ) {
      k = (*ieface)[b];
      eo[k] += (*utmp[1])[k] * (*normals)[j][b] * (*bmass)[b];
    }
  }

  // -= D^{T,j} [ eta d u^i (D_j u_i) ] terms:
  for ( auto j=0; j<nxy; j++ ) { 
   *utmp[1] = 0.0;
    for ( auto i=0; i<nxy; i++ ) {
       grid_->deriv(*u[i], j+1, *utmp[0], *utmp[2]);
       utmp[2]->pointProd(*u[i]);
       *utmp[1] += *utmp[2];
    }
    // Point-multiply by eta before taking 'divergence':
    utmp[1]->pointProd(d);
    utmp[1]->pointProd(*eta_);
    grid_->wderiv(*utmp[1], j+1, TRUE, *utmp[0], *utmp[2]);
    eo -= *utmp[2];

    // Do the surface terms for jth component of normal:
    for ( auto b=0; usebdy && b<ieface->size(); b++ ) {
      k = (*ieface)[b];
      eo[k] += (*utmp[1])[k] * (*normals)[j][b] * (*bmass)[b];
    }
  }

  // Compute dilitation term:
  //   -= D^{T,j} (lambda d (Div u) delta_ij):
  //   ... First, compute Div u:
  // (NOTE: we'll use MTK to compute Div u eventually):

  // eo = [ eta  d u_i (Del_i u_j + Del_j u_i)],j 
  //    + [ lambda d u_i (Div u delta_ij) ],j 

  grid_->deriv(*u[0]  , 1, *utmp[0], *utmp[1]); // store Div in utmp[1]]
  for ( auto j=1; j<nxy; j++ ) { 
    grid_->deriv(*u[j], j+1, *utmp[0], *utmp[2]); 
    *utmp[1] += *utmp[2];
  }

  utmp[1]->pointProd(d);
  utmp[1]->pointProd(*lambda_);

#if defined(DO_COMPRESS_MODES_ONLY)
  for ( auto i=0; i<utmp[1]->size(); i++ ) {
    isgn      = sgn<Ftype>((*utmp[1])[i]);
    fsgn      = static_cast<Ftype>(isgn);
    tfact_[i] = isgn == 0 ? 1.0 : 0.5*(1.0-fsgn);
  }
#endif

  // Now compute
  //  -= D^{T,j} [lambda d u^i (Div u) delta_ij]:
  for ( auto j=0; j<nxy; j++ ) { 
    u[j]->pointProd(*utmp[1],*utmp[2]); 
    grid_->wderiv(*utmp[2], j+1, TRUE, *utmp[0], *utmp[3]); 
#if defined(DO_COMPRESS_MODES_ONLY)
    utmp[3]->pointProd(tfact_);
#endif
    eo -= *utmp[3];

    // Do the surface terms for jth component of normal:
    for ( auto b=0; usebdy && b<ieface->size(); b++ ) {
      k = (*ieface)[b];
    #if defined(DO_COMPRESS_MODES_ONLY)
      eo[k] += (*utmp[2])[k]*tfact_[k] * (*normals)[j][b] * (*bmass)[b];
    #else
      eo[k] += (*utmp[2])[k] * (*normals)[j][b] * (*bmass)[b];
    #endif
    }
  }

} // end of method energy_update_full_cons


//**********************************************************************************
//**********************************************************************************
// METHOD : mom_update_reduced
// DESC   : Compute application of this operator to input momentum vector:
//            so_i = d [nu u_i,j ],j
//          where mu = nu d, and nu is the kinematic viscosity. This 
//          operator uses 'collocation' discretization only.
// ARGS   : d   : density
//          u   : input vector field
//          idir: which momentum component we're computing for
//          utmp: array of tmp arrays. At least 3 required.
//          so  : output (result) vector component, idir
//             
// RETURNS:  none
//**********************************************************************************
template<typename TypePack>
void GStressEnOp<TypePack>::mom_update_reduced(StateComp &d, State &u, GINT idir, State &utmp, StateComp &so) 
{

  assert( utmp.size() >= 4
       && "Insufficient temp space specified");

  GBOOL      usebdy = grid_->usebdydata();
  GINT       nxy = grid_->gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;

  assert( idir > 0 && idir <= nxy );


  // so_idir = d D^{j} [nu Dj u_idir) ]:
  // Below, i = idir:

  so = 0.0;

  // Do D^{j} [nu d (D_j u_i) ] terms:
  for ( auto j=0; j<nxy; j++ ) { 
    grid_->deriv(*u[idir-1], j+1, *utmp[0], *utmp[1]);
    // Point-multiply by nu before taking 'divergence':
//  utmp[1]->pointProd(d);
    utmp[1]->pointProd(*nu_);
    grid_->deriv(*utmp[1], j+1, *utmp[0], *utmp[2]);
    so += *utmp[2];
  }
  so.pointProd(d);


} // end of method mom_update_reduced


//**********************************************************************************
//**********************************************************************************
// METHOD : energy_update_reduced
// DESC   : Compute application of this operator to input energy:
//            eo = [ eta  u_i  Del_j u_i],j 
//          where
//            eta  = d eta, and eta is the kinematic viscosity. This
//          operator uses 'collocation' discretization only.
// ARGS   : d   : density
//          u   : input vector field
//          utmp: array of tmp arrays. At least 4 required.
//          eo  : output (result) vector component, idir
//             
// RETURNS:  none
//**********************************************************************************
template<typename TypePack>
void GStressEnOp<TypePack>::energy_update_reduced(StateComp &d, State &u, State &utmp, StateComp &eo) 
{

  assert( utmp.size() >= 4
       && "Insufficient temp space specified");

  GBOOL      usebdy = grid_->usebdydata();
  GINT       nxy = grid_->gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;

  // eo = d D^{j} [ eta u^i Dj u_i] 

  eo = 0.0;

  // = D^{j} [ eta d u^i (D_j u_i) ] terms:
  for ( auto j=0; j<nxy; j++ ) { 
   *utmp[1] = 0.0;
    for ( auto i=0; i<nxy; i++ ) {
       grid_->deriv(*u[i], j+1, *utmp[0], *utmp[2]);
       utmp[2]->pointProd(*u[i]);
       *utmp[1] += *utmp[2];
    }
    // Point-multiply by eta before taking 'divergence':
//  utmp[1]->pointProd(d);
    utmp[1]->pointProd(*eta_);
    grid_->deriv(*utmp[1], j+1, *utmp[0], *utmp[2]);
    eo += *utmp[2];
  }
  eo.pointProd(d);

} // end of method energy_update_reduced


