//==================================================================================
// Module       : gstress.hpp
// Date         : 09/05/20 (DLR)
// Description  : Represents the SEM discretization of the full viscous
//                stress-energy operator. The effect of the viscous stress in the 
//                momentum eqution is
//                    F_i = [2  mu s_{ij}],j + (zeta Div u delta_{ij}),j,
//                where
//                    s_{ij} = (u_j,i + u_i,j)/2 - 1/2d Div u delta_{ij}, and
//                d is the problem dimension. The viscous stress-energy for the 
//                energy equation is
//                    [2 kappa u_i s_{ij}],j - [lambda u_i Div u delta_{ij}],j
//                where u_i is the velocity, and mu, the (shear) viscosity, zeta is
//                the 'bulk' viscosity. Strictly speaking, kappa=mu, and lambda=zeta,
//                but we allow these to be set independently for now. Repeated
//                indices are summed here.  mu, zeta, kappa, lambda, may vary
//                in space or be constant. Currently, the so-called Stokes 
//                approximation is used by default s.t.
//                      (zeta -  mu/d) = -2/3 mu, and
//                      (lambda -  kappa/d ) = -2/3 kappa.
//               
//                For the energy, this operator is nonlinear, 
//                so it should not derive from GLinOp. 
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
mu_                  (NULLPTR),
zeta_                (NULLPTR),
kappa_               (NULLPTR),
lambda_              (NULLPTR)
{
  assert(grid_->ntype().multiplicity(0) == GE_MAX-1 
        && "Only a single element type allowed on grid");
  
  if ( traits_.indep_diss) {
    assert(traits_.mu.size() > 0 && traits_.kappa.size()  > 0);
    if  ( traits_.Stokes_hyp ) {
      traits_.zeta.resize(traits_.mu.size());
      traits_.lambda.resize(traits_.kappa.size());
      traits_.zeta  = traits_.mu   ; traits_.zeta   *= -2.0/3.0;
      traits_.lambda= traits_.kappa; traits_.lambda *= -2.0/3.0;
    } else {
      assert(traits_.zeta.size() > 0 && traits_.lambda.size()  > 0);
      traits_.zeta   -= (traits_.mu * (1.0/GDIM));
      traits_.lambda -= (traits_.kappa * (1.0/GDIM));
    }
    mu_     = &traits_.mu;
    zeta_   = &traits_.zeta;
    kappa_  = &traits_.kappa;
    lambda_ = &traits_.lambda;
  } else { // kappa, lambda depend only on mu, zeta:
    assert(traits_.mu.size() > 0 );
    if  ( traits_.Stokes_hyp ) {
      traits_.zeta.resize(traits_.mu.size());
      traits_.zeta = traits_.mu; traits_.zeta *= -2.0/3.0;
    } else {
      assert(traits_.zeta.size() > 0 );
      traits_.zeta   -= (traits_.mu * (1.0/GDIM));
    }
    mu_     = &traits_.mu;
    zeta_   = &traits_.zeta;
    kappa_  = &traits_.mu;
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
// DESC   : Compute application of this operator to input momentum vector:
//            so_i = [ mu (u_j,i + u_i,j) ] + zeta Div u delta_ij],j
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

  assert( utmp.size() >= 4
       && "Insufficient temp space specified");

  GBOOL      usebdy = grid_->usebdydata();
  GINT       nxy = grid_->gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;
  GINT                       isgn;
  GSIZET                     k;
  Ftype                      fsgn;
  GTVector<GSIZET>          *igbdy   = &grid_->igbdy() ;
  GTVector<GTVector<Ftype>> *normals = &grid_->bdyNormals();
  StateComp                 *mass    =  grid_->massop().data();
  StateComp                 *bmass   = &grid_->bdyMass();

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
    utmp[1]->pointProd(*mu_);
    grid_->wderiv(*utmp[1]  , j+1, TRUE, *utmp[0], *utmp[2]);
    so -= *utmp[2];

#if defined(DO_BDY)
    // Compute bdy terms for this component, j:
    for ( auto b=0; usebdy && b<igbdy->size(); b++ ) {
      k = (*igbdy)[b];
      so[k] += (*utmp[1])[k] * (*normals)[j][b] * (*bmass)[b];
    }
#endif
  }

  // Do -D^{T,j} [mu d (D_j u_i) ] terms:
  for ( auto j=0; j<nxy; j++ ) { 
    grid_->deriv(*u[idir-1], j+1, *utmp[0], *utmp[1]);
    // Point-multiply by mu before taking 'divergence':
    utmp[1]->pointProd(d);
    utmp[1]->pointProd(*mu_);
    grid_->wderiv(*utmp[1]  , j+1, TRUE, *utmp[0], *utmp[2]);
    so -= *utmp[2];

#if defined(DO_BDY)
    // Compute surface terms for this component, j:
    for ( auto b=0; usebdy && b<igbdy->size(); b++ ) {
      k = (*igbdy)[b];
      so[k] += (*utmp[1])[k] * (*normals)[j][b] * (*bmass)[b];
    }
#endif
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
 
#if defined(DO_BDY)
  // Compute surface terms for
  //  Integral zeta d (Div u) delta_ij.n^j dV:
  // Use kernel above, for i=idir:
  for ( auto b=0; usebdy && b<igbdy->size(); b++ ) {
    k = (*igbdy)[b];
  #if defined(DO_COMPRESS_MODES_ONLY)
    so[k] += (*utmp[1])[k]*tfact_[k] * (*normals)[idir-1][b] * (*bmass)[b];
  #else
    so[k] += (*utmp[1])[k] * (*normals)[idir-1][b] * (*bmass)[b];
  #endif
  }
#endif 



} // end of method apply (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : apply (2)
// DESC   : Compute application of this operator to input energy:
//            eo = [ kappa  u_i (Del_i u_j + Del_j u_i)],j 
//               + [ lambda u_i (Div u delta_ij) ],j 
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

  assert( utmp.size() >= 4
       && "Insufficient temp space specified");

  GBOOL      usebdy = grid_->usebdydata();
  GINT       nxy = grid_->gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;
  Ftype                      isgn;
  GSIZET                     k;
  Ftype                      fsgn;
  GTVector<GSIZET>          *igbdy   = &grid_->igbdy() ;
  GTVector<GTVector<Ftype>> *normals = &grid_->bdyNormals();
  StateComp                 *mass    =  grid_->massop().data();
  StateComp                 *bmass   = &grid_->bdyMass();

#if defined(DO_COMPRESS_MODES_ONLY)
  tfact_.resizem(u[0]->size());
#endif

  // eo -= D^{T,j} [ kappa d  u^i [D_i u_j + Dj u_i) 
  //    + lambda d u^i (Dk u^k) delta_ij ]
  //    + surface terms:

  // - D^{T,j} [ kappa d u^i (D_i u_j) ] terms:
  eo = 0.0;
  for ( auto j=0; j<nxy; j++ ) { 
   *utmp[1] = 0.0;
    for ( auto i=0; i<nxy; i++ ) {
       grid_->deriv(*u[j], i+1, *utmp[0], *utmp[2]);
       utmp[2]->pointProd(*u[i]);
      *utmp[1] += *utmp[2];
    }
    // Point-multiply by kappa before taking 'divergence':
    utmp[1]->pointProd(d);
    utmp[1]->pointProd(*kappa_);
    grid_->wderiv(*utmp[1], j+1, TRUE, *utmp[0], *utmp[2]);
    eo -= *utmp[2];

#if defined(DO_BDY)
    // Do the surface terms for jth component of normal:
    for ( auto b=0; usebdy && b<igbdy->size(); b++ ) {
      k = (*igbdy)[b];
      eo[k] += (*utmp[1])[k] * (*normals)[j][b] * (*bmass)[b];
    }
#endif
  }

  // -= D^{T,j} [ kappa d u^i (D_j u_i) ] terms:
  for ( auto j=0; j<nxy; j++ ) { 
   *utmp[1] = 0.0;
    for ( auto i=0; i<nxy; i++ ) {
       grid_->deriv(*u[i], j+1, *utmp[0], *utmp[2]);
       utmp[2]->pointProd(*u[i]);
       *utmp[1] += *utmp[2];
    }
    // Point-multiply by kappa before taking 'divergence':
    utmp[1]->pointProd(d);
    utmp[1]->pointProd(*kappa_);
    grid_->wderiv(*utmp[1], j+1, TRUE, *utmp[0], *utmp[2]);
    eo -= *utmp[2];

#if defined(DO_BDY)
    // Do the surface terms for jth component of normal:
    for ( auto b=0; usebdy && b<igbdy->size(); b++ ) {
      k = (*igbdy)[b];
      eo[k] += (*utmp[1])[k] * (*normals)[j][b] * (*bmass)[b];
    }
#endif
  }

  // Compute dilitation term:
  //   -= D^{T,j} (lambda d (Div u) delta_ij):
  //   ... First, compute Div u:
  // (NOTE: we'll use MTK to compute Div u eventually):

  // eo = [ kappa  d u_i (Del_i u_j + Del_j u_i)],j 
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

#if defined(DO_BDY)
    // Do the surface terms for jth component of normal:
    for ( auto b=0; usebdy && b<igbdy->size(); b++ ) {
      k = (*igbdy)[b];
    #if defined(DO_COMPRESS_MODES_ONLY)
      eo[k] += (*utmp[2])[k]*tfact_[k] * (*normals)[j][b] * (*bmass)[b];
    #else
      eo[k] += (*utmp[2])[k] * (*normals)[j][b] * (*bmass)[b];
    #endif
    }
#endif
  }

} // end of method apply (2)


