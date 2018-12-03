//**********************************************************************************
//**********************************************************************************
// METHOD : compute_grefderivs
// DESC   : Compute tensor product derivs of specified field, u, in ref space
//          for grid, using grid object to determine which to compute. Compute:
//            du = [ I_X_I_X_Dx
//                   I_X_Dy_X_I
//                   Dz_X_I_X_I].
//     
//          where Dx, Dy, Dz are 1d derivative objects from basis functions     
// ARGS   : 
//          grid   : GGrid object 
//          u      : input field whose derivative we want, allocated globally 
//                   (e.g., for all elements).
//          du     : vector of length 2 or 3 containing the derivatives.
//                   If using GE_REGULAR in 2D, we only need to vector
//                   elements; else we need 3. These should be allocated globally.
//          dotrans: flag telling us to tak transpose of deriv operators (TRUE) or
//                   not (FALSE).
//             
// RETURNS:  none
//**********************************************************************************
void GHelmholtz::compute_grefderivs( GGrid &grid, GTVector<GFTYPE> &u,
                                GTVector<GTVector<GFTYPE>*> &du, GBOOL dotrans)
{
  assert((grid.itype(GE_REGULAR)   .size() > 0 && du.size() >= GDIM)
       ||(grid.itype(GE_DEFORMED)  .size() > 0 && du.size() >= GDIM)
       && "Insufficient number of derivatives specified");
  

  GSIZET                       ibeg, iend; // beg, end indices for global array
  GTVector<GSIZET>             N(GDIM);
  GTVector<GTMatrix<GFTYPE>*>  Di(GDIM);   // element-based 1d derivative operators
  GElemList                   *gelems = &grid.elems();

#if defined(_G_IS2D)

  for ( GSIZET e=0; e<grid.elems().size(); e++ ) {
    ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
    u.range(ibeg, iend); // restrict global vecs to local range
    for ( GSIZET k=0; k<GDIM; k++ ) du[k]->range(ibeg, iend);
    for ( GSIZET k=0; k<GDIM; k++ ) N[k]= (*gelems)[e]->size(k);
    Di[0] = (*gelems)[e]->gbasis(0)->getDerivMatrix (dotrans);
    Di[1] = (*gelems)[e]->gbasis(1)->getDerivMatrix(!dotrans);
    GMTK::I2_X_D1(*Di[0], u, N[0], N[1], *du[0]); 
    GMTK::D2_X_I1(*Di[1], u, N[0], N[1], *du[1]); 
  }
  u.range(0, grid.ndof()); // restrict global vec to local range
  for ( GSIZET k=0; k<GDIM+1; k++ ) du[k]->range(0, grid.ndof()-1);

#elif defined((_G_IS3D)

  for ( GSIZET e=0; e<grid.elems().size(); e++ ) {
    ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
    u.range(ibeg, iend); // restrict global vecs to local range
    for ( GSIZET k=0; k<GDIM+1; k++ ) du[k]->range(ibeg, iend);
    for ( GSIZET k=0; k<GDIM  ; k++ ) N[k]= (*gelems)[e]->size(k);
    Di[0] = (*gelems)[e]->gbasis(0)->getDerivMatrix (dotrans); 
    Di[1] = (*gelems)[e]->gbasis(1)->getDerivMatrix(!dotrans); 
    Di[2] = (*gelems)[e]->gbasis(1)->getDerivMatrix(!dotrans); 
    GMTK::I3_X_I2_X_D1(*Di[0], u, N[0], N[1], N[2], *du[0]); 
    GMTK::I3_X_D2_X_I1(*Di[1], u, N[0], N[1], N[2], *du[1]); 
    GMTK::D3_X_I2_X_I1(*Di[2], u, N[0], N[1], N[2], *du[2]); 
  }
  u.range(0, grid.ndof()); // reset global vec to globalrange
  for ( GSIZET k=0; k<GDIM+1; k++ ) du[k]->range(0, grid.ndof()-1);

#endif

} // end of method compute_grefderivs


//**********************************************************************************
//**********************************************************************************
// METHOD : compute_grefderivsW
// DESC   : Compute tensor product derivs of specified field, u, in ref space 
//          for grid using grid object to determine which to compute, and 
//          include weights.
//          Compute:
//            du = (Mz_X_My_Mx) [ I_X_I_X_Dx
//                                I_X_Dy_X_I
//                                Dz_X_I_X_I].
//     
//          where Dx, Dy, Dz are 1d derivative objects from basis functions, and
//          Mi are the (diagonal) 1d-weights (or 'mass functions'). This can be 
//          re-written as
//            du = [ Mz_X_My_X_(Mx Dx)
//                   Mz_X_(My Dy)_X_Mx
//                  (Mz Dz)_X_My_X_Mx],
//           with comparable expressions for 2d.
// ARGS   : 
//          grid   : GGrid object 
//          u      : input field whose derivative we want, allocated globally 
//                   (e.g., for all elements).
//          du     : vector of length 2 or 3 containing the derivatives.
//                   If using GE_REGULAR in 2D, we only need to vector
//                   elements; else we need 3. These should be allocated globally.
//          dotrans: flag telling us to tak transpose of deriv operators (TRUE) or
//                   not (FALSE).
//             
// RETURNS:  none
//**********************************************************************************
void GHelmholtz::compute_grefderivsW(GGrid &grid, GTVector<GFTYPE> &u,
                                 GTVector<GTVector<GFTYPE>*> &du, GBOOL dotrans)
{
  assert((grid.itype(GE_REGULAR)   .size() > 0 && du.size() >= GDIM)
       ||(grid.itype(GE_DEFORMED)  .size() > 0 && du.size() >= GDIM)
       && "Insufficient number of derivatives specified");
  

  GSIZET                       ibeg, iend; // beg, end indices for global array
  GTVector<GSIZET>             N(GDIM);
  GTVector<GTVector<GFTYPE>*>  W(GDIM);    // element weights
  GTVector<GTMatrix<GFTYPE>*>  Di(GDIM);   // element-based 1d derivative operators
  GElemList                   *gelems = &grid.elems();

#if defined(_G_IS2D)

  for ( GSIZET e=0; e<grid.elems().size(); e++ ) {
    ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
    u.range(ibeg, iend); // restrict global vecs to local range
    for ( GSIZET k=0; k<GDIM; k++ ) du[k]->range(ibeg, iend);
    for ( GSIZET k=0; k<GDIM; k++ ) N[k]= (*gelems)[e]->size(k);
    for ( GSIZET k=0; k<GDIM; k++ ) W[k]= (*gelems)[e]->gbasis(k)->getWeights();
    Di[0] = (*gelems)[e]->gbasis(0)->getDerivMatrix (dotrans);
    Di[1] = (*gelems)[e]->gbasis(1)->getDerivMatrix(!dotrans);
    GMTK::Dg2_X_D1(*Di[0], *W[1], u, tmp1_, *du[0]); 
    GMTK::D2_X_Dg1(*W[0], *Di[1], u, tmp1_, *du[1]); 
  }
  u.range(0, grid.ndof()); // restrict global vec to local range
  for ( GSIZET k=0; k<GDIM+1; k++ ) du[k]->range(0, grid.ndof()-1);

#elif defined((_G_IS3D)

  for ( GSIZET e=0; e<grid.elems().size(); e++ ) {
    ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
    u.range(ibeg, iend); // restrict global vecs to local range
    for ( GSIZET k=0; k<GDIM+1; k++ ) du[k]->range(ibeg, iend);
    for ( GSIZET k=0; k<GDIM; k++ ) {
      N[k]= (*gelems)[e]->size(k);
      W[k]= (*gelems)[e]->gbasis(k)->getWeights();
    }
    Di[0] = (*gelems)[e]->gbasis(0)->getDerivMatrix (dotrans); 
    Di[1] = (*gelems)[e]->gbasis(1)->getDerivMatrix(!dotrans); 
    Di[2] = (*gelems)[e]->gbasis(1)->getDerivMatrix(!dotrans); 
    GMTK::Dg3_X_Dg2_X_D1(*Di[0], *W [1], *W [2], u, N[0], N[1], N[2], tmp1_, *du[0]); 
    GMTK::Dg3_X_D2_X_Dg1(*W [0], *Di[1], *W [2], u, N[0], N[1], N[2], tmp1_, *du[1]); 
    GMTK::D3_X_Dg2_X_Dg1(*W [0], *W [1], *Di[2], u, N[0], N[1], N[2], tmp1_, *du[2]); 
  }
  u.range(0, grid.ndof()); // reset global vec to globalrange
  for ( GSIZET k=0; k<GDIM+1; k++ ) du[k]->range(0, grid.ndof()-1);

#endif

} // end of method compute_grefderivsW


//**********************************************************************************
//**********************************************************************************
// METHOD : compute_grefdiv
// DESC   : Compute tensor product 'divergence' of input fields in ref space
//          for grid:
//             Div u = [I_X_I_X_Dx     |u1|
//                      I_X_Dy_X_I   . |u2|
//                      Dz_X_I_X_I]    |u3|
//          or
//     
//             Div u = [I_X_I_X_DxT     |u1|
//                      I_X_DyT_X_I   . |u2|
//                      DzT_X_I_X_I]    |u3|
//          where Dx(T), Dy(T), Dz(T) are 1d derivative objects from basis functions     
// ARGS   : 
//          grid   : GGrid object 
//          u      : input vector field whose divergence we want, allocated globally 
//                   (e.g., for all elements). Must have GDIM components, unless
//                   we're using an embedded grid, when GDIM=2, when it should have
//                   3 components.
//          divu   : scalar result
//          dotrans: flag telling us to tak transpose of deriv operators (TRUE) or
//                   not (FALSE).
//             
// RETURNS:  none
//**********************************************************************************
void GHelmholtz::compute_grefdiv(GGrid &grid, GTVector<GTVector<GFTYPE>*> &u, 
                             GTVector<GFTYPE> &divu, GBOOL dotrans)
{

  GBOOL                        bembedded;
  GSIZET                       ibeg, iend; // beg, end indices for global array
  GTVector<GSIZET>             N(GDIM);
  GTVector<GTMatrix<GFTYPE>*>  Di(GDIM);   // element-based 1d derivative operators
  GElemList                   *gelems = &grid.elems();

  bembedded=grid.itype(GE_DEFORMED).size() > 0 && GDIM == 2;
  assert(( (bembedded && u.size()==3) 
        || (!bembedded&& u.size()==GDIM) )
       && "Insufficient number of vector field components provided");

  divu = 0.0;
#if defined(_G_IS2D)

  for ( GSIZET e=0; e<grid.elems().size(); e++ ) {
    // restrict global vecs to local range
    ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
    divu.range(ibeg,iend); 
    for ( GSIZET k=0; k<u.size(); k++ ) u[k]->range(ibeg, iend); 
    Di[0] = (*gelems)[e]->gbasis(0)->getDerivMatrix (dotrans);
    Di[1] = (*gelems)[e]->gbasis(1)->getDerivMatrix(!dotrans);
    for ( GSIZET k=0; k<GDIM; k++ ) N[k]= (*gelems)[e]->size(k);
    tmp1_.resizem((*gelems)[e]->nnodes());
    GMTK::I2_X_D1(*Di[0], *u[0], N[0], N[1], tmp1_); // D1 u1
    divu += tmp1_;
    GMTK::D2_X_I1(*Di[1], *u[1], N[0], N[1], tmp1_); // D2 u2
    divu += tmp1_;
    if ( bembedded ) divu += *u[2]; // D3 acts as I
  }

#elif defined((_G_IS3D)

  for ( GSIZET e=0; e<grid.elems().size(); e++ ) {
    ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
    divu.range(ibeg,iend); 
    for ( GSIZET k=0; k<u.size(); k++ ) u[k]->range(ibeg, iend); 
    for ( GSIZET k=0; k<GDIM; k++ ) N[k]= (*gelems)[e]->size(k);
    tmp1_.resizem((*gelems)[e]->nnodes());
    Di[0] = (*gelems)[e]->gbasis(0)->getDerivMatrix (dotrans); 
    Di[1] = (*gelems)[e]->gbasis(1)->getDerivMatrix(!dotrans); 
    Di[2] = (*gelems)[e]->gbasis(1)->getDerivMatrix(!dotrans); 

    GMTK::I3_X_I2_X_D1(*Di[0], *u[0], N[0], N[1], N[2], tmp1_); // D1 u1
    divu += tmp1_;
    GMTK::I3_X_D2_X_I1(*Di[1], *u[1], N[0], N[1], N[2], tmp1_); // D2 u2
    divu += tmp1_;
    GMTK::D3_X_I2_X_I1(*Di[2], *u[2], N[0], N[1], N[2], tmp1_); // D3 u3
    divu += tmp1_;
  }
  ibeg = 0; iend = grid.ndof()-1;
  divu.range(ibeg,iend); 
  for ( GSIZET k=0; k<u.size(); k++ ) u[k]->range(ibeg, iend); 

#endif

} // end, method compute_grefdiv


