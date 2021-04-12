//==================================================================================
// Module       : gprojection_filter.hpp
// Date         : 4/10/21 (DLR)
// Description  : Computes a projectioon filter for stabilization.
//                Taken from Deville, Fischer & Mund "High-Order
//                Methods for Incompressible Flow"
//                 
//                Define interpolation matrices,
//                    I_N^M(i,j) = h_N,j(xi_M,i)
//                where h_N is the Lagrange interpolating polynomial
//                of order N, evaluated at nodes, x_M,i from the Mth 
//                polynomial. Then define
//                    P_N^M = I_M^N I_N^M.
//                The filter is then defined in 1d as
//                    F = alpha Pi_N^M + 1-alpha) I_N^N
//                where
//                    I_N^N 
//                is the Nth-order identify matrix. Filter, F, is then 
//                applied in tensor product form.
// Copyright    : Copyright 2021. Colorado State University. All rights reserved.
// Derived From : FilterBase
//==================================================================================

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Default constructor
// ARGS   : grid    : Grid object
//          ifilter : starting mode for filtering
//          mufilter: filter factor (by which to truncatte)
// RETURNS: none
//**********************************************************************************
template<typename TypePack>
GProjectionFilter<TypePack>::GProjectionFilter(Traits &traits, Grid &grid)
:
bInit_               (FALSE),
traits_              (traits),
grid_                (&grid)
{

//assert(grid_->ntype().multiplicity(0) == GE_MAX-1 
//      && "Only a single element type allowed on grid");
  assert(grid_->ispconst() ); // order must not vary 
  assert(traits_.pdelta.size() >= GDIM && traits_.alpha.size() >= GDIM);
} // end of constructor method (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor method 
// DESC   : 
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<typename TypePack>
GProjectionFilter<TypePack>::~GProjectionFilter()
{
} // end, destructor


//**********************************************************************************
//**********************************************************************************
// METHOD : apply_impl
// DESC   : Compute application of this filter to input vector
//           
// ARGS   : t   : Time
//          u   : input vector field
//          utmp: array of tmp arrays; not used here
//          uo  : output (result) vector
//             
// RETURNS:  none
//**********************************************************************************
template<typename TypePack>
void GProjectionFilter<TypePack>::apply_impl(const Time &t, State &u, State &utmp, State &uo) 
{

  GINT             is, nstate;
  GSIZET           ibeg, iend; // beg, end indices in global array
  GTVector<Ftype>  tmp;
  typename TypePack::GElemList       *gelems=&grid_->elems();

  if ( !bInit_ ) init();

  nstate = traits_.istate.size() == 0 ? u.size() 
         : traits_.istate.size();
  for ( auto j=0; j<nstate; j++ ) { // over required states
    is = traits_.istate.size() == 0 ? j
       : traits_.istate[j];

    for ( auto e=0; e<gelems->size(); e++ ) {
      ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
      u [is]->range(ibeg, iend); // restrict global vecs to local range
      uo[is]->range(ibeg, iend); 
#if defined(_G_IS2D)
      GMTK::D2_X_D1<Ftype>(F_[0], FT_[1], *u[is], tmp_, *uo[is]);
#elif defined(_G_IS3D)
      GMTK::D3_X_D2_X_D1<Ftype>(F_[0], FT_[1], FT_[2],  *u[is], tmp_, *uo[is]);
#endif

    }
    u [is]->range_reset(); 
    uo[is]->range_reset(); 
  }

} // end of method apply_impl


//**********************************************************************************
//**********************************************************************************
// METHOD : apply_impl
// DESC   : In-place application of this filter to input vector
//           
// ARGS   : t   : Time
//          u   : input vector field
//          utmp: array of tmp arrays
//             
// RETURNS:  none
//**********************************************************************************
template<typename TypePack>
void GProjectionFilter<TypePack>::apply_impl(const Time &t, State &u, State &utmp) 
{
  assert( utmp.size() >= u.size()
       && "Insufficient temp space provided");

  GINT  is, nstate;
  State unew(u.size());

  nstate = traits_.istate.size() == 0 ? u.size() 
         : traits_.istate.size();

  for ( auto j=0; j<u.size(); j++ ) unew[j] = utmp[u.size()+j];
  apply_impl(t, u, utmp, unew);                

  // Deep copy filtered components back to u:
  for ( auto j=0; j<nstate; j++ ) {
     is = traits_.istate.size() == 0 ? j
        : traits_.istate[j];
    *u[is] = *unew[is];
  }

} // end of method apply_impl


//**********************************************************************************
//**********************************************************************************
// METHOD : init
// DESC   : Initilize operators
// ARGS   : none.
// RETURNS:  none
//**********************************************************************************
template<typename TypePack>
void GProjectionFilter<TypePack>::init()
{
  GINT             nnodes;
  Ftype            a, b, xi, xf0, xN;
  GTVector<GINT>   Nhi(GDIM), Nlow(GDIM);
  GTVector<Ftype>  xihi, xilow;
  GTVector<GNBasis<GCTYPE,Ftype>*> *bhi;
  GTVector<GLLBasis<GCTYPE,Ftype>> blow(GDIM); 
  GTMatrix<Ftype>  Id, Ihi, Ilow, M;
  typename TypePack::GElemList    *gelems=&grid_->elems();

  // For now, let's assume this filter only works
  // when order is constant among elements.

  F_.resize(GDIM);
  FT_.resize(GDIM);

  // First, compute I_M^N I_N^M;
  bhi   = &(*gelems)[0]->gbasis();
  Nhi   = (*gelems)[0]->size();
  for ( auto j=0; j<GDIM; j++ ) { // allocate matrices
    // Limit new p to be in [pold-1, pold/2]:
    assert(traits_.pdelta[j] >= 1 && traits_.pdelta[j] < (Nhi[j]-1)/2); 
    (*bhi)[j]->getXiNodes(xihi);
    Nlow[j]   = Nhi[j] - traits_.pdelta[j];
    blow[j]   .resize(Nlow[j]-1);      // create low order bases
    blow[j]   .getXiNodes(xilow);
    Ilow      .resize(Nlow[j],Nhi[j]); // interp to low order basis
    Ihi       .resize(Nhi[j],Nlow[j]); // interp to high order basis
    F_  [j]   .resize(Nhi[j],Nhi[j]);  // 1d filter matrices
    FT_ [j]   .resize(Nhi[j],Nhi[j]);  // 1d filter matrix transposes
    M         .resize(Nhi[j],Nhi[j]);  // tmp matrix
    Id        .resize(Nhi[j],Nhi[j]); Id.createIdentity();
    (*bhi)[j]->evalBasis(xilow,Ilow); // create Ilow
    blow[j].evalBasis(xihi,Ihi);     // create Ihi

    // Compute 1d filter matrices: F = alpha Ihi Ilow + (1-alpha) I;
    M        = Ihi * Ilow;
    F_  [j]  = M * traits_.alpha[j] ;
    F_  [j] += ( Id * (1.0-traits_.alpha[j]) );
    F_  [j]  .transpose(FT_[j]);

  }

  bInit_ = TRUE;

} // end of method init

