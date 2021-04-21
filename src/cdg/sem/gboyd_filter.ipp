//==================================================================================
// Module       : gboyd_filter.hpp
// Date         : 9/14/20 (DLR)
// Description  : Computes the Boyd filter to diminish aliasing errors.
//                Taken from Giraldo & Rosemont 2004, MWR:132 133:
//                    u <-- F u
//                where
//                    F = L Lambda L^-1; s.t.
//                and 
//                    Lambda = 1 if i< ifilter
//                             mu [(i-ifilter)/(N - ifilter)]^2, i>= ifilter.
//                L is the Legendre transform matrix:
//                    L = | P_0(xi0), P_1(xi0) ... P_i(xi0)-P_{i-2)(xi0) ... |
//                        | P_0(xi1), P_1(xi1) ... P_i(xi1)-P_{i-2)(xi1) ... |
//                        |  ...                                             |.
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
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
GBoydFilter<TypePack>::GBoydFilter(Traits &traits, Grid &grid)
:
bInit_               (FALSE),
traits_              (traits),
grid_                (&grid)
{
  assert(grid_->ntype().multiplicity(0) == GE_MAX-1 
        && "Only a single element type allowed on grid");
  assert(traits_.pdelta.size() >= GDIM && traits_.strength.size() >= GDIM);
  auto emax = *max_element(std::begin(traits_.strength), std::end(traits_.strength));
  auto emin = *min_element(std::begin(traits_.strength), std::end(traits_.strength));
  assert(emin >= 0 && emax <= 1);

} // end of constructor method (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor method 
// DESC   : 
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<typename TypePack>
GBoydFilter<TypePack>::~GBoydFilter()
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
void GBoydFilter<TypePack>::apply_impl(const Time &t, State &u, State &utmp, State &uo) 
{

  GINT             is, nstate;
  GSIZET           ibeg, iend; // beg, end indices in global array
  GTMatrix<Ftype> *F[GDIM];
  typename TypePack::GElemList       *gelems=&grid_->elems();

  if ( !bInit_ ) init();

  nstate = traits_.istate.size() == 0 ? u.size() 
         : traits_.istate.size();
  for ( auto j=0; j<nstate; j++ ) {
    is = traits_.istate.size() == 0 ? j
       : traits_.istate[j];
    for ( auto e=0; e<gelems->size(); e++ ) {
      ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
      u[is]->range(ibeg, iend); // restrict global vecs to local range
      uo[is]->range(ibeg, iend); 
      F [0] = (*gelems)[e]->gbasis(0)->getFilterMat();
      F [1] = (*gelems)[e]->gbasis(1)->getFilterMat(TRUE);
#if defined(_G_IS2D)
      GMTK::D2_X_D1<GFTYPE>(*F[0], *F[1], *u[is], tmp_, *uo[is]);
#elif defined(_G_IS3D)
      F [2] = (*gelems)[e]->gbasis(2)->getFilterMat(TRUE);
      GMTK::D3_X_D2_X_D1<GFTYPE>(*F[0], *F[1], *F[2], *u[is], tmp_, *uo[is]);
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
void GBoydFilter<TypePack>::apply_impl(const Time &t, State &u, State &utmp) 
{

  assert( utmp.size() >= 2*u.size()
       && "Insufficient temp space provided");

  State unew(u.size());
  
  for ( auto j=0; j<u.size(); j++ ) unew[j] = utmp[u.size()+j];
  apply_impl(t, u, utmp, unew); 

  for ( auto j=0; j<u.size(); j++ ) *u[j] = *unew[j];

} // end of method apply_impl


//**********************************************************************************
//**********************************************************************************
// METHOD : init
// DESC   : Initilize operators
// ARGS   : none.
// RETURNS:  none
//**********************************************************************************
template<typename TypePack>
void GBoydFilter<TypePack>::init()
{

  GINT             ifilter, nnodes;
  GTVector<GNBasis<GCTYPE,Ftype>*>
                   ipool;
  GTMatrix<Ftype> *F, *FT, *iL, *L, Lambda;
  GTMatrix<Ftype> tmp;
  typename TypePack::GElemList       *gelems=&grid_->elems();

  // Build the filter matrix, F, and store within 
  // each basis object for later use:
  //   F = L Lambda L^-1; s.t.
  // u <-- F u
  // where
  //   Lambda = 1 if i< ifilter
  //            mu [(i-ifilter)/(N - ifilter)]^2, i>= ifilter

  for ( auto e=0; e<gelems->size(); e++ ) {
    for ( auto k=0; k<GDIM; k++ ) {
//    (*gelems)[e]->gbasis(k)->computeLegTransform(ifilter_); 
      F  = (*gelems)[e]->gbasis(k)->getFilterMat(); // storage for filter
      FT = (*gelems)[e]->gbasis(k)->getFilterMat(TRUE); // transpose 
      nnodes = (*gelems)[e]->gbasis(k)->getOrder()+1;
      if ( ipool.contains((*gelems)[e]->gbasis(k)) ) continue;
      ipool.push_back((*gelems)[e]->gbasis(k));
      ifilter = nnodes - traits_.pdelta[k] - 1;
      assert(ifilter > 0 && ifilter < nnodes);
      Lambda.resize(nnodes,nnodes); Lambda = 0.0;
      tmp   .resize(nnodes,nnodes);
      
      L      = (*gelems)[e]->gbasis(k)->getLegTransform();
      iL     = (*gelems)[e]->gbasis(k)->getiLegTransform();

      for ( auto i=0; i<nnodes; i++ ) { // build weight matrix, Lambda
        Lambda(i,i) = 1.0;
        if ( i >= ifilter ) {
//cout << " ..................... i=" << i << " ifilter=" << ifilter << endl;
#if 1
          Lambda(i,i) = traits_.strength[k] 
                      * ( 1.0 -  pow( fabs( (Ftype)(i-ifilter) / ( (Ftype)(nnodes-ifilter) ) ), 1.0) ) 
                      + 1 - traits_.strength[k];
#else
          Lambda(i,i) = traits_.strength[k] 
                      * pow( fabs( (Ftype)(i-ifilter) / ( (Ftype)(nnodes-ifilter) ) ), 2.0); 
#endif
        }
      } // end, node/mode loop 
      tmp = Lambda * (*iL);
     *F   = (*L) * tmp;
      F   ->transpose(*FT);
//cout << "F=" << *F << endl;
    } // end, k-loop
  } // end, element loop

  bInit_ = TRUE;

} // end of method init


