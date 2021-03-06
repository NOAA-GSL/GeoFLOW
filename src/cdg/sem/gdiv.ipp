//==================================================================================
// Module       : gdiv.ipp
// Date         : 09/05/20 (DLR)
// Description  : Represents the SEM discretization of the divergence
//                operator,
//                      Div(rho \vec{v})
//                where rho is a scalar field, and  \vec{v} is
//                a vector field. This isn't the strictly conservative
//                form that uses Gauss' theorem to add bdy fluxes; 
//                it is volume-integrated.
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
template<typename Types>
GDivOp<Types>::GDivOp(Traits &traits, Grid &grid)
:
traits_       (traits),
grid_         (&grid),
massop_       (&grid.massop())
{
  assert(grid_->ntype().multiplicity(0) == GE_MAX-1 
        && "Only a single element type allowed on grid");

  gadvect_ = new GAdvect<Types>(*grid_);
} // end of constructor method (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor method 
// DESC   : 
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<typename Types>
GDivOp<Types>::~GDivOp()
{
  if ( gadvect_ != NULLPTR ) delete gadvect_;
} // end, destructor


//**********************************************************************************
//**********************************************************************************
// METHOD : apply (1)
// DESC   : Compute application of this now nonlinear operator to scalar
//          field:
//            div = Div (d \vec{v})
//          Remember, normally, this discretization would be multiplied by
//          -1 to represent this operation. We do not apply this sign here.
// ARGS   : d   : scalar field
//          u   : input vector field
//          utmp: array of tmp arrays. 3 arrays required.
//          div : output (result) 
//          ivec: if d is a vector component, specifies which component. 
//                Default is -1 (meaning, a scalar).
//             
// RETURNS:  none
//**********************************************************************************
template<typename Types>
void GDivOp<Types>::apply(StateComp &d, State &u, State &utmp, StateComp &div, GINT ivec) 
{

  assert( utmp.size() >= 3
       && "Insufficient temp space specified");

  GBOOL      usebdy = grid_->usebdydata();
  GINT       nxy = grid_->gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;
  GSIZET                     k;
  GTVector<GSIZET>          *ieface  = &grid_->gieface();
  GTVector<GTVector<Ftype>> *normals = &grid_->faceNormals(); 
  StateComp                 *mass    =  grid_->massop().data();
  StateComp                 *bmass   = &grid_->faceMass();

  typename Grid::BinnedBdyIndex *igb = &grid_->igbdy_binned();


  if ( !traits_.docollocation ) {
    // div = -D^{T,j} ( d u_j MJ )
    //     + bdy terms:
    div  = 0.0;
    for ( auto j=0; j<nxy; j++ ) { 
       grid_->dealias(d, *u[j], *utmp[1]);
       grid_->wderiv(*utmp[1], j+1, TRUE, *utmp[0], *utmp[2]);
#if 0 // defined(GEOFLOW_USE_NEUMANN_HACK)
if ( ivec == -1 || ivec == 2 ) {
GMTK::zero<Ftype>(div,(*igb)[1][GBDY_0FLUX]);
GMTK::zero<Ftype>(div,(*igb)[3][GBDY_0FLUX]);
}
else if ( ivec == -1 || ivec == 1 ) {
GMTK::zero<Ftype>(*utmp[2],(*igb)[0][GBDY_0FLUX]);
GMTK::zero<Ftype>(*utmp[2],(*igb)[2][GBDY_0FLUX]);
}
#endif
       div -= *utmp[2];

      // Global bdy terms:
      //  Note: utmp[1] should contain effect of bdy conditions,
      //        which are imposed on entry
      for ( auto b=0; usebdy && b<ieface->size(); b++ ) {
        k = (*ieface)[b];
        div[k] += (*utmp[1])[k] * (*normals)[j][b] * (*bmass)[b];
      }
    }
  }
  else {
#if 1
    // Do collocation form:
    // div = MJ d/dx_j ( d u_j ):
//  d.pointProd(*u[0], *utmp[1]);
    grid_->dealias(d, *u[0], *utmp[1]);
    grid_->deriv(*utmp[1], 1, *utmp[0], div);
#if 0 // defined(GEOFLOW_USE_NEUMANN_HACK)
if ( ivec == -1 || ivec == 2 ) {
GMTK::zero<Ftype>(div,(*igb)[1][GBDY_0FLUX]);
GMTK::zero<Ftype>(div,(*igb)[3][GBDY_0FLUX]);
}
#endif
    for ( auto j=1; j<nxy; j++ ) { 
       grid_->dealias(d, *u[j], *utmp[1]);
       grid_->deriv(*utmp[1], j+1, *utmp[0], *utmp[2]);
#if 0 //defined(GEOFLOW_USE_NEUMANN_HACK)
if ( ivec == -1 || ivec == 1 ) {
GMTK::zero<Ftype>(*utmp[2],(*igb)[0][GBDY_0FLUX]);
GMTK::zero<Ftype>(*utmp[2],(*igb)[2][GBDY_0FLUX]);
}
#endif
       div += *utmp[2];
    }
    div *= *(massop_->data());

#else
    // Do collocation form: Div (dv) = p Div v + v.Grad p:

    gadvect_->apply(d, u, utmp, div, ivec); // v.Grad d
    *utmp[2] = 0.0;
    for ( auto j=0; j<nxy; j++ ) {  // d Div v
       grid_->deriv(*u[j], j+1, *utmp[0], *utmp[1]);
       *utmp[2] += *utmp[1];
    }
    utmp[2]->pointProd(d);
    utmp[2]->pointProd(*(massop_->data()));
    div += *utmp[2];
#endif

  }


} // end of method apply (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : apply (2)
// DESC   : Compute application of this operator to pure vector 
//          field:
//            div = Div (\vec{v})
//          Remember, normally, this discretization would be multiplied by
//          -1 to represent this operation. We do not apply this sign here.
// ARGS   : u   : input vector field
//          utmp: array of tmp arrays. 2 arrays required.
//          div : output (result) 
//          ivec: if d is a vector component, specifies which component. 
//                Default is -1 (meaning, a complete scalar).
//             
// RETURNS:  none
//**********************************************************************************
template<typename Types>
void GDivOp<Types>::apply(State &u, State &utmp, StateComp &div, GINT ivec) 
{

  assert( utmp.size() >= 2
       && "Insufficient temp space specified");

  GBOOL      usebdy = grid_->usebdydata();
  GINT       nxy = grid_->gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;
  GSIZET     k;
  GTVector<GSIZET>          *ieface  = &grid_->gieface() ;
  GTVector<GTVector<Ftype>> *normals = &grid_->faceNormals(); 
  StateComp                 *mass    =  grid_->massop().data();
  StateComp                 *bmass   = &grid_->faceMass();


  if ( !traits_.docollocation ) {

    // div = -D^{T,j} ( u_j ) 
    //     + bdy. terms:
    div = 0.0;
    for ( auto j=0; j<nxy; j++ ) { 
       grid_->wderiv(*u[j], j+1, TRUE, *utmp[0], *utmp[1]);
       div -= *utmp[1];

      for ( auto b=0; usebdy && b<ieface->size(); b++ ) {
        k = (*ieface)[b];
//cout << "GDiv:: apply(2): k=" << k << " b=" << b << " nx=" << (*normals)[0][b] << " ny=" << (*normals)[1][b] << " bmass=" << (*bmass)[b] << endl;
        div[k] += (*u[j])[k] * (*normals)[j][b] * (*bmass)[b];
      }
    }
  }
  else {

    // Do collocation form:
    grid_->deriv(*u[0], 1, *utmp[0], div);
    for ( auto j=1; j<nxy; j++ ) { 
      grid_->deriv(*u[j], j+1, *utmp[0], *utmp[1]);
      div += *utmp[1];
    }
    div.pointProd(*mass);

  }

} // end of method apply (2)


