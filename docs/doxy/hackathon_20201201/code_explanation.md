# Hackathon Source Code Explanation

## Relevant Files

### [GeoFLOW/src/cdg/grid/ggrid.cpp](../../../src/cdg/grid/ggrid.cpp)

#### compute_grefderiv
Top level entry point for derivative calculations is the function
**compute_grefderiv**.  This function contains two branches depending on
which approach to calculating the derivatives we would like to take. 

* grefderiv_varp - When the order of polynomial can "vary" between elements 
* grefderiv_constp - When each element has a "contant" polynomial order

They both call the same lower level matrix multiplication routines but
the constant polynomial simplification (and restriction) enables a single 
large matrix multiply to occur.  The variable polynomial requires numerous
smaller matrix multiply operations which might lend itself to a batch
matrix multiply operation.  

```cpp
void GGrid::compute_grefderiv(GTVector<GFTYPE> &u, GTVector<GFTYPE> &etmp,
                              GINT idir, GBOOL dotrans, GTVector<GFTYPE> &du)
{
  switch ( gderivtype_ ) {
    case GDV_VARP:
      grefderiv_varp(u, etmp, idir, dotrans, du);
      break;
    case GDV_CONSTP:
      grefderiv_constp (u, etmp, idir, dotrans, du);
      break;
    default:
      assert(false);
  }
}
```

#### grefderiv_constp
The constant polynomial routines contain both 2D and 3D code but for simplicity only 
a portion of the 2D code will be displayed and examined here. As can be seen, the code 
branches depending on the direction of the derivative (idir) and a single matrix multiply 
is performed (GMTK::I2_X_D1 or GMTK::D2_X_I1) for the large matrix enabled by the constant 
polynomial order of the elements.  

```cpp
void GGrid::grefderiv_constp(GTVector<GFTYPE> &u, GTVector<GFTYPE> &etmp,
                            GINT idir, GBOOL dotrans, GTVector<GFTYPE> &du)
{
  GSIZET               ibeg, iend, Ne,  NN; // beg, end indices for global array
  GTVector<GSIZET>     N(GDIM);
  GTMatrix<GFTYPE>    *Di;         // element-based 1d derivative operators
  GElemList           *gelems = &this->elems();


  for ( auto k=0; k<GDIM; k++ ) N[k]= (*gelems)[0]->size(k);
  Ne = gelems->size();

  switch (idir) {
  case 1:
    Di = (*gelems)[0]->gbasis(0)->getDerivMatrix (dotrans);
    GMTK::I2_X_D1(*Di, u, N[0], N[1], Ne,  du); 
    break;
  case 2:
    Di = (*gelems)[0]->gbasis(1)->getDerivMatrix(!dotrans);
    GMTK::D2_X_I1(*Di, u, N[0], N[1], Ne, du); 
    break;
  default:
    assert(FALSE && "Invalid coordinate direction");
  }
  
  ... More Stuff ...
}
```

#### grefderiv_varp
The variable polynomial derivatives are more complex since each element can 
contain a different number of solutions. This requires a loop over the elements 
and each matrix multiply to be performed individually. It is important to notice 
both derivative calculations (constant and variable) use the same matrix multipy 
routines (GMTK::I2_X_D1 or GMTK::D2_X_I1) but the number of inputs can be different 
and the range over the solution vectors is modified to reflect the individual 
element sizes.

```cpp
void GGrid::grefderiv_varp(GTVector<GFTYPE> &u, GTVector<GFTYPE> &etmp,
                             GINT idir, GBOOL dotrans, GTVector<GFTYPE> &du)
{
  GSIZET               ibeg, iend; // beg, end indices for global array
  GTVector<GSIZET>     N(GDIM);
  GTMatrix<GFTYPE>    *Di;         // element-based 1d derivative operators
  GElemList           *gelems = &this->elems();


#if defined(_G_IS2D)
  switch (idir) {
  case 1:
    for ( auto e=0; e<gelems->size(); e++ ) {
      ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
      u.range(ibeg, iend); // restrict global vecs to local range
      du.range(ibeg, iend);
      for ( auto k=0; k<GDIM; k++ ) N[k]= (*gelems)[e]->size(k);
      Di = (*gelems)[e]->gbasis(0)->getDerivMatrix (dotrans);
      GMTK::I2_X_D1(*Di, u, N[0], N[1], du); 
    }
    break;
  case 2:
    for ( auto e=0; e<gelems->size(); e++ ) {
      ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
      u.range(ibeg, iend); // restrict global vecs to local range
      du.range(ibeg, iend);
      for ( auto k=0; k<GDIM; k++ ) N[k]= (*gelems)[e]->size(k);
      Di = (*gelems)[e]->gbasis(1)->getDerivMatrix(!dotrans);
      GMTK::D2_X_I1(*Di, u, N[0], N[1], du); 
    }
    break;
  case 3:
    assert( GDIM == 3
         && "Only GDIM reference derivatives");
    du = 0.0; //u;
    break;
  default:
    assert(FALSE && "Invalid coordinate direction");
  }
  
  ... More Stuff ...
}
```

### [GeoFLOW/src/cdg/blas/gmtk.cpp](../../../src/cdg/blas/gmtk.ipp)

#### I2\_X\_D1 (5 & 6 Arguments Nearly Same)
The matrix multiply routines forward the call to BLAS like routine to perform the operations 
required.  Since a templated function such as I2_X_D1 can take matricies and arrays of 
any data type we must proved 3 branches to select the proper call to interface with a 
library such as BLAS which are not templated.

```cpp
template<typename T>
void I2_X_D1(GTMatrix<T> &D1,
             GTVector<T> &u, GSIZET N1, GSIZET N2, GTVector<T> &y)
{
  GSIZET ND1, ND2;

  ND1 = D1.size(1);
  ND2 = D1.size(2);

  // Compute y = I2_X_D1 u:
  if      ( std::is_same<T,GFLOAT>::value ) {
    fmxm((GFLOAT*)y.data(), (GFLOAT*)(D1.data().data()), &ND1, &ND2, (GFLOAT*)(u.data()), &N1, &N2, &szMatCache_);
  }
  else if ( std::is_same<T,GDOUBLE>::value ) {
    dmxm((GDOUBLE*)y.data(), (GDOUBLE*)(D1.data().data()), &ND1, &ND2, (GDOUBLE*)(u.data()), &N1, &N2, &szMatCache_);
  }
  else if ( std::is_same<T,GQUAD>::value ) {
    qmxm((GQUAD*)y.data(), (GQUAD*)(D1.data().data()), &ND1, &ND2, (GQUAD*)(u.data()), &N1, &N2, &szMatCache_);
  }
  else {
    assert(FALSE);
  }

} // end of method I2_X_D1 (1)
```

#### D2\_X\_I1 (5 Argument)

```cpp
template<typename T>
void D2_X_I1(GTMatrix<T> &D2T, 
              GTVector<T> &u, GSIZET N1, GSIZET N2, GTVector<T> &y)
{
  GSIZET N21, N22;

  N21 = D2T.size(1);
  N22 = D2T.size(2);

  // Compute y = I2_X_D1 u = u * D2T:
  if      ( std::is_same<T,GFLOAT>::value ) {
    fmxm((GFLOAT*)y.data(), (GFLOAT*)u.data(), &N1, &N2, (GFLOAT*)D2T.data().data(), &N21, &N22, &szMatCache_);
  }
  else if ( std::is_same<T,GDOUBLE>::value ) {
    dmxm((GDOUBLE*)y.data(), (GDOUBLE*)u.data(), &N1, &N2, (GDOUBLE*)D2T.data().data(), &N21, &N22, &szMatCache_);
  }
  else if ( std::is_same<T,GQUAD>::value ) {
    qmxm((GQUAD*)y.data(), (GQUAD*)u.data(), &N1, &N2, (GQUAD*)D2T.data().data(), &N21, &N22, &szMatCache_);
  }
  else {
    assert(FALSE);
  }

} // end of method D2_X_I1 (1)
```


#### D2\_X\_I1 (6 Argument)

```cpp
template<typename T>
void D2_X_I1(GTMatrix<T> &D2T, 
              GTVector<T> &u, GSIZET N1, GSIZET N2, GSIZET Ne, GTVector<T> &y)
{
  GSIZET N21, N22, Nu;

  N21 = D2T.size(1);
  N22 = D2T.size(2);
  Nu  = N1 * N2;

  // Compute y = I2_X_D1 u = u * D2T:
  if      ( std::is_same<T,GFLOAT>::value ) {
    for ( auto i=0; i<Ne; i++ ) {
      fmxm((GFLOAT*)y.data()+i*Nu, (GFLOAT*)u.data()+i*Nu, &N1, &Nu, (GFLOAT*)D2T.data().data(), &N21, &N22, &szMatCache_);
    }
  }
  else if ( std::is_same<T,GDOUBLE>::value ) {
    for ( auto i=0; i<Ne; i++ ) {
      dmxm((GDOUBLE*)y.data()+i*Nu, (GDOUBLE*)u.data()+i*Nu, &N1, &Nu, (GDOUBLE*)D2T.data().data(), &N21, &N22, &szMatCache_);
    }
  }
  else if ( std::is_same<T,GQUAD>::value ) {
    for ( auto i=0; i<Ne; i++ ) {
      qmxm((GQUAD*)y.data()+i*Nu, (GQUAD*)u.data()+i*Nu, &N1, &Nu, (GQUAD*)D2T.data().data(), &N21, &N22, &szMatCache_);
    }
  }
  else {
    assert(FALSE);
  }
} // end of method D2_X_I1 (2)
```


 