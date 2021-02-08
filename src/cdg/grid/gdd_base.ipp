//==================================================================================
// Module       : gdd_base.ipp
// Date         : 8/28/18 (DLR)
// Description  : Forms virtual abstract base class for domain 
//                decomposition objects. Default behaviour is to
//                do a simple partition that divides the number of 
//                elememt representations by the number of MPI tasks.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#include "gcomm.hpp"
#include "gdd_base.hpp"

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Instantiate with radius, refinement level, directional basis, and
//          domain decompoisition object. This generates a 2d grid.
// ARGS   : radius: radius
//          ilevel: refinement level. Level 0 is just the icos generator 
//          b     : vector of basis pointers, of size at least GDIM. 
// RETURNS: none
//**********************************************************************************
template<typename T>
GDD_base<T>::GDD_base(GINT nprocs)
:
nprocs_  (nprocs)
{
  GEOFLOW_TRACE();
} // end of constructor method (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : Copy constructor method 
// DESC   : 
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<typename T>
GDD_base<T>::GDD_base(const GDD_base &obj)
{
  GEOFLOW_TRACE();
  nprocs_ = obj.nprocs_;
} // end, destructor


//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor method 
// DESC   : 
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<typename T>
GDD_base<T>::~GDD_base()
{
  GEOFLOW_TRACE();
} // end, destructor


//**********************************************************************************
//**********************************************************************************
// METHOD : doDD (1)
// DESC   : Do simple distribution among tasks by dividing element
//          representations among tasks as evenly as possible.
//          
// ARGS   : x    : coords (e.g. centroids) representing position of 
//                 each element. All tasks see the same array, so it's 'global'
//          irank: MPI task whose elements are requested
//          iret : indirection indices into x that give the elements to 
//                 be 'ownded' by task irank. Size will be set here.
// RETURNS: number of elements belonging to rank irank.
//**********************************************************************************
template<typename T>
GSIZET GDD_base<T>::doDD(const GTVector<GTVector<T>>&x, GINT irank, GTVector<GINT> &iret)
{
  GEOFLOW_TRACE();
  GString serr = "GDD_base<T>::doDD (1): ";
  assert(irank >=0 && irank < nprocs_ && "Invalid rank");

  GSIZET nxp  = x[0].size() / nprocs_;
  GSIZET nrem = x[0].size() % nprocs_;
  GSIZET ibeg, iend;

  // Divide any remainder among tasks:

  if ( irank < nrem ) {
    ibeg = irank * (nxp + 1);
    iend = ibeg + nxp;
    iret.resize(nxp+1);
  }
  else {
    ibeg = nrem * (nxp + 1) + (irank - nrem) * nxp;
    iend = ibeg + nxp - 1;
    iret.resize(nxp);
  }
  assert((iend-ibeg+1) == iret.size() && "Incompatible number of elements");
  for ( GSIZET j=0; j<iret.size(); j++ ) iret[j] = ibeg + j;
  
  return iret.size();
  
} // end of method doDD (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : doDD (2)
// DESC   : Do simple distribution among tasks by dividing element
//          representations among tasks as evenly as possible.
//          
// ARGS   : x    : points (e.g. centroids) representing position of each element
//          irank: MPI task whose elements are requested
//          iret : indirection indices into x that give the elements to 
//                 be 'ownded' by task irank.
// RETURNS: number of elements belonging to rank irank.
//**********************************************************************************
template<typename T>
GSIZET GDD_base<T>::doDD(const GTVector<GTPoint<T>> &x, GINT irank, GTVector<GINT> &iret)
{
  GEOFLOW_TRACE();
  GString serr = "GDD_base<T>::doDD (2): ";
  assert(irank >=0 && irank < nprocs_ && "Invalid rank");

  GSIZET nxp  = x.size() / nprocs_;
  GSIZET nrem = x.size() % nprocs_;
  GSIZET ibeg, iend;


  // Divide any remainder among tasks:
  
  if ( irank < nrem ) {
    ibeg = irank * (nxp + 1);
    iend = ibeg + nxp;
    iret.resize(nxp+1);
  }
  else {
    ibeg = nrem * (nxp + 1) + (irank - nrem) * nxp;
    iend = ibeg + nxp - 1;
    iret.resize(nxp);
  }
  assert((iend-ibeg+1) == iret.size() && "Incompatible number of elements");
  for ( GSIZET j=0; j<iret.size(); j++ ) iret[j] = ibeg + j;
  
  return iret.size();
  
} // end of method doDD (2)


