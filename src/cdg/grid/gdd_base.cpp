//==================================================================================
// Module       : gdd_default.cpp
// Date         : 8/28/18 (DLR)
// Description  : Forms virtual abstract base class for domain 
//                decomposition objects. Default behaviour is to
//                do a simple partition that divides the number of 
//                elememt representations by the number of MPI tasks.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
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
GDD_base::GDD_base(GINT nprocs)
:
nprocs_  (nprocs)
{
} // end of constructor method (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : Copy constructor method 
// DESC   : 
// ARGS   : none
// RETURNS: none
//**********************************************************************************
GDD_base::GDD_base(const GDD_base &obj)
{
  nprocs_ = obj.nprocs_;
} // end, destructor


//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor method 
// DESC   : 
// ARGS   : none
// RETURNS: none
//**********************************************************************************
GDD_base::~GDD_base()
{
} // end, destructor


//**********************************************************************************
//**********************************************************************************
// METHOD : doDD (1)
// DESC   : Do simple distribution among tasks by dividing element
//          representations among tasks as evenly as possible.
//          
// ARGS   : x    : coords (e.g. centroids) representing position of each element
//          irank: MPI task whose elements are requested
//          iret : indirection indices into x that give the elements to 
//                 be 'ownded' by task irank. Size will be set here.
// RETURNS: number of elements belonging to rank irank.
//**********************************************************************************
GSIZET GDD_base::doDD(const GTVector<GTVector<GFTYPE>>&x, GINT irank, GTVector<GINT> &iret)
{
  GString serr = "GDD_base::doDD (1): ";
  assert(irank >=0 && irank < nprocs_ && "Invalid rank");

  GINT nxp  = x[0].size() / nprocs_;
  GINT nrem = x[0].size() % nprocs_;

  iret.resize(irank==nprocs_-1 ? nxp+nrem : nxp);

  // Place any remainder in final task:
  for ( GSIZET j=0; j<iret.size(); j++ ) iret[j] = irank*nxp + j;
  
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
GSIZET GDD_base::doDD(const GTVector<GTPoint<GFTYPE>> &x, GINT irank, GTVector<GINT> &iret)
{
  GString serr = "GDD_base::doDD (2): ";
  assert(irank >=0 && irank < nprocs_ && "Invalid rank");

  GINT nxp  = x.size() / nprocs_;
  GINT nrem = x.size() % nprocs_;

  iret.resize(irank==nprocs_-1 ? nxp+nrem : nxp);

  // Place any remainder in final task:
  for ( GSIZET j=0; j<iret.size(); j++ ) iret[j] = irank*nxp + j;
  
  return iret.size();
  
} // end of method doDD (2)


