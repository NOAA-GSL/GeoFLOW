//==================================================================================
// Module       : gutils.cpp
// Date         : 1/31/19 (DLR)
// Description  : GeoFLOW utilities namespace
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#include "gutils.hpp"
#include "gcomm.hpp"

namespace geoflow
{

//**********************************************************************************
//**********************************************************************************
// METHOD : str2bdytype
// DESC   : Convert string to GBdyType
// ARGS   : stype: string type
// RETURNS: GBdyType
//**********************************************************************************
GBdyType str2bdytype(const GString &stype)
{
  GString s0;
  for ( auto j=0; j<GBDY_NONE; j++ ) {
    s0 = sGBdyType[j];
    if ( stype.compare(s0) == 0 ) return static_cast<GBdyType>(j);
  }
  assert(FALSE && "Invalid boundary type");
} // end, str2bdytype



//**********************************************************************************
//**********************************************************************************
// METHOD : str2comptype
// DESC   : Convert string to GStateCompType
// ARGS   : stype: string type
// RETURNS: GStateCompType
//**********************************************************************************
GStateCompType str2comptype(const GString &stype)
{
  GString s0;
  for ( auto j=0; j<GBDY_NONE; j++ ) {
    s0 = sGStateCompType[j];
    if ( stype.compare(s0) == 0 ) return static_cast<GStateCompType>(j);
  }
  assert(FALSE && "Invalid state component type");
} // end, str2comptype


//**********************************************************************************
//**********************************************************************************
// METHOD : file_empty
// DESC   : Determine if specified file exists or is empty
// ARGS   : sfile: file name
// RETURNS: GBOOL
//**********************************************************************************
GBOOL file_empty(GString sfile)
{

  GBOOL         bret=FALSE;
  std::ifstream itst;

  GComm::Synch();
  itst.open(sfile);
  bret = itst.peek() == std::ofstream::traits_type::eof();
  itst.close();

  return bret;

} // end, method file_empty



} // end, namespace

