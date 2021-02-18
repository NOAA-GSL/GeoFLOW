//==================================================================================
// Module       : gutils.cpp
// Date         : 1/31/19 (DLR)
// Description  : GeoFLOW utilities namespace
// Copyright    : Copyright 2019. Colorado State University. All rights reserved.
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
  for ( auto j=0; j<GBDY_MAX; j++ ) {
    s0 = sGBdyType[j];
    if ( stype.compare(s0) == 0 ) return static_cast<GBdyType>(j);
  }
  assert(FALSE && "Invalid boundary type");
  return static_cast<GBdyType>(0);
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
  return static_cast<GStateCompType>(0);
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


//**********************************************************************************
//**********************************************************************************
// METHOD : get_bdy_block
// DESC   : Get bdy config block info
// ARGS   : sptree : prop tree for the block
//          stblock: stBdyBlock containing return info
// RETURNS: GBOOL
//**********************************************************************************
void get_bdy_block(const geoflow::tbox::PropertyTree &sptree, stBdyBlock &stblock)
{
  GBOOL                bvalreq; // Dirichlet value vec required?
  GString              sconf;
  std::vector<GString> svec;
  std::vector<std::vector<GINT>> 
                       ivecvec;
  std::vector<std::vector<GFTYPE>> 
                       fvecvec;
  
  // Clear out structure data:
  stblock.istate.clear();
  stblock.value .clear();
  stblock.tbdy  .clear();
  stblock.config_method.clear();
  stblock.inflow_method.clear();

  // Get bdy block data:
  
#if 0
  if ( !sptree.isArray<GString>("base_type") ) {
    cout << "GUtils::get_bdy_block: base_type array is missing" << endl;
    assert(FALSE); 
  }
  if ( !sptree.isArray<GINT>("istate") ) {
    cout << "GUtils::get_bdy_block: istate array is missing" << endl;
    assert(FALSE); 
  }
#endif
  svec = sptree.getArray<GString>("base_type");
  bvalreq = FALSE;
  for ( auto j=0; j<svec.size(); j++ ) {
    bvalreq = bvalreq || svec[j] == "GBDY_DIRICHLET";
  }
  if ( bvalreq && !sptree.keyExists("value") )  {
    cout << "GUtils::get_bdy_block: Dirichlet 2D 'value' array is missing" << endl;
    assert(FALSE); 
  }

  ivecvec = sptree.getArray2D<GINT>("istate");
  if ( bvalreq && !sptree.isArray2D<GFTYPE>("value") ) {
    cout << "GUtils::get_bdy_block: Dirichlet 'value' array is not a 2D array" << endl;
    assert(FALSE); 
  }

  if ( bvalreq ) {
    fvecvec = sptree.getArray2D<GFTYPE>("value");
  }
  else { // make arrays of len 0:
    fvecvec.resize(svec.size());
  }
  
  if ( ivecvec.size() != svec.size() ) {
    cout << "GUtils::get_bdy_block: istate and base_type array sizes not equal" << endl;
    assert(FALSE); 
  }
  if ( fvecvec.size() != svec.size() ) {
    cout << "GUtils::get_bdy_block: value and base_type array sizes not equal" << endl;
    assert(FALSE); 
  }

  stblock.tbdy.resize(svec.size());
  for ( auto j=0; j<svec.size(); j++ ) {
    stblock.tbdy[j] = geoflow::str2bdytype(svec[j]);
  }

  stblock.istate.resize(svec.size());
  stblock.value .resize(svec.size());
  for ( auto j=0; j<svec.size(); j++ ) {
    stblock.istate[j].resize(ivecvec[j].size());
    stblock.value [j].resize(fvecvec[j].size());
    stblock.istate[j] = ivecvec[j];
    stblock.value [j] = fvecvec[j];
  }


} // end, method get_bdy_block



} // end, namespace geoflow

