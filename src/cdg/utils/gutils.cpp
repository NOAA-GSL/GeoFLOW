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
// METHOD : bdy_block_conform_per
// DESC   : Determine if specified bdy block conforms to requirements 
//          for specifying PERIODIC bcs
//          
// ARGS   : sptree : prop tree for the block
// RETURNS: 
//          0: if bdy condition is ok, but it does not represent a PERIODIC 
//             boundary;
//          1: if sptree represents a valid specification;
//          2: if there's an attempt to represent a PERIODIC bdy,
//             but it's invalid
//          
//**********************************************************************************
GINT bdy_block_conform_per(const geoflow::tbox::PropertyTree &sptree)
{
  GBOOL                battempt;
  GString              bdyclass;
  std::vector<GString> stypes;

  bdyclass = sptree.getValue<GString>("base_class"); // required



  // Get base bdy condition type:
  stypes = sptree.getArray<GString>("base_type"); // is a vector

  battempt = std::find(stypes.begin(), stypes.end(), "GBDY_PERIODIC") != stypes.end();
  if ( !battempt ) return 0;
  
  if ( "uniform" != bdyclass 
   ||  "GBDY_PERIODIC" != stypes[0] 
   ||   stypes.size() != 1 ) return 2;

  return 1;

} // end, method bdy_block_conform_per


//**********************************************************************************
//**********************************************************************************
// METHOD : get_bdy_block
// DESC   : Get bdy config block info
// ARGS   : ptree  : main prop tree 
//          sbdy   : bdy block name
//          ibc    : which bdy condition to retrieve from within 
//                   block; may continue to retrieve until return is FALSE
//          stblock: stBdyBlock containing return info
// RETURNS: TRUE if stblock contains valid data; else FALSE
//**********************************************************************************
GBOOL get_bdy_block(const geoflow::tbox::PropertyTree &ptree, GString &sbdy, GINT ibc, stBdyBlock &stblock)
{
  GBOOL                bvalreq; // Dirichlet value vec required?
  GINT                 nbc, nstate;
  GString              bdyclass, ss;
  std::vector<GBOOL>   bvec;
  std::vector<GINT>    ivec;
  std::vector<GFTYPE>  fvec;
  std::vector<GString> stypes, svec;
  std::vector<std::vector<GINT>> 
                       ivecvec;
  std::vector<std::vector<GFTYPE>> 
                       fvecvec;
  PropertyTree         sptree = ptree.getPropertyTree(sbdy);
  
  // Clear out structure data; initialize:
  stblock.istate   .clear();
  stblock.value    .clear();
  stblock.farfield .clear();
  stblock.falloff  .clear();
  stblock.diffusion.clear();
  stblock.bdyclass .clear();
  stblock.smethod  .clear();

  stblock.use_init    = FALSE;
  stblock.bdyid       = -1;
  stblock.tbdy        = GBDY_NONE;
  stblock.xstart      = 0.0;
  stblock.xmax        = 0.0;
  stblock.bdyid       = -1;
  stblock.tbdy        = GBDY_NONE;
  stblock.bdyclass    = "";
  stblock.smethod     = "";

  // Get bdy block data:
  
  stblock.bdyclass = sptree.getValue<GString>("base_class"); // required

  // Get base bdy condition type:
  stypes = sptree.getArray<GString>("base_type"); // is a vector
  nbc = stypes.size(); // number of bcs specified

  // If ibc out of range, return:
  if ( ibc < 0 || ibc >= nbc ) return FALSE; 
  stblock.tbdy = geoflow::str2bdytype(svec[ibc]); // set bdy type id

  // Get state ids to operate on:
  ivecvec = sptree.getArray2D<GINT>("istate");
  if ( ivecvec.size() != nbc ) {
    cout << "GUtils::get_bdy_block: A 2d JSON vector of size(base_type) must specify vector of state ids for each base_type" << endl;
      cout << "GUtils::get_bdy_block: so 2d JSON vector of size(base_type) must specify vector of istate ids for each bc entry in 'base_type'" << endl;
    assert(FALSE); 
  }
  stblock.istate.resize(ivecvec[ibc].size()); 
  stblock.istate = ivecvec[ibc];
  nstate = stblock.istate.size();
  
  // If DIRICHLET bdy, retrieve vector of values for state ids:
  if ( svec[ibc] == "GBDY_DIRICHLET" ) {
    fvecvec = sptree.getArray2D<GFTYPE>("value");
    if ( fvecvec.size() != nbc ) {
      cout << "GUtils::get_bdy_block: DIRICHLET bc is specified; 2d JSON vector of size(base_type) must specify DIRICHLET values for each istate ('[]' is valid for non-DIRICHLET entries)" << endl;
      assert(FALSE); 
    }
    stblock.value.resize(nstate);
    stblock.value = fvecvec[ibc];
  }
  
  // If INFLOW bdy, retrieve method name, or other data:
  if ( svec[ibc] == "GBDY_INFLOW" ) {
    bvec = sptree.getArray<GBOOL>("use_init");
    if ( bvec.size() != nbc ) {
      cout << "GUtils::get_bdy_block: INFLOW bc is specified; a vector of size(base_type) must specify Boolean 'use_init' flags for each bc entry in 'base_type' " << endl;
      assert(FALSE); 
    }
    stblock.use_init= bvec[ibc];

    if ( !stblock.use_init ) { // get user method if required
      svec = sptree.getArray<GString>("method");
      if ( svec.size() != nbc ) {
        cout << "GUtils::get_bdy_block: INFLOW bc is specified; a vector of size(base_type) must specify inflow methods for each bc entry in 'base_type' (\"\" is valid for non-INFLOW entries)" << endl;
        assert(FALSE); 
      }
      stblock.smethod = svec[ibc];
    }
  }

  // If SPONGE bdy, retrieve required data:
  if ( svec[ibc] == "GBDY_SPONGE" ) { 
    fvecvec = sptree.getArray2D<GFTYPE>("farfield");
    if ( fvecvec.size() != nbc ) {
      cout << "GUtils::get_bdy_block: SPONGE bc is specified; a vector of size(base_type) must specify farfield values for each state group for each bc entry in 'base_type'   ('[]' is valid for non-SPONGE entries)" << endl;
      assert(FALSE); 
    }
    stblock.farfield.resize(nstate);
    stblock.farfield  = fvecvec[ibc];

    fvecvec = sptree.getArray2D<GFTYPE>("falloff");
    if ( fvecvec.size() != nbc ) {
      cout << "GUtils::get_bdy_block: SPONGE bc is specified; a vector of size(base_type) must specify falloff rates for each state for each bc entry in 'base_type'   ('[]' is valid for non-SPONGE entries)" << endl;
      assert(FALSE); 
    }
    stblock.falloff.resize(nstate);
    stblock.falloff = fvecvec[ibc];

    fvecvec = sptree.getArray2D<GFTYPE>("diffusion");
    if ( fvecvec.size() != nbc ) {
      cout << "GUtils::get_bdy_block: SPONGE bc is specified; a vector of size(base_type) must specify diffusion factors for each state for each bc entry in 'base_type'   ('[]' is valid for non-SPONGE entries)" << endl;
      assert(FALSE); 
    }
    stblock.falloff.resize(nstate);
    stblock.falloff = fvecvec[ibc];

    ivec = sptree.getArray<GINT>("idir");
    if ( fvecvec.size() != nbc ) {
      cout << "GUtils::get_bdy_block: SPONGE bc is specified; a vector of size(base_type) must specify sponge surface direction idir for each bc entry in 'base_type'" << endl;
      assert(FALSE); 
    }
    stblock.idir = ivec[ibc];

    fvec = sptree.getArray<GFTYPE>("xstart");
    if ( fvecvec.size() != nbc ) {
      cout << "GUtils::get_bdy_block: SPONGE bc is specified; a vector of size(base_type) must specify start positions in direction idir for each bc entry in 'base_type'" << endl;
      assert(FALSE); 
    }
    stblock.xstart = fvec[ibc];

    fvec = sptree.getArray<GFTYPE>("xmax");
    if ( fvecvec.size() != nbc ) {
      cout << "GUtils::get_bdy_block: SPONGE bc is specified; a vector of size(base_type) must specify max positions in direction idir for each bc entry in 'base_type'" << endl;
      assert(FALSE); 
    }
    stblock.xmax = fvec[ibc];

  }

  return TRUE;

} // end, method get_bdy_block



} // end, namespace geoflow

