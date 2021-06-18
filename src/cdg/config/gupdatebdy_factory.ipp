//==================================================================================
// Module       : gupdatebdy_factory
// Date         : 7/11/19 (DLR)
// Description  : GeoFLOW bdy update object initialization factory
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================


//**********************************************************************************
//**********************************************************************************
// METHOD : build
// DESC   : Build bdy update object
// ARGS   : ptree   : main property tree
//          sbdy    : bdy condition name
//          grid    : grid object
//          bcblock : stBdyBlock structure
//          ibdy    : bdy indirection indices into computational volume
//          igbdy_start:
//                    where ibdy start in global bdy index list
//                   
// RETURNS: none.
//**********************************************************************************
template<typename Types>
typename GUpdateBdyFactory<Types>::BdyBasePtr
GUpdateBdyFactory<Types>::build(const PropertyTree& ptree, const GString &sbdy, Grid &grid, stBdyBlock &bcblock, GTVector<GSIZET> &ibdy, GTVector<GUINT> &dbdy, GSIZET igbdy_start)
{
  GBOOL              bret = FALSE;
  BdyBasePtr         base_ptr;


  if ( "none"         == sbdy
    || ""             == sbdy ) {
    using UpdateImpl = NullUpdateBdy<Types>;

    // Allocate update implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl());

    // Set back to base type
    base_ptr = update_impl;
    return base_ptr;
  }

  if ( !ptree.isPropertyTree(sbdy) ) {
    cout << "GUpdateBdyFactory<Types>::build: PropertyTree " << sbdy << " not found" << endl;
    assert(FALSE);
  }

  base_ptr = GUpdateBdyFactory<Types>::get_bdy_class(ptree, grid, bcblock, ibdy, dbdy, igbdy_start);

  return base_ptr;

} // end, init method build


//**********************************************************************************
//**********************************************************************************
// METHOD : get_inflow_callback
// DESC   : Gets CallbackPtr corresponding to sname arg for inflow conditions.
//          the function are gottne from the ginflow_user.* struct
//          collection of methods, which the user may modify.
// ARGS   : sname : inflow function name
//          id    : canonical bdy id
// RETURNS: CallbackPtr for callback function
//**********************************************************************************
template<typename Types>
typename GUpdateBdyFactory<Types>::CallbackPtr
GUpdateBdyFactory<Types>::get_inflow_callback(const GString& sname, const GINT id)
{
  GBOOL              bret = FALSE;
  CallbackPtr        callback;


  if      ( ""     == sname 
   ||       "none" == sname ) {
    assert(FALSE); // Must have a named method
  }
  if      ( "myinflow"     == sname ) {
    callback = 

         [](EqnBasePtr &eqn,
            Grid       &grid,
            Time       &time,
            const GINT id,
            State      &u,
            State      &utmp,
            State      &ub)-> GBOOL{return GInflowUser<Types>::myinflow(eqn, grid, time, id, u, utmp, ub);}; 

  }
  else {
    assert(FALSE && "Specified inflow bdy update method unknown");
  }

  return callback;

} // end, init method get_inflow_callback


//**********************************************************************************
//**********************************************************************************
// METHOD : get_bdy_class
// DESC   : get bdy class for specified bdy type
// ARGS   : ptree   : main prop tree
//          grid    : grid object
//          bcblock : stBdyBlock structure
//          ibdy    : indirection indices into computational volume, 
//                    representing the bdy nodes this method applies to
//          dbdy    : bdy node descriptors
//          igbdy_start:
//                    where ibdy start in global bdy index list
//                   
// RETURNS: none.
//**********************************************************************************
template<typename Types>
typename GUpdateBdyFactory<Types>::BdyBasePtr
GUpdateBdyFactory<Types>::get_bdy_class(const PropertyTree &ptree, Grid &grid, stBdyBlock &bcblock, GTVector<GSIZET> &ibdy, GTVector<GUINT> &dbdy, GSIZET igbdy_start)
{
  GINT               nstate;
  GBdyType           bdytype = bcblock.tbdy;
  GLONG              iloc;
  BdyBasePtr         base_ptr;
  

  nstate = bcblock.istate.size();
  

  if       ( GBDY_DIRICHLET == bdytype ) {
    using UpdateImpl = GDirichletBdy<Types>;
    typename GDirichletBdy<Types>::Traits traits;

    traits.istate.resize(bcblock.istate.size()) ; 
    traits.ibdyvol.resize(ibdy.size()); 
  
    traits.bdyid       = bcblock.bdyid;
    traits.istate      = bcblock.istate;
    traits.ibdyvol     = ibdy;
    traits.value .resize(bcblock.value.size())  ; 
    traits.value       = bcblock.value;;

    // Allocate update implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl(traits));

    // Set back to base type
    base_ptr = update_impl;
    return base_ptr;
  }
  else if ( GBDY_INFLOW == bdytype ) {
    using UpdateImpl = GInflowBdy<Types>;
    typename GInflowBdy<Types>::Traits traits;

    traits.istate.resize(bcblock.istate.size()) ; 
    traits.ibdyvol.resize(ibdy.size()); 
  
    traits.bdyid        = bcblock.bdyid;
    traits.istate       = bcblock.istate;
    traits.ibdyvol      = ibdy;
    traits.use_init     = bcblock.use_init;
    traits.compute_once = bcblock.compute_once;
    if ( traits.use_init ) {
      traits.smethod = ptree.getValue<GString>("initstate_block");
    }
    else {
      traits.callback = GUpdateBdyFactory<Types>::get_inflow_callback(bcblock.smethod, bcblock.bdyid);
    }
    traits.ptree = ptree;
    // Allocate update implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl(traits));

    // Set back to base type
    base_ptr = update_impl;
    return base_ptr;
  }
#if 0
  else if ( GBDY_NOSLIP == bdytype ) {
    using UpdateImpl = GNoSlipBdy<Types>;
    typename GNoSlipBdy<Types>::Traits traits;

    // Allocate update implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl(traits));

    // Set back to base type
    base_ptr = update_impl;
    return base_ptr;
  }
#endif
  else if ( GBDY_0FLUX == bdytype ) {
    using UpdateImpl = G0FluxBdy<Types>;
    typename G0FluxBdy<Types>::Traits traits;

    traits.istate.resize(bcblock.istate.size()) ; 
    traits.ibdyvol.resize(ibdy.size()); 
  
    traits.bdyid       = bcblock.bdyid;
    traits.istate      = bcblock.istate;
    traits.ibdyvol     = ibdy;
    traits.ibdydsc     = dbdy;

    traits.ibdyloc.resize(traits.ibdyvol.size());
    // Find index of ibdyvol in global bdy vector:
    for ( auto j=0; j<traits.ibdyloc.size(); j++ ) {
      iloc = j + igbdy_start;
      assert(iloc >= 0);
      traits.ibdyloc[j] = iloc;
    }
    
    // Allocate update implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl(traits));

    // Set back to base type
    base_ptr = update_impl;
    return base_ptr;
  }
  else if ( GBDY_OUTFLOW == bdytype ) {
    assert(FALSE); // not available yet
    using UpdateImpl = GOutflowBdy<Types>;
    typename GOutflowBdy<Types>::Traits traits;

    traits.istate.resize(bcblock.istate.size()) ; 
    traits.ibdyvol.resize(ibdy.size()); 
  
    traits.bdyid       = bcblock.bdyid;
    traits.istate      = bcblock.istate;
    traits.ibdyvol     = ibdy;
    
    // Allocate update implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl(traits));

    // Set back to base type
    base_ptr = update_impl;
    return base_ptr;
  }
  else if ( GBDY_SPONGE == bdytype ) {
    using UpdateImpl = GSpongeBdy<Types>;
    typename GSpongeBdy<Types>::Traits traits;

    traits.istate.resize(bcblock.istate.size()) ; 
    traits.ibdyvol.resize(ibdy.size()); 
  
    traits.bdyid       = bcblock.bdyid;
    traits.istate      = bcblock.istate;
    traits.ibdyvol     = ibdy;
    traits.idir        = bcblock.idir;
    traits.xstart      = bcblock.xstart;
//  traits.xend        = bcblock.xmax;

    traits.farfield.resize(bcblock.farfield.size());
    traits.falloff .resize(bcblock.falloff .size());
    traits.exponent.resize(bcblock.exponent.size());

    traits.farfield = bcblock.farfield;
    traits.falloff  = bcblock.falloff;
    traits.exponent = bcblock.exponent;
    
    // Allocate update implementation
    std::shared_ptr<UpdateImpl> update_impl(new UpdateImpl(traits));

    // Set back to base type
    base_ptr = update_impl;
    return base_ptr;
  }
  else {
    assert(FALSE && "Specified bdy update method unknown");
  }

  return base_ptr;

} // end, init method get_bdy_class


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
template<typename Types>
GINT GUpdateBdyFactory<Types>::bdy_block_conform_per(const geoflow::tbox::PropertyTree &sptree)
{
  GBOOL                battempt;
  GString              bdyclass;
  std::vector<GString> stypes;

  bdyclass = sptree.getValue<GString>("bdy_class"); // required



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
template<typename Types>
GBOOL GUpdateBdyFactory<Types>::get_bdy_block(const geoflow::tbox::PropertyTree &ptree, GString &sbdy, GINT ibc, stBdyBlock &stblock)
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
  stblock.exponent .clear();
  stblock.bdyclass .clear();
  stblock.smethod  .clear();

  stblock.use_init    = FALSE;
  stblock.idir        = 0;
  stblock.bdyid       = -1;
  stblock.tbdy        = GBDY_NONE;
  stblock.xstart      = 0.0;
  stblock.bdyid       = -1;
  stblock.tbdy        = GBDY_NONE;
  stblock.bdyclass    = "";
  stblock.smethod     = "";

  // Get bdy block data:
  
  stblock.bdyclass = sptree.getValue<GString>("bdy_class"); // required

  // Get base bdy condition type:
  stypes = sptree.getArray<GString>("base_type"); // is a vector
  nbc = stypes.size(); // number of bcs specified

  // If ibc out of range, return:
  if ( ibc < 0 || ibc >= nbc ) return FALSE; 
  stblock.tbdy = geoflow::str2bdytype(stypes[ibc]); // set bdy type id

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
  if ( stypes[ibc] == "GBDY_DIRICHLET" ) {
    fvecvec = sptree.getArray2D<GFTYPE>("value");
    if ( fvecvec.size() != nbc ) {
      cout << "GUtils::get_bdy_block: DIRICHLET bc is specified; 2d JSON vector of size(base_type) must specify DIRICHLET values for each istate ('[]' is valid for non-DIRICHLET entries)" << endl;
      assert(FALSE); 
    }
    stblock.value.resize(nstate);
    stblock.value = fvecvec[ibc];
  }
  
  // If INFLOW bdy, retrieve method name, or other data:
  else if ( stypes[ibc] == "GBDY_INFLOW" ) {
    bvec = sptree.getArray<GBOOL>("use_init");
    if ( bvec.size() != nbc ) {
      cout << "GUtils::get_bdy_block: INFLOW bc is specified; a vector of size(base_type) must specify Boolean 'use_init' flags for each bc entry in 'base_type' " << endl;
      assert(FALSE); 
    }
    stblock.use_init= bvec[ibc];

    bvec = sptree.getArray<GBOOL>("compute_once");
    if ( bvec.size() != nbc ) {
      cout << "GUtils::get_bdy_block: INFLOW bc is specified; a vector of size(base_type) must specify Boolean 'compute_once' flags for each bc entry in 'base_type' " << endl;
      assert(FALSE); 
    }
    stblock.compute_once = bvec[ibc];

    if ( !stblock.use_init ) { // get user method if required
      svec = sptree.getArray<GString>("method");
      if ( svec.size() != nbc ) {
        cout << "GUtils::get_bdy_block: INFLOW bc is specified; a vector of size(base_type) must specify inflow methods for each bc entry in 'base_type' (\"\" is valid for non-INFLOW entries)" << endl;
        assert(FALSE); 
      }
      stblock.smethod = svec[ibc];
    }
  } 

  // Note: the only data required for 0FLUX bdys
  //       is set above already

  // If SPONGE bdy, retrieve required data:
  else if ( stypes[ibc] == "GBDY_SPONGE" ) { 
    fvecvec = sptree.getArray2D<GFTYPE>("farfield");
    if ( fvecvec.size() != nbc ) {
      cout << "GUtils::get_bdy_block: SPONGE bc is specified; a vector of size(base_type) must specify farfield values for each state group for each bc entry in 'base_type'   ('[]' is valid for non-SPONGE entries)" << endl;
      assert(FALSE); 
    }
    stblock.farfield.resize(nstate);
    stblock.farfield  = fvecvec[ibc];

    fvecvec = sptree.getArray2D<GFTYPE>("falloff_rate");
    if ( fvecvec.size() != nbc ) {
      cout << "GUtils::get_bdy_block: SPONGE bc is specified; a vector of size(base_type) must specify falloff rates for each state for each bc entry in 'base_type'   ('[]' is valid for non-SPONGE entries)" << endl;
      assert(FALSE); 
    }
    stblock.falloff.resize(nstate);
    stblock.falloff = fvecvec[ibc];

    fvecvec = sptree.getArray2D<GFTYPE>("exponent");
    if ( fvecvec.size() != nbc ) {
      cout << "GUtils::get_bdy_block: SPONGE bc is specified; a vector of size(base_type) must specify diffusion exponents for each state for each bc entry in 'base_type'   ('[]' is valid for non-SPONGE entries)" << endl;
      assert(FALSE); 
    }
    stblock.exponent.resize(nstate);
    stblock.exponent= fvecvec[ibc];

    if ( !sptree.keyExists("idir" ) ) {
      cout << "GUtils::get_bdy_block: SPONGE bc is specified; a constant integer is required that specifies sponge surface direction, idir, for all fields in 'istate'" << endl;
      assert(FALSE); 
    }
    stblock.idir = sptree.getValue<GINT>("idir");

    if ( !sptree.keyExists("xstart" ) ) {
      cout << "GUtils::get_bdy_block: SPONGE bc is specified; a coordinate value in direction idir must be provided specifying start positions in direction idir for all fields in 'istate'" << endl;
      assert(FALSE); 
    }
    stblock.xstart = sptree.getValue<GFTYPE>("xstart");

#if 0
    fvec = sptree.getArray<GFTYPE>("xmax");
    if ( fvecvec.size() != nbc ) {
      cout << "GUtils::get_bdy_block: SPONGE bc is specified; a vector of size(base_type) must specify max positions in direction idir for each bc entry in 'base_type'" << endl;
      assert(FALSE); 
    }
    stblock.xmax = fvec[ibc];

#endif
  }

  return TRUE;

} // end, method get_bdy_block



