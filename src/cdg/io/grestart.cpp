//**********************************************************************************
//**********************************************************************************
// METHOD : grestart
// DESC   : Read restart files for entire state, based on prop tree configuration.
//          NOTE: grid is handled in ggrid_factory currently.
// ARGS   : ptree : main property tree
//          tindex: time index to restart
//          u     : state object
// RETURNS: none
//**********************************************************************************
void grestart(const geoflow::tbox::PropertyTree& ptree, GSIZET tindex, GTVector<GFTYPE> &u)
{
  GINT                 myrank = GCOMM::WorldRank(comm_);
  GINT                 dim, ivers;
  GElemType            gtype;
  GTMatrix<GINT>       porder;
  char                 sfname[2048];
  GFTYPE               time;
  GString              def = ".";
  GString              sfname;
  GString              sin;
  GTVector<GString>    stnames;
  std::vector<GString> stdobslist;
  std::vector   <GINT> stdilist;
  geoflow::tbox::PropertyTree& inputptree;
  stdilist   = ptree.getArray<GINT>("state_indices");
  stdobslist = ptree.getArray<GString>("state_names");
  stnames    = stdobslist;

  if ( gobslist.contains("posixio_observer") ) { 
    inputptree = ptree.getPropertyTree("posixio_observer");
    sin         = inputptree.getValue<GString>("indirectory",def);
    for ( j=0; j<stdilist.size(); j++ ) { // Retrieve all grid state vectors
      sfname      = sin + "/" + stnames[stdilist[j]];
      printf(fname, "%s/%s.%05d.out", sin.c_str(), sfname.c_str(),  myrank);
      sfname.assign(fname);
      gio_read(sfname, comm_, u[stdlist[j]]);
    }
  }
  else {
    assert( FALSE && "Configuration does not exist for grid files");
  }

} // end, grestart method
