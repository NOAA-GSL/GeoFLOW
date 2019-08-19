
namespace gspecb {


//**********************************************************************************
//**********************************************************************************
// METHOD : impl_uniform
// DESC   : Do bdy specification assuming that base_type bdy condition
//          in property tree is to be set for all indices
// ARGS   : sptree : specification prop tree
//          grid   : grid
//          id     : may serve as canonical bdy id
//          ibdy   : indirection array into state indicating global bdy
//          tbdy   : array of size ibdy.size giving bdy condition type, returned
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_uniform(const PropertyTree &sptree, GGrid &grid, const GINT id, GTVector<GSIZET> &ibdy, GTVector<GBdyType> &tbdy)
{
  GBdyType btype = geoflow::str2bdytype(sptree.getValue<GString>("base_type", "GBDY_NONE"));

  tbdy.resizem(ibdy.size());
  tbdy = btype;

  return TRUE;

} // end, method impl_uniform


} // end, gspecb namespace
