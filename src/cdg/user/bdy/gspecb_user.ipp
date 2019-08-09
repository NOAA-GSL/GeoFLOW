
namespace gspecb {


//**********************************************************************************
//**********************************************************************************
// METHOD : impl_mybdyspec
// DESC   : Do bdy specification using my new method
// ARGS   : ptree: main prop tree
//          grid : grid
//          ibdy  : indirection array into state indicating global bdy
//          tbdy  : array of size ibdy.size giving bdy condition type, returned
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_mybdyspec(const PropteryTree &ptree, GGrid &grid, const GINT id, GTVector<GSIZET> &ibdy, GTVector<GBdyType> &tbdy)
{

  /* 
     Fill in here; change function name 
     if desired. Add (unique) function name to 
     gspecb_factory.ipp.
  */

  return FALSE;

} // end, method impl_mybdyspec


} // end, gupdateb namespace
