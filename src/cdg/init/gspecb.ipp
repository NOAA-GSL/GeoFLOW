
namespace gspecb {


//**********************************************************************************
//**********************************************************************************
// METHOD : impl_mybdyspec
// DESC   : Do bdy specification assuming that base_type bdy condition
//          in property tree is to be set for all indices
// ARGS   : ptree: specification prop tree
//          grid : grid
//          ibdy  : indirection array into state indicating global bdy
//          tbdy  : array of size ibdy.size giving bdy condition type, returned
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_uniform(const PropteryTree &ptree, GGrid &grid, GTVector<GSIZET> &ibdy, GTVector<GBdyType> &tbdy)
{
  GBdyType btype = ptree.value<GBdyType>("base_type");

  tbdy = btype;

  return TRUE;

} // end, method impl_uniform


} // end, gspecb namespace
