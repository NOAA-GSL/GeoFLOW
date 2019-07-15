#include "ginitf.hpp"


namespace ginitf {


//**********************************************************************************
//**********************************************************************************
// METHOD : initf_impl_rand
// DESC   : Initialize for random force function
// ARGS   : ftree: force prop tree
//          t    : time
//          u    : current state
//          ub   : bdy vectors (one for each state element)
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL initf_impl_rand(const PropertyTree &ftree, const Time &t, State &u, State &ub)
{

  assert(FALSE);

  for ( auto j=0; j<uf.size(); j++ ) {
    if ( uf[j] != NULLPTR ) *uf[j] = 0.0;
  }

  return TRUE;

} // end of method initf_rand



} // end, namespace ginitf
