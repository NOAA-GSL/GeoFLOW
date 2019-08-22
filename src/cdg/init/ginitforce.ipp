
namespace ginitforce {




//**********************************************************************************
//**********************************************************************************
// METHOD : impl_myforceinit
// DESC   : Initialize force method 
// ARGS   : ptree: main prop tree
//          grid : grid
//          t    : time
//          utmp : tmp arrays
//          ub   : bdy vectors (one for each state element)
//          uf   : current state, overwritten here
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_myforceinit(const PropertyTree &ptree, GGrid &grid, Time &time, State &utmp, State &ub, State &uf)
{

  Time             tt = t;
  GString          serr = "impl_myforceinit: ";


  return FALSE;

} // end, impl_myforceinit




} // end, ginitforce namespace
