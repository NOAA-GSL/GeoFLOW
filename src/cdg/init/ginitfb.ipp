
namespace ginitfb {


//**********************************************************************************
//**********************************************************************************
// METHOD : impl_rand
// DESC   : Inititialize velocity with Gaussian-randomized values
// ARGS   : vtree  : initial condition property tree
//          grid   : grid object
//          peqn   : pointer to EqnBase
//          time   : initialization time
//          utmp   : tmp arrays
//          ub     : boundary state (also initialized here)
//          uf     : state to be initialized.
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_rand(const PropertyTree &vtree, GGrid &grid, EqnBasePtr &peqn,  Time &time, State &utmp, State &ub, State &uf)
{

  return FALSE;

} // end, method impl_rand


} // end, ginitfb  namespace
