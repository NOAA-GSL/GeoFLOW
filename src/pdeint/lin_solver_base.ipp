/*
 * lin_solver_base.ipp
 *
 *  Created on: Mar. 16, 2020
 *      Author: D. Rosenberg
 */


namespace geoflow {
namespace pdeint {


template<typename Types>
typename LinSolverBase<Types>::GNormType
LinSolverBase<Types>::str2normtype(std::string snorm)
{
  int j0 = static_cast<int>(GCG_NORM_INF);
  static const constexpr char* const
    sGNormType[] = {"GCG_NORM_INF","GCG_NORM_EUC","GCG_NORM_L2","GCG_NORM_L1","GCG_NORM_NONE"};
  std::string s0;

  for ( auto j=j0; j<GCG_NORM_NONE; j++ ) {
    s0 = sGNormType[j];
    if ( snorm.compare(s0) == 0 ) return static_cast<GNormType>(j);
  }
  assert(FALSE && "Invalid norm type");
  return static_cast<GNormType>(GCG_NORM_INF);
}


} // namespace pdeint
} // namespace geoflow

