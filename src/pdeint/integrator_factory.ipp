/*
 * integrator_factory.ipp
 *
 *  Created on: Nov 29, 2018
 *      Author: bryan.flynt
 */

#include <limits>

#include "integrator_factory.hpp"

namespace geoflow {
namespace pdeint {

template<typename EquationType>
typename IntegratorFactory<EquationType>::IntegratorPtr
IntegratorFactory<EquationType>::build(const tbox::PropertyTree& ptree, const EqnBasePtr& eqn, const ObsBasePtr& obs){

	// Get the integrator traits
	typename Integrator<Equation>::Traits traits;
	traits.cfl_min = ptree.getValue("cfl_min", std::numeric_limits<Value>::lowest() );
	traits.cfl_max = ptree.getValue("cfl_max", std::numeric_limits<Value>::max() );
	traits.dt_min  = ptree.getValue("dt_min",  std::numeric_limits<Time>::lowest() );
	traits.dt_max  = ptree.getValue("dt_max",  std::numeric_limits<Time>::max() );

	// Allocate Integrator Implementation
	IntegratorPtr integrator_ptr(new Integrator<Equation>(eqn,obs,traits));

	// Set any parameters we need to set
	// NA

	// Return
	return integrator_ptr;
}


} // namespace pdeint
} // namespace geoflow
