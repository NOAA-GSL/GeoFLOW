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
IntegratorFactory<EquationType>::build(const tbox::PropertyTree& ptree, const EqnBasePtr& eqn, const ObsBasePtr& obs,
                                                                               Grid&       grid){

	// Get the integrator traits
	typename Integrator<Equation>::Traits traits;
        std::string stype;
	stype             = ptree.getValue("integ_type", Integrator<EquationType>::INTEG_CYCLE );
        if ( stype == "cycle" ) traits.integ_type = Integrator<EquationType>::INTEG_CYCLE;
        if ( stype == "time"  ) traits.integ_type = Integrator<EquationType>::INTEG_TIME;
	traits.cycle_end  = ptree.getValue("cycle_end", std::numeric_limits<Value>::max() );
	traits.dt         = ptree.getValue("dt",  std::numeric_limits<Time>::lowest() );
	traits.time_max  = ptree.getValue("time_end",  std::numeric_limits<Time>::max() );

	// Allocate Integrator Implementation
	IntegratorPtr integrator_ptr(new Integrator<Equation>(eqn, obs, grid, traits));

	// Set any parameters we need to set
	// NA

	// Return
	return integrator_ptr;
}


} // namespace pdeint
} // namespace geoflow
