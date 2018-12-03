/*
 * observer_factory.ipp
 *
 *  Created on: Nov 29, 2018
 *      Author: bryan.flynt
 */

//#include "observer_factory.hpp"

#include <string>

#include "tbox/error_handler.hpp"
#include "odeint/null_observer.hpp"

namespace geoflow {
namespace odeint {


template<typename ET>
typename ObserverFactory<ET>::ObsBasePtr
ObserverFactory<ET>::build(const tbox::PropertyTree& ptree){


	// Set the default observer type
	const std::string default_observer = "none";

	// Get the type of observer
	const std::string observer_name = ptree.getValue("observer", default_observer);

	// Create the Stepper and cast to base type
	ObsBasePtr base_ptr;
	if( "none" == observer_name ){
		using ObsImpl = NullObserver<Equation>;

		// Allocate Stepper Implementation
		std::shared_ptr<ObsImpl> obs_impl(new ObsImpl());

		// Set any parameters we need to set
		// NA

		// Set back to base type
		base_ptr = obs_impl;
	}
	else {
		EH_ERROR("Requested observer not found: " << observer_name);
	}

	return base_ptr;
}


} // namespace odeint
} // namespace geoflow

