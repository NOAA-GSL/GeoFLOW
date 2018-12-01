/*
 * stepper_factory.ipp
 *
 *  Created on: Nov 29, 2018
 *      Author: bryan.flynt
 */

//#include "stepper_factory.hpp"

#include "odeint/euler_stepper.hpp"
#include "tbox/error_handler.hpp"

namespace geoflow {
namespace odeint {


template<typename EquationType>
typename StepperFactory<EquationType>::StepBasePtr
StepperFactory<EquationType>::build(const tbox::PropertyTree& ptree, const EqnBasePtr& eqn){

	// Set the default stepper type
	const std::string default_stepper = "bdf1";

	// Get the type of stepper
	const std::string stepper_name = ptree.getValue("time stepper", default_stepper);

	// Create the Stepper and cast to base type
	StepBasePtr base_ptr;
	if( "bdf1" == stepper_name ){
		using StpImpl = EulerStepper<Equation>;

		// Allocate Stepper Implementation
		std::shared_ptr<StpImpl> stp_impl(new StpImpl(eqn));

		// Set any parameters we need to set
		// NA

		// Set back to base type
		base_ptr = stp_impl;
	}
	else {
		EH_ERROR("Requested time stepper not found: " << stepper_name);
	}

	return base_ptr;
}


} // namespace odeint
} // namespace geoflow




