/*
 * stepper_factory.hpp
 *
 *  Created on: Nov 29, 2018
 *      Author: bryan.flynt
 */

#ifndef SRC_ODEINT_STEPPER_FACTORY_HPP_
#define SRC_ODEINT_STEPPER_FACTORY_HPP_


#include "odeint/stepper_base.hpp"
#include "tbox/property_tree.hpp"

namespace geoflow {
namespace odeint {


template<typename EquationType>
struct StepperFactory {
	using Equation     = EquationType;
	using StepBase     = StepperBase<Equation>;
	using StepBasePtr  = std::shared_ptr<StepBase>;
	using EqnBasePtr   = typename StepBase::EquationPtr;


	static StepBasePtr build(const tbox::PropertyTree& ptree, const EqnBasePtr& eqn);

};

} // namespace odeint
} // namespace geoflow

#include "stepper_factory.ipp"


#endif /* SRC_ODEINT_STEPPER_FACTORY_HPP_ */
