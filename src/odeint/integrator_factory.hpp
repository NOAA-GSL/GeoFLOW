/*
 * integrator_factory.hpp
 *
 *  Created on: Nov 29, 2018
 *      Author: bryan.flynt
 */

#ifndef SRC_ODEINT_INTEGRATOR_FACTORY_HPP_
#define SRC_ODEINT_INTEGRATOR_FACTORY_HPP_


#include "odeint/integrator.hpp"
#include "tbox/property_tree.hpp"

namespace geoflow {
namespace odeint {


template<typename EquationType>
struct IntegratorFactory {
	using Equation      = EquationType;
	using StepBase      = StepperBase<Equation>;
	using ObsBase       = ObserverBase<Equation>;
	using StepBasePtr   = std::shared_ptr<StepBase>;
	using ObsBasePtr    = std::shared_ptr<ObsBase>;
	using IntegratorPtr = std::shared_ptr<Integrator<Equation>>;
	using Value         = typename Equation::Value;
	using Time          = typename Equation::Time;

	static IntegratorPtr build(const tbox::PropertyTree& ptree, const StepBasePtr& stp, const ObsBasePtr& obs);

};

} // namespace odeint
} // namespace geoflow

#include "integrator_factory.ipp"

#endif /* SRC_ODEINT_INTEGRATOR_FACTORY_HPP_ */
