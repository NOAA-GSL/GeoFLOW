/*
 * integrator_factory.hpp
 *
 *  Created on: Nov 29, 2018
 *      Author: bryan.flynt
 */

#ifndef SRC_PDEINT_INTEGRATOR_FACTORY_HPP_
#define SRC_PDEINT_INTEGRATOR_FACTORY_HPP_


#include "pdeint/mixer_base.hpp"
#include "pdeint/integrator.hpp"
#include "tbox/property_tree.hpp"

namespace geoflow {
namespace pdeint {


template<typename EquationType>
struct IntegratorFactory {
	using Equation      = EquationType;
//      using EqnBase       = typename Equation::Base;
	using EqnBase       = EquationBase<EquationType>;
	using ObsBase       = ObserverBase<Equation>;
	using EqnBasePtr    = std::shared_ptr<EqnBase>;
	using MixerBasePtr  = std::shared_ptr<MixerBase<Equation>>;
	using ObsBasePtr    = std::shared_ptr<std::vector<std::shared_ptr<ObsBase>>>;
	using IntegratorPtr = std::shared_ptr<Integrator<Equation>>;
	using Grid          = typename Equation::Grid;
	using Value         = typename Equation::Value;
	using Time          = typename Equation::Time;

	static IntegratorPtr build(const tbox::PropertyTree& ptree, const EqnBasePtr& eqn, const MixerBasePtr& mixer, const ObsBasePtr& obs, Grid& grid);

};

} // namespace pdeint
} // namespace geoflow

#include "integrator_factory.ipp"

#endif /* SRC_PDEINT_INTEGRATOR_FACTORY_HPP_ */