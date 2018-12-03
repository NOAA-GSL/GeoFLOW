/*
 * observer_factory.hpp
 *
 *  Created on: Nov 29, 2018
 *      Author: bryan.flynt
 */

#ifndef SRC_ODEINT_OBSERVER_FACTORY_HPP_
#define SRC_ODEINT_OBSERVER_FACTORY_HPP_

#include <memory>

#include "odeint/observer_base.hpp"
#include "tbox/property_tree.hpp"

namespace geoflow {
namespace odeint {


template<typename EquationType>
struct ObserverFactory {
	using Equation    = EquationType;
	using ObsBase     = ObserverBase<Equation>;
	using ObsBasePtr  = std::shared_ptr<ObsBase>;

	static ObsBasePtr build(const tbox::PropertyTree& ptree);

};

} // namespace odeint
} // namespace geoflow

#include "odeint/observer_factory.ipp"

#endif /* SRC_ODEINT_OBSERVER_FACTORY_HPP_ */
