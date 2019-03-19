/*
 * observer_factory.hpp
 *
 *  Created on: Nov 29, 2018
 *      Author: bryan.flynt
 */

#ifndef SRC_PDEINT_OBSERVER_FACTORY_HPP_
#define SRC_PDEINT_OBSERVER_FACTORY_HPP_

#include <memory>

#include "pdeint/observer_base.hpp"
#include "tbox/property_tree.hpp"

namespace geoflow {
namespace pdeint {


template<typename EquationType>
struct ObserverFactory {
	using Equation    = EquationType;
	using ObsBase     = ObserverBase<Equation>;
	using ObsBasePtr  = std::shared_ptr<ObsBase>;
	using Grid        = EquationType::Grid;

	static ObsBasePtr build(const tbox::PropertyTree& ptree, Grid& grid);

};

} // namespace pdeint
} // namespace geoflow

#include "pdeint/observer_factory.ipp"

#endif /* SRC_PDEINT_OBSERVER_FACTORY_HPP_ */
