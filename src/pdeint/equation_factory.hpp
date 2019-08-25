/*
 * equation_factory.hpp
 *
 *  Created on: July 8, 2019 
 *      Author: bryan.flynt, d.rosenberg
 */

#ifndef SRC_PDEINT_EQUATION_FACTORY_HPP_
#define SRC_PDEINT_EQUATION_FACTORY_HPP_

#include <memory>
#include <string>
#include "tbox/error_handler.hpp"
#include "pdeint/equation_base.hpp"
#include "tbox/property_tree.hpp"


namespace geoflow {
namespace pdeint {


template<typename EquationType>
struct EquationFactory {
	using Equation    = EquationType;
	using EqnBase     = EquationBase<Equation>;
	using EqnBasePtr  = std::shared_ptr<EqnBase>;
        using State       = typename EquationType::State;
	using Grid        = typename EquationType::Grid;

	static EqnBasePtr build(const tbox::PropertyTree& ptree, Grid& grid, State& utmp);

};

} // namespace pdeint
} // namespace geoflow

#include "pdeint/equation_factory.ipp"

#endif /* SRC_PDEINT_EQUATION_FACTORY_HPP_ */
