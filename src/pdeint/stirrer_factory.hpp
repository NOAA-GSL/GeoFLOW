/*
 * stir_factory.hpp
 *
 *  Created on: MAr 27, 2019 
 *      Author: bryan.flynt, d.rosenberg
 */

#ifndef SRC_PDEINT_STIRRER_FACTORY_HPP_
#define SRC_PDEINT_STIRRER_FACTORY_HPP_

#include <memory>
#include <string>
#include "tbox/error_handler.hpp"
#include "pdeint/null_stirrer.hpp"
#include "pdeint/stirrer_base.hpp"
#include "tbox/property_tree.hpp"


namespace geoflow {
namespace pdeint {


template<typename EquationType>
struct StirrerFactory {
	using Equation    = EquationType;
	using StirBase    = StirrerBase<Equation>;
	using StirBasePtr = std::shared_ptr<StirBase>;
	using Grid        = typename EquationType::Grid;

	static StirBasePtr build(const tbox::PropertyTree& ptree, Grid& grid);

};

} // namespace pdeint
} // namespace geoflow

#include "pdeint/stirrer_factory.ipp"

#endif /* SRC_PDEINT_STIRRER_FACTORY_HPP_ */
