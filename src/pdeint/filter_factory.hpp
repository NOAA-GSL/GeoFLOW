/*
 * filter_factory.hpp
 *
 *  Created on: Sep. 14, 2020
 *      Author: d.rosenberg
 */

#ifndef SRC_PDEINT_FILTER_FACTORY_HPP_
#define SRC_PDEINT_FILTER_FACTORY_HPP_

#include <memory>
#include <string>
#include "tbox/error_handler.hpp"
#include "tbox/property_tree.hpp"
#include "pdeint/filter_base.hpp"
#include "pdeint/null_filter.hpp"
#include "gboyd_filter.hpp"
#include "gprojection_filter.hpp"


namespace geoflow {
namespace pdeint {


template<typename TypePack>
struct FilterFactory {
	using Types         = TypePack;
	using FiltBase      = FilterBase<Types>;
	using FilterBasePtr = std::shared_ptr<FiltBase>;
	using Grid          = typename TypePack::Grid;

	static FilterBasePtr build(const tbox::PropertyTree& ptree, std::string filtername,  Grid& grid);

};

} // namespace pdeint
} // namespace geoflow

#include "pdeint/filter_factory.ipp"

#endif /* SRC_PDEINT_FILTER_FACTORY_HPP_ */
