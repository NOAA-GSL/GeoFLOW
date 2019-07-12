/*
 * mixer_factory.hpp
 *
 *  Created on: Mar 27, 2019 
 *      Author: bryan.flynt, d.rosenberg
 */

#ifndef SRC_PDEINT_MIXER_FACTORY_HPP_
#define SRC_PDEINT_MIXER_FACTORY_HPP_

#include <memory>
#include <string>
#include "tbox/error_handler.hpp"
#include "pdeint/null_mixer.hpp"
#include "pdeint/mixer_base.hpp"
#include "tbox/property_tree.hpp"


namespace geoflow {
namespace pdeint {


template<typename EquationType>
struct MixerFactory {
	using Equation    = EquationType;
	using MixBase     = MixerBase<Equation>;
	using MixBasePtr  = std::shared_ptr<MixBase>;
	using Grid        = typename EquationType::Grid;

	static MixBasePtr build(const tbox::PropertyTree& ptree, Grid& grid);

};

} // namespace pdeint
} // namespace geoflow

#include "pdeint/mixer_factory.ipp"

#endif /* SRC_PDEINT_MIXER_FACTORY_HPP_ */
