/*
 * integrator_factory.hpp
 *
 *  Created on: Nov 29, 2018
 *      Author: bryan.flynt
 */

#ifndef SRC_PDEINT_INTEGRATOR_FACTORY_HPP_
#define SRC_PDEINT_INTEGRATOR_FACTORY_HPP_


#include "pdeint/equation_base.hpp"
#include "pdeint/mixer_base.hpp"
#include "pdeint/observer_base.hpp"
#include "pdeint/null_observer.hpp"
#include "pdeint/integrator.hpp"
#include "tbox/property_tree.hpp"

namespace geoflow {
namespace pdeint {


template<typename TypePack>
struct IntegratorFactory {
        using Types        = TypePack;
        using EqnBase      = EquationBase<TypePack>;
        using EqnBasePtr   = std::shared_ptr<EqnBase>;
        using State        = typename Types::State;
        using StateInfo    = typename Types::StateInfo;
        using Grid         = typename Types::Grid;
        using Ftype        = typename Types::Ftype;
        using Derivative   = typename Types::Derivative;
        using Time         = typename Types::Time;
        using Jacobian     = typename Types::Jacobian;
        using Size         = typename Types::Size;
        using MixerBasePtr = std::shared_ptr<MixerBase<Types>>;
        using ObsBasePtr   = std::shared_ptr<std::vector<std::shared_ptr<ObserverBase<Types>>>>;
	using IntegratorPtr = std::shared_ptr<Integrator<Types>>;

	static IntegratorPtr build(const tbox::PropertyTree& ptree, const EqnBasePtr& eqn, const MixerBasePtr& mixer, const ObsBasePtr& obs, Grid& grid);

};

} // namespace pdeint
} // namespace geoflow

#include "integrator_factory.ipp"

#endif /* SRC_PDEINT_INTEGRATOR_FACTORY_HPP_ */
