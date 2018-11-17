/*
 * equation_base.hpp
 *
 *  Created on: Nov 16, 2018
 *      Author: bflynt
 */

#ifndef SRC_ODEINT_EQUATION_BASE_HPP_
#define SRC_ODEINT_EQUATION_BASE_HPP_


namespace geoflow {
namespace odeint {



/**
 * \brief
 * Base class to provided an interface strategy for all
 * equation types.
 *
 * \details
 * The base class provides the strategy for all equation
 * implementations.
 */
template<typename TypePack>
class EquationBase {

public:
	using State      = typename TypePack::State;
	using Value      = typename TypePack::Value;
	using Derivative = typename TypePack::Derivative;
	using Time       = typename TypePack::Time;
	using Jacobian   = typename TypePack::Jacobian;

	EquationBase() = default;
	EquationBase(const EquationBase& eb) = default;
	virtual ~EquationBase() = default;
	EquationBase& operator=(const EquationBase& eb) = default;

	void dt(const State& u, Time& dt){
		this->dt_impl(u,dt);
	}

	void dudt(const State& u, Derivative& dudt, const Time& t){
		this->dudt_impl(u,dudt,t);
	}

	void dfdu(const State& u, Jacobian& dfdu, const Time& t){
		this->dfdu_impl(u,dfdu,t);
	}



protected:

	virtual void dt_impl(const State& u, Time& dt) = 0;

	virtual void dudt_impl(const State& u, Derivative& dudt, const Time& t) = 0;

	virtual void dfdu_impl(const State& u, Jacobian& dfdu, const Time& t) = 0;
};


} // namespace odeint
} // namespace geoflow


#endif /* SRC_ODEINT_EQUATION_BASE_HPP_ */
