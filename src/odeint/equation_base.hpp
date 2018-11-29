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
	using Size       = typename TypePack::Size;

	EquationBase() = default;
	EquationBase(const EquationBase& eb) = default;
	virtual ~EquationBase() = default;
	EquationBase& operator=(const EquationBase& eb) = default;


	bool has_dt() const{
		return this->has_dt_impl();
	}

	/** Time difference corresponding to a CFL of one.
	 *
	 */
	void dt(const Time& t, State& u, Time& dt){
		this->dt_impl(t,u,dt);
	}

	void dudt(const Time& t, State& u, Derivative& dudt){
		this->dudt_impl(t,u,dudt);
	}

	void dfdu(const Time& t, State& u, Jacobian& dfdu){
		this->dfdu_impl(t,u,dfdu);
	}



protected:

	virtual bool has_dt_impl() const = 0;

	virtual void dt_impl(const Time& t, State& u, Time& dt) = 0;

	virtual void dudt_impl(const Time& t, State& u, Derivative& dudt) = 0;

	virtual void dfdu_impl(const Time& t, State& u, Jacobian& dfdu) = 0;
};


} // namespace odeint
} // namespace geoflow


#endif /* SRC_ODEINT_EQUATION_BASE_HPP_ */
