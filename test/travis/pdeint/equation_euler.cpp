/*
 * odeint.cpp
 *
 *  Created on: Nov 16, 2018
 *      Author: bflynt
 */



#include <iostream>
#include <memory>
#include <valarray>

#include "pdeint/equation_base.hpp"

using namespace geoflow::pdeint;


template<
typename StateType,
typename ValueType = double,
typename DerivType = StateType,
typename TimeType  = ValueType,
typename JacoType  = std::nullptr_t,
typename SizeType  = std::size_t
>
struct EquationTypes {
	using State      = StateType;
	using Value      = ValueType;
	using Derivative = DerivType;
	using Time       = TimeType;
	using Jacobian   = JacoType;
	using Size       = SizeType;
};




template<typename TypePack>
struct HarmonicOscillator : public EquationBase<TypePack> {
	using Interface  = EquationBase<TypePack>;
	using Base       = Interface;
	using State      = typename Interface::State;
	using Value      = typename Interface::Value;
	using Derivative = typename Interface::Derivative;
	using Time       = typename Interface::Time;
	using Jacobian   = typename Interface::Jacobian;

	static constexpr Value m_gam = 0.15;


	virtual ~HarmonicOscillator() = default;

protected:

	bool has_dt_impl() const{
		return false;
	}

	void dt_impl(const Time& t, State& u, Time& dt){
		dt = 1.0;
	}


	void step_impl(const Time& t, const Time& dt, State& u){
		Derivative dudt(u.size());
		this->dudt_impl_(t,u,dudt);
		u += dt * dudt;
	}


	void step_impl(const Time& t, const State& uin, const Time& dt, State& uout){
		uout = uin;
		this->step_impl(t,dt,uout);
	}


	void dudt_impl_(const Time& t, State& u, Derivative& dudt){
		dudt[0] = +u[1];
		dudt[1] = -u[0] - m_gam*u[1];
	}


};





int main(){

	using MyTypes = EquationTypes<std::valarray<double>>;
	using EqnBase = EquationBase<MyTypes>;
	using EqnImpl = HarmonicOscillator<MyTypes>;

	std::shared_ptr<EqnImpl> eqn_impl(new EqnImpl());
	std::shared_ptr<EqnBase> eqn_base = eqn_impl;


	const int N = 10;
	typename MyTypes::State u(2);
	typename MyTypes::State uerr(2);
	typename MyTypes::Time  t  = 0;
	typename MyTypes::Time  dt = 0.1;

	for(int i = 0; i < N; ++i){
		eqn_base->step(t,dt,u);
		t += dt;
	}


	return 0;
}
