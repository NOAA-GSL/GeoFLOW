/*
 * odeint.cpp
 *
 *  Created on: Nov 16, 2018
 *      Author: bflynt
 */



#include <iostream>
#include <memory>
#include <valarray>

#include "odeint/equation_base.hpp"
#include "odeint/euler_stepper.hpp"
#include "odeint/stepper_base.hpp"

using namespace geoflow::odeint;


template<
typename StateType,
typename ValueType = double,
typename DerivType = StateType,
typename TimeType  = ValueType,
typename JacoType  = std::nullptr_t
>
struct EquationTypes {
	using State      = StateType;
	using Value      = ValueType;
	using Derivative = DerivType;
	using Time       = TimeType;
	using Jacobian   = JacoType;
};




template<typename TypePack>
struct HarmonicOscillator : public EquationBase<TypePack> {
	using Interface  = EquationBase<TypePack>;
	using State      = typename Interface::State;
	using Value      = typename Interface::Value;
	using Derivative = typename Interface::Derivative;
	using Time       = typename Interface::Time;
	using Jacobian   = typename Interface::Jacobian;

	static constexpr Value m_gam = 0.15;


	virtual ~HarmonicOscillator() = default;

protected:

	void dt_impl(const State& u, Time& dt){
		dt = 0.1;
	}

	void dudt_impl(const State& u, Derivative& dudt, const Time& t){
		dudt[0] = +u[1];
		dudt[1] = -u[0] - m_gam*u[1];
	}

	void dfdu_impl(const State& u, Jacobian& dfdu, const Time& t){
		(void*)(&dfdu);
	}

};





int main(){

	using MyTypes = EquationTypes<std::valarray<double>>;
	using EqnBase = EquationBase<MyTypes>;
	using EqnImpl = HarmonicOscillator<MyTypes>;
	using StpBase = StepperBase<EqnBase>;
	using StpImpl = EulerStepper<EqnBase>;

	std::shared_ptr<EqnBase> sys(new EqnImpl());
	std::shared_ptr<StpBase> stepper(new StpImpl());

	const int N = 10;
	typename MyTypes::State u(2);
	typename MyTypes::Time  t  = 0;
	typename MyTypes::Time  dt = 0.1;

	for(int i = 0; i < N; ++i){
		stepper->step(sys,u,t,dt);
		t += dt;
	}

	std::valarray<std::valarray<double>> a = {{1,2},{3,4}};
	std::valarray<std::valarray<double>> b = {{1,2},{3,4}};
	b = std::valarray<double>(3.1415,2) * a;


	return 0;
}
