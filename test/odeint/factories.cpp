/*
 * factories.cpp
 *
 *  Created on: Nov 29, 2018
 *      Author: bryan.flynt
 */



#include <iostream>
#include <memory>
#include <valarray>

#include "tbox/pio.hpp"
#include "tbox/property_tree.hpp"

#include "odeint/equation_base.hpp"
#include "odeint/integrator_factory.hpp"
#include "odeint/observer_factory.hpp"
#include "odeint/stepper_factory.hpp"

using namespace geoflow::tbox;
using namespace geoflow::odeint;


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

	void dudt_impl(const Time& t, State& u, Derivative& dudt){
		dudt[0] = +u[1];
		dudt[1] = -u[0] - m_gam*u[1];
	}

	void dfdu_impl(const Time& t, State& u, Jacobian& dfdu){
		//(void*)(&dfdu);
	}

};





int main(){

	using MyTypes = EquationTypes<std::valarray<double>>;
	using EqnBase = EquationBase<MyTypes>;
	using EqnImpl = HarmonicOscillator<MyTypes>;

	pio::initialize(0);

	std::shared_ptr<EqnImpl> eqn_impl(new EqnImpl());
	std::shared_ptr<EqnBase> eqn_base = eqn_impl;

	// Load the Property Tree
	PropertyTree ptree;
	ptree.load_string("{ \"time stepper\":\"bdf1\",\"observer\":\"none\",\"dt_max\":1.0}");

	// Create the Stepper
	auto stp_ptr = StepperFactory<EqnBase>::build(ptree,eqn_base);

	// Create the Observer
	auto obs_ptr = ObserverFactory<EqnBase>::build(ptree);

	// Create the Integrator
	auto integrator_ptr = IntegratorFactory<EqnBase>::build(ptree,stp_ptr,obs_ptr);

	const int N = 10;
	typename MyTypes::State u(2);
	typename MyTypes::State uerr(2);
	typename MyTypes::Time  t = 0;
	typename MyTypes::Time  t0 = 0;
	typename MyTypes::Time  dt = 0.1;

	integrator_ptr->steps(t0,dt,N,u,t);

	return 0;
}


