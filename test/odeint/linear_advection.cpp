/*
 * wave.cpp
 *
 *  Created on: Nov 18, 2018
 *      Author: bflynt
 */

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <vector>
#include <valarray>

#include "odeint/equation_base.hpp"
#include "odeint/stepper_base.hpp"
#include "odeint/euler_stepper.hpp"

using namespace geoflow::odeint;

/**
 * Helper type to define all types.
 *
 * Used to package all types into a single
 * parameter and help reduce the number of
 * template variables a user has to define.
 */
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



/**
 * Linear Advection Equation in 1-D
 *
 *  pu       pu
 * ---- + a ---- = 0
 *  pt       px
 */
template<typename TypePack>
struct LinearAdvection : public EquationBase<TypePack> {
	using Interface  = EquationBase<TypePack>;
	using State      = typename Interface::State;
	using Value      = typename Interface::Value;
	using Derivative = typename Interface::Derivative;
	using Time       = typename Interface::Time;
	using Jacobian   = typename Interface::Jacobian;

	LinearAdvection() = delete;
	LinearAdvection(const LinearAdvection& la) = default;
	virtual ~LinearAdvection() = default;
	LinearAdvection& operator=(const LinearAdvection& la) = default;

	LinearAdvection(const Value& speed, const std::vector<Value>& grid) :
		wave_speed_(speed), x_(grid){
	}

	void apply_bc(State& u, const Time& dt){
		u[0] = u[u.size()-2];
		u[u.size()-1] = u[1];
	}

protected:

	void dt_impl(State& u, Time& dt){
		using std::min;
		const Value cfl = 1.0;
		Value dxmin = 9999999;
		for(int i = 1; i < x_.size(); ++i){
			dxmin = min(dxmin, (x_[i] - x_[i-1]));
		}
		dt = cfl * dxmin/wave_speed_;
	}

	/**
	 * Calculate the dudt derivative using state u
	 *
	 * Method to calculate the derivative of statue u with respect to time t.
	 * - Resizes dudt to be same size as u
	 * - Applies BC's to state u
	 * - Calculates dudt using a 1st order up winding
	 * - Applies BC's to dudt (needed ?)
	 */
	void dudt_impl(State& u, Derivative& dudt, const Time& t){
		dudt.resize(u.size(),0.0);
		using std::abs;
		//this->apply_bc(u,t); // Update BC's before calculating dudt
		for(int i = 1; i < x_.size()-1; ++i){
			const Value ws_pls = 0.5*(wave_speed_+abs(wave_speed_));
			const Value ws_mns = 0.5*(wave_speed_-abs(wave_speed_));
			Value dudt_pls = ws_pls*(u[i]-u[i-1])/(x_[i]-x_[i-1]);
			Value dudt_mns = ws_mns*(u[i+1]-u[i])/(x_[i+1]-x_[i]);
			dudt[i] = -1.0 * (dudt_pls + dudt_mns);
		}
		this->apply_bc(dudt,t); // Not always valid
	}

	void dfdu_impl(State& u, Jacobian& dfdu, const Time& t){

	}

private:
	std::vector<Value> x_;
	Value              wave_speed_;
};


/**
 * Static functions to create 1-D grid.
 */
struct HelperFunctions {

	static std::vector<double> createGrid(double X0, double X1, int N){
		std::vector<double> X(N+3); // 1 layer of ghost cells +2
		double dx = (X1-X0)/N;
		for(int i = 0; i < N+3; ++i){
			X[i] = X0 + ((i-1) * dx);
		}
		return X;
	}

	template<typename T>
	static void initConditions(T& u){
		const int sz = u.size();
		for(auto& val : u){
			val = 0;
		}
		for(int i = 1; i < sz/2+1; ++i){
			u[i] = (i-1) * 1.0/static_cast<double>((sz-2)/2);
		}
	}

};





int main(){

	// Get types used for equations and solver
	using MyTypes = EquationTypes<std::valarray<double>>; // Define types used
	using EqnBase = EquationBase<MyTypes>;                // Equation Base Type
	using EqnImpl = LinearAdvection<MyTypes>;             // Equation Implementation Type
	using StpBase = StepperBase<EqnBase>;                 // Stepper Base Type
	using StpImpl = EulerStepper<EqnBase>;                // Stepper Implementation Type

	//using IntBase = IntegratorBase<StpBase>;              // Integrator Base Type
	//using IntImpl = StepIntegrator<StpBase>;              // Integrator Implementation Type


	// Create the Grid
	// - We create simple 1D array
	// - Application could read from disk or create using other methods
	// - Probably should return as a shared pointer to a grid base type
	int N = 20; // Grid dimensions
	auto grid = HelperFunctions::createGrid(0,1,N);

	// Pass the Grid and other equation parameters to Equation Constructor
	double wave_speed = +1.0;
	std::shared_ptr<EqnImpl> eqn_impl(new EqnImpl(wave_speed, grid));
	std::shared_ptr<EqnBase> eqn_base = eqn_impl;

	// Create the Stepper Implementation
	std::shared_ptr<StpImpl> stp_impl(new StpImpl(eqn_base));
	std::shared_ptr<StpBase> stp_base = stp_impl;

	// Create the Integrator Implementation
	//std::shared_ptr<IntImpl> int_impl(new IntImpl(eqn_base,stp_base));
	//std::shared_ptr<IntBase> int_base = int_impl;

	const int MaxSteps = 21; // <-- 21 is one full lap
	typename MyTypes::State u(grid.size());
	typename MyTypes::Time  t  = 0;
	typename MyTypes::Time  dt = grid[1] - grid[0];

	HelperFunctions::initConditions(u);
	auto u_init = u;

	//
	// Complete one full lap around grid
	//
	//Integrator<StpBase>::time(stp_base, u, t0, t1, dt, obs);
	//Integrator<StpBase>::steps(stp_base, u, t, dt, MaxSteps, obs);
	//Integrator<StpBase>::list(stp_base, u, tlist, obs);
	for(int i = 0; i < MaxSteps; ++i){
		stp_base->step(u,t,dt);
		t += dt;
	}


	double inf_err = std::abs(u_init-u).max();
	std::cout << "Max Error = " << inf_err << std::endl;
	assert(inf_err < 1.0e-13);

	return 0;
}


