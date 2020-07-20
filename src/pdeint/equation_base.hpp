/*
 * equation_base.hpp
 *
 *  Created on: Nov 16, 2018
 *      Author: bflynt
 */

#ifndef SRC_PDEINT_EQUATION_BASE_HPP_
#define SRC_PDEINT_EQUATION_BASE_HPP_

#include <functional>


namespace geoflow {
namespace pdeint {

/**
 * Base class to provided an interface for all equation types.
 *
 * The base class provides the interface using the strategy pattern
 * for all equation implementations.
 */
template<typename TypePack>
class EquationBase {

public:
	using Types      = TypePack;
	using State      = typename Types::State;
        using StateInfo  = typename Types::StateInfo;
	using Grid       = typename Types::Grid;
	using Value      = typename Types::Value;
	using Derivative = typename Types::Derivative;
	using Time       = typename Types::Time;
	using Jacobian   = typename Types::Jacobian;
	using CompDesc   = typename Types::CompDesc;
	using Size       = typename Types::Size;


/*
        struct GStateInfo {
          int         sttype  = 1;        // state type index (grid=0 or state=1)
          int         gtype   = 0;        // check src/cdg/include/gtypes.h
          size_t      index   = 0;        // time index
          size_t      nelems  = 0;        // num elems
          size_t      cycle   = 0;        // continuous time cycle
          Value       time    = 0.0;      // state time
          std::vector<std::string>
                      svars;              // names of state members
          CompDesc    icomptype;          // encoding of state component types    
          GTMatrix<GINT>
                      porder;             // if ivers=0, is 1 X GDIM; else nelems X GDIM;
          std::string idir    = ".";      // input directory
          std::string odir    = ".";      // output directory
        };
        static_assert(std::is_same<StateInfo,GStateInfo>::value,
               "StateInfo is of incorrect type");
*/


	EquationBase() { update_bdy_callback_ = nullptr; }
	EquationBase(const EquationBase& eb) = default;
	virtual ~EquationBase() = default;
	EquationBase& operator=(const EquationBase& eb) = default;

	/**
	 * Return if the equation system can determine a time step.
	 *
	 * Query if the equation system can return a time step value
	 * which is equivalent to a CFL of one.
	 *
	 * @return Boolean value
	 */
	bool has_dt() const{
		return this->has_dt_impl();
	}

	/**
	 * Time difference corresponding to a CFL of one.
	 *
	 * Calculate the time step size which equates to a CFL of one
	 * using the modeled equation.
	 *
	 * @param[in]     t  Current time of state u before taking step
	 * @param[in,out] u  Is the state of the system of equations
	 * @param[out]    dt Size of time step corresponding to CFL of one
	 */
	void dt(const Time& t, State& u, Time& dt){
		this->dt_impl(t,u,dt);
	}

	/**
	 * Apply Boundary Condition.
	 *
	 * Apply the provided ub boundary condition to the state u.
	 *
	 * @param[in]     t  Current time of state u before taking step
	 * @param[in,out] u  Is the state of the system of equations
	 * @param[in]     ub Boundary condition values
	 */
	void apply_bc(const Time &t, State &u, State &ub){
		this->apply_bc_impl(t,u,ub);
	}


	/** Take one step from time t to time t+dt.
	 *
	 * Take exactly one time step from t with current solution u to
	 * time t + dt and return new solution within u.
	 *
	 * \param[in]     t  Current time of state u before taking step
	 * \param[in]     ub Boundary conditions for u
	 * \param[in,out] u  State of the system of equations
	 * \param[in,out] uf Forcing tendency
	 * \param[in,out] ub Boundary vector
	 * \param[in]     dt Size of time step to take
	 */
	void step(const Time& t, State& u, State& uf, State& ub, const Time &dt){
		this->step_impl(t,u,uf,ub,dt);
	}

	/** Take one step from time t to time t + dt.
	 *
	 * Take exactly one time step from t with current solution uin to
	 * time t + dt and return new solution within uout.
	 *
	 * \param[in]     t     Current time of state u before taking step
	 * \param[in]     uin   State of the system of equations at time t
	 * \param[in,out] uf    Forcing tendency
	 * \param[in,out] ub    Boundary conditions for u
	 * \param[in]     dt    Size of time step to take
	 * \param[out]    uout  New state of the system of equations at t + dt
	 */
	void step(const Time& t, const State& uin, State& uf, State& ub, const Time& dt, State& uout){
		this->step_impl(t,uin,uf,ub,dt,uout);
	}

	/** Return StateInfo data
         * 
         */
        StateInfo& stateinfo() {
                return stateinfo_;
        }


	/** Set boundary update callback function
	 *
	 *
	 * \param[in]     fcn   bdy update function
	 */
	virtual void set_bdy_update_callback(std::function<void(const Time& t, State& u, State& ub)> fcn){
		update_bdy_callback_ = fcn;
	}

protected:

	/**
	 * Must be provided by implementation
	 */
	virtual bool has_dt_impl() const = 0;

	/**
	 * Must be provided by implementation
	 */
	virtual void dt_impl(const Time& t, State& u, Time& dt) = 0;

	/**
	 * Must be provided by implementation
	 */
	virtual void apply_bc_impl(const Time &t, State &u, State &ub) = 0;

	/**
	 * Must be provided by implementation
	 */
	virtual void step_impl(const Time& t, State& u, State& uf,  State& ub, const Time& dt) = 0;

	/**
	 * Must be provided by implementation
	 */
	virtual void step_impl(const Time& t, const State& uin, State& uf, State& ub, const Time& dt, State& uout) = 0;
  
        std::function<void(const Time &t, State &u, State &ub)> 
                 update_bdy_callback_;

        StateInfo stateinfo_;

};


} // namespace pdeint
} // namespace geoflow


#endif /* SRC_PDEINT_EQUATION_BASE_HPP_ */
