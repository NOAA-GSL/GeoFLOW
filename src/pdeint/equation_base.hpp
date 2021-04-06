/*
 * equation_base.hpp
 *
 *  Created on: Nov 16, 2018
 *      Author: bflynt
 */

#ifndef SRC_PDEINT_EQUATION_BASE_HPP_
#define SRC_PDEINT_EQUATION_BASE_HPP_

#include <functional>
#include "filter_base.hpp"


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
	using StateComp  = typename Types::StateComp;
        using StateInfo  = typename Types::StateInfo;
        using Mass       = typename Types::Mass;
	using Grid       = typename Types::Grid;
	using Ftype      = typename Types::Ftype;
	using Derivative = typename Types::Derivative;
	using Time       = typename Types::Time;
	using Jacobian   = typename Types::Jacobian;
	using CompDesc   = typename Types::CompDesc;
	using Size       = typename Types::Size;
        using FilterBasePtr = std::shared_ptr<FilterBase<Types>>;
	using FilterList    = std::vector<FilterBasePtr>;


/*
        struct GStateInfo {
          int         sttype  = 1;        // state type index (grid=0 or state=1)
          int         gtype   = 0;        // check src/cdg/include/gtypes.h
          size_t      index   = 0;        // time index
          size_t      nelems  = 0;        // num elems
          size_t      cycle   = 0;        // continuous time cycle
          Ftype       time    = 0.0;      // state time
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


	EquationBase() : update_bdy_callback_ (nullptr) 
                       {}
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
	 * Apply the boundary condition to the state u.
	 *
	 * @param[in]     t  Current time of state u before taking step
	 * @param[in,out] u  Is the state of the system of equations
	 */
	void apply_bc(const Time &t, State &u){
		this->apply_bc_impl(t,u);
	}


	/** Take one step from time t to time t+dt.
	 *
	 * Take exactly one time step from t with current solution u to
	 * time t + dt and return new solution within u.
	 *
	 * \param[in]     t  Current time of state u before taking step
	 * \param[in,out] u  State of the system of equations
	 * \param[in,out] uf Forcing tendency
	 * \param[in]     dt Size of time step to take
	 */
	void step(const Time& t, State& u, State& uf,  const Time &dt){
		this->step_impl(t,u,uf,dt);
	}

	/** Take one step from time t to time t + dt.
	 *
	 * Take exactly one time step from t with current solution uin to
	 * time t + dt and return new solution within uout.
	 *
	 * \param[in]     t     Current time of state u before taking step
	 * \param[in]     uin   State of the system of equations at time t
	 * \param[in,out] uf    Forcing tendency
	 * \param[in]     dt    Size of time step to take
	 * \param[out]    uout  New state of the system of equations at t + dt
	 */
	void step(const Time& t, const State& uin, State& uf, const Time& dt, State& uout){
		this->step_impl(t,uin,uf,dt,uout);
	}

	/** Initialize solver
	 *
	 * \param[in]     u     state vectors
	 * \param[in]     utmp  tmp state vectors
	 */
	void init(State &u, State& utmp){
		this->init_impl(u, utmp);
	}

	/** Compute derived quantity 
	 *
	 * \param[in]     u     state vector
	 * \param[in]     sop   operation/quantity to compute
	 * \param[in]     utmp  tmp state vectors
	 * \param[out]    uout  resultant quantity
	 * \param[out]    iuout indices pointing to quantity(-ies) in uout
	 */
	void compute_derived(const State& u, const std::string sop, State& utmp, State& uout, std::vector<int>& iuout){
                 this->compute_derived_impl(u, sop, utmp, uout, iuout);
        }

  	/** Return state size 
         * 
         */
        int state_size(){
		return this->state_size_impl();
        }

  	/** Return solve size 
         * 
         */
        int solve_size(){
		return this->solve_size_impl();
        }

  	/** Return tmp size 
         * 
         */
        int tmp_size(){;
		return this->tmp_size_impl();
        }

  	/** Return state component ids that are forced
         * 
         */
        std::vector<int>& iforced(){
		return this->iforced_impl();
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
	virtual void set_bdy_update_callback(std::function<void(const Time& t, State& u)> fcn){
		update_bdy_callback_ = fcn;
	}


	/** Set filter list: 
         *    one filter for each evolved state comp
	 *    
	 *
	 * \param[in]     pFilter
	 */
	virtual void set_filter_list(FilterList &pFilter) {
		filterList_.resize(pFilter.size());;
		filterList_ = pFilter;
	}

	/** Get filter list: 
	 *    
	 *
	 * \param[out]     FilterList &
	 */
	virtual FilterList& get_filter_list() {
		return filterList_;
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
	virtual void apply_bc_impl(const Time &t, State &u) = 0;

	/**
	 * Must be provided by implementation
	 */
	virtual void step_impl(const Time& t, State& u, State& uf,  const Time& dt) = 0;

	/**
	 * Must be provided by implementation
	 */
	virtual void step_impl(const Time& t, const State& uin, State& uf, const Time& dt, State& uout) = 0;
  
	/**
	 * Must be provided by implementation
	 */
	virtual void init_impl(State& u, State& utmp) = 0;
  
	/**
	 * Must be provided by implementation
	 */
	virtual void compute_derived_impl(const State& u, const std::string sop, State& utmp, State& uout, std::vector<int>& iuout) = 0;
  
	/**
	 * Must be provided by implementation
	 */
	virtual int state_size_impl() = 0;
  
	/**
	 * Must be provided by implementation
	 */
	virtual int solve_size_impl() = 0;
  
	/**
	 * Must be provided by implementation
	 */
	virtual int tmp_size_impl() = 0;
  
	/**
	 * Must be provided by implementation
	 */
	virtual std::vector<int> &iforced_impl() = 0;
  
        std::function<void(const Time &t, State &u)> 
                 update_bdy_callback_;

        StateInfo  stateinfo_;
        FilterList filterList_;

};


} // namespace pdeint
} // namespace geoflow


#endif /* SRC_PDEINT_EQUATION_BASE_HPP_ */
