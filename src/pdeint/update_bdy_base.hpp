/*
 * update_bdy_base.hpp
 *
 *  Created on: June 25, 2020 
 *      Author: d.rosenberg
 */

#ifndef SRC_PDEINT_UPDATEBDYBASE_HPP_
#define SRC_PDEINT_UPDATEBDYBASE_HPP_


#include <memory>
#include <vector>
#include "equation_base.hpp"

namespace geoflow {
namespace pdeint {

/**
 * Base to handle specification and time updates of bdy information 
 * (or possibly state due to bdy) 
 *
 */
template<typename TypePack>
class UpdateBdyBase {

public:

        using Types      = TypePack;
        using State      = typename Types::State;
        using StateComp  = typename Types::StateComp;
        using StateInfo  = typename Types::StateInfo;
        using Grid       = typename Types::Grid;
        using Mass       = typename Types::Mass;
        using Ftype      = typename Types::Ftype;
        using Derivative = typename Types::Derivative;
        using Time       = typename Types::Time;
        using CompDesc   = typename Types::CompDesc;
        using Jacobian   = typename Types::Jacobian;
	using IBdyVol    = typename Types::IBdyVol;
	using TBdyVol    = typename Types::TBdyVol;
        using Size       = typename Types::Size;
        using EqnBase    = typename Types::EqnBase;
        using EqnBasePtr = typename Types::EqnBasePtr;

#if 0
        struct Traits {
          bool   compute_once=FALSE;
          bool       use_init=FALSE;
          int                  idir;    // // canonical coord direction definining surfaces 
          int                 bdyid;    // bdy id
          Ftype               xstart;   // number defining sponge surface start
          std::vector<Ftype>  value;    // Diriclet value for each istate 
          std::vector<int>    istate;   // state indices to operate on
          std::vector<size_t> ibdyvol;  // indir. inidices into comput volume
          std::vector<size_t> ibdyloc;  // inidices of ubdyvol in global bdy array, grid->igbdy
          std::vector<unsigned int> 
                              ibdysrc;  // bdy descriptor
          std::vector<Ftype>  farfield; // far-field solution for each istate
          std::vector<Ftype>  falloff;  // fall-off rate for solution
          std::vector<Ftype>  exponent; // decay exponents 
          std::vector<size_t> isponge;  // contains grid indices of sponge layer points; 
          std::string         smethod;  // init method if use_init==TRUE
          std::string         sconfig;  // ptree config method

          std::function<bool(const PropertyTree& ptree,
                              std::string &sconfig,
                              EqnBasePtr  &eqn,
                              Grid        &grid,
                              Time        &time,
                              const int    id,
                              State       &utmp,
                              State       &u,
                              State       &ub)> callback = NULLPTR;
          geoflow::tbox::PropertyTree ptree;

        }
#endif
      
	UpdateBdyBase() = default;
	UpdateBdyBase(const UpdateBdyBase& I) = default;
	~UpdateBdyBase() = default;
	UpdateBdyBase& operator=(const UpdateBdyBase& I) = default;

	/**
	 * Update bdy conditions with state at t
	 *
	 * @param[in,out] eqn    : equation base (not the shared pointer)
	 * @param[in,out] grid   : initial time at start, and final time
	 * @param[in,out] time   : current time
	 * @param[in,out] utmp   : tmp arrays
	 * @param[in,out] u      : current state array
	 */
	bool update (EqnBasePtr &eqn,
                     Grid       &grid, 
                     Time       &time, 
                     State      &utmp, 
                     State      &u){ 
                        return this->update_impl(eqn, grid, time, utmp, u);
                     }

	/**
	 * Get list of state ids that bdy condition acts upon
	 *
	 */
	std::vector<int> 
             get_istate(){
               return get_istate_impl();
             }

#if 0
        /**
         * Get traits.
         *
         */
        Traits &get_traits() {return traits_;}
#endif


protected:

#if 0
        Traits       traits_;
#endif
	virtual bool update_impl (
                     EqnBasePtr &eqn,
                     Grid       &grid, 
                     Time       &time, 
                     State      &utmp, 
                     State      &u) = 0; 

	virtual std::vector<int>&
             get_istate_impl() = 0;
};


} // namespace pdeint
} // namespace geoflow



#endif /* SRC_PDEINT_UPDATEBDYBASE_HPP_ */
