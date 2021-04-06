/*
 * integrator_factory.ipp
 *
 *  Created on: Nov 29, 2018
 *      Author: bryan.flynt
 */


namespace geoflow {
namespace pdeint {

template<typename ET>
typename IntegratorFactory<ET>::IntegratorPtr
IntegratorFactory<ET>::build(const tbox::PropertyTree& ptree, const EqnBasePtr& eqn, const MixerBasePtr& mixer, const ObsBasePtr& obs,
                                                                               Grid&       grid){

	// Get the integrator traits
        PropertyTree subtree;
	typename Integrator<ET>::Traits traits;
        std::string stype;

        subtree = ptree.getPropertyTree  ("time_integration");

        // Get integrator cadence type: integrate by time or by cycle:
	stype             = subtree.getValue<std::string>("integ_type", "cycle" );
        if ( "cycle" == stype ) traits.integ_type = Integrator<ET>::INTEG_CYCLE;
        if ( "time"  == stype ) traits.integ_type = Integrator<ET>::INTEG_TIME;
	traits.cycle_end  = subtree.getValue<size_t>("cycle_end", 1);
	traits.dt         = subtree.getValue<Ftype>("dt",  std::numeric_limits<Time>::lowest() );
	traits.time_end   = subtree.getValue<Ftype>("time_end",  std::numeric_limits<Time>::max() );

	// Allocate Integrator Implementation
	IntegratorPtr integrator_ptr(new Integrator<ET>(eqn, mixer, obs, grid, traits));

	// Set any parameters we need to set
	// NA

	// Return
	return integrator_ptr;
}


} // namespace pdeint
} // namespace geoflow
