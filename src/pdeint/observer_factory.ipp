/*
 * observer_factory.ipp
 *
 *  Created on: Nov 29, 2018
 *      Author: bryan.flynt
 */


namespace geoflow {
namespace pdeint {


template<typename ET>
typename ObserverFactory<ET>::ObsBasePtr
ObserverFactory<ET>::build(const tbox::PropertyTree& ptree, EqnBasePtr& equation, Grid& grid){


	// Set the default observer type
	const std::string default_observer = "none";

	// Get the type of observer
	const std::string observer_name = ptree.getValue("observer_name", default_observer);

        // Set traits from prop tree:
        typename ObserverBase<ET>::Traits traits;

        // Get whether 'observation' cadence should be by cycle or time:
        std::string stype = ptree.getValue<std::string>("cadence_type","none");
        if      ( "cycle" == stype )  traits.itype = ObserverBase<ET>::OBS_CYCLE;
        else if ( "time"  == stype )  traits.itype = ObserverBase<ET>::OBS_TIME;
        else EH_ERROR("Invalid observer type specified");

        std::vector<std::string> defst_names = {"u1","u2","u3","u4","u5","u6","u7","u8","u9"};
        std::vector<std::string> defgr_names = {"xgrid","ygrid","zgrid"};;
        std::vector<int> def_ids;
        traits.state_index   = ptree.getArray<int>        ("state_index",def_ids);    // state ids to 'observe' [0, 1, 2...]
        traits.state_names   = ptree.getArray<std::string>("state_names",defst_names);  // state names 
        traits.grid_names    = ptree.getArray<std::string>("grid_names",defgr_names);   // grid comp names 
        traits.cycle_interval= ptree.getValue<size_t>     ("cycle_interval", 10);       // cadence for cycle type
        traits.time_interval = ptree.getValue<double>     ("time_interval", 1.0);       // cadence for time type
        traits.freq_fact     = ptree.getValue<double>     ("interval_freq_fact", 1.0);  // freq factor relative to, say restart
        traits.idir          = ptree.getValue<std::string>("indirectory",".");          // input directory
        traits.odir          = ptree.getValue<std::string>("outdirectory",".");         // outputdirectory
        traits.start_ocycle  = ptree.getValue<size_t>     ("start_ocycle",0);           // starting output cycle 
//      traits.start_cycle   = ptree.getValue<size_t>     ("start_cycle",0);            // start evol cycle
        traits.start_time    = ptree.getValue<double>     ("start_time",0.0);           // start evol time
     
	// Create the observer and cast to base type
	ObsBasePtr base_ptr;
	if( "none" == observer_name ){
		using ObsImpl = NullObserver<ET>;

		// Allocate observer Implementation
		std::shared_ptr<ObsImpl> obs_impl(new ObsImpl(equation, grid));

		// Set any parameters we need to set
		// NA

		// Set back to base type
		base_ptr = obs_impl;
	}
        else if( "posixio_observer" == observer_name ) {
		using ObsImpl = GPosixIOObserver<ET>;

                traits.itag1  = ptree.getValue <GINT>("time_field_width",6);  
                traits.itag2  = ptree.getValue <GINT>("task_field_width",5);  
                traits.itag3  = ptree.getValue <GINT>("filename_size",2048);  

		// Allocate observer Implementation
		std::shared_ptr<ObsImpl> obs_impl(new ObsImpl(equation, grid, traits));

		// Set back to base type
		base_ptr = obs_impl;
        }
        else if( "global_diag_basic" == observer_name ) {
		using ObsImpl = GGlobalDiag_basic<Equation>;

		// Allocate observer Implementation
		std::shared_ptr<ObsImpl> obs_impl(new ObsImpl(equation, grid, traits));

		// Set back to base type
		base_ptr = obs_impl;
        }
	else {
		EH_ERROR("Requested observer not found: " << observer_name);
	}

	return base_ptr;
}


} // namespace pdeint
} // namespace geoflow

