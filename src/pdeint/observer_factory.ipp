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
ObserverFactory<ET>::build(const tbox::PropertyTree& ptree, const std::string obsname, EqnBasePtr& equation, Grid& grid){

        // Note: the entry ptree is the main property tree!


	// Get the prop tree specified by obsname:
        PropertyTree obstree = ptree.getPropertyTree(obsname);

	// Set the default observer type
	const std::string default_observer = "none";

	// Get the type of observer
	const std::string observer_name = obstree.getValue("observer_name", default_observer);

        // Set traits from prop tree:
        typename ObserverBase<ET>::Traits traits;
    
        // Set default derived quantity info:
        PropertyTree dqtree;
        std::vector<int>         defi;
        std::vector<std::string> defq;
        std::vector<std::string> dqnames;
        

        // Get whether 'observation' cadence should be by cycle or time:
        std::string stype = obstree.getValue<std::string>("cadence_type","none");
        if      ( "cycle" == stype )  traits.itype = ObserverBase<ET>::OBS_CYCLE;
        else if ( "time"  == stype )  traits.itype = ObserverBase<ET>::OBS_TIME;
        else if ( "none" != stype ) {
          cout << "ObserverFactory<ET>::build: stype=" << stype << endl;
          EH_ERROR("Invalid observer type specified");
        }

        std::vector<std::string> defst_names = {"u1","u2","u3","u4","u5","u6","u7","u8","u9"};
        std::vector<std::string> defgr_names = {"xgrid","ygrid","zgrid"};;
        std::vector<int> def_ids;
        traits.treat_as_1d   = obstree.getValue<bool>       ("treat_as_1d",false);        // treat-as-1d flag
        traits.state_index   = obstree.getArray<int>        ("state_index",defi);         // state ids to 'observe' [0, 1, 2...]
        traits.state_names   = obstree.getArray<std::string>("state_names",defst_names);  // state names 
        traits.grid_names    = obstree.getArray<std::string>("grid_names",defgr_names);   // grid comp names 
        traits.cycle_interval= obstree.getValue<size_t>     ("cycle_interval", 10);       // cadence for cycle type
        traits.time_interval = obstree.getValue<double>     ("time_interval", 1.0);       // cadence for time type
        traits.freq_fact     = obstree.getValue<double>     ("interval_freq_fact", 1.0);  // freq factor relative to, say restart
        traits.idir          = obstree.getValue<std::string>("indirectory",".");          // input directory
        traits.odir          = obstree.getValue<std::string>("outdirectory",".");         // outputdirectory
        traits.start_ocycle  = obstree.getValue<size_t>     ("start_ocycle",0);           // starting output cycle 
//      traits.start_cycle   = obstree.getValue<size_t>     ("start_cycle",0);            // start evol cycle
        traits.start_time    = obstree.getValue<double>     ("start_time",0.0);           // start evol time
     
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

                traits.itag1  = obstree.getValue <GINT>("time_field_width",6);  
                traits.itag2  = obstree.getValue <GINT>("task_field_width",5);  
                traits.itag3  = obstree.getValue <GINT>("filename_size",2048);  

                // Fill derived quantities strutures, if any:
                dqnames       = obstree.getArray<std::string> ("derived_quantities",defq);  // list of derived quantities to output
                traits.derived_quantities.resize(dqnames.size());
                for ( auto j=0; j<dqnames.size(); j++ ) {
                  dqtree = ptree.getPropertyTree(dqnames[j]); // get prop tree for named quantity
                  traits.derived_quantities[j].icomponents = dqtree.getArray<int>("state_index");
                  traits.derived_quantities[j].snames      = dqtree.getArray<std::string>("names"  ,defq);
                  traits.derived_quantities[j].smath_op    = dqtree.getValue<std::string>("mathop" ,"");
                }
                

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

