/*
 * equation_factory.ipp
 *
 *  Created on: July 8, 2019 
 *      Author: bryan.flynt, d.rosenberg
 */

namespace geoflow {
namespace pdeint {


template<typename ET>
typename EquationFactory<ET>::EqnBasePtr
EquationFactory<ET>::build(const tbox::PropertyTree& ptree, Grid& grid){

        int itmp;

	// Set the default eqution type
	const std::string default_equation = "none";

	// Get the equation name
	std::string equation_name = ptree.getValue("pde_name", default_equation);

        // Set traits from prop tree:
        typename GBurgers<ET>::Traits btraits;
        typename GMConv  <ET>::Traits ctraits;

	// Set the default state components to force:
	std::vector<int> comps, default_comps;
	std::vector<GFTYPE> dstd;

        if ( !ptree.isPropertyTree(equation_name) ) {
          cout << "EquationFactory::build: PropertyTree does not exist: " << equation_name << endl;
        }
        PropertyTree eqn_ptree = ptree.getPropertyTree(equation_name);
        PropertyTree stp_ptree = ptree.getPropertyTree("stepper_props");
        PropertyTree dis_ptree = ptree.getPropertyTree("dissipation_traits");

     
	// Create the equation and cast to base type
	EqnBasePtr base_ptr;
	if( "pde_burgers" == equation_name ){
		using EqnImpl = GBurgers<ET>;

                btraits.doheat    = eqn_ptree.getValue<bool>  ("doheat",false);
                btraits.bpureadv  = eqn_ptree.getValue<bool>  ("bpureadv",false);
                btraits.bconserved= eqn_ptree.getValue<bool>  ("bconserved",false);
                btraits.bforced   = eqn_ptree.getValue<bool>  ("use_forcing",false);
                btraits.variabledt= stp_ptree.getValue<bool>  ("variable_dt",false);
                btraits.nsolve    = btraits.doheat || btraits.bpureadv ? 1 : GDIM; 
                btraits.nstate    =  btraits.nsolve
                                  + (btraits.bpureadv ? GDIM : 0);
                btraits.courant   = stp_ptree.getValue<double>("courant",0.5);
                btraits.itorder   = stp_ptree.getValue<int>   ("time_deriv_order",4);
                btraits.nstage    = stp_ptree.getValue<int>   ("nstage",4);
                btraits.bSSP      = stp_ptree.getValue<int>   ("stab_preserving",false);
                btraits.inorder   = stp_ptree.getValue<int>   ("extrap_order",2);
                btraits.ssteptype = stp_ptree.getValue<std::string>
                                                             ("stepping_method","GSTEPPER_EXRK");
                btraits.nu        = dis_ptree.getValue<double>("nu");
                for ( auto i=0; i<GDIM; i++ ) default_comps.push_back(i);
                comps            = eqn_ptree.getArray<int>   ("forcing_comp",default_comps);
                btraits.iforced.resize(comps.size());
                btraits.iforced   = comps; // traits.iforced may be a different d.structure

		// Allocate equation Implementation
		std::shared_ptr<EqnImpl> eqn_impl(new EqnImpl(grid, btraits));

		// Set back to base type
		base_ptr = eqn_impl;
	}
	else if( "pde_mconv" == equation_name ) {
		using EqnImpl = GMConv<ET>;

                int nc = grid.gtype() == GE_2DEMBEDDED ? 3 : GDIM;
                ctraits.docoriolis  = eqn_ptree.getValue<bool>  ("docoriolis");
                ctraits.dodry       = eqn_ptree.getValue<bool>  ("dodry");
                ctraits.dofallout   = eqn_ptree.getValue<bool>  ("dofallout");
                ctraits.dograv      = eqn_ptree.getValue<bool>  ("dogravity");
                ctraits.usebase     = eqn_ptree.getValue<bool>  ("usebase_state");
                ctraits.divopcolloc = eqn_ptree.getValue<bool>  ("divopcolloc");
                ctraits.usebdydata  = eqn_ptree.getValue<bool>  ("usebdydat",true);
                ctraits.Stokeshyp   = eqn_ptree.getValue<bool>  ("Stokeshyp");
                ctraits.bindepdiss  = eqn_ptree.getValue<bool>  ("bindepdiss",false);
                ctraits.nbase       = ctraits.usebase ? 2 : 0;
                ctraits.nlsector    = eqn_ptree.getValue<bool>  ("nliq",0);
                ctraits.nisector    = eqn_ptree.getValue<bool>  ("nice",0);
                ctraits.nsolve      = nc + 2                 // mom + denTot + energy_den
                                    + (!ctraits.dodry ? 1 : 0) // vapor
                                    + ( ctraits.dofallout ? ctraits.nlsector 
                                                        + ctraits.nisector : 0); // q_i
                ctraits.nstate      =  ctraits.nsolve
                                    + (ctraits.dofallout ? ctraits.nlsector 
                                                       + ctraits.nisector : 0)
                                    + (ctraits.usebase ? 2 : 0);
                ctraits.bconserved  = eqn_ptree.getValue<bool>  ("bconserved",false);
                ctraits.bforced     = eqn_ptree.getValue<bool>  ("use_forcing",false);
                ctraits.Ts_base     = eqn_ptree.getValue<double>("T_surf"); // K
                ctraits.P0_base     = eqn_ptree.getValue<double>("P0"); // millibar = hPa
                ctraits.P0_base     *= 100.0; // convert from mb to Pa
                ctraits.variabledt  = stp_ptree.getValue<bool>  ("variable_dt",false);
                ctraits.bvarvterm   = stp_ptree.getValue<bool>  ("variable_term_vel",false);
                ctraits.itorder     = stp_ptree.getValue<int>   ("time_deriv_order",4);
                ctraits.nstage      = stp_ptree.getValue<int>   ("nstage",4);
                ctraits.bSSP        = stp_ptree.getValue<int>   ("stab_preserving",false);
                ctraits.inorder     = stp_ptree.getValue<int>   ("extrap_order",2);
                ctraits.courant     = stp_ptree.getValue<double>("courant",0.5);
                ctraits.ssteptype   = stp_ptree.getValue<std::string>
                                                             ("stepping_method","GSTEPPER_EXRK");
                ctraits.nu          = dis_ptree.getValue<double>("nu");
                ctraits.kappa       = dis_ptree.getValue<double>("kappa");
                ctraits.zeta        = dis_ptree.getValue<double>("zeta");
                ctraits.lambda      = dis_ptree.getValue<double>("lambda");
                for ( auto i=0; i<GDIM; i++ ) default_comps.push_back(i);
                comps               = eqn_ptree.getArray<int>   ("forcing_comp",default_comps);
                ctraits.iforced.resize(comps.size());
                ctraits.iforced     = comps; // traits.iforced may be a different d.structure
                if ( ctraits.docoriolis ) {
                  dstd              = eqn_ptree.getArray<GFTYPE>("omega");
                } 
                else {
                  dstd.resize(0);
                }
                ctraits.omega.resize(dstd.size());
                ctraits.omega       = dstd; 

		// Allocate equation Implementation
		std::shared_ptr<EqnImpl> eqn_impl(new EqnImpl(grid, ctraits));

                // Configure filter list:
                int nsolve = eqn_impl->solve_size();
                eqn_impl->get_filter_list().resize(nsolve);
                EquationFactory<ET>::config_filters(ptree, equation_name, grid, eqn_impl->get_filter_list());
                
		// Set back to base type
		base_ptr = eqn_impl;
	}
	else {
		EH_ERROR("Requested equation not found: " << equation_name);
	}

	return base_ptr;
}


template<typename ET>
void
EquationFactory<ET>::config_filters(const PropertyTree& ptree, const std::string equation_name,  Grid& grid, FilterList& filter_list) {

        // Config blks for each filter are contained as sub-trees
        // within the equation-configuration block, so...
        // ...Get equation_ptree:
        if ( !ptree.isPropertyTree(equation_name) ) {
          cout << "EquationFactory::config_filters: PropertyTree does not exist: " << equation_name << endl;
        }
        PropertyTree eqn_ptree = ptree.getPropertyTree(equation_name);
        PropertyTree fptree;

        std::vector<std::string> vsdef, fblocklist;


        // On entry, filter_list should be sized to the number
        // of solved-for state-components and into NULL:
        for ( auto j=0; j<filter_list.size(); j++ ) filter_list[j] = nullptr;

        // Cycle over each filter config block, and configure
        // required filter:
        std::vector<int> icomp; // state comp ids being filtered
        FilterBasePtr pfilter;
        fblocklist = eqn_ptree.getArray<std::string>("filter_list", vsdef);
        for ( auto j=0; j<fblocklist.size(); j++ ) { // cycle over filters
          if ( !eqn_ptree.isPropertyTree(fblocklist[j]) ) {
            cout << "EquationFactory::config_filters: PropertyTree doesn't exist: " << fblocklist[j] << endl;
            exit(1);
          }
          fptree  = eqn_ptree.getPropertyTree(fblocklist[j]);
          icomp   = fptree.getArray<int>("state_index");
          pfilter = FilterFactory<ET>::build(eqn_ptree, fblocklist[j], grid);

          // Cycle over specified state indices, and set 
          // filter_list for those indices to this filter:
          for ( auto i=0; i<icomp.size(); i++ ) {
            if ( icomp[i] < 0 || icomp[i] >= filter_list.size() ) {
              cout << "EquationFactory::config_filters: state index " << icomp[i] << " invalid in PropertTree: " << fblocklist[j] << endl;
              exit(1);
            }
            filter_list[icomp[i]] = pfilter;
          }
        }


}




} // namespace pdeint
} // namespace geoflow

