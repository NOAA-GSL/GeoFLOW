//==================================================================================
// Module       : gburgers.hpp
// Date         : 10/18/18 (DLR)
// Description  : Object defining a multidimensional Burgers (advection-diffusion) 
//                PDE. 
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#if !defined(_GBURGERS_HPP)
#define _GBURGERS_HPP

#include "gtypes.h"
#include <functional>
#include "gtvector.hpp"
#include "gdd_base.hpp"
#include "ggrid.hpp"

class GBurgers 
{

public:
                            GBurgers(GGrid &grid, GTVector<GFTYPE>> &u, GTVector<GFTYPE> &v, GTVector<GFTYPE> &tmp);
                            GBurgers(GGrid &grid, GTVector<GFTYPE>> &u, GTVector<GFTYPE> &v, GTVector<GFTYPE> &w, GTVector<GFTYPE>> &tmp);
                           ~GBurgers();

        GINT                step();                                      // take a time step
        void                set_bdy_callback(
                            std::function<void(GGrid &)> &callback);     // set bdy-set callback


friend  std::ostream&       operator<<(std::ostream&, GBurgers &);       // Output stream operator
 

protected:
         void               init2d();                                       // initialize base icosahedron for 2d grid
         void               init3d();                                       // initialize for 3d grid
       

private:

std::function<void(GGrid&)>
                       *bdycallback_ ; // callback object+method to set bdy conditions

};

#endif
