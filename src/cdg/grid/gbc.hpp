//==================================================================================
// Module       : gbc
// Date         : 2/6/19 (DLR)
// Description  : GeoFLOW CDG boundary condition object
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#if !defined(_GBC_HPP)
#define _GBC_HPP 

#include "gtvector.hpp"
#include "ggrid.hpp"



class GBC 
{
public:
                             GBC() = delete;
                             GBC(GGrid &grid);
                             GBC(const GBC &gbc) = default;
                            ~GBC();

      void                   update_bdy_state(GFTYPE &t, GTVector<GTVector<GFTYPE>*> &u,        
                             GTVector<GTVector<GFTYPE>*> &ub);                 // entry function
      void                   set_update_callback(
                             std::function<void(const GFTYPE &t, GTVector<GTVector<GFTYPE>*> &u,
                             GTVector<GTVector<GFTYPE>*> &ub)> callback)
                             {update_dirichlet_callback_ = callback; 
                              bupdatebc_ = TRUE;}                              // set bdy-update callback


private:
        void                  do_reflective_bcs(GGrid &grid, GFTYPE t, 
                              GTVector<GTVector<GFTYPE>*> &u, 
                              GTVector<GTVector<GFTYPE>*> &ub); 

        GBOOL                  bupdatebc_;                           // bdy update callback set?
        GGrid                 *grid_;                                // element list
        std::function<void(const GFTYPE &t, GTVector<GTVector<GFTYPE>*> &u, 
        GTVector<GTVector<GFTYPE>*> &ub)>
                               update_dirichlet_callback_;           // dirichlet bdy update callback function

};

#endif
