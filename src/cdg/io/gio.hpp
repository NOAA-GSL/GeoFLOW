//==================================================================================
// Module       : gio.hpp
// Date         : 1/20/20 (DLR)
// Description  : GIO object encapsulating methods for POSIX and collective IO
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : IOBase.
//==================================================================================
#if !defined(_GIO_HPP)
#define _GIO_HPP

#include "gtvector.hpp"
#include "ggrid.hpp"
#include "pdeint/io_base.hpp"
#include "tbox/property_tree.hpp"
#include "tbox/mpixx.hpp"

using namespace geoflow::pdeint;
using namespace std;

template<typename TypePack>
class GIO : public IOBase<TypePack>
{

public:
        enum GIOType        {GIO_POSIX=0, GIO_COLL}; // POSIX or collective
        using Interface   = EquationBase<TypePack>;
        using IOBaseType  = IOBase<TypePack>;
        using IOBasePtr   = std::shared_ptr<IOBaseType>;
        using State       = typename Interface::State;
        using Grid        = typename Interface::Grid;
        using Value       = typename Interface::Value;
        using Time        = typename Interface::Time;
        using Size        = typename Interface::Size;
        using StateInfo   = typename Interface::StateInfo; 

        struct Traits {
          GIOType     io_type = GIO_COLL; // default to collective IO
          GINT        ivers   = 0;        // IO version tag
          GBOOL       multivar= false;    // multiple vars in file (only of COLL types)?
          GINT        wtime   = 6;        // time-field width
          GINT        wtask   = 5;        // task-field width (only for POSIX types)
          GINT        wfile   = 2048;     // file name max
          GINT        dim     = GDIM;     // problem dimension
          GString     idir    = ".";      // input directory
          GString     odir    = ".";      // output directory
        }; 

        static_assert(std::is_same<Value,GFTYPE>::value,
               "Value is of incorrect type");
        static_assert(std::is_same<State,GTVector<GTVector<GFTYPE>*>>::value,
               "State is of incorrect type");
        static_assert(std::is_same<StateInfo,GStateIOTraits>::value,
               "StateInfo is of incorrect type");

                           GIO() = delete;
                           GIO(Grid &grid, Traits &traits, GC_COMM comm);
                          ~GIO(); 
                           GIO(const GIO &a) = default;
                           GIO &operator=(const GIO &bu) = default;

        void               write_state_impl(std::string filename, StateInfo &info, State &u);
        void               read_state_impl(std::string filename, StateInfo &info, State &u);
        Traits            &get_traits() { return traits_; }

private:
// Private methods:
        void               init(const Time &t, const State &u);
//      void               read_state_posix (StateInfo &info,       State  &u);
//      void               read_state_coll  (StateInfo &info,       State  &u);
        GSIZET             write_posix(GString filename, StateInfo &info, const GTVector<Value> &u);
        GSIZET             read_posix (GString filename, StateInfo &info,       GTVector<Value> &u);
        GSIZET             write_coll (GString filename, StateInfo &info, const State           &u);
        GSIZET             read_coll  (GString filename, StateInfo &info,       State           &u);
        GSIZET             read_header(GString filename, StateInfo &info, Traits &traits);
        GSIZET             write_header(GFPTR, StateInfo &info, Traits &traits);
        GSIZET             sz_header(StateInfo &info, Traits &traits);
        void               resize(GINT n);


// Private data:
        GBOOL              bInit_;      // is initialized?
        GINT               nfname_;
        GSIZET             nbheader_;   // # bytes in header
        GSIZET             nbfield_;    // # bytes in each state component
        GSIZET             nbgdof_;     // total # global dof bytes
        #if defined(_G_USE_MPI)
          MPI_Datatype     mpi_state_type_;
          MPI_Offset       state_disp_;
          MPI_Aint         state_extent_;
        #endif
        GC_COMM            comm_;
        std::stringstream  svarname_;
        GString            fname_;
        char              *cfname_;
        Traits             traits_;
        Grid              *grid_;
        GTVector<GINT>     state_index_;// list of state indices to print
        GTVector<GString>  state_names_;// list of names of state files
        GTVector<GString>  grid_names_ ;// list of names of grid comp files

};

#include "gio.ipp"

#endif
