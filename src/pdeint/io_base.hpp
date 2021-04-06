/*
 * io_base.hpp
 *
 *  Created on: Jan 21, 2020 
 *      Author: d.rosenberg
 */

#ifndef SRC_PDEINT_IOBASE_HPP_
#define SRC_PDEINT_IOBASE_HPP_


#include <memory>
#include <vector>

namespace geoflow {
namespace pdeint {

/**
 * IO base to carry out state IO 
 *
 */
template<typename TypePack>
class IOBase {

public:
        enum GIOType        {GIO_POSIX=0, GIO_COLL}; // POSIX or collective
        using Types        = TypePack;
	using State        = typename Types::State;
	using StateInfo    = typename Types::StateInfo; // May contain time, time index, var name etc
	using Grid         = typename Types::Grid;
	using Ftype        = typename Types::Ftype;
        using Time         = typename Types::Time;
	using Size         = typename Types::Size;

        struct Traits {
          GIOType      io_type = GIO_COLL; // default to collective IO
          int          ivers   = 0;        // IO version tag
          bool         multivar= false;    // multiple vars in file (only of COLL types)?
          bool         prgrid  = false;    // flag to print grid
          int          wtime   = 6;        // time-field width
          int          wtask   = 5;        // task-field width (only for POSIX types)
          int          wfile   = 2048;     // file name max
          int          dim     = GDIM;     // problem dimension
          std::string  idir         ;      // input directory
          std::string  odir         ;      // output directory
        };

      
	IOBase() = delete;

	/**
	 * Constructor to initialize everything needed to do IO
	 *
	 */
	IOBase(Grid& grid, Traits& traits):
          grid_ (&grid), traits_ (traits) {}
	IOBase(const IOBase& I) = default;
	~IOBase() = default;
	IOBase& operator=(const IOBase& I) = default;

	/**
	 * Write state at t
	 *
	 * @param[in,out] t  Initial time at start, and final time
	 * @param[in,out] u  Current state values
	 */
	void write_state( std::string  filename,
                          StateInfo&   info,
                          const State&       u){
               return this->write_state_impl(filename, info, u);
             }
        /**
	 * Read state 
	 *
	 * @param[in    ] filename: filename if reading in aggregate
	 * @param[in,out] info    : StateInfo structure
	 * @param[in,out] u       : state values
	 * @param[in    ] bstate  : flag to set to FALSE if not reading state
	 */
	void read_state ( std::string  filename,
                          StateInfo&   info,
                          State&       u,
                          bool         bstate=true){
               return this->read_state_impl(filename, info, u, bstate);
             }
        /**
	 * Read state info
	 *
	 * @param[in    ] filename: filename 
	 * @param[in,out] info    : StateInfo structure
	 */
	void read_state_info ( std::string  filename,
                          StateInfo&   info ){
               return this->read_state_info_impl(filename, info);
             }
        /**
	 * Get traits
	 *
	 */
	Traits &get_traits( ) { 
               return traits_;
             }

protected:
        virtual void write_state_impl(std::string  filename,
                                      StateInfo&   info,
                                      const State& u) = 0;
        virtual void read_state_impl (std::string  filename,
                                      StateInfo&   info,
                                      State&       u,
                                      bool         bstate) = 0;
        virtual void read_state_info_impl (std::string  filename,
                                           StateInfo&   info) = 0;
        Grid   *grid_;
        Traits  traits_;
};


} // namespace pdeint
} // namespace geoflow



#endif /* SRC_PDEINT_IOBASE_HPP_ */
