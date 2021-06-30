
#ifndef TIMER_H_
#define TIMER_H_



#include <chrono>
#include <cstdint>
#include <ostream>
#include <type_traits>

#include <sys/times.h>
#include <unistd.h>

namespace geoflow {
namespace tbox {


/**
 * @brief
 * Detect the highest resolution clock between the system and steady
 */
using maxres_sys_or_steady =
		std::conditional<
		std::chrono::steady_clock::period::den <= std::chrono::system_clock::period::den,
		std::chrono::steady_clock, std::chrono::system_clock>::type;

/**
 * @brief
 * Detect the highest resolution non-sleeping clock
 */
using maxres_nonsleep_clock =
		std::conditional<
		std::chrono::high_resolution_clock::is_steady,
		std::chrono::high_resolution_clock, maxres_sys_or_steady>::type;



namespace detail {

/**
 * Get the tick factor
 *
 * The tick_factor calculated as the number of
 * nanoseconds per system tick.
 */
inline
std::chrono::nanoseconds::rep tick_factor(){
	std::int64_t factor = ::sysconf( _SC_CLK_TCK );
	factor = static_cast<std::int64_t>(1000000000l) / factor;
    return factor;
}

} /* namespace detail */


/// Real/Wall CPU Clock
/**
 */
class process_real_cpu_clock {
public:
	typedef std::chrono::nanoseconds                        duration;
	typedef duration::rep                                   rep;
	typedef duration::period                                period;
	typedef std::chrono::time_point<process_real_cpu_clock> time_point;
	static constexpr bool is_steady =                       true;

	static inline time_point now() noexcept {
	    tms tm;
	    clock_t c = ::times( &tm );
	    return time_point(std::chrono::nanoseconds(c*detail::tick_factor()));
	}
};

/// User CPU Clock
/**
 */
class process_user_cpu_clock {
public:
	typedef std::chrono::nanoseconds                        duration;
	typedef duration::rep                                   rep;
	typedef duration::period                                period;
	typedef std::chrono::time_point<process_user_cpu_clock> time_point;
	static constexpr bool is_steady =                       true;

	static inline time_point now() noexcept {
	    tms tm;
	    clock_t c = ::times( &tm );
	    return time_point(std::chrono::nanoseconds((tm.tms_utime + tm.tms_cutime)*detail::tick_factor()));
	}
};

/// System CPU Clock
/**
 */
class process_system_cpu_clock {
public:
	typedef std::chrono::nanoseconds                          duration;
	typedef duration::rep                                     rep;
	typedef duration::period                                  period;
	typedef std::chrono::time_point<process_system_cpu_clock> time_point;
	static constexpr bool is_steady =                         true;

	static inline time_point now() noexcept {
	    tms tm;
	    clock_t c = ::times( &tm );
	    return time_point(std::chrono::nanoseconds((tm.tms_stime + tm.tms_cstime)*detail::tick_factor()));
	}
};


/// Data structure to hold real, user, system times
/**
 * Simple structure to hold real, user, system times
 * and provide basic add, subtract operations.
 *
 * TODO:
 * Create a "process_all_cpu_clock" which will calculate
 * and return a cpu_times structure with calling ::times( &tm )
 * multiple times.  Seems easy enough but getting the types
 * such as time_point, rep, duration, correct is a mystery
 * at this point.
 */
struct cpu_times {
	using size_type = std::int64_t;

	size_type wall;
    size_type user;
    size_type system;

    cpu_times& operator+=(cpu_times const& other) noexcept {
    	wall   += other.wall;
    	user   += other.user;
    	system += other.system;
    	return *this;
    }

    cpu_times& operator-=(cpu_times const& other) noexcept {
    	wall   -= other.wall;
    	user   -= other.user;
    	system -= other.system;
    	return *this;
    }

    cpu_times operator+(cpu_times const& other) const noexcept {
    	cpu_times ans(*this);
    	ans += other;
    	return ans;
    }

    cpu_times operator-(cpu_times const& other) const noexcept {
    	cpu_times ans(*this);
    	ans -= other;
    	return ans;
    }
};


/// Get current CPU times
/**
 * TODO:
 * Create a "process_all_cpu_clock" which will calculate
 * and return a cpu_times structure with calling ::times( &tm )
 * multiple times.  Seems easy enough but getting the types
 * such as time_point, rep, duration, correct is a mystery
 * at this point.
 */
inline
void get_cpu_times(cpu_times& current) noexcept {
	//current.wall   = xstd::chrono::process_real_cpu_clock::now();
	//current.user   = xstd::chrono::process_user_cpu_clock::now();
	//current.system = xstd::chrono::process_system_cpu_clock::now();
    tms tm;
    const clock_t c = ::times( &tm );
    const auto tick = detail::tick_factor();
	current.wall   = std::chrono::nanoseconds(c*tick).count();
	current.user   = std::chrono::nanoseconds((tm.tms_utime + tm.tms_cutime)*tick).count();
	current.system = std::chrono::nanoseconds((tm.tms_stime + tm.tms_cstime)*tick).count();
}


/// Write to ostream the current cpu_times
/**
 * Write an elapsed cpu_times value to the provided ostream. Writing
 * a value directly returned from get_cpu_times() will print but the
 * values will be those since some time in the past. To display
 * meaningful times the cpu_times passed in need to be a difference
 * such as the one returned from xstd::chrono::Timer::elapsed().
 *
 * \param current Current cpu_times to write
 * \param os Output Stream to write to
 * \param precision Number of decimal place to display for seconds
 */
inline
void show_time(cpu_times const& current, std::ostream& os, const short precision = 2) noexcept {
	const double sec = 1.0e9;
	const double wall_sec  = static_cast<double>(current.wall)   / sec;
	const double user_sec  = static_cast<double>(current.user)   / sec;
	const double syst_sec  = static_cast<double>(current.system) / sec;
	const double total_sec = static_cast<double>(current.system + current.user) / sec;

	//os.setf(std::ios_base::fixed, std::ios_base::floatfield);
	os.setf(std::ios_base::scientific, std::ios_base::floatfield);
	os.precision(precision);
	os << wall_sec  << "s wall, ";
	os << user_sec  << "s user + ";
	os << syst_sec  << "s system = ";
	os << total_sec << "s CPU ";

    os.setf(std::ios_base::fixed, std::ios_base::floatfield);
	os.precision(1);
	const auto perc_tol = static_cast<double>(0.001);
	if( (wall_sec > perc_tol) && (total_sec > perc_tol) ){
		os << "(" << (100.0*total_sec/wall_sec) << "%)\n";
	}
	else {
		os << "(n/a)\n";
	}
}



/// Timer class for calculating the time of things
/**
 * Timer class which calculates the real, user and system
 * times between start/stop calls.
 */
class Timer {

public:

	Timer() {
		this->start();
	}

	/**
	 * Return if the timer is currently stopped
	 */
	bool is_stopped() const noexcept {
		return is_stopped_;
	}

	/**
	 * Start Timer from Zero
	 */
	void start() noexcept {
		is_stopped_ = false;
		get_cpu_times(times_);
	}

	/**
	 * Stop Timer
	 */
	void stop() noexcept {
		if(this->is_stopped()) {
			return;
		}

		cpu_times current;
		get_cpu_times(current);
		times_ = (current - times_);
		is_stopped_ = true;
	}

	/**
	 * Resume Timer from current count
	 */
	void resume() noexcept {
		if(this->is_stopped()) {
			cpu_times current(times_);
			this->start();
			times_ -= current;
		}
	}

	/**
	 * Returns current timings but does not stop timer
	 */
	cpu_times elapsed() const noexcept {
		if(this->is_stopped()) {
			return times_;
		}
		cpu_times current;
		get_cpu_times(current);
		current -= times_;
		return current;
	}

	/**
	 * Display the current elapsed timings
	 */
	void display(std::ostream& os) noexcept {
		show_time(this->elapsed(),os);
	}

protected:
	cpu_times  times_;
	bool       is_stopped_;
};





} // namespace tbox
} // namespace geoflow





#endif // TIMER_H_