/*
 * configure.hpp
 *
 *  Created on: Oct 23, 2018
 *      Author: bflynt
 */
#ifndef GEOFLOW_CONFIG_HPP_IN_
#define GEOFLOW_CONFIG_HPP_IN_


//====================================================
//                    Debug Flags
//====================================================

/** Turn on static assertions.
 */
#define GEOFLOW_USE_STATIC_ASSERTIONS

/** Turn on runtime assertions.
 */
#define GEOFLOW_USE_RUNTIME_ASSERTIONS

/** Turn on traceback when asserts.
 */
#define GEOFLOW_USE_TRACEBACK

/** Enable program tracer output.
 */
#define GEOFLOW_USE_TRACER

/** Insure empty loops are not optimized away
 */
#define GEOFLOW_INSURE_LOOP

/** Insure unused variables do not report as errors
 */
#define GEOFLOW_INSURE_USE

//====================================================
//                    Clock Settings
//====================================================

// Clock mechanism
//#define USE_POSIX_TIME
#define USE_BOOST_TIME
//#define USE_C11_TIME
//#define USE_C_TIME


#endif /* GEOFLOW_CONFIG_HPP_IN_ */
