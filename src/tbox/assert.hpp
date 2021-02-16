/*
 * Filename:	assert.hpp
 * Author:      bflynt
 * Created:		Jun 11, 2016
 * Copyright:   2016, Bryan Flynt, All rights reserved.
 */
#ifndef GEOFLOW_DEBUG_ASSERT_H_
#define GEOFLOW_DEBUG_ASSERT_H_

#include "configure.hpp"

/*
 * Macros for:
 * ASSERT - Assert condition evaluates to true at
 * runtime, otherwise exit program
 *
 * ASSERT_MSG - Assert condition evaluates to true
 * at runtime, otherwise exit program and report message
 *
 * VERIFY - Verify condition evaluates to true at
 * runtime, otherwise exit program.  Condition statements
 * are still evaluated even when ASSERT checking is turned
 * off
 *
 * VERIFY_MSG - Verify condition evaluates to true at
 * runtime, otherwise exit program and report message.
 * Condition statements are still evaluated even when
 * ASSERT checking is turned off.
 *
 * STATIC_ASSERT - Assert condition evaluates to true at
 * compile time.
 *
 * STATIC_ASSERT_MSG - Assert condition evaluates to true at
 * compile time, print message if not.
 *
 * NULL_USE - Wrap any used variable with macro to disable
 * compiler warning about unused variables
 */



// ------------------------------------------------------------------- //
//                       NULL Statement Definition                     //
// ------------------------------------------------------------------- //

/*!
 * The Intel compiler will not do the loop if there isn't a
 * statement inside so make a do nothing statement to insure it
 * does the loop.  Maybe it optimizes the loop away?
 */
#ifdef GEOFLOW_INSURE_LOOP
#define NULL_STATEMENT if (0) int nullstatement = 0
#else
#define NULL_STATEMENT
#endif

// ------------------------------------------------------------------- //
//                           NULL Use Definition                       //
// ------------------------------------------------------------------- //

/**
 * A null use of a variable, use to avoid GNU compiler
 * warnings about unused variables.
 */
#ifdef GEOFLOW_INSURE_USE
#define NULL_USE(EXP) do { \
       if(0) {char *temp = (char *)&variable; temp++;} \
    } while (0);
#else
#define NULL_USE(EXP)
#endif


//
// Force a Hang
// Infinite Loop to force a hang to attach debugger
//
#ifdef GEOFLOW_ASSERT_HANG
	#include "tbox/pio.hpp"
    #define MANUAL_FORCE_HANG								  	        \
		do {														    \
			pio::perr << std::endl;                                     \
			pio::perr << "******************************" << std::endl; \
			pio::perr << "********* FORCED HANG ********" << std::endl; \
			pio::perr << "******************************" << std::endl; \
			pio::perr << std::endl;                                     \
			while( true ) {                                             \
				if (0) int nullstatement = 0;                           \
			}                                 			                \
		} while (0);
#else
	#define MANUAL_FORCE_HANG
#endif

//
// Trigger a Core Dump
// Dumps core but keeps running
//
#ifdef GEOFLOW_ASSERT_CORE
	#include <sys/types.h>
    #include <unistd.h>
	#include <sys/wait.h>
	#include <cstdlib>
    #define MANUAL_CORE_DUMP								            \
		do {														    \
			pio::perr << std::endl;                                     \
			pio::perr << "******* Dumping Core ******" << std::endl;    \
			pio::perr << std::endl;                                     \
            int st = 0;        			                                \
			pid_t p = fork();                                           \
            if (!p) {                                                   \
				signal(SIGABRT, SIG_DFL);                               \
				std::abort();                                           \
			} else {                                                    \
				waitpid(p, &st, 0);                                     \
			}                                                           \
		} while (0);
#else
	#define MANUAL_CORE_DUMP
#endif

//
// Manual Traceback
// Note: Deallocating memory using free() was
// colliding with other names of free so we deallocate using
// realloc() with 0 size;
//
#ifdef GEOFLOW_ASSERT_STACK
	#include "tbox/pio.hpp" // pio::
	#include <execinfo.h>   // backtrace
    #define MANUAL_TRACEBACK								  	        \
		do {														    \
			void* callstack[128];                             			\
			int frames = backtrace(callstack, 128);           			\
			char** strs = backtrace_symbols(callstack,frames); 			\
			pio::perr << "********* TraceBack ********" << std::endl;  	\
			for (int i=0; i<frames; ++i) {                  		    \
				pio::perr << strs[i] << std::endl;           			\
			}                                               		    \
			free(strs);                                  			    \
		} while (0);
#else
	#define MANUAL_TRACEBACK
#endif

//
// Static Asserts
//
#ifdef GEOFLOW_USE_STATIC_ASSERTIONS
	#define STATIC_ASSERT(...) static_assert(__VA_ARGS__)
	#define STATIC_ASSERT_MSG(...) static_assert(__VA_ARGS__)
#else
	#define STATIC_ASSERT(...)
	#define STATIC_ASSERT_MSG(...)
#endif

//
// RunTime Asserts
//
#ifdef GEOFLOW_USE_RUNTIME_ASSERTIONS
	#include "tbox/error_handler.hpp" // EH::abort();
	#include "tbox/pio.hpp"   // pio::
	#define ASSERT(EXP)                                                      \
		do {                                                                 \
			using namespace geoflow::tbox;                                   \
			if (!(EXP)) {                                                    \
				pio::perr << std::endl;                                      \
				pio::perr << "***** Failed assertion *****" << std::endl;    \
				pio::perr << "Failed expression: " << # EXP << std::endl;    \
				pio::perr << "File: " << __FILE__ << std::endl;              \
				pio::perr << "Line: " << __LINE__ << std::endl;              \
				pio::perr << std::endl;                                      \
				MANUAL_TRACEBACK									         \
				MANUAL_CORE_DUMP                                             \
				MANUAL_FORCE_HANG                                            \
				EH::abort();                                                 \
			}                                                                \
		} while (0);

	#define ASSERT_MSG(EXP,MSG)                                              \
		do {                                                                 \
			using namespace geoflow::tbox;                                   \
			if (!(EXP)) {                                                    \
				pio::perr << std::endl;                                      \
				pio::perr << "***** Failed assertion *****" << std::endl;    \
				pio::perr << "Failed expression: " << # EXP << std::endl;    \
				pio::perr << "File: " << __FILE__ << std::endl;              \
				pio::perr << "Line: " << __LINE__ << std::endl;              \
				pio::perr << "Message: " << MSG << std::endl;                \
				pio::perr << std::endl;                                      \
				MANUAL_TRACEBACK									         \
				MANUAL_CORE_DUMP                                             \
				MANUAL_FORCE_HANG                                            \
				EH::abort();                                                 \
			}                                                                \
		} while (0);

	#define VERIFY(EXP)         ASSERT(EXP)
	#define VERIFY_MSG(EXP,MSG) ASSERT_MSG(EXP,MSG)
#else
	#define ASSERT(EXP)
	#define ASSERT_MSG(EXP,MSG)
	#define VERIFY(EXP)         do{(void)(EXP);}while(0)
	#define VERIFY_MSG(EXP,MSG) do{(void)(EXP);}while(0)
#endif


#endif /* GEOFLOW_DEBUG_ASSERT_H_ */
