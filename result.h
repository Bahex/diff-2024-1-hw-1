#ifndef _RESULT_H_
#define _RESULT_H_

#include <stdbool.h>
#include <stdint.h>

#include "panic.h"

/**
 * @file result.h
 * Provides an alternative to exceptions, whose semantics are slightly better
 * than exceptions
 *
 * C does not support exceptions. Many C functions instead return a status code
 * reporting what happened. Alternatively, many C libraries allow a callback to
 * be logged in order to report a detailed description of any errors that
 * occurred. However, the programming language known as Rust deliberately left
 * out exceptions; instead they opted for return values to report errors, with a
 * `RESULT` type. This type contains the actual value if there is one, or some
 * error otherwise.
 *
 * This is the approach that this file defines.
 */

/**
 * @brief A macro which gives the value of a `RESULT` if present, else returns
 *        the result.
 * @details Do not use this twice on the same line.
 *
 *          Example Usage:
 *          @code
 *          RESULT(int) someResultFromFn = ...;
 *
 *          int value;
 *
 *          RESULT_M_TRY(int, value, someResultFromFn);
 *          @endcode
 *
 * @param type The type of the `RESULT(type)`
 * @param dest The location to write the new result into
 * @param result The old result
 * @return nothing. (wrapped in `do-while`)
 *
 * @see RESULT_M_TRY_DECL
 */
#define RESULT_M_TRY(type, dest, result)                          \
	do {                                                          \
		RESULT(type) _result_##type##_try_at_##__LINE__ = result; \
		if (_result_##type##_try_at_##__LINE__.err) {             \
			return _result_##type##_try_at_##__LINE__;            \
		}                                                         \
		dest = _result_##type##_try_at_##__LINE__.value;          \
	} while (0)

/**
 * @brief A macro which gives the value of a `RESULT` if present, else returns
 *        the result.
 * @details Do not use this twice on the same line.
 *
 *          Example Usage:
 *          @code
 *          RESULT(int) myFn(int someParam) {
 *              ....
 *          }
 *
 *          RESULT(int) function() {
 *              RESULT_M_TRY_DECL(int, val, myFn(1));
 *              // equivalent to:
 *              RESULT(int) r = myFn(1);
 *              if (!RESULT_IS_OK(int)(r)) return r;
 *              int val = RESULT_UNWRAP(r);
 *          }
 *          @endcode
 *
 * @param type The type of the `RESULT(type)`
 * @param dest The location to write the new result into
 * @param result The old result
 * @return nothing. (wrapped in `do-while`)
 *
 * @see RESULT_M_TRY
 */
#define RESULT_M_TRY_DECL(type, dest, result) \
	type dest;                                \
	RESULT_M_TRY(type, dest, result)

/**
 * @brief Template for a RESULT class (similar to the one that the language Rust
 *        has). This is the typename.
 * @details This macro expands to a unique name per type, so `RESULT(int)` is
 *          roughly equivalent to the C++ `RESULT<int>`. Note that
 *          `RESULT(int*)` will not compile; if you want to use something along
 *          those lines, you must typedef it first.
 *
 *          For the error strings, you must use string literals or something
 * with storage duration, so the pointer doesn't become invalidated.
 *
 * @param  The template parameter
 * @return `Result_ type`, so `RESULT(int)` is `Result_int` Each macro of this
 * sort follows that convention, and you can rely on it
 */
#define RESULT(type) Result_##type

/* Class functions */
/**
 * @brief Creates a result which holds a value of `type`
 * @details Indicates that no error occurred, and that the computed value from
 *          the function can be retrieved via `RESULT_UNWRAP`.
 *
 * @param type value
 * @return Returns a `RESULT(type)` containing a `type`
 */
#define RESULT_OK(type) Result_##type##_Ok

/**
 * @brief Creates a result which indicates that an error occurred.
 * @details Indicates that an error occurred in the function. Similar to an
 *          exception in other languages.
 *
 * @param const char * error message
 * @return Returns a `RESULT(type)` that does not contain a `type`
 */
#define RESULT_ERR(type) Result_##type##_Err

/* methods */
/**
 * @brief Tests if the `RESULT` actually contains a `type` (else it was an
 * error)
 *
 * @return `true` if calling `RESULT_UNWRAP` is safe, else `false`
 */
#define RESULT_IS_OK(type) Result_##type##_is_ok

/**
 * @brief Tests if the `RESULT` doesn't contain a `type` (else it was ok)
 *
 * @return `true` if calling `RESULT_UNWRAP_ERR` is safe, else `false`
 */
#define RESULT_IS_ERR(type) Result_##type##_is_err

/**
 * @brief Unwraps the `RESULT` into a `type`, panicking if impossible
 * @details If `RESULT_IS_OK`, this returns the `type` value of the `RESULT`.
 *          Otherwise, this panics.
 *
 * @return `type` that is contained in the `RESULT`, or could `exit()`
 *
 * @see RESULT_EXPECT
 */
#define RESULT_UNWRAP(type) Result_##type##_unwrap

/**
 * @brief Unwraps the `RESULT` into a `type`, panicking with the error message
 *        if impossible
 * @details Similar to `RESULT_UNWRAP`, except this has a specialized error
 *          message
 *
 * @param const char * the error message on failure
 * @return `type` contained in `RESULT`, else terminates the program
 *
 * @see RESULT_UNWRAP
 */
#define RESULT_EXPECT(type) Result_##type##_expect

/**
 * @brief Unwraps the `RESULT` into a `type`, returning the supplied alternative
 *        otherwise
 * @details If `RESULT_IS_OK`, this returns the `type` value of the `RESULT`.
 *          Otherwise, this returns the alternative supplied in the function
 *          arguments.
 *
 * @param type or_else The value returned if the unwrap failed
 * @return `type` contained in `RESULT`, or the alternative
 *
 * @see RESULT_UNWRAP
 * @see RESULT_EXPECT
 */
#define RESULT_UNWRAP_OR(type) Result_##type##_unwrap_or

/**
 * @brief Unwraps the `RESULT` into an error string, panicking if impossible
 * @details If `RESULT_IS_ERR`, this returns error message of the `RESULT`.
 *          Otherwise, this panics.
 *
 * @param const char * the error string
 * @return error string in the RESULT, or could terminate
 *
 * @see RESULT_EXPECT
 */
#define RESULT_UNWRAP_ERR(type) Result_##type##_unwrap_err

/**
 * @brief Declares a RESULT type for use.
 * @details This declares the struct and functions for a given type. Note that
 *          the functions still need to be defined. The values inside the struct
 *          are implementation details, and may change at any time; do not use
 *          them.
 *
 * @param  template parameter
 *
 * @see RESULT_DEFINE
 */
#define RESULT_DECLARE(type)                                      \
	typedef struct RESULT(type) {                                 \
		type value;                                               \
		/* string literals only */                                \
		const char *err;                                          \
	} RESULT(type);                                               \
	RESULT(type) RESULT_OK(type)(type value);                     \
	RESULT(type) RESULT_ERR(type)(const char *err);               \
	bool RESULT_IS_OK(type)(const RESULT(type) *);                \
	bool RESULT_IS_ERR(type)(const RESULT(type) *);               \
	type RESULT_UNWRAP(type)(const RESULT(type) *);               \
	type RESULT_EXPECT(type)(const RESULT(type) *, const char *); \
	type RESULT_UNWRAP_OR(type)(const RESULT(type) *, type);      \
	const char *RESULT_UNWRAP_ERR(type)(const RESULT(type) *);

/**
 * @brief Defines all the functions for the given `RESULT(type)`
 * @details Defines each and every function necessary to use the `RESULT`.
 *
 * @param  template parameter
 *
 * @see RESULT_DECLARE
 */
#define RESULT_DEFINE(type)                                                \
	RESULT(type) RESULT_OK(type)(type value) {                             \
		RESULT(type) res = {.value = value, .err = NULL};                  \
		return res;                                                        \
	}                                                                      \
	RESULT(type) RESULT_ERR(type)(const char *err) {                       \
		RESULT(type) res = {.err = err};                                   \
		return res;                                                        \
	}                                                                      \
	bool RESULT_IS_OK(type)(const RESULT(type) * res) {                    \
		return res->err == NULL;                                           \
	}                                                                      \
	bool RESULT_IS_ERR(type)(const RESULT(type) * res) {                   \
		return res->err != NULL;                                           \
	}                                                                      \
	type RESULT_UNWRAP(type)(const RESULT(type) * res) {                   \
		if (RESULT_IS_ERR(type)(res))                                      \
			panicf("Attempted to unwrap empty Result of type " #type       \
				   ". Instead had error: %s",                              \
				   res->err);                                              \
		return res->value;                                                 \
	}                                                                      \
	type RESULT_EXPECT(type)(const RESULT(type) * res,                     \
							 const char *message_on_err) {                 \
		if (RESULT_IS_ERR(type)(res))                                      \
			panic(message_on_err);                                         \
		return res->value;                                                 \
	}                                                                      \
	type RESULT_UNWRAP_OR(type)(const RESULT(type) * res, type else_val) { \
		if (RESULT_IS_ERR(type)(res))                                      \
			return else_val;                                               \
		return res->value;                                                 \
	}                                                                      \
	const char *RESULT_UNWRAP_ERR(type)(const RESULT(type) * res) {        \
		if (RESULT_IS_OK(type)(res))                                       \
			panic("Result was not an error; type: " #type);                \
		return res->err;                                                   \
	}

#endif /* _RESULT_H_ */
