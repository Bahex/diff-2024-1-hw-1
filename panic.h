#ifndef _PANIC_H_
#define _PANIC_H_

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * @brief panics with the error message, terminating the program
 * @details This function never returns. It writes `str` to `stderr` and exits
 *          with an error code of `1`.
 *
 * @param str The error message displayed on termination
 */
#define panic(STR)          \
	do {                    \
		fputs(STR, stderr); \
		exit(1);            \
	} while (0)

/**
 * @brief panics with the error message, terminating the program
 * @details This function never returns. It writes to `stderr` and exits
 *          with an error code of `1`. This function calls `printf` to
 *          allow message formatting
 *
 * @param str The error message displayed on termination
 */
#define panicf(FMT, ...)                   \
	do {                                   \
		fprintf(stderr, FMT, __VA_ARGS__); \
		exit(1);                           \
	} while (0)

#endif /* _PANIC_H_ */
