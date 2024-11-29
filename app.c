#include <stdio.h>

#include "project.h"
#include "result.h"
#include "panic.h"

#define READ_BUFFER_SIZE 20

RESULT_DECLARE(int)
RESULT_DECLARE(Polynomial_p)
RESULT_DEFINE(int)
RESULT_DEFINE(Polynomial_p)

typedef enum {
	M_bisection,
	M_regula_falsi,
	M_newton_rhapson,
	M_matrix_inverse,
	M_gauss_elimination,
	M_gauss_seidel,
	M_derivative,
	M_integration_simpson,
	M_integration_trapezoidal,
	M_newton_interpolation,
	M_exit,
	M_invalid,
} MenuItem_t;

#define read_T(T) read_##T
#define read_T_define(T, RF)                                      \
	RESULT(T) read_T(T)() {                                       \
		char buffer[READ_BUFFER_SIZE];                            \
		fputs("(" #T ") > ", stdout);                             \
		if (fgets(buffer, READ_BUFFER_SIZE, stdin) == NULL) {     \
			panic("Can't read from standard input, exiting.");    \
		}                                                         \
		T read;                                                   \
		if (sscanf(buffer, RF, &read) == 1) {                     \
			return RESULT_OK(T)(read);                            \
		}                                                         \
		return RESULT_ERR(T)("Couldn't parse " #T " from input"); \
	}

/* read_loop_T does not return FAILURE ever, either SUCCESS or ABORT */
#define read_loop_T(T) read_loop_##T
#define read_loop_T_define(T, RF)                  \
	T read_loop_T(T)() {                           \
		RESULT(T) result;                          \
		do {                                       \
			result = read_T(T)();                  \
			if (RESULT_IS_ERR(T)(&result)) {       \
				puts("Invalid input, try again."); \
			}                                      \
		} while (RESULT_IS_ERR(T)(&result));       \
		return RESULT_UNWRAP(T)(&result);          \
	}

read_T_define(int, "%d")
read_T_define(number, "%lf")

read_loop_T_define(int, "%d")
read_loop_T_define(number, "%lf")

Polynomial_p readPolynomial();
Matrix_p readMatrix();
Matrix_p readMatrixSquare();

void menu();

int main() {
	menu();

	return 0;
}

void printBar() {
	puts("------------------------------");
}

void menu() {
	MenuItem_t choice = M_invalid;
	while (choice != M_exit) {
		puts("------------------------------\n"
			 " [ 1] Bisection\n"
			 " [ 2] Regula Falsi\n"
			 " [ 3] Newton Rhapson\n"
			 " [ 4] Matrix Inverse (Not Implemented)\n"
			 " [ 5] Gauss Elimination\n"
			 " [ 6] Gauss Seidel\n"
			 " [ 7] Derivative\n"
			 " [ 8] Integration Simpson\n"
			 " [ 9] Integration Trapezoidal\n"
			 " [10] Newton Interpolation (Not Implemented)\n\n"
			 " [11] Exit\n"
			 "------------------------------\n");

		Result_int result = read_T(int)();
		if (Result_int_is_err(&result)) {
			choice = M_invalid;
		} else {
			const int choiceID = Result_int_unwrap(&result);
			if (choiceID >= 1 && choiceID <= M_invalid) {
				choice = (MenuItem_t)(choiceID - 1);
			} else {
				choice = M_invalid;
			}
		}

		switch (choice) {
			case M_bisection: {
				puts("Bisection");
				Polynomial_t *p = readPolynomial();
				puts("Left bracket:");
				const number l = read_loop_T(number)();
				puts("Right bracket:");
				const number r = read_loop_T(number)();
				puts("epsilon:");
				const number epsilon = read_loop_T(number)();
				printBar();
				fputs("Polynomial: ", stdout);
				Polynomial_print(p);
				const Result_number result = Polynomial_bisection(p, l, r, epsilon);
				Polynomial_free(&p);
				if (Result_number_is_ok(&result)) {
					printf("\nRoot: %f\n", Result_number_unwrap(&result));
				} else {
					printf("\nError: %s", Result_number_unwrap_err(&result));
				}
			} break;
			case M_regula_falsi: {
				puts("Regula Falsi");
				Polynomial_t *p = readPolynomial();
				puts("Left bracket:");
				const number l = read_loop_T(number)();
				puts("Right bracket:");
				const number r = read_loop_T(number)();
				puts("epsilon:");
				const number epsilon = read_loop_T(number)();
				printBar();
				fputs("Polynomial: ", stdout);
				Polynomial_print(p);
				const Result_number result = Polynomial_regula_falsi(p, l, r, epsilon);
				Polynomial_free(&p);
				if (Result_number_is_ok(&result)) {
					printf("\nRoot: %f\n", Result_number_unwrap(&result));
				} else {
					printf("\nError: %s", Result_number_unwrap_err(&result));
				}
			} break;
			case M_newton_rhapson: {
				puts("Newton Rhapson");
				Polynomial_t *p = readPolynomial();
				puts("Initial X value:");
				const number x = read_loop_T(number)();
				puts("Derivative difference (h):");
				const number h = read_loop_T(number)();
				puts("epsilon:");
				const number epsilon = read_loop_T(number)();
				printBar();
				fputs("Polynomial: ", stdout);
				Polynomial_print(p);
				const Result_number result = Polynomial_newton_raphson(p, x, h, epsilon);
				Polynomial_free(&p);
				if (Result_number_is_ok(&result)) {
					printf("\nRoot: %f\n", Result_number_unwrap(&result));
				} else {
					printf("\nError: %s", Result_number_unwrap_err(&result));
				}
			} break;
			case M_matrix_inverse: {
				puts("Error: Not Implemented");
			} break;
			case M_gauss_elimination: {
				puts("Gauss Elimination");
				Matrix_p m = readMatrix();
				printBar();
				puts("Input:");
				Matrix_print(m, " |");
				Matrix_gauss_elimination(m);
				puts("Output:");
				Matrix_print(m, " |");
				Matrix_free(&m);
			} break;
			case M_gauss_seidel: {
				puts("Gauss Seidel");
				Matrix_p m = readMatrix();
				puts("epsilon:");
				const number epsilon = read_loop_T(number)();
				printBar();
				puts("Input:");
				Matrix_print(m, " |");
				Result_Matrix_p result = Matrix_gauss_seidel(m, epsilon);
				Matrix_free(&m);
				puts("Output:");
				if (Result_Matrix_p_is_ok(&result)) {
					Matrix_p x = Result_Matrix_p_unwrap(&result);
					Matrix_print(x, " |");
					Matrix_free(&x);
				} else {
					printf("Error: %s", Result_Matrix_p_unwrap_err(&result));
				}
			} break;
			case M_derivative: {
				puts("Derivative");
				Polynomial_t *p = readPolynomial();
				puts("X value:");
				const number x = read_loop_T(number)();
				puts("Derivative difference (h):");
				const number h = read_loop_T(number)();
				printBar();
				fputs("Polynomial: ", stdout);
				Polynomial_print(p);
				putchar('\n');
				printf("Backwards difference derivative: %f\n", Polynomial_derivative_backward(p, x, h));
				printf("  Central difference derivative: %f\n", Polynomial_derivative_central(p, x, h));
				printf("  Forward difference derivative: %f\n", Polynomial_derivative_forward(p, x, h));
				Polynomial_free(&p);
			} break;
			case M_integration_simpson: {
				puts("Integration Simpson");
				Polynomial_t *p = readPolynomial();
				puts("a value:");
				const number a = read_loop_T(number)();
				puts("b value:");
				const number b = read_loop_T(number)();
				puts("Steps (n):");
				const number n = read_loop_T(number)();
				printBar();
				fputs("Polynomial: ", stdout);
				Polynomial_print(p);
				putchar('\n');
				printf("Result: %f\n", Polynomial_simpsons_rule(p, a, b, n));
				Polynomial_free(&p);
			} break;
			case M_integration_trapezoidal: {
				puts("Integration Trapezoidal");
				Polynomial_t *p = readPolynomial();
				puts("a value:");
				const number a = read_loop_T(number)();
				puts("b value:");
				const number b = read_loop_T(number)();
				puts("Steps (n):");
				const number n = read_loop_T(number)();
				printBar();
				fputs("Polynomial: ", stdout);
				Polynomial_print(p);
				putchar('\n');
				printf("Result: %f\n", Polynomial_trapezoidal(p, a, b, n));
				Polynomial_free(&p);
			} break;
			case M_newton_interpolation: {
				puts("Error: Not Implemented");
			} break;
			case M_exit: {
				puts("Exiting");
				return;
			} break;
			case M_invalid: {
				puts("Invalid input");
			} break;
		}
	}
}

Polynomial_p readPolynomial() {
	int degree;

	puts("Degree of the polynomial:");
	do {
		degree = read_loop_T(int)();
		if (degree < 0) { puts("Degree can't be negative."); }
	} while (degree < 0);

	Polynomial_t *p = Polynomial_new(degree);

	int i;
	for (i = degree; i >= 0; i--) {
		number read;
		printf("Enter the coefficient of %d degree term:\n", i);
		read = read_loop_T(number)();
		*Polynomial_coef_ptr(p, i) = read;
	}

	return p;
}

Matrix_p readMatrixRect(const int square) {
	int N;
	int M;

	puts(square ? "Enter N for the NxN matrix:"
				: "Enter N for the NxM matrix:");
	do {
		N = read_loop_T(int)();
		if (N < 1) { puts("N can't be less than 1."); }
	} while(N < 1);

	if (square) {
		M = N;
	} else {
		puts("Enter M for the NxM matrix:");
		do {
			M = read_loop_T(int)();
			if (M < 1) { puts("M can't be less than 1."); }
		} while(M < 1);
	}

	Matrix_p matrix = Matrix_new(N, M);

	int i, j;
	for (i = 0; i < N; i++) {
		for (j = 0; j < M; j++) {
			number read;
			printf("Enter the value of element %d, %d of matrix:\n", i + 1, j + 1);
			read = read_loop_T(number)();
			matrix->data[i][j] = read;
		}
	}

	return matrix;
}

Matrix_p readMatrix() {
	return readMatrixRect(0);
}

Matrix_p readMatrixSquare() {
	return readMatrixRect(1);
}
