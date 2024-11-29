#ifndef _PROJECT_H_
#define _PROJECT_H_

#include "result.h"

/* Preferred type for numeric operations
 * Can be any of these:
 *  - `float`
 *  - `double`
 *  - `long double`
 * `typedef` can changed with minimal modification to rest of the code, only
 * required modification being format strings. */
typedef double number;

RESULT_DECLARE(number)

typedef struct matrix Matrix_t;
typedef struct matrix *Matrix_p;
struct matrix {
	int row;
	int column;
	number *data[];
};

RESULT_DECLARE(Matrix_p)

void Matrix_free(Matrix_p *matrix);
Matrix_p Matrix_new(int row, int column);
Matrix_p Matrix_new_square(int N);
Matrix_p Matrix_new_identity(int N);
Matrix_p Matrix_transposed(Matrix_t const *orig);
void Matrix_row_swap(Matrix_p matrix, int r1, int r2);
void Matrix_row_multiply(Matrix_p matrix, int row, number factor);
void Matrix_row_multiply_add(Matrix_p matrix, int from_row, int to_row, number factor);
Result_Matrix_p Matrix_multiply(Matrix_t const *A, Matrix_t const *B);
void Matrix_gauss_elimination(Matrix_p matrix);
Result_Matrix_p Matrix_gauss_seidel(Matrix_p matrix, number epsilon);
number Matrix_determinant(Matrix_t const *matrix);
Matrix_p Matrix_inverse(Matrix_t const *matrix);
void _matrix_print(number const *matrix, int row, int column, char const *prefix);
void Matrix_print(Matrix_t const *matrix, char const *prefix);

typedef struct polynomial Polynomial_t;
typedef struct polynomial *Polynomial_p;
struct polynomial {
	int degree;
	number coefficents[];
};

Polynomial_t *Polynomial_new(int degree);
void Polynomial_free(Polynomial_t **p);
number *Polynomial_coef_ptr(Polynomial_t *p, int degree);
void Polynomial_print(Polynomial_t *p);
number Polynomial_calc(Polynomial_t const *p, number x);
number Polynomial_derivative_central(Polynomial_t const *p, number x, number h);
number Polynomial_derivative_backward(Polynomial_t const *p, number x, number h);
number Polynomial_derivative_forward(Polynomial_t const *p, number x, number h);
Result_number Polynomial_bisection(Polynomial_t const * p, number l, number r, number  epsilon);
Result_number Polynomial_regula_falsi(Polynomial_t const * p, number l, number r, number  epsilon);
Result_number Polynomial_newton_raphson(Polynomial_t const * p, number x, number  h, number  epsilon);
number Polynomial_simpsons_rule(Polynomial_t const *p, number a, number b, int n);
number Polynomial_trapezoidal(Polynomial_t const *p, number a, number b, int n);

#endif /* !_PROJECT_H_ */
