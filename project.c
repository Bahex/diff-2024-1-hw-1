/**
 * Bahadır Cüneyd Birtürk    19011028
 * Numerical Analysis Project
 *
 * Algorithms To Implement
 * √ 1. Bisection yöntemi
 * √ 2. Regula-Falsi yöntemi
 * √ 3. Newton-Raphson yöntemi
 *
 *   4. NxN’lik matrisin tersi
 * √ 5. Gauss Eleminasyon yöntemi
 * √ 6. Gauss Seidal yöntemi
 *
 * √ 7. Sayısal Türev (merkezi, ileri ve geri fark)
 * √ 8. Simpson yöntemi
 * √ 9. Trapez yöntemi
 *   10. Değişken dönüşümsüz Gregory Newton Enterpolasyonu
 *
 * Needed Data Structures
 * √ Matrices and related operations
 * √ Polynomials
 *
 * Code style:
 *   1. Use const unless needed otherwise
 *   2. Use "Type_method" style names for functions
 *   3. There is no borrow checking in C, use different types to indicate owned
 *      and borrowed values. Wrapper structs might be useful:
 *        OwnedType { Type *ptr; }
 *   4. MAYBE tagged unions for results.
 */

#include <stdio.h>
#include <stdlib.h>

#include "project.h"

RESULT_DEFINE(number)
RESULT_DEFINE(Matrix_p)

const number ALMOST_ZERO = 1.0 / (1l << 60);

int number_sign(const number n) {
	return n < 0.0 ? -1 : 1;
}

number number_abs(const number n) {
	if (n < 0.0) {
		return -n;
	}
	return n;
}

Matrix_p Matrix_new(int const row, int const column) {
	size_t const header_size = sizeof(Matrix_t);
	size_t const ptr_size = sizeof(number *) * row;
	size_t const matrix_size = sizeof(number) * row * column;

	Matrix_p const matrix = calloc(1, header_size + ptr_size + matrix_size);

	matrix->row = row;
	matrix->column = column;
	number *const matrix_data = (number *)(matrix->data + row);

	int i;
	for (i = 0; i < row; i++) {
		matrix->data[i] = matrix_data + column * i;
	}

	return matrix;
}

void Matrix_free(Matrix_p *matrix) {
	free(*matrix);
	*matrix = NULL;
}

Matrix_p Matrix_new_square(const int N) {
	return Matrix_new(N, N);
}

Matrix_p Matrix_new_identity(const int N) {
	Matrix_p matrix = Matrix_new_square(N);
	if (matrix == NULL) {
		return NULL;
	}

	int i;
	for (i = 0; i < N; i++) {
		matrix->data[i][i] = 1.0;
	}

	return matrix;
}

Matrix_p Matrix_transposed(Matrix_t const *const orig) {
	Matrix_p transposed = Matrix_new(orig->column, orig->row);

	int i, j;
	for (i = 0; i < orig->row; i++) {
		for (j = 0; j < orig->column; j++) {
			transposed->data[j][i] = orig->data[i][j];
		}
	}

	return transposed;
}

void Matrix_row_swap(Matrix_p const matrix, int const r1, int const r2) {
	if (r1 >= matrix->row || r2 >= matrix->row || r1 == r2) {
		return;
	}
	number temp;
	int i;
	for (i = 0; i < matrix->column; i++) {
		temp = matrix->data[r1][i];
		matrix->data[r1][i] = matrix->data[r2][i];
		matrix->data[r2][i] = temp;
	}
}

void Matrix_row_multiply(Matrix_p const matrix, int const row, number factor) {
	int i;
	for (i = 0; i < matrix->column; i++) {
		matrix->data[row][i] *= factor;
	}
}

void Matrix_row_multiply_add(Matrix_p const matrix, int const from_row, int const to_row, number factor) {
	int i;
	for (i = 0; i < matrix->column; i++) {
		matrix->data[to_row][i] += matrix->data[from_row][i] * factor;
	}
}

void Matrix_multiply_into(Matrix_t const *const A, Matrix_t const *const B, Matrix_p const R) {
	int i, j, k;
	for (i = 0; i < R->row; i++) {
		for (j = 0; j < R->column; j++) {
			for (k = 0; k < A->column; k++) {
				R->data[i][j] += A->data[i][k] * B->data[k][j];
			}
		}
	}
}

/**
 * Multiplies two matrices in the given order.
 * @param A
 * @param B
 * @return `A*B` if `A.column == B.row` else `NULL`
 */
Result_Matrix_p Matrix_multiply(Matrix_t const *const A, Matrix_t const *const B) {
	if (A->column != B->row) {
		return Result_Matrix_p_Err("Dimensions of matrices don't match.\n");
	}
	Matrix_p const rt = Matrix_new(A->row, B->column);
	Matrix_multiply_into(A, B, rt);
	return Result_Matrix_p_Ok(rt);
}

/* TODO:
 * - Should take in the augmented matrix
 * - Should return an array or one dimensional matrix of roots/solutions
 * - Decide if it should just be Gauss or Gauss-Jordan and fix accordingly */
void Matrix_gauss_elimination(Matrix_p const matrix) {
	int lead = 0;
	int r, i;
	for (r = 0; r < matrix->row; r++) {
		if (matrix->column <= lead) {
			return;
		}
		i = r;
		while (matrix->data[i][lead] == 0) {
			i++;
			if (matrix->row == i) {
				i = r;
				lead++;
				if (matrix->column <= lead) {
					return;
				}
			}
		}
		if (i != r) {
			Matrix_row_swap(matrix, i, r);
		}
		Matrix_row_multiply(matrix, r, 1.0 / matrix->data[r][lead]);
		for (i = 0; i < matrix->row; i++) {
			if (i != r) {
				Matrix_row_multiply_add(matrix, r, i, -matrix->data[i][lead]);
			}
		} 
		lead++;
	}
}

/**
 * Solves a square system of linear equations Ax = b
 * @param matrix Augmented matrix `A | b`
 * @return unkowns x
 */
Result_Matrix_p Matrix_gauss_seidel(Matrix_p const matrix, number const epsilon) {
	/* NOTE: Allocate `x` */
	Matrix_p x = Matrix_new(matrix->row, 1);

	/* NOTE: Allocate `diffs` */
	number *diffs = calloc(x->row, sizeof(number));

	int within_epsilon;
	int is_first_loop = 1;

	int i, j;
	do {
		/* Reset at the beginning of the loop */
		within_epsilon = 0;

		for (i = 0; i < x->row; i++) {
			/* Load the value */
			number X_i = x->data[i][0];

			X_i = matrix->data[i][matrix->column - 1];

			for (j = 0; j < i; j++) {
				X_i -= matrix->data[i][j] * x->data[j][0];
			}
			/* Need to skip the `j == i` condtition
			 * Could have been done with an `if` inside the loop instead of
			 * having two separete loops, but that adds one extra comparision to
			 * all iterations */
			for (j = i + 1; j < x->row; j++) {
				X_i -= matrix->data[i][j] * x->data[j][0];
			}

			/* Check if diogonal element is zero */
			if (matrix->data[i][i] == 0.0) {
				/* NOTE: Free `x` */
				Matrix_free(&x);
				/* NOTE: Free `diffs` */
				free(diffs);
				return Result_Matrix_p_Err("Zero on diogonal.\n");
			}
			X_i /= matrix->data[i][i];

			const number diff = number_abs(x->data[i][0] - X_i);

			if (!is_first_loop && diff > diffs[i]) {
				/* NOTE: Free `x` */
				Matrix_free(&x);
				/* NOTE: Free `diffs` */
				free(diffs);
				return Result_Matrix_p_Err("Result diverges.\n");
			}

			if (diff <= epsilon) {
				within_epsilon++;
			}

			/* Store the diff */
			diffs[i] = diff;
			/* Write the value back */
			x->data[i][0] = X_i;
		}
		is_first_loop = 0;
	} while (within_epsilon < x->row);

	/* NOTE: Free `diffs` */
	free(diffs);

	/* NOTE: Return `x` */
	return Result_Matrix_p_Ok(x);
}

/* TODO: After determinant  */
Matrix_p Matrix_inverse(Matrix_t const *const matrix) {
	(void) matrix;
	return NULL;
}

/**
 * Prints an array of `m_elem`s as a matrix. The "matrix" must be continous in memory
 */
void _matrix_print(number const *const matrix, int const row, int const column, char const *const prefix) {
	int i, j;
	for (i = 0; i < row; i++) {
		fputs(prefix ? prefix : "", stdout);
		for (j = 0; j < column; j++) {
			printf(" %5.5f", matrix[i * column + j]);
		}
		putchar('\n');
	}
}

void Matrix_print(Matrix_t const *const matrix, char const *const prefix) {
	_matrix_print(&matrix->data[0][0], matrix->row, matrix->column, prefix);
}

Polynomial_t *Polynomial_new(int const degree) {
	size_t const header_size = sizeof(Polynomial_t);
	size_t const array_size = sizeof(number) * (degree + 1);

	if (degree < 0) {
		return NULL;
	}

	Polynomial_t *const new = calloc(1, header_size + array_size);
	if (new == NULL) {
		return  NULL;
	}

	new->degree = degree;
	return new;
}

void Polynomial_free(Polynomial_t **p) {
	free(*p);
	*p = NULL;
}

number *Polynomial_coef_ptr(Polynomial_t *const p, int degree) {
	if (p == NULL || degree > p->degree || degree < 0) {
		return NULL;
	}

	return &p->coefficents[p->degree - degree];
}

void Polynomial_print(Polynomial_t *const p) {
	int i;
	for (i = 0; i <= p->degree; i++) {
		const int degree = p->degree - i;
		const number coef = p->coefficents[i];

		if (coef != 0) {
			printf("%+.2lf", coef);
			if (degree == 1) {
				printf("*x ");
			} else if (degree == 0) {
				printf(" ");
			} else {
				printf("*x^%d ", degree);
			}
		}
	}
}

number Polynomial_calc(Polynomial_t const *const p, number const x) {
	int i;
	number result = 0;
	for (i = 0; i <= p->degree; i++) {
		result *= x;
		result += p->coefficents[i];
	}

	return result;
}

number Polynomial_derivative_backward(Polynomial_t const *const p, number const x, number const h) {
	return (Polynomial_calc(p, x) - Polynomial_calc(p, x - h)) / h;
}

number Polynomial_derivative_forward(Polynomial_t const *const p, number const x, number const h) {
	return (Polynomial_calc(p, x + h) - Polynomial_calc(p, x)) / h;
}

number Polynomial_derivative_central(Polynomial_t const *const p, number const x, number const h) {
	return (Polynomial_calc(p, x + h / 2.0) - Polynomial_calc(p, x - h / 2.0)) / h;
}

/* Behavior is undefined for any case other than `F(l) * F(r) < 0` */
Result_number Polynomial_bisection(Polynomial_t const *const p, number l, number r, number const epsilon) {
	number mid;

	/* Early returns in case one the initial values is a root */
	number l_y = Polynomial_calc(p, l);
	number r_y = Polynomial_calc(p, r);
	const int l_sign = number_sign(l_y);

	if (number_abs(l_y) <= epsilon) {
		return Result_number_Ok(l);
	}
	if (number_abs(r_y) <= epsilon) {
		return Result_number_Ok(r);
	}

	while (1) {
		mid = (l + r) / 2.0;

		number m_y = Polynomial_calc(p, mid);
		const int m_sign = number_sign(m_y);
		m_y = number_abs(m_y);

		if (m_y <= epsilon) {
			return Result_number_Ok(mid);
		}

		if (number_abs(r - l) <= ALMOST_ZERO) {
			return Result_number_Err("Bracket ends met on non root value.\n");
		}

		if (m_sign == l_sign) {
			l = mid;
		} else {
			r = mid;
		}
	}
}

/* Behavior is undefined for any case other than `F(l) * F(r) < 0` */
Result_number Polynomial_regula_falsi(Polynomial_t const *const p, number l, number r, number const epsilon) {
	number c_y, l_y, r_y;

	l_y = Polynomial_calc(p, l);
	r_y = Polynomial_calc(p, r);

	/* Early returns in case one the initial values is a root */
	if (number_abs(l_y) <= epsilon) {
		return Result_number_Ok(l);
	}
	if (number_abs(r_y) <= epsilon) {
		return Result_number_Ok(r);
	}

	/* Initial value used to compare the current iteration's value with the last
	 * one. A `is_first_loop` variable could also have been used, but this is
	 * cheaper in most cases as there is no additional check in the loop. */
	number last_y = number_abs(l_y) > number_abs(r_y) ? l_y : r_y;

	while (1) {
		l_y = Polynomial_calc(p, l);
		r_y = Polynomial_calc(p, r);
		const int l_sign = number_sign(l_y);
		number c = (l * r_y - r * l_y) / (r_y - l_y);

		/* if `r_y` and `l_y` are equal, that results in a division by zero
		 * which can make `c` `inf` or `nan`. In that case it's preferable to
		 * use bisection. `c != c` checks for `nan`, `c * 0 != c * 0` checks for
		 * `inf`, as `inf * 0` is `nan` */
		if (c != c || c * 0 != c * 0) {
			c = (l + r) / 2.0;
		}

		c_y = Polynomial_calc(p, c);
		const int c_sign = number_sign(c_y);
		c_y = number_abs(c_y);

		if (c_y <= epsilon) {
			return Result_number_Ok(c);
		}

		/* Can't trust comparing floats for equality */
		if (number_abs(r - l) <= ALMOST_ZERO) {
			return Result_number_Err("Bracket ends met on non root value.\n");
		}
		/* Exit if a root hasn't been found and there is no change on the value
		 * from the last iteration */
		if (number_abs(last_y - c_y) <= ALMOST_ZERO) {
			return Result_number_Err("Converged on non root value.\n");
		}

		if (c_sign == l_sign) {
			l = c;
		} else {
			r = c;
		}
		last_y = c_y;
	}
}

Result_number Polynomial_newton_raphson(Polynomial_t const *const p, number x, number const h, number const epsilon) {
	number F_x, F_x_abs, Fd_x, last_diff;
	int is_first_loop = 1;

	last_diff = number_abs(Polynomial_calc(p, x));
	while (1) {
		F_x = Polynomial_calc(p, x);
		F_x_abs = number_abs(F_x);

		if (F_x_abs <= epsilon) {
			return Result_number_Ok(x);
		}
		if (F_x_abs > last_diff) {
			return Result_number_Err("Diverging, getting further away from zero.\n");
		}
		if (is_first_loop) {
			is_first_loop = 0;
		} else if ((last_diff - F_x_abs) <= ALMOST_ZERO) {
			return Result_number_Err("Converged on non root value.\n");
		}

		last_diff = F_x_abs;
		Fd_x = Polynomial_derivative_central(p, x, h);
		if (Fd_x == 0.0) {
			return Result_number_Err("Derivative is zero.\n");
		}

		x = x - (F_x / Fd_x);
	}
}

number Polynomial_simpsons_rule(Polynomial_t const *const p, const number a, const number b, const int n) {
	const number h = (b - a) / n;
	number rt;
	int i;

	rt = 0.0;
	rt += Polynomial_calc(p, a);

	for (i = 1; i < n - 1; i += 2) {
		rt += Polynomial_calc(p, a + (h * (number)i)) * 4.0;
		rt += Polynomial_calc(p, a + (h * (number)(i + 1))) * 2.0;
	}
	rt += Polynomial_calc(p, a + (h * (number)(n-1))) * 4.0;

	rt += Polynomial_calc(p, b);

	rt *= h / 3.0;

	return rt;
}

number Polynomial_trapezoidal(Polynomial_t const *const p, const number a, const number b, const int n) {
	const number h = (b - a) / n;
	number rt;
	int i;

	rt = 0.0;
	rt += Polynomial_calc(p, a);
	rt += Polynomial_calc(p, b);
	rt *= 0.5;

	for (i = 1; i < n; i++) {
		rt += Polynomial_calc(p, a + (h * (number)i));
	}

	rt *= h;

	return rt;
}
