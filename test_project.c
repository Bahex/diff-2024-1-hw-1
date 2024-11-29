#include <stdio.h>
#include <string.h>

#include "acutest.h"
#include "project.h"

#define TEST_FAILED (acutest_cond_failed_)
#define TEST_DUMP_MATRIX(P, M, R, C)               \
	do {                                           \
		if (!TEST_FAILED)                          \
			break;                                 \
		puts("    " P);                            \
		_matrix_print(&(M)[0][0], R, C, "     |"); \
	} while (0)


#define is_array_equal(A1, A2, LEN) \
	(memcmp((void *)(A1), (void *)(A2), sizeof(*(A1)) * (LEN)) == 0)

#define is_matrix_equal(M1, M2, R, C) \
	is_array_equal(&(M1)[0][0], &(M2)[0][0], ((R) * (C)))

#define array_copy(DEST, SRC, LEN) \
	(memcpy((void *)(DEST), (void *)(SRC), sizeof(*(DEST)) * (LEN)))

#define ARRAY_SIZE(A) (sizeof(A) / sizeof((A)[0]))

#define POLYNOMIAL_NEW(NAME, DEGREE, ARRAY)                       \
	do {                                                          \
		Polynomial_free(&NAME);                                   \
		NAME = Polynomial_new(DEGREE);                            \
		array_copy(&NAME->coefficents[0], (ARRAY), (DEGREE) + 1); \
	} while (0)

void Matrix_load_from(Matrix_p const matrix, number const *const array) {
	const size_t size = sizeof(number) * matrix->row * matrix->column;
	memcpy(&matrix->data[0][0], array, size);
}

Matrix_p Matrix_new_from(int const row, int const column, number const *const array) {
	Matrix_p const matrix = Matrix_new(row, column);
	Matrix_load_from(matrix, array);
	return matrix;
}

number f_abs(number const f) {
	return f >= 0.0 ? f : -f;
}

int f_eq(number const f_1, number const f_2, number const epsilon) {
	return f_abs(f_1 - f_2) <= epsilon;
}

int _is_f_array_approximate(number const *const fa_1, number const *const fa_2, size_t const len, number const epsilon) {
	size_t i;
	for (i = 0; i < len; i++) {
		if (!f_eq(fa_1[i], fa_2[i], epsilon)) {
			return 0;
		}
	}
	return 1;
}

void test_tool_matrix_equal() {
	const number M1[2][2] = {
		{1.0, 2.0},
		{-1.0, 3.0},
	};
	const number M2[2][2] = {
		{1.0, 2.0},
		{-1.0, 3.0},
	};
	const number M3[2][2] = {
		{1.0, 2.0},
		{1.0, 3.0},
	};
	TEST_CASE("Equal matrix test");
	TEST_CHECK(is_matrix_equal(M1, M2, 2, 2));

	TEST_CASE("Non-Equal matrix test");
	TEST_CHECK(!is_matrix_equal(M2, M3, 2, 2));
}

void test_tool_matrix_load_new_from() {
	const number zero[2][2] = {0};
	const number src[2][2] = {
		{1.0, 2.0},
		{-1.0, 3.0},
	};

	Matrix_p m = Matrix_new(2, 2);

	TEST_CHECK(is_matrix_equal(zero, m->data, 2, 2));

	Matrix_load_from(m, &src[0][0]);
	TEST_CHECK(!is_matrix_equal(zero, m->data, 2, 2));
	TEST_CHECK(is_matrix_equal(src, m->data, 2, 2));

	Matrix_free(&m);

	m = Matrix_new_from(2, 2, &src[0][0]);
	TEST_CHECK(!is_matrix_equal(zero, m->data, 2, 2));
	TEST_CHECK(is_matrix_equal(src, m->data, 2, 2));

	Matrix_free(&m);
}

void test_tool_float() {
	TEST_CASE("Float absolute");
	TEST_CHECK(f_abs(1.0) == 1.0);
	TEST_CHECK(f_abs(-1.0) == 1.0);
	TEST_CHECK(f_abs(0.0) == 0.0);
	TEST_CHECK(f_abs(-2.0001) == 2.0001);

	TEST_CASE("Float approximate");
	/* (1 << 16) == 2^16 */
	number epsilon = 1.0 / (1l << 60);
	TEST_CHECK(f_eq(0.0, epsilon, epsilon * 2));
	TEST_CHECK(!f_eq(0.0, epsilon * 2, epsilon));
}

void test_Matrix_new_free() {
	const number expected[2][3] = {{0.0}};

	TEST_CASE("Matrix_new");
	Matrix_p matrix = Matrix_new(2, 3);
	TEST_ASSERT_(matrix != NULL, "Memory allocation failed.");

	TEST_CHECK(matrix->row == 2 && matrix->column == 3);
	TEST_CHECK(is_matrix_equal(expected, matrix->data, 2, 3));

	TEST_CASE("Matrix_free");
	Matrix_free(&matrix);
	TEST_CHECK(matrix == NULL);
}

void test_Matrix_new_identity() {
	const number expected[3][3] = {
		{1.0, 0.0, 0.0},
		{0.0, 1.0, 0.0},
		{0.0, 0.0, 1.0},
	};
	Matrix_p matrix = Matrix_new_identity(3);
	TEST_CHECK(is_matrix_equal(expected, matrix->data, 3, 3));

	Matrix_free(&matrix);
}

void test_Matrix_transposed() {
	const number orig[2][3] = {
		{0.0, 1.0, 2.0},
		{3.0, 4.0, 5.0},
	};
	const number expected[3][2] = {
		{0.0, 3.0},
		{1.0, 4.0},
		{2.0, 5.0},
	};

	Matrix_p matrix = Matrix_new_from(2, 3, &orig[0][0]);

	Matrix_p tranposed = Matrix_transposed(matrix);
	TEST_CHECK(tranposed->row == 3 && tranposed->column == 2);
	TEST_CHECK(is_matrix_equal(expected, tranposed->data, 3, 2));

	Matrix_free(&matrix);
	Matrix_free(&tranposed);
}

void test_Matrix_elemantary_row_operations() {
	const number orig[3][4] = {
		{0.0, 1.0, 2.0, 0.0},
		{3.0, 4.0, 5.0, 0.0},
		{0.0, 0.0, 0.0, 0.0},
	};
	const number swapped[3][4] = {
		{0.0, 1.0, 2.0, 0.0},
		{0.0, 0.0, 0.0, 0.0},
		{3.0, 4.0, 5.0, 0.0},
	};
	const number multiplied[3][4] = {
		{0.0, 2.5, 5.0, 0.0},
		{0.0, 0.0, 0.0, 0.0},
		{3.0, 4.0, 5.0, 0.0},
	};
	const number multiply_add[3][4] = {
		{0.0, 2.5, 5.0, 0.0},
		{0.0, 5.0, 10.0, 0.0},
		{3.0, 4.0, 5.0, 0.0},
	};

	Matrix_p matrix = Matrix_new_from(3, 4, &orig[0][0]);

	TEST_CASE("Swap rows 1 and 2");
	Matrix_row_swap(matrix, 1, 2);
	TEST_CHECK(is_matrix_equal(swapped, matrix->data, 3, 4));
	TEST_DUMP_MATRIX("Expected:", swapped, 3, 4);
	TEST_DUMP_MATRIX("Produced:", matrix->data, 3, 4);
	if (TEST_FAILED) {
		Matrix_load_from(matrix, &swapped[0][0]);
	}

	TEST_CASE("Multiply row 0 by 2.5");
	Matrix_row_multiply(matrix, 0, 2.5);
	TEST_CHECK(is_matrix_equal(multiplied, matrix->data, 3, 4));
	TEST_DUMP_MATRIX("Expected:", multiplied, 3, 4);
	TEST_DUMP_MATRIX("Produced:", matrix->data, 3, 4);
	if (TEST_FAILED) {
		Matrix_load_from(matrix, &multiplied[0][0]);
	}

	TEST_CASE("Multiply row 0 by 2 and add to row 1");
	Matrix_row_multiply_add(matrix, 0, 1, 2);
	TEST_CHECK(is_matrix_equal(multiply_add, matrix->data, 3, 4));
	TEST_DUMP_MATRIX("Expected:", multiply_add, 3, 4);
	TEST_DUMP_MATRIX("Produced:", matrix->data, 3, 4);

	Matrix_free(&matrix);
}

void test_Matrix_multiply() {
	const number epsilon = 1.0 / (1l << 62);

	const number a_A[2][3] = {
		{1.0, 2.0, 3.0},
		{4.0, 5.0, 6.0},
	};
	const number a_B[3][2] = {
		{10.0, 11.0},
		{20.0, 21.0},
		{30.0, 31.0},
	};
	const number expected[2][2] = {
		{140.0, 146.0},
		{320.0, 335.0},
	};
	Matrix_p A = Matrix_new_from(2, 3, &a_A[0][0]);
	Matrix_p B = Matrix_new_from(3, 2, &a_B[0][0]);

	Result_Matrix_p result = Matrix_multiply(A, B);
	TEST_CHECK(Result_Matrix_p_is_ok(&result));

	Matrix_p product = Result_Matrix_p_unwrap(&result);
	TEST_CHECK(product->row == 2 && product->column == 2);
	TEST_CHECK(is_matrix_equal(expected, product->data, 2, 2));
	TEST_CHECK(_is_f_array_approximate(&expected[0][0], &product->data[0][0], 2 * 2, epsilon));

	Matrix_free(&A);
	Matrix_free(&B);
	Matrix_free(&product);
}

void test_Matrix_gauss_elimination() {
	const number epsilon = 1.0 / (1l << 62);
	const number pair_1[2][2][3] = {
		{
			{1.0, 3.0, -1.0},
			{0.0, 1.0, 7.0},
		},
		{
			{1.0, 0.0, -22.0},
			{0.0, 1.0, 7.0},
		},
	};
	const number pair_2[2][3][4] = {
		{
			{2.0, 1.0, -1.0, 8.0},
			{-3.0, -1.0, 2.0, -11.0},
			{-2.0, 1.0, 2.0, -3.0},
		},
		{
			{1.0, 0.0, 0.0, 2.0},
			{0.0, 1.0, 0.0, 3.0},
			{0.0, 0.0, 1.0, -1.0},
		},
	};

	struct matrix_compare {
		number const *orig, *expected;
		int const row, column;
	} data[2] = {
		{&pair_1[0][0][0], &pair_1[1][0][0], 2, 3},
		{&pair_2[0][0][0], &pair_2[1][0][0], 3, 4},
	};

	int i;
	int length = ARRAY_SIZE(data);
	for (i = 0; i < length; i++) {
		number const *orig = data[i].orig, *expected = data[i].expected;
		int const row = data[i].row, column = data[i].column;

		TEST_CASE_("Gauss-Jordan Elimanation Case %d", i + 1);

		Matrix_p m = Matrix_new(row, column);
		array_copy(&m->data[0][0], orig, row * column);

		Matrix_gauss_elimination(m);

		TEST_CHECK(_is_f_array_approximate(expected, &m->data[0][0], row * column, epsilon));
		TEST_DUMP_MATRIX("Expected:", &expected, row, column);
		TEST_DUMP_MATRIX("Produced:", m->data, row, column);

		Matrix_free(&m);
	}
}

void test_Matrix_gauss_seidel() {
	const number epsilon = 0.0001;

	TEST_CASE("Convergent matrix");
	Matrix_p matrix = Matrix_new_from(2, 3, (number[]){16.0, 3.0, 11.0, 7.0, -11.0, 13.0});
	const number expected[2][1] = {
		{0.8122},
		{-0.6650}
	};

	Result_Matrix_p result = Matrix_gauss_seidel(matrix, epsilon);
	TEST_CHECK(Result_Matrix_p_is_ok(&result));

	Matrix_p solution = Result_Matrix_p_unwrap(&result);
	TEST_CHECK(_is_f_array_approximate(&solution->data[0][0], &expected[0][0], 2, epsilon));
	TEST_DUMP_MATRIX("Expected:", expected, 2, 1);
	TEST_DUMP_MATRIX("Produced:", solution->data, 2, 1);

	Matrix_free(&matrix);
	Matrix_free(&solution);

	TEST_CASE("Divergent matrix");
	matrix = Matrix_new_from(2, 3, (number[]){2.0, 3.0, 11.0, 5.0, 7.0, 13.0});
	result = Matrix_gauss_seidel(matrix, epsilon);
	TEST_CHECK(Result_Matrix_p_is_err(&result));

	Matrix_free(&matrix);
}

void test_Matrix_determinant() {
	TEST_CHECK_(0, "Not implemented!");
}

void test_Matrix_inverse() {
	TEST_CHECK_(0, "Not implemented!");
}

void test_polynomial_new() {
	const double expected[4] = {0.0};

	TEST_CASE("Polynomial new");
	Polynomial_t *p = Polynomial_new(3);
	TEST_ASSERT_(p != NULL, "Memory allocation failed.");

	TEST_CHECK(p->degree == 3);
	TEST_CHECK(is_array_equal(expected, p->coefficents, 4));
	TEST_CHECK(_is_f_array_approximate(expected, p->coefficents, 4, 0.0));

	TEST_CASE("Polynomial_free");
	Polynomial_free(&p);
	TEST_CHECK(p == NULL);
}

void test_polynomial_coef_ptr() {
	Polynomial_t *p = Polynomial_new(3);
	p->coefficents[0] = 2;
	p->coefficents[1] = 0;
	p->coefficents[2] = 3;
	p->coefficents[3] = 5;

	TEST_CHECK(Polynomial_coef_ptr(p, -1) == NULL);
	TEST_CHECK(Polynomial_coef_ptr(p, 4) == NULL);

	TEST_CHECK(*Polynomial_coef_ptr(p, 3) == 2);
	TEST_CHECK(*Polynomial_coef_ptr(p, 2) == 0);
	TEST_CHECK(*Polynomial_coef_ptr(p, 1) == 3);
	TEST_CHECK(*Polynomial_coef_ptr(p, 0) == 5);

	Polynomial_free(&p);
}

void test_polynomial_calc() {
	Polynomial_t *p = NULL;
	POLYNOMIAL_NEW(p, 3, ((number[]){2, -1, 3, 5}));

	struct {
		number input;
		number expected;
	} const cases[] = {
		{1, 9},
		{-1, -1},
		{2, 23},
		{-2, -21},
		{0, 5},
		{10, 1935},
	};

	int i;
	for (i = 0; i < (int)ARRAY_SIZE(cases); i++) {
		const number input = cases[i].input;
		const number expected = cases[i].expected;

		TEST_CASE_("Case %d: f(x) = 2x^3 - x^2 + 3x + 5, x = %.1f", i, input);

		const number produced = Polynomial_calc(p, input);

		TEST_CHECK(produced == expected);
		TEST_MSG("Expected: %f", expected);
		TEST_MSG("Produced: %f", produced);
	}

	Polynomial_free(&p);
}

/* Forward and backward differantial derivatives are less accurate than their
 * central counterpart, to the point where this unit test would fail at bigger
 * `h` or smaller `epsilon values`. */
void test_polynomial_derivative() {
	const number h = 1.0 / (1l << 26);
	const number epsilon = 1.0 / (1l << 62);
	Polynomial_t *p = NULL;

	struct {
		int degree;
		number *arr;
		number x;
		number expected;
	} const cases[] = {
		{0, (number[]){5}, 10, 0.0},
		{1, (number[]){3, 5}, 10, 3.0},
		{1, (number[]){3, 5}, 20, 3.0},
		{2, (number[]){3, 5, 1}, 1, 11.0},
		{2, (number[]){3, 5, 1}, 2, 17.0},
	};

	struct {
		number (*fn)(const struct polynomial *, number, number);
		char *name;
	} const fnname[] = {
		{Polynomial_derivative_backward, "Backward differance derivative"},
		{Polynomial_derivative_central, "Central differance derivative"},
		{Polynomial_derivative_forward, "Forward differance derivative"},
	};

	int i, j;
	for (i = 0; i < (int)ARRAY_SIZE(cases); i++) {
		TEST_CASE_("Case %d: %d degree polynomial derivative at x=%.1f",
				   i,
				   cases[i].degree,
				   cases[i].x);
		POLYNOMIAL_NEW(p, cases[i].degree, cases[i].arr);

		for (j = 0; j < (int)ARRAY_SIZE(fnname); j++) {
			const number produced = fnname[j].fn(p, cases[i].x, h);

			TEST_CHECK_(f_eq(produced, cases[i].expected, epsilon), "%s", fnname[j].name);
			TEST_MSG("Expected: %.20f", cases[i].expected);
			TEST_MSG("Produced: %.20f", produced);
		}
	}

	Polynomial_free(&p);
}

void test_polynomial_bisection() {
	const number epsilon = 1.0 / (1l << 62);
	Polynomial_t *p = NULL;

	struct {
		int degree;
		number *arr;
		number x_1;
		number x_2;
		number expected;
	} const cases[] = {
		{1, (number[]){4, -5}, 0, 10, 5.0 / 4.0},
		{1, (number[]){4, 5}, -10, 0, -5.0 / 4.0},
		{3, (number[]){1, 1, -9, -9}, -10, 10, 3.0},
		/* {2, (number[]){3, 5, 1}, 1, 3, 0.0}, */
		/* {2, (number[]){3, 5, 1}, 2, 5, 0.0}, */
	};

	int i;
	for (i = 0; i < (int)ARRAY_SIZE(cases); i++) {
		TEST_CASE_("Case %d: %d degree polynomial bisection in [%.2f, %.2f]",
				   i,
				   cases[i].degree,
				   cases[i].x_1,
				   cases[i].x_2);
		POLYNOMIAL_NEW(p, cases[i].degree, cases[i].arr);

		const Result_number result = Polynomial_bisection(p, cases[i].x_1, cases[i].x_2, epsilon);
		TEST_CHECK(Result_number_is_ok(&result));

		const number produced = Result_number_unwrap(&result);

		TEST_CHECK(f_eq(produced, cases[i].expected, epsilon));
		TEST_MSG("Expected: %.20f", cases[i].expected);
		TEST_MSG("Produced: %.20f", produced);
	}
	Polynomial_free(&p);
}

/* TODO: Incomplete */
void test_polynomial_regula_falsi() {
	const number epsilon = 1.0 / (1l << 62);
	Polynomial_t *p = NULL;

	struct {
		int degree;
		number *arr;
		number x_1;
		number x_2;
		number expected;
	} const cases[] = {
		{1, (number[]){4, -5}, 0, 10, 5.0 / 4.0},
		{1, (number[]){4, 5}, -10, 5, -5.0 / 4.0},
		{3, (number[]){1, 1, -9, -9}, -10, 10, -1.0},
		/* {2, (number[]){3, 5, 1}, 1, 3, 0.0}, */
		/* {2, (number[]){3, 5, 1}, 2, 5, 0.0}, */
	};

	int i;
	for (i = 0; i < (int)ARRAY_SIZE(cases); i++) {
		TEST_CASE_("Case %d: %d degree polynomial regula-falsi in [%.2f, %.2f]",
				   i,
				   cases[i].degree,
				   cases[i].x_1,
				   cases[i].x_2);
		POLYNOMIAL_NEW(p, cases[i].degree, cases[i].arr);

		const Result_number result = Polynomial_regula_falsi(p, cases[i].x_1, cases[i].x_2, epsilon);
		TEST_CHECK(Result_number_is_ok(&result));

		const number produced = Result_number_unwrap(&result);

		TEST_CHECK(f_eq(produced, cases[i].expected, epsilon));
		TEST_MSG("Expected: %.20f", cases[i].expected);
		TEST_MSG("Produced: %.20f", produced);
	}
	Polynomial_free(&p);
}

void test_polynomial_newton_raphson() {
	const number epsilon = 1.0 / (1l << 16);
	Polynomial_t *p = NULL;

	struct {
		int degree;
		number *arr;
		number init;
		number expected;
	} const cases[] = {
		{1, (number[]){4, -5}, 10, 5.0 / 4.0},
		{1, (number[]){4, 5}, 5, -5.0 / 4.0},
		{2, (number[]){3, 5, 1}, 1, -0.23241},
		{2, (number[]){3, 5, 1}, -2, -1.43425854},
	};

	int i;
	for (i = 0; i < (int)ARRAY_SIZE(cases); i++) {
		TEST_CASE_("Case %d: %d degree newton raphson in %.2f",
				   i,
				   cases[i].degree,
				   cases[i].init);
		POLYNOMIAL_NEW(p, cases[i].degree, cases[i].arr);

		const Result_number result = Polynomial_newton_raphson(p, cases[i].init, epsilon, epsilon);
		TEST_CHECK(Result_number_is_ok(&result));

		const number produced = Result_number_unwrap(&result);

		TEST_CHECK(f_eq(produced, cases[i].expected, epsilon));
		TEST_MSG("Expected: %.20f", cases[i].expected);
		TEST_MSG("Produced: %.20f", produced);
	}
	Polynomial_free(&p);
}

TEST_LIST = {
	{"test_tool_matrix_equal", test_tool_matrix_equal},
	{"test_tool_matrix_load_new_from", test_tool_matrix_load_new_from},
	{"test_tool_float", test_tool_float},
	{"Matrix_new and Matrix_free", test_Matrix_new_free},
	{"Matrix_new_identity", test_Matrix_new_identity},
	{"Matrix_transposed", test_Matrix_transposed},
	{"Matrix elemantary row operations", test_Matrix_elemantary_row_operations},
	{"test_Matrix_multiply", test_Matrix_multiply},
	{"Matrix_gauss_elimination", test_Matrix_gauss_elimination},
	{"Matrix_gauss_seidel", test_Matrix_gauss_seidel},
	// {"Matrix_determinant", test_Matrix_determinant},
	// {"Matrix_inverse", test_Matrix_inverse},
	{"Polynomial new", test_polynomial_new},
	{"Polynomial coefficent getter", test_polynomial_coef_ptr},
	{"Polynomial calc", test_polynomial_calc},
	{"Polynomial derivative", test_polynomial_derivative},
	{"Polynomial bisection", test_polynomial_bisection},
	{"Polynomial regula falsi", test_polynomial_regula_falsi},
	{"Polynomial newton raphson", test_polynomial_newton_raphson},
	{NULL, NULL} /* zeroed record marking the end of the list */
};
