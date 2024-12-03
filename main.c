#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define STB_IMAGE_IMPLEMENTATION
#include "external/stb_image.h"

#include "sensible/common.h"
#include "sensible/str.h"
#include "sensible/slice/f64.h"
#include "sensible/util.h"

#define T f64
#include "sensible/result_T.h"

const f64 DERIVATIVE_H = 1.0 / (1l << 10);
const f64 TRAINING_RATIO = 0.8;

#define CLASS_A str("1")
#define CLASS_B str("0")

static f64
vector_dot_product_unchecked(usize len, const f64 a[len], const f64 b[len]) {
	f64 r = 0;
	for (usize i = 0; i < len; i++) {
		r += a[i] * b[i];
	}
	return r;
}

static Result_f64 vector_dot_product(f64s a, f64s b) {
	if (a.len != b.len) {
		return Result_f64_Err("vectors arent the same length");
	}
	if (a.data == NULL || b.data == NULL) {
		return Result_f64_Err("one or both of the vectors are NULL");
	}
	f64 r = vector_dot_product_unchecked(a.len, a.data, b.data);
	return Result_f64_Ok(r);
}

static bool is_newline(u8 c) {
	return c == '\n';
}

static bool is_tab(u8 c) {
	return c == '\t';
}

strbuf scratch_str_buffer = {0};

typedef struct Example {
	f64_buf data;
	int label;
} Example;

#define T Example
#include "sensible/slice_T.h"
#undef T

static Example load_example(str file, str label) {
	strbuf_clone_from_str(&scratch_str_buffer, file);
	char *filename = strbuf_as_cstr(&scratch_str_buffer);

	int x, y, n;
	unsigned char *data = stbi_load(filename, &x, &y, &n, 1);

	usize dims = (usize)x * (usize)y + 1;
	f64_buf buf = f64_buf_with_capacity(dims);
	buf.data[0] = 1;
	for (buf.len = 1; buf.len < dims; buf.len++) {
		buf.data[buf.len] = ((f64)data[buf.len]) / (f64)255.0;
	}
	stbi_image_free(data);

	Example r = {0};
	r.data = buf;
	if (str_equals(label, CLASS_A)) {
		r.label = 1;
	} else if (str_equals(label, CLASS_B)) {
		r.label = -1;
	} else {
		panicf("Unexpected label encountered: \"" fstr "\"", pstr(label));
	}

	return r;
}

static f64 mean_sq_error(Examples examples, f64s weights) {
	f64 err = 0.0;
	for (usize i = 0; i < examples.len; i++) {
		const Example e = examples.data[i];
		Result_f64 result_y = vector_dot_product(e.data.slice, weights);
		f64 y = tanh(Result_f64_unwrap(&result_y));
		err += pow(y - (f64)e.label, 2.0);
	}
	err /= (f64)examples.len;
	return err;
}

f64_buf gradient(Examples examples, f64s weights) {
	f64_buf grad = f64_buf_with_capacity(weights.len);
	for (usize i = 0; i < weights.len; i++) {
		const f64 w_i = weights.data[i];

		weights.data[i] = w_i + DERIVATIVE_H / 2.0;
		f64 l = mean_sq_error(examples, weights);
		weights.data[i] = w_i - DERIVATIVE_H / 2.0;
		f64 r = mean_sq_error(examples, weights);

		weights.data[i] = w_i;

		f64 pdv = (l - r) / DERIVATIVE_H;
		grad.data[grad.len] = pdv;
		grad.len = i + 1;
	}
	return grad;
}

static f64 accuracy(Examples examples, f64s weights) {
	usize correct = 0;
	for (usize i = 0; i < examples.len; i++) {
		const Example e = examples.data[i];
		Result_f64 result_y = vector_dot_product(e.data.slice, weights);
		f64 y = tanh(Result_f64_unwrap(&result_y));
		if (y >= 0) {
			if (e.label == 1) {
				correct++;
			}
		} else {
			if (e.label == -1) {
				correct++;
			}
		}
	}
	return (f64)correct / (f64)examples.len;
}

typedef struct LogState {
	usize epoch;
	usize update;
	usize time_start;
	Examples training_set;
	Examples test_set;
	f64s weights;
} LogState;

static void log_header(void) {
	printf( //
		"epoch\t"
		"update\t"
		"time\t"
		"training_loss\t"
		"testing_loss\t"
		"training_accuracy\t"
		"testing_accuracy\t"
		"weights"
		"\n"
	);
}

static void log_weights(f64s weights) {
	for (usize i = 0; i < weights.len; i++) {
		printf("%.17g,", weights.data[i]);
	}
}

static void log_state(LogState *s) {
	usize epoch = s->epoch;
	usize update = s->update;
	f64 train_loss = mean_sq_error(s->training_set, s->weights);
	f64 train_acc = accuracy(s->training_set, s->weights);
	f64 test_loss = mean_sq_error(s->test_set, s->weights);
	f64 test_acc = accuracy(s->test_set, s->weights);
	f64 time = (f64)((usize)clock() - s->time_start) / (f64)CLOCKS_PER_SEC;

	printf(
		"%zu\t%zu\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t",
		epoch,
		update,
		time,
		train_loss,
		test_loss,
		train_acc,
		test_acc
	);

	log_weights(s->weights);
	putchar('\n');
}

static void gd_step(Examples examples, f64s *weights, f64 learning_rate) {
	f64_buf grad = gradient(examples, *weights);
	for (usize i = 0; i < weights->len; i++) {
		weights->data[i] -= learning_rate * grad.data[i];
	}
}

static void
gd_epoch(Examples examples, f64s *weights, f64 learning_rate, LogState *lgs) {
	lgs->update = 0;
	lgs->weights = *weights;
	log_state(lgs);
	gd_step(examples, weights, learning_rate);
}

static void sgd_epoch(
	Examples examples,
	f64s *weights,
	f64 learning_rate,
	usize batch_size,
	LogState *lgs
) {
	usize updates = examples.len / batch_size;
	usize remaining = examples.len % batch_size;

	for (usize i = 0; i < updates; i++) {
		lgs->update = i;
		lgs->weights = *weights;
		log_state(lgs);
		const usize start = i * batch_size;
		const usize end = start + batch_size;
		gd_step(Examples_slice(examples, start, end), weights, learning_rate);
	}
	if (remaining > 0) {
		lgs->update += 1;
		lgs->weights = *weights;
		log_state(lgs);
		const usize start = examples.len - remaining;
		const usize end = examples.len;
		gd_step(Examples_slice(examples, start, end), weights, learning_rate);
	}
}

typedef struct AdamParams {
	f64 stepsize;
	f64 beta_1, beta_2;
	f64 beta_1_t, beta_2_t;
	f64s parameters;
	f64_buf m, v;
	f64 eps;
} AdamParams;

static AdamParams
adam_new(f64 stepsize, f64 beta_1, f64 beta_2, f64 eps, f64s initial_params) {
	AdamParams p = {0};
	p.stepsize = stepsize;
	p.beta_1 = beta_1;
	p.beta_2 = beta_2;
	p.beta_1_t = beta_1;
	p.beta_2_t = beta_2;
	p.eps = eps;
	p.parameters = initial_params;
	p.m = f64_buf_with_capacity(p.parameters.len);
	p.v = f64_buf_with_capacity(p.parameters.len);
	for (usize i = 0; i < p.parameters.len; i++) {
		p.m.data[i] = 0.0;
		p.v.data[i] = 0.0;
	}
	p.m.len = p.parameters.len;
	p.v.len = p.parameters.len;
	return p;
}

static void adam_step(AdamParams *p, Examples examples) {
	p->beta_1_t *= p->beta_1;
	p->beta_2_t *= p->beta_2;

	f64_buf grad = gradient(examples, p->parameters);

	f64_buf *m = &p->m;
	for (usize i = 0; i < m->len; i++) {
		m->data[i] *= p->beta_1;
		m->data[i] += (1 - p->beta_1) * grad.data[i];
	}

	f64_buf *v = &p->v;
	for (usize i = 0; i < v->len; i++) {
		v->data[i] *= p->beta_2;
		v->data[i] += (1 - p->beta_2) * (grad.data[i] * grad.data[i]);
	}

	f64s *params = &p->parameters;
	for (usize i = 0; i < params->len; i++) {
		f64 m_hat_i = m->data[i] / (1 - p->beta_1_t);
		f64 v_hat_i = v->data[i] / (1 - p->beta_2_t);
		params->data[i] -= p->stepsize * m_hat_i / (sqrt(v_hat_i) + p->eps);
	}
}

static void
adam_epoch(AdamParams *p, Examples examples, usize batch_size, LogState *lgs) {
	usize updates = examples.len / batch_size;
	usize remaining = examples.len % batch_size;

	for (usize i = 0; i < updates; i++) {
		lgs->update = i;
		lgs->weights = p->parameters;
		log_state(lgs);
		const usize start = i * batch_size;
		const usize end = start + batch_size;
		adam_step(p, Examples_slice(examples, start, end));
	}
	if (remaining > 0) {
		lgs->update += 1;
		lgs->weights = p->parameters;
		log_state(lgs);
		const usize start = examples.len - remaining;
		const usize end = examples.len;
		adam_step(p, Examples_slice(examples, start, end));
	}
}

f64_buf initial_weights(usize dims, f64 seed) {
	f64_buf weights = f64_buf_with_capacity(dims);
	srand((unsigned int)(seed * 1000));
	for (usize i = 0; i < dims; i++) {
		weights.data[i] = (((f64)rand() / (f64)RAND_MAX) * 2 - 1) * seed;
	}
	weights.len = dims;
	return weights;
}

#ifndef ACUTEST_H

int main(int argc, char *argv[]) {
	char *arg_program = argv[0];

	if (argc != 5) {
		panicf(
			"Usage: %s <examples.tsv> <algo: gd/sgd/adam> <epochs>\n",
			arg_program
		);
	}

	char *arg_examples_file = argv[1];

	enum algo {
		GD,
		SGD,
		ADAM,
	};
	enum algo algo;

	{
		str arg_algo = str_from_cstr(argv[2]);
		if (str_equals(arg_algo, str("gd"))) {
			algo = GD;
		} else if (str_equals(arg_algo, str("sgd"))) {
			algo = SGD;
		} else if (str_equals(arg_algo, str("adam"))) {
			algo = ADAM;
		} else {
			panicf(
				"Unexpected algorithm name: \"" fstr
				"\".\nMust be one of: gd, sgd, adam",
				pstr(arg_algo)
			);
		}
	}

	usize epochs;
	{
		char *arg_epochs = argv[3];
		char *end;
		epochs = strtoul(arg_epochs, &end, 10);
		if (end == arg_epochs || *end != '\0') {
			panicf("epochs must be an integer. received \"%s\"", arg_epochs);
		}
		if (epochs < 1) {
			panic("epochs must be at least 1");
		}
	}

	f64 seed;
	{
		char *arg_seed = argv[4];
		char *end;
		seed = strtod(arg_seed, &end);
		if (end == arg_seed || *end != '\0') {
			panicf(
				"seed must be an floating point number. received \"%s\"",
				arg_seed
			);
		}
	}

	strbuf file = {0};
	{
		Result_strbuf r_file = read_file(arg_examples_file);
		file = Result_strbuf_unwrap(&r_file);
	}

	Examples training_set = {0};
	Examples test_set = {0};
	{
		Example_buf examples = Example_buf_with_capacity(0);

		for (
			str line, rest = file.slice;
			(line = str_sep_fn(&rest, is_newline)).data != NULL;
			/* no-stepper */
		) {
			str line_trimmed = str_trim(line);

			str file_name = str_sep_fn(&line_trimmed, is_tab);
			str label = str_sep_fn(&line_trimmed, is_tab);

			Example e = load_example(file_name, label);
			Example_buf_push(&examples, e);
		}

		training_set = Examples_slice(
			examples.slice, 0, (usize)((f64)examples.len * TRAINING_RATIO)
		);
		test_set =
			Examples_slice(examples.slice, training_set.len, examples.len);
	}

	usize dims = training_set.data[0].data.len;
	f64_buf weights = initial_weights(dims, seed);

	AdamParams adam_params = adam_new(0.001, 0.9, 0.999, 1e-8, weights.slice);
	LogState lgs = {
		.epoch = 0,
		.update = 0,
		.time_start = (usize)clock(),
		.training_set = training_set,
		.test_set = test_set,
		.weights = weights.slice
	};

	log_header();

	for (usize i = 0; i < epochs; i++) {
		lgs.epoch = i;

		switch (algo) {
		case GD:
			// Full/Batch GD
			gd_epoch(training_set, &weights.slice, 0.001, &lgs);
			break;
		case SGD:
			// Stochastic GD
			sgd_epoch(training_set, &weights.slice, 0.001, 1, &lgs);
			break;
		case ADAM:
			// ADAM
			adam_epoch(&adam_params, training_set, 1, &lgs);
			break;
		}
	}

	lgs.epoch += 1;
	lgs.update = 0;
	lgs.weights = weights.slice;
	log_state(&lgs);

	return EXIT_SUCCESS;
}

#endif /* ifndef ACUTEST_H */
