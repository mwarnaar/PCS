#include <sys/time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "compute.h"
#include "omp.h"

#ifdef GEN_PICTURES
static void do_draw(const struct parameters *p, size_t key, size_t height, size_t width, double (*restrict g)[height][width])
{
	size_t i, j;

	begin_picture(key, width - 2, height - 2, p->io_tmin, p->io_tmax);

	for (i = 1; i < height - 1; ++i)
		for (j = 1; j < width - 1; ++j)
			draw_point(j - 1, i - 1, (*g)[i][j]);

	end_picture();
}
#endif

static void do_copy(size_t height, size_t width, double (*restrict g)[height][width])
{
	size_t i;

	/* copy left and right column to opposite border */
	for (i = 0; i < height; ++i) {
		(*g)[i][width - 1] = (*g)[i][1];
		(*g)[i][0] = (*g)[i][width - 2];
	}
}

static void fill_report(const struct parameters *p, struct results *r, size_t height, size_t width, double (*a)[height][width], double maxdiff, double iter, struct timeval *before)
{
	/* compute min/max/avg */
	double tmin = INFINITY, tmax = -INFINITY;
	double sum = 0.0;
	struct timeval after;

	for (size_t i = 1; i < height - 1; ++i) {
		for (size_t j = 1; j < width - 1; ++j) {
			double v = (*a)[i][j];
			sum += v;
			if (tmin > v) tmin = v;
			if (tmax < v) tmax = v;
		}
	}

	r->niter = iter;
	r->maxdiff = maxdiff;
	r->tmin = tmin;
	r->tmax = tmax;
	r->tavg = sum / (p->N * p->M);

	gettimeofday(&after, NULL);
	r->time = (double)(after.tv_sec - before->tv_sec) +
		(double)(after.tv_usec - before->tv_usec) / 1e6;
}

void do_compute(const struct parameters* p, struct results *r)
{
	size_t i, j;

	/* alias input parameters */
	const double (*restrict tinit)[p->N][p->M] = (const double (*)[p->N][p->M])p->tinit;
	const double (*restrict cinit)[p->N][p->M] = (const double (*)[p->N][p->M])p->conductivity;

	/* allocate grid data */
	const size_t height = p->N + 2;
	const size_t width = p->M + 2;
	double (*restrict g1)[height][width] = malloc(height * width * sizeof(double));
	double (*restrict g2)[height][width] = malloc(height * width * sizeof(double));

	/* allocate halo for conductivities */
	double (*restrict c)[height][width] = malloc(height * width * sizeof(double));

	struct timeval before;
	gettimeofday(&before, NULL);

	static const double c_cdir = 0.25 * M_SQRT2 / (M_SQRT2 + 1.0);
	static const double c_cdiag = 0.25 / (M_SQRT2 + 1.0);

	/* set the correct amount of threads to be used by omp */
	omp_set_num_threads(p->nthreads);

	/* set initial temperatures and conductivities */
	for (i = 1; i < height - 1; ++i) {
		for (j = 1; j < width - 1; ++j) {
			(*g1)[i][j] = (*tinit)[i - 1][j - 1];
			(*c)[i][j] = (*cinit)[i - 1][j - 1];
		}
	}

	/* smear outermost row to border */
	for (j = 1; j < width - 1; ++j) {
		(*g1)[0][j] = (*g2)[0][j] = (*g1)[1][j];
		(*g1)[height - 1][j] = (*g2)[height - 1][j] = (*g1)[height - 2][j];
	}

	printf("Amount of threads available in single: %d", omp_get_num_threads());
	/* compute */
	size_t iter;
	double maxdiff = 0.0;
	double (*restrict source)[height][width] = g2;
	double (*restrict destination)[height][width] = g1;

	for (iter = 1; iter <= p->maxiter; ++iter)
	{
#ifdef GEN_PICTURES
		do_draw(p, iter, height, width, source);
#endif
		/* swap source and destination */
		{ void *tmp = source; source = destination; destination = tmp; }

		/* initialize halo on source */
		do_copy(height, width, source);

		/* compute */
		maxdiff = 0.0;

#pragma omp parallel for private(i, j)
		for (i = 1; i < height - 1; ++i) {
			for (j = 1; j < width - 1; ++j) {
				double width = (*c)[i][j];
				double restw = 1.0 - width;

				(*destination)[i][j] = width * (*source)[i][j] +

					((*source)[i+1][j  ] + (*source)[i-1][j  ] +
					 (*source)[i  ][j+1] + (*source)[i  ][j-1]) * (restw * c_cdir) +

					((*source)[i-1][j-1] + (*source)[i-1][j+1] +
					 (*source)[i+1][j-1] + (*source)[i+1][j+1]) * (restw * c_cdiag);

				double diff = fabs((*source)[i][j] - (*destination)[i][j]);

#pragma omp critical
				{
					if (diff > maxdiff) {
						maxdiff = diff;
					}
				}
			}
		}

		/* check for convergence */
		if (maxdiff < p->threshold)
		{ ++iter; break; }

		/* conditional reporting */
		if (iter % p->period == 0) {
			fill_report(p, r, height, width, destination, maxdiff, iter, &before);
			report_results(p, r);
		}
	}

	/* report at end in all cases */
	fill_report(p, r, height, width, destination, maxdiff, iter - 1, &before);

	free(c);
	free(g2);
	free(g1);
}

