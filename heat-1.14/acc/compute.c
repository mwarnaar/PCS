#include <sys/time.h>
#include <math.h>
#include <stdlib.h>
#include <openacc.h>
#include "compute.h"

#ifdef GEN_PICTURES
static void do_draw(const struct parameters *p,
		size_t key, size_t h, size_t w, size_t len,
		double* g)
{
	begin_picture(key, w-2, h-2, p->io_tmin, p->io_tmax);
	size_t i, j;
	for (i = 1; i < h-1; ++i)
		for (j = 1; j < w-1; ++j)
			draw_point(j-1, i-1, (*g)[i*w+j]);
	end_picture();
}
#endif

static void do_copy(size_t h, size_t w, size_t len, double* g)
{
	size_t i;

	/* copy left and right column to opposite border */
	for (i = 0; i < h; ++i) {
		g[i*w + w-1] = g[i*w+1];
		g[i*w + 0] = g[i*w + w-2];
	}
}

static void fill_report(const struct parameters *p, struct results *r,
		size_t h, size_t w, size_t len,
		double* a,
		double maxdiff, 
		double iter,
		struct timeval *before)
{
	/* compute min/max/avg */
	double tmin = 9000000000000, tmax = -9000000000000;
	double sum = 0.0;
	struct timeval after;

	for (size_t i = 1; i < h - 1; ++i)
		for (size_t j = 1; j < w - 1; ++j) 
		{
			double v = a[i*w + j];
			sum += v;
			if (tmin > v) tmin = v;
			if (tmax < v) tmax = v;
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
	//const double (*restrict tinit)[p->N][p->M] = (const double (*)[p->N][p->M])p->tinit;
	//const double (*restrict cinit)[p->N][p->M] = (const double (*)[p->N][p->M])p->conductivity;

	/* allocate grid data */
	const size_t h = p->N + 2;
	const size_t w = p->M + 2;
	double (*restrict g1)[h * w] = malloc(h * w * sizeof(double));
	double (*restrict g2)[h * w] = malloc(h * w * sizeof(double));

	double *restrict tmp = (double*) malloc(h * w * sizeof(double));

	/* allocate halo for conductivities */
	double *restrict c = (double*) malloc(h * w * sizeof(double));

	struct timeval before;
	gettimeofday(&before, NULL);

	static const double c_cdir = 0.25 * M_SQRT2 / (M_SQRT2 + 1.0);
	static const double c_cdiag = 0.25 / (M_SQRT2 + 1.0);

	int ndev = acc_get_num_devices(acc_device_nvidia);
	printf("Num NVIDIA (%d): %d\n", acc_device_nvidia, ndev);

	/* set initial temperatures and conductivities */
	for (i = 1; i < h - 1; ++i)
		for (j = 1; j < w - 1; ++j) 
		{
			//(*g1)[i][j] = (*tinit)[i-1][j-1];
			//(*c)[i][j] = (*cinit)[i-1][j-1];
			(*g1)[i*w + j] = p->tinit[(p->N * (i-1)) + j-1];
			c[i*w + j] = p->conductivity[(p->N * (i-1)) + j-1];
		}

	/* smear outermost row to border */
	for (j = 1; j < w-1; ++j) {
		(*g1)[0*w + j] = (*g2)[0*w + j] = (*g1)[1*w + j];
		(*g1)[(h-1)*w + j] = (*g2)[(h-1)*w + j] = (*g1)[(h-2)*w + j];
	}

	/* compute */
	int len = h*w;
	size_t iter;
	double maxdiff = 0.0;
	double *dst = (double *) g1;
	double *src = (double *) g2;
#pragma acc data copyin(src[0:len], c[0:len], dst[0:len], tmp[0:len], maxdiff, h, w)
	{
		const size_t height = h;
		const size_t width = w;
		for (iter = 1; iter <= p->maxiter; ++iter) {
#ifdef GEN_PICTURES
			do_draw(p, iter, h, w, len, src);
#endif

			/* initialize halo on source */
			do_copy(h, w, len, src);

			/* compute */
#pragma acc region
			{
#pragma acc loop gang worker collapse(2) reduction(max:maxdiff) independent
				for (int i = 1; i < height - 1; ++i) {
					double temp_maxdiff = 0.0;
					for (int j = 1; j < width - 1; ++j) {
						double c_factor = c[i*width + j];
						double restw = 1.0 - c_factor;

						dst[i*width + j] = c_factor * src[i*width + j] + 

							(src[(i+1)*width + j] + src[(i-1)*width + j] + 
							 src[i*width + j + 1] + src[i*width + j - 1]) * (restw * c_cdir) +

							(src[(i-1)*width + j - 1] + src[(i-1)*width + j + 1] + 
							 src[(i+1)*width + j - 1] + src[(i+1)*width + j + 1]) * (restw * c_cdiag);

						double diff = fabs(src[i*width + j] - dst[i*width + j]);
						if (diff > temp_maxdiff)
							temp_maxdiff = diff;
					}
					if (temp_maxdiff > maxdiff)
						maxdiff = temp_maxdiff;
				}

			/* swap source and destination */
#pragma acc loop gang independent	
				for (int i = 0; i < len; ++i) {
					tmp[i] = src[i];
					src[i] = dst[i];
					dst[i] = tmp[i];			
				}
			}
#pragma acc update host(maxdiff)

			/* check for convergence */
			if (maxdiff < p->threshold) { 
/*#pragma acc update host(dst[0:len])*/
				++iter; 
				break; 
			}

			/* conditional reporting */
			if (iter % p->period == 0) {
#pragma acc update host(maxdiff, dst[0:len])
				fill_report(p, r, height, width, len, dst, maxdiff, iter, &before);
				report_results(p, r);
			}
		}

		/* report at end in all cases */
		fill_report(p, r, height, width, len, dst, maxdiff, iter-1, &before);
	/* end of data block */
	}
	free(c);
	free(g1);
	free(g2);
	free(tmp);
}
