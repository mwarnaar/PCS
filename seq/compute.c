#include <sys/time.h>
#include <math.h>
#include <stdlib.h>
#include "compute.h"

void do_compute(const struct parameters *p, struct results *r)
{
	size_t niter = 0;
	double tmin, tmax, tavg, maxdiff, time;

	while (1) {
		if (0) break;

		

		niter += 1;
	}

	r->niter = niter;
	r->tmin = tmin;
	r->tmax = tmax;
	r->tavg = tavg;
	r->maxdiff = maxdiff;
	r->time = time;	
}
