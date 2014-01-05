use common;
use Time;
use BlockDist;

proc do_compute(p : params) 
{
	/* Results structure */
	var r : results;

	/* Timer */
	var t : Timer;

	/* Start timer */
	t.start();

	/* Allocate grid data */
	var height : int = p.N + 2;
	var width  : int = p.M + 2;
	var srcMatrix: [0..height-1, 0..width-1] real;
	var dstMatrix: [0..height-1, 0..width-1] real;

	/* Allocate halo for conductivities */
	var cndMatrix: [0..height-1, 0..width-1] real;

	var c_cdir  : real = 0.25*sqrt(2)/(sqrt(2)+1.0);
	var c_cdiag : real = 0.25/(sqrt(2)+1.0);

	/* Set initial temperatures and conductivities */ 
	for i in 1..height-2 {
		for j in 1..width-2 {
			dstMatrix[i, j] = p.tinit[i, j];
			cndMatrix[i, j] = p.tcond[i, j];
		}
	}

	/* Smear outermost row to border */
	for j in 1..width-1 {
		dstMatrix[0, j] = dstMatrix[1, j];
		srcMatrix[0, j] = dstMatrix[1, j];
		dstMatrix[height-1, j] = dstMatrix[height-2, j];
		srcMatrix[height-1, j] = dstMatrix[height-2, j];
	}
	var maxdiff : real = INFINITY;
	/* Compute */
	const CylinderWithBorders: domain(2) dmapped Block(boundingBox={0..p.N+1, 0..p.M+1}) = {0..p.N+1, 0..p.M+1};
	var iteration : int = 0;
	while((iteration < p.maxiter) && (maxdiff > p.threshold)) {

		/* Copy dst to src */
		for i in 1..height-2 {
			for j in 1..width-2 {
				srcMatrix[i, j] = dstMatrix[i, j];
			}
		}

		/* Copy left and right column to opposite border */
		for i in 0..height-1 {
			srcMatrix[i, width-1] = srcMatrix[i, 1];
			srcMatrix[i, 0] = srcMatrix[i, width-2];
		}

		maxdiff = 0.0;

		forall (i, j) in CylinderWithBorders do {
			var width : real = cndMatrix[i, j];
			var restw : real = 1.0 - width;

			dstMatrix[i, j] = width * srcMatrix[i, j] +

				(srcMatrix[i+1, j] + srcMatrix[i-1, j] +
				 srcMatrix[i, j+1] + srcMatrix[i, j-1]) * (restw * c_cdir) +
				(srcMatrix[i-1, j-1] + srcMatrix[i-1, j+1] +
				 srcMatrix[i+1, j-1] + srcMatrix[i+1, j+1]) * (restw * c_cdiag);

			/* Take absolute value of all differences in the two matrices */
			var diff : real = srcMatrix[i, j] - dstMatrix[i, j];
			if (diff < 0) {
				diff = -diff;
			}
			if (diff > maxdiff) {
				maxdiff = diff;
			}

		}
		/* Conditional reporting */
		if (iteration % p.period == 0) {
			r = fill_report(r, dstMatrix, height, width, t, iteration, maxdiff);
			report_results(p, r);
		}
		iteration = iteration + 1;
	}
	/* Stop timer, set final timing */
	r = fill_report(r, dstMatrix, height, width, t, iteration, maxdiff);
	t.stop();
	return r;

}

proc fill_report(report, matrix, height, width, timer, itr, maxdiff)
{
	var tmin: real =  INFINITY;
	var tmax: real = -INFINITY;
	var sum:  real = 0.0;

	for i in 1..height-2 {
		for j in 1..width-2 {
			var v: real = matrix[i, j];
			sum += v;
			if (v < tmin) {
				tmin = v;
			}
			if (v > tmax) {
				tmax = v;
			}
		}	
	}

	report.niter = itr;
	report.maxdiff = maxdiff;
	report.tmin = tmin;
	report.tmax = tmax;
	report.tavg = sum / (width * height);
	report.time = timer.elapsed();

	return report;
}
