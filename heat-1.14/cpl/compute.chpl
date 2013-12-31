use common;
use Time;

proc do_compute(p : params) 
{
	/* Results structure */
	var r : results;

	/* Timer */
	var t : Timer;

	/* Start timer */
	t.start();

	var height : int = p.N + 2;
	var width : int = p.M + 2;

	var srcMatrix: [0..height-1, 0..width-1] real;
	var dstMatrix: [0..height-1, 0..width-1] real;
	var cndMatrix: [0..height-1, 0..width-1] real;

	var c_cdir  : real = 0.25*sqrt(2)/(sqrt(2)+1.0);
	var c_cdiag : real = 0.25/(sqrt(2)+1.0);

	/* Set initial temperatures and conductivities */ 
	for h in 1..height-1 {
		for w in 1..width-1 {
			srcMatrix(h, w) = p.tinit(h-1, w-1);
			cndMatrix(h, w) = p.tcond(h-1, w-1);
		}
	}

	/* Smear outermost row to border */
	for w in 1..width-1 {
		dstMatrix(0, w) = srcMatrix(1, w);
		srcMatrix(0, w) = dstMatrix(0, w);
		dstMatrix(height-1, w) = srcMatrix(height-2, w);
		srcMatrix(height-1, w) = dstMatrix(height-1, w);
	}

	/* Real default value is 0.0 */
	var maxdiff : real;

	/* Compute */
	for i in 1..p.maxiter {

		/* Copy left and righit column to opposite border */
		for i in 0..height {
			srcMatrix[i, width-1] = srcMatrix[i,1];
			srcMatrix[i, 0] = srcMatrix[i, width-2];
		}

		maxdiff = 0.0;

		for i in 1..height-1 {
			for j in 1..width-1 {
				var width : real = cndMatrix[i, j];
				var restw : real = 1.0 - width;

				dstMatrix[i, j] = width * srcMatrix[i, j] +
					(srcMatrix[i+1,j] + srcMatrix[i-1, j] +
					 srcMatrix[i,j] + srcMatrix[i, j-1]) * (restw * c_cdir) +

					(srcMatrix[i-1, j-1] + srcMatrix[i-1, j+1] +
					 srcMatrix[i+1, j-1] + srcMatrix[i+1, j+1]) * (restw*c_cdiag);

				/* Take absolute value of all differences in the two matrices */
				var diff: real = srcMatrix[i, j] - dstMatrix[i, j];
				if (diff < 0) {
					diff = -diff;
				}
				if (diff > maxdiff) {
					maxdiff = diff;
				}
			}
		}

		/* Check for convergence */
		if (maxdiff < p.threshold) {
			i = i+1;
			break;
		}

		/* Conditional reporting */
		if (i % p.period == 0) {
			/* Fill report */
			var tmin: real =  9e10;
			var tmax: real = -9e10;
			var sum:  real;

			for i in 1..height {
				for j in 1..width {
					var v: real = dstMatrix[i, j];
					sum += v;
					if (tmin > v) {
						tmin = v;
					}
					if (tmax < v) {
						tmax = v;
					}
				}	
			}

			r.niter = i;
			r.maxdiff = maxdiff;
			r.tmin = tmin;
			r.tmax = tmax;
			r.tavg = sum / (width * height);
			r.time = t.elapsed();
			report_results(p, r);
		}
		
		/* Copy dst to src */
		for i in 1..height {
			for j in 1..width {
				srcMatrix = dstMatrix[i, j];
			}
		}
	}
	/* Stop timer, set final timing */
	t.stop();
	r.time = t.elapsed();
	return r;
}
