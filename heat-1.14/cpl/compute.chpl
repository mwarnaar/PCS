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
	for i in 1..height-1 {
		for j in 1..width-1 {
			dstMatrix[i, j] = p.tinit[i-1, j-1];
			cndMatrix[i, j] = p.tcond[i-1, j-1];
		}
	}

	/* Smear outermost row to border */
	for j in 1..width-1 {
		dstMatrix[0, j] = dstMatrix[1, j];
		srcMatrix[0, j] = dstMatrix[1, j];
		dstMatrix[height-1, j] = dstMatrix[height-2, j];
		srcMatrix[height-1, j] = dstMatrix[height-2, j];
	}

	writeln("dst[0, 1] " + srcMatrix[0, 1]);
	var maxdiff : real = 0.0;
	/* Compute */
	
	for itr in 1..p.maxiter {

		/* Copy dst to src */
		for i in 1..height-1 {
			for j in 1..width-1 {
				srcMatrix[i, j] = dstMatrix[i, j];
			}
		}

		/* Copy left and right column to opposite border */
		for i in 0..height {
			srcMatrix[i, width-1] = srcMatrix[i, 1];
			srcMatrix[i, 0] = srcMatrix[i, width-2];
		}

		maxdiff = 0.0;

		for i in 1..height-1 {
			for j in 1..width-1 {
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
		}
		
		/* Check for convergence */
		if (maxdiff < p.threshold) {
			itr = itr+1;
			break;
		}

		/* Conditional reporting */
		if (itr % p.period == 0) {
			var tmin: real =  INFINITY;
			var tmax: real = -INFINITY;
			var sum:  real = 0.0;

			for i in 1..height-1 {
				for j in 1..width-1 {
					var v: real = dstMatrix[i, j];
					sum += v;
					if (v < tmin) {
						tmin = v;
					}
					if (v > tmax) {
						tmax = v;
					}
				}	
			}
			
			r.niter = itr;
			r.maxdiff = maxdiff;
			r.tmin = tmin;
			r.tmax = tmax;
			r.tavg = sum / (width * height);
			r.time = t.elapsed();
			report_results(p, r);
		}
		
	}
	/* Stop timer, set final timing */
	t.stop();
	r.time = t.elapsed();
	return r;
}
