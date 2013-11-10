#include <sys/time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "compute.h"

void do_compute(const struct parameters *p, struct results *r)
{
	int i, j, k, row, column;
	double tmin, tmax, tavg, maxdiff, maxdiff_temp, time, total_time;
	struct timeval begin, end, total_begin, total_end;

	double top_constants[p->M];
	double bottom_constants[p->M];

	for (i = 0; i < p->M; i++) {
		top_constants[i] = p->tinit[i];
		bottom_constants[i] = p->tinit[(p->N - 1) * p->M + i];
	}

	double direct_neighbour = (sqrt(2) / (sqrt(2) + 1)) / 4;
	double diagonal_neighbour = (1 / (sqrt(2) + 1)) / 4;

	double neighbours[3][3];

	double temperature_neighbours;

	gettimeofday(&total_begin, NULL);

	double temperatures_old[p->N][p->M];
	double temperatures_new[p->N][p->M];
	double conductivity[p->N][p->M];

	memcpy(temperatures_old, p->tinit, p->N * p->M * sizeof(double));
	memcpy(conductivity, p->conductivity, p->N * p->M * sizeof(double));

	printf("Initial temps:\n");
	for (row = 0; row < p->N; row++) {
		for (column = 0; column < p->M; column++) {
			printf("Temperature [%d][%d] = %f\n", row, column, temperatures_old[row][column]);
		}
	}

	printf("Initial conductivity\n");
	for (row = 0; row < p->N; row++) {
		for (column = 0; column < p->M; column++) {
			printf("Conductivity [%d][%d] = %f\n", row, column, conductivity[row][column]);
		}
	}

	for (i = 0; i < p->maxiter; i++) {
		gettimeofday(&begin, NULL);
		/*printf("Iteration %d\n", i + 1);*/
		tavg = 0;
		tmax = -99999;
		tmin = 99999;
		maxdiff = 0;

		// Compute
		for (row = 0; row < p->N; row++) {
			for (column = 0; column < p->M; column++) {
				temperature_neighbours = 0;

				neighbours[3][3];

				/* TOP ROW */

				if (row == 0 && column == 0) {
					/*
					   o | o | o
					   o | s | i
					   o | i | i
					   */
					   neighbours[0][0] = top_constants[p->M - 1];
					   neighbours[0][1] = top_constants[column];
					   neighbours[0][2] = top_constants[column + 1];
					}

					if (row == 0 && column > 0 && column < (p->M - 1)) {
					/*
					o | o | o
					i | s | i
					i | i | i
					 */
					neighbours[0][0] = top_constants[column - 1];
					neighbours[0][1] = top_constants[column];
					neighbours[0][2] = top_constants[column + 1];
				}

				if (row == 0 && column == (p->M - 1)) {
					/*
					o | o | o
					i | s | o
					i | i | o
					 */
					neighbours[0][0] = top_constants[column - 1];
					neighbours[0][1] = top_constants[column];
					neighbours[0][2] = top_constants[0];
				}

				/* BOTTOM ROW */

				if (row == (p->N - 1) && column == 0) {
					/*
					   o | i | i
					   o | s | i
					   o | o | o
					   */
					   neighbours[2][0] = bottom_constants[p->M - 1];
					   neighbours[2][1] = bottom_constants[column];
					   neighbours[2][2] = bottom_constants[column + 1];
					}

					if (row == (p->N - 1) && column > 0 && column < (p->M - 1)) {
					/*
					i | i | i
					i | s | i
					o | o | o
					 */
					neighbours[2][0] = bottom_constants[column - 1];
					neighbours[2][1] = bottom_constants[column];
					neighbours[2][3] = bottom_constants[column + 1];
				}

				if (row == (p->N - 1) && column == (p->M - 1)) {
					/*
					i | i | o
					i | s | o
					o | o | o
					 */
					neighbours[2][0] = bottom_constants[column - 1];
					neighbours[2][1] = bottom_constants[column];
					neighbours[2][2] = bottom_constants[0];
				}

				/*
 				io | i | io
				io | s | io
				io | i | io
				*/
				for (j = 0; j < 3; j++) {
					for (k = 0; k < 3; k++) {
						if (neighbours[j][k]) {
							/*printf("Neighbours at [%d, %d] (neighbour %d) is not empty\n", row, column, j);*/
							continue;
						}

						if (column == 1) {
							// If we have no direct left neighbours we take the right-most ones (cyclic)
							neighbours[0][0] = temperatures_old[row - 1][p->M - 1];
							neighbours[1][0] = temperatures_old[row][p->M - 1];
							neighbours[2][0] = temperatures_old[row + 1][p->M - 1];
						} else {
							// Otherwise we take the direct left neighbours
							neighbours[0][0] = temperatures_old[row - 1][column - 1];
							neighbours[1][0] = temperatures_old[row][column - 1];
							neighbours[2][0] = temperatures_old[row + 1][column - 1];
						}

						if (column == p->M) {
							// If we have no direct right neighbours we take the left-most ones (cyclic)
							neighbours[0][2] = temperatures_old[row - 1][0];
							neighbours[1][2] = temperatures_old[row][0];
							neighbours[2][2] = temperatures_old[row + 1][0];
						} else {
							// Otherwise we take the direct right neighbours
							neighbours[0][2] = temperatures_old[row - 1][column + 1];
							neighbours[1][2] = temperatures_old[row][column + 1];
							neighbours[2][2] = temperatures_old[row + 1][column + 1];
						}

						neighbours[1][1] = 0; /* We are not self involved in the sum of the neighbours */

						neighbours[0][1] = temperatures_old[row - 1][column];
						neighbours[2][1] = temperatures_old[row + 1][column];
					}
				}

				/* END CALCULATING NEIGHBOURS */

				/* START CALCULATING NEIGHBOURS AVERAGE TEMPERATURE */

				temperature_neighbours += (neighbours[0][0] * diagonal_neighbour);
				temperature_neighbours += (neighbours[0][1] * direct_neighbour);
				temperature_neighbours += (neighbours[0][2] * diagonal_neighbour);

				temperature_neighbours += (neighbours[1][0] * direct_neighbour);
				temperature_neighbours += (neighbours[1][2] * direct_neighbour);

				temperature_neighbours += (neighbours[2][0] * diagonal_neighbour);
				temperature_neighbours += (neighbours[2][1] * direct_neighbour);
				temperature_neighbours += (neighbours[2][2] * diagonal_neighbour);

				// New temperature = (conductivity * old temperature) + ((1 - conductivity) * neighbours temperature)
				temperatures_new[row][column] = (conductivity[row][column] * temperatures_old[row][column]) + ((1 - conductivity[row][column]) * temperature_neighbours);
			}
		}

		/* Calculate results */
		for (row = 0; row < p->M; row++) {
			for (column = 0; column < p->N; column++) {
				tavg += temperatures_new[row][column];
				if (temperatures_new[row][column] > tmax) {
					tmax = temperatures_new[row][column];
				}

				if (temperatures_new[row][column] < tmin) {
					tmin = temperatures_new[row][column];
				}

				maxdiff_temp = temperatures_new[row][column] - temperatures_old[row][column];
				if (maxdiff_temp < 0.0) {
					maxdiff_temp = -maxdiff_temp;
				}

				if (maxdiff_temp > maxdiff) {
					maxdiff = temperatures_new[row][column] - temperatures_old[row][column];
				}
			}
		}

		tavg = tavg / (p->M * p->N);

		memcpy(temperatures_old, temperatures_new, p->N * p->M * sizeof(double));

		gettimeofday(&end, NULL);
		time = (end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0);
		r->niter = i;
		r->tmin = tmin;
		r->tmax = tmax;
		r->tavg = tavg;
		r->maxdiff = maxdiff;
		r->time = time;

		if (i % p->period == 0) {
			report_results(p, r);
		}
	}
	gettimeofday(&total_end, NULL);
	total_time = (total_end.tv_sec - total_begin.tv_sec) + ((total_end.tv_usec - total_begin.tv_usec)/1000000.0);
	printf("Total runtime: %f seconds\n", total_time);
}
