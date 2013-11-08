#include <sys/time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "compute.h"

void do_compute(const struct parameters *p, struct results *r)
{
	int i, j, row, column;
	double tmin, tmax, tavg, maxdiff, time;

	double top_constants[p->M];
	double bottom_constants[p->M];

	for (i = 0; i < p->M; i++) {
		top_constants[i] = p->tinit[i];
		bottom_constants[i] = p->tinit[(p->N - 1) * p->M + i];
	}

	double direct_neighbour = (sqrt(2) / (sqrt(2) + 1));
	double diagonal_neighbour = (1 / (sqrt(2) + 1));

	double *temperatures_old, *temperatures_new, *conductivity;
	double *neighbours;

	double temperature_neighbours;

	temperatures_old = malloc(p->N * p->M * sizeof(double));
	temperatures_new = malloc(p->N * p->M * sizeof(double));
	conductivity = malloc(p->N * p->M * sizeof(double));
	neighbours = malloc(3 * 3 * sizeof(double));

	memcpy(temperatures_old, p->tinit, p->N * p->M * sizeof(double));
	memcpy(conductivity, p->conductivity, p->N * p->M * sizeof(double));

	for (i = 0; i < p->maxiter; i++) {
		printf("Iteration %d\n", i + 1);

		// Compute
		for (row = 1; row <= p->N; row++) {
			for (column = 1; column <= p->M; column++) {
				/* printf("Temperature [%d][%d] = %f\n", row, column, p->tinit[row * p->M + column]); */
				/* printf("Conductivity [%d][%d] = %f\n", row, column, p->conductivity[row * p->M + column]); */

				temperature_neighbours = 0;

				/* TOP ROW */

				if (row == 1 && column == 1) {
					/*
					   o | o | o
					   o | s | i
					   o | i | i
					   */
					neighbours[0] = top_constants[p->M - 1];
					neighbours[1] = top_constants[column - 1];
					neighbours[2] = top_constants[column];
				}

				if (row == 1 && column > 1 && column < p->M) {
					/*
					o | o | o
					i | s | i
					i | i | i
					 */
					neighbours[0] = top_constants[column - 2];
					neighbours[1] = top_constants[column - 1];
					neighbours[2] = top_constants[column];
				}

				if (row == 1 && column == p->M) {
					/*
					o | o | o
					i | s | o
					i | i | o
					 */
					neighbours[0] = top_constants[column - 2];
					neighbours[1] = top_constants[column - 1];
					neighbours[2] = top_constants[0];
				}

				/* BOTTOM ROW */

				if (row == p->N && column == 1) {
					/*
					   o | i | i
					   o | s | i
					   o | o | o
					   */
					neighbours[6] = bottom_constants[p->M - 1];
					neighbours[7] = bottom_constants[column - 1];
					neighbours[8] = bottom_constants[column];
				}

				if (row == p->N && column > 1 && column < p->M) {
					/*
					i | i | i
					i | s | i
					o | o | o
					 */
					neighbours[6] = bottom_constants[column - 2];
					neighbours[7] = bottom_constants[column - 1];
					neighbours[8] = bottom_constants[column];
				}

				if (row == p->N && column == p->M) {
					/*
					i | i | o
					i | s | o
					o | o | o
					 */
					neighbours[6] = bottom_constants[column - 2];
					neighbours[7] = bottom_constants[column - 1];
					neighbours[8] = bottom_constants[0];
				}

				/*
 				io | i | io
				io | s | io
				io | i | io
				*/
				for (j = 0; j < 9; j++) {
					if (neighbours[j]) continue;

					if (column == 1) {
						// If we have no direct left neighbours we take the right-most ones (cyclic)
						neighbours[0] = temperatures_old[(row - 2) * p->M - 1];
						neighbours[3] = temperatures_old[(row - 1) * p->M - 1];
						neighbours[6] = temperatures_old[row * p->M - 1];
					} else {
						// Otherwise we take the direct left neighbours
						neighbours[0] = temperatures_old[(row - 2) * p->M  + column - 2];
						neighbours[3] = temperatures_old[(row - 1) * p->M  + column - 2];
						neighbours[6] = temperatures_old[row * p->M  + column - 2];
					}

					if (column == p->M) {
						// If we have no direct right neighbours we take the left-most ones (cyclic)
						neighbours[2] = temperatures_old[(row - 2) * p->M];
						neighbours[5] = temperatures_old[(row - 1) * p->M];
						neighbours[8] = temperatures_old[row * p->M];
					} else {
						// Otherwise we take the direct right neighbours
						neighbours[2] = temperatures_old[(row - 2) * p->M  + column];
						neighbours[5] = temperatures_old[(row - 1) * p->M  + column];
						neighbours[8] = temperatures_old[row * p->M  + column];
					}

					neighbours[4] = 0; /* We are not self involved in the sum of the neighbours */

					neighbours[1] = temperatures_old[(row - 2) * p->M  + column - 1];
					neighbours[7] = temperatures_old[row * p->M  + column - 1];
				}

				/* END CALCULATING NEIGHBOURS */

				/* START CALCULATING NEIGHBOURS AVERAGE TEMPERATURE */

				temperature_neighbours += (neigbhours[0] * diagonal_neighbour);
				temperature_neighbours += (neigbhours[1] * direct_neighbour);
				temperature_neighbours += (neigbhours[2] * diagonal_neighbour);
				temperature_neighbours += (neigbhours[3] * direct_neighbour);

				temperature_neighbours += (neigbhours[5] * direct_neighbour);
				temperature_neighbours += (neigbhours[6] * diagonal_neighbour);
				temperature_neighbours += (neigbhours[7] * direct_neighbour);
				temperature_neighbours += (neigbhours[8] * diagonal_neighbour);

				// New temperature = (conductivity * old temperature) + ((1 - conductivity) * neighbours temperature)
				temperatures_new[(row - 1) * p->M + column - 1] = (conductivity[(row - 1) * p->M + column - 1] * temperatures_old[(row - 1) * p->M + column - 1]) + ((1 - conductivity[(row - 1) * p->M + column - 1]) * temperature_neighbours);
			}
		}

		r->niter = i;
		r->tmin = tmin;
		r->tmax = tmax;
		r->tavg = tavg;
		r->maxdiff = maxdiff;
		r->time = time;

		if (i % p->period == 0) {
			printf("Report %d\n", i);
			report_results(p, r);
		}
	}
}
