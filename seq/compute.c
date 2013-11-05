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


					if (row == 1 && column == 1) {
						/*
						   o | o | o	
						   o | s | i
						   o | i | i
						   */
						neighbours[0] = top_constants[p->M - 1];
						neighbours[1] = top_constants[column - 1];
						neighbours[2] = top_constants[1];
						neighbours[3] = temperatures_old[(row - 1) * p->M];
						neighbours[7] = temperatures_old[row * p->M];
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
						neighbours[5] = temperatures_old[(row - 1) * p->M];
						neighbours[8] = temperatures_old[row * p->M];
					}

					if ( 
				} 	

				

				/*
 				i | i | i
				i | s | i
				i | i | i
				*/
				for (j = 0; j < 9; j++) {
					if (neighbours[j]) continue;

					neighbours[0] = temperatures_old[(row - 2) * p->M  + column - 1];		
					neighbours[1] = temperatures_old[(row - 2) * p->M  + column];		
					neighbours[2] = temperatures_old[(row - 2) * p->M  + column + 1];		
					neighbours[3] = temperatures_old[(row - 1) * p->M  + column - 1];		
					neighbours[4] = 0; /* We are not self involved in the sum of the neighbours */
					neighbours[5] = temperatures_old[(row - 2) * p->M  + column + 1];		
					neighbours[6] = temperatures_old[row * p->M  + column - 1];		
					neighbours[7] = temperatures_old[row * p->M  + column];		
					neighbours[8] = temperatures_old[row * p->M  + column + 1];		
				}
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
