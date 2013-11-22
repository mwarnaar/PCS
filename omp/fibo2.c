#include<stdio.h>
#include<stdlib.h>
#include<sys/time.h>
#include<omp.h>
#include<vector.h>
#include<vector.c>

int fib(int n);

int fib(int n)
{
	int i, j;

	if (n < 2) {
		return n;
	} else {
		i = fib(n - 1);
		j = fib(n - 2);

		return i + j;
	}
}

int main ( int argc, char *argv[] )
{
	if (argc != 2) {
		printf("usage: %s num_threads\n", argv[0]);
		return 0;
	}

	int i;
	int vector[45];
	double time_spent;
	struct timeval before, after;

	for (i = 0; i < 45; ++i) {
		// vector[i] = i + 1;
		vector[i] = rand() % 45 + 1;
	}		

	gettimeofday(&before, NULL);

	omp_set_num_threads(atoi(argv[1]));

	#pragma omp parallel for schedule(guided)
	for (i = 0; i < sizeof(vector) / sizeof(int); i++) {
		int n;
		
		n = vector[i];
		fib(n);
	}

	gettimeofday(&after, NULL);
	time_spent = (double)(after.tv_sec - before.tv_sec) + (double)(after.tv_usec - before.tv_usec) / 1e6;

	printf("FINISHED: timespent = %f, threads = %d\n", time_spent, atoi(argv[1]));

	return 0;
}
