#include<stdio.h>
#include<stdlib.h>
#include<sys/time.h>
#include<omp.h>

int fib(int n);

int fib(int n)
{
	int i, j;

	if (n < 2) {
		return n;
	} else {
#pragma omp task shared(i) firstprivate(n)
		{
			i = fib(n - 1);
		}

#pragma omp task shared(j) firstprivate(n)
		{
			j = fib(n - 2);
		}

#pragma omp taskwait

		return i + j;
	}
}

int main ( int argc, char *argv[] )
{
	if (argc != 3) {
		printf("usage: %s number num_threads", argv[0]);
		return 0;
	}

	int n, f;
	double time_spent;
	struct timeval before, after;

	gettimeofday(&before, NULL);
	n = atoi(argv[1]);

	omp_set_num_threads(atoi(argv[2]));

#pragma omp parallel shared(n, f)
	{
#pragma omp single
		f = fib(n);
	}

	gettimeofday(&after, NULL);
	time_spent = (double)(after.tv_sec - before.tv_sec) +
		(double)(after.tv_usec - before.tv_usec) / 1e6;

	printf("\nfib(%d) = %d, time spent = %f, threads = %d\n", n, f, time_spent, atoi(argv[2]));

	return 0;
}
