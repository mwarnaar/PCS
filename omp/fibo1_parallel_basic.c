#include<stdio.h>
#include<stdlib.h>
#include<time.h>

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

int main()
{
	int n, f;
	clock_t begin, end;
	double time_spent;

	printf("Enter number:");
	scanf("%d", &n);

	begin = clock();

	#pragma omp parallel shared(n, f)
	{
		#pragma omp single
		f = fib(n);
	}

	end = clock();
	time_spent = (double) (end - begin) / CLOCKS_PER_SEC;

	printf("\nfib(%d) = %d, time spent = %f\n", n, f, time_spent);
	return 0;
}
