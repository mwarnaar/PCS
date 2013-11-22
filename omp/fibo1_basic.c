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
		i = fib(n - 1);
		j = fib(n - 2);

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

	f = fib(n);

	end = clock();
	time_spent = (double) (end - begin) / CLOCKS_PER_SEC;

	printf("\nfib(%d) = %d, time spent = %f\n", n, f, time_spent);
	return 0;
}
