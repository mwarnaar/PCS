#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define A(i,j) a[i*n + j]
#define B(i,j) b[i*n + j]
#define X(i,j) x[i*n + j]


#ifdef MATMUL_ACC
#include <openacc.h>
#endif

#ifdef MATMUL_OMP
#include <omp.h>
#endif


#ifdef MATMUL_SEQ
void matmul(float *a, float *b, float *x, int len){

	/* Perform sequential multiplication */
	int sum = 0;
	
	for (int i = 0; i < len; i++) {

		for (int j = 0; j < len; j++) {

			for (int k = 0; k < len; k++) {
				sum += a[i * len + k] * b[k * len + j];  
			}

			x[i * len + j] = sum;
			sum = 0;
		}

	}  
}
#endif

#ifdef MATMUL_ACC
void matmul(float *a, float *b, float *x, int len){
	/* Perform parallel GPU threads computation */
	int lensq = len * len;
#pragma acc parallel copyin(a[0:lensq], b[0:lensq], len) copyout(x[0:lensq]) num_gangs(480)
	{
#pragma acc loop gang
		for (int i = 0; i < len; i++) {
#pragma acc loop worker
			for (int j = 0; j < len; j++) {
				int sum = 0;
				for (int k = 0; k < len; k++) {
					sum += a[i * len + k] * b[k * len + j];  
				}

				x[i * len + j] = sum;
			}
		}  
	}
}
#endif

#ifdef MATMUL_OMP
void matmul(float *a, float *b, float *x, int len){
	int i;
	
	/* Perform parallel CPU threads multiplication */
	int sum = 0;
#pragma omp parallel for private(i)
	for (int i = 0; i < len; i++) {

		for (int j = 0; j < len; j++) {

			for (int k = 0; k < len; k++) {
				sum += a[i * len + k] * b[k * len + j];  
			}

			x[i * len + j] = sum;
			sum = 0;
		}

	}  
}

#endif

void info(){
#ifdef MATMUL_SEQ
	printf("Running matmul sequential\n");
#endif

#ifdef MATMUL_SSE
	printf("Running matmul using vectorization (sse)\n");
#endif

#ifdef MATMUL_ACC
	printf("Running matmul using OpenACC\n");
#endif

#ifdef MATMUL_OMP
#pragma omp parallel
	printf("Running matmul using OpenMP (Thread %d of %d threads)\n", 
			omp_get_thread_num(),
			omp_get_num_threads());
#endif
}

void do_compute(int n){
	float *a = NULL, *b = NULL, *x = NULL;
	struct timeval before, after;

	size_t memsize =  3*(n*n*sizeof(float));
  	printf("Try to allocate %llu * %llu * %llu = %llu bytes (%llu MB) of memory.\n", 3, n * n, sizeof(float), memsize, memsize/(1024*1024));

	a = malloc(n * n * sizeof(float));
	b = malloc(n * n * sizeof(float));
	x = malloc(n * n * sizeof(float));

	/* Did we have sufficient memory? */
	if(a == NULL || b == NULL || x == NULL){
		printf("Memory Alloc error...\n");
		exit(-1);
	}

	printf("\nInitialize...\n");
	
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			A(i, j) = i + j;
			B(i, j) = i - j;
			X(i, j) = 0;
		} 
	}

	printf("Ready to compute...\n");

	gettimeofday(&before, NULL);
	matmul(a, b, x, n);
	gettimeofday(&after, NULL);

	double time = (double)(after.tv_sec - before.tv_sec) +
		(double)(after.tv_usec - before.tv_usec) / 1e6;

	printf("Computed in %f seconds!\n", time);

	free(a);
	free(b);
	free(x);
}

int main(int argc, char **argv){

	int n = 100;

	if(argc > 1){
		n = atoi(argv[1]);
	}

#ifdef MATMUL_OMP
	int t = 1;

	if(argc > 2){
		t = atoi(argv[2]);
	}
	omp_set_num_threads(t);
#endif

	info();

	printf("Multiply %dx%d matrices.\n", n, n);

	do_compute(n);

	return 0;

}

