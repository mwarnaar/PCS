#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
	double *cylinder;
	int x, y, dimension_x, dimension_y;
	
	dimension_x = 5;
	dimension_y = 5;

	cylinder = (double *) malloc(dimension_x * dimension_y * sizeof(double));

	for (x = 0; x < dimension_x; x++) {
		for (y = 0; y < dimension_y; y++) {
			cylinder[x * dimension_x + y] = x + y;
		}
	}

	for (x = 0; x < dimension_x; x++) {
		for (y = 0; y < dimension_y; y++) {
			printf("Cylinder[%d][%d] = %f\n", x, y, cylinder[x * dimension_x + y]);	
		}
	}

	while (1) {
		
	}

	free(cylinder);

	return 0;
}
