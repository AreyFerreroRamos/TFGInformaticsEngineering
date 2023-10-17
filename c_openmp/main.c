#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>

#include <assert.h>
#include <omp.h>

#define SIZE_X 2048
#define SIZE_Y 4096


double nestedness(bool M[][SIZE_X]) {
	unsigned long row_interactions = 0;
	unsigned long col_interactions = 0;
	unsigned long row_count[SIZE_Y] = {0};
	unsigned long col_count[SIZE_X] = {0};
	unsigned long row_isa = 0;
	unsigned long col_isa = 0;
	for (int y = 0; y < SIZE_Y; y++) {
		for (int x = 0; x < SIZE_X; x++) {
			if (M[y][x]) {
				row_count[y]++;
			}
		}
	}
	for (int x = 0; x < SIZE_X; x++) {
		for (int y = 0; y < SIZE_Y; y++) {
			col_count[x] += M[y][x];
		}
	}
	for (int i = 0; i < SIZE_Y-1; i++) {
		for (int ii = i+1; ii < SIZE_Y; ii++) {
			for (int x = 0; x < SIZE_X; x++) {
				if (M[i][x] & M[ii][x]) {
					row_interactions++;
				}
			}
			row_isa += row_count[i] < row_count[ii] ? row_count[i] : row_count[ii];
		}
	}
	for (int j = 0; j < SIZE_X-1; j++) {
		for (int jj = j+1; jj < SIZE_X; jj++) {
			for (int y = 0; y < SIZE_Y; y++) {
				if (M[y][j] & M[y][jj]) {
					col_interactions++;
				}
			}
			col_isa += col_count[j] < col_count[jj] ? col_count[j] : col_count[jj];
		}
	}
	return (double) (row_interactions + col_interactions) / (double) (row_isa + col_isa);
}


void calc_counts(bool M[][SIZE_X], unsigned int row_count[SIZE_Y], unsigned int col_count[SIZE_X]) {
	for (int y = 0; y < SIZE_Y; y++) {
		for (int x = 0; x < SIZE_X; x++) {
			if (M[y][x]) {
				row_count[y]++;
				col_count[x]++;
			}
		}
	}
}


double nestedness_aux(bool M[][SIZE_X], unsigned int row_count[SIZE_Y], unsigned int col_count[SIZE_X]) {
	unsigned long row_interactions = 0;
	unsigned long col_interactions = 0;
	unsigned long row_isa = 0;
	unsigned long col_isa = 0;

	#pragma omp parallel for default(none) shared(M, row_count) reduction(+:row_interactions,row_isa) schedule(static)
	for (int i = 0; i < SIZE_Y-1; i++) {
		for (int ii = i+1; ii < SIZE_Y; ii++) {
			for (int x = 0; x < SIZE_X; x++) {
				if (M[i][x] && M[ii][x]) {
					row_interactions++;
				}
			}
			row_isa += (row_count[i] < row_count[ii]) ? row_count[i] : row_count[ii];
		}
	}
	#pragma omp parallel for default(none) shared(M, col_count) reduction(+:col_interactions,col_isa) schedule(static)
	for (int j = 0; j < SIZE_X-1; j++) {
		for (int jj = j+1; jj < SIZE_X; jj++) {
			for (int y = 0; y < SIZE_Y; y++) {
				if (M[y][j] && M[y][jj]) {
					col_interactions++;
				}
			}
			col_isa += (col_count[j] < col_count[jj]) ? col_count[j] : col_count[jj];
		}
	}
	double interactions = row_interactions + col_interactions;
	double isa = row_isa + col_isa;
	return interactions / isa;
}

/// FIXME: Incorrect results
double single(bool M[][SIZE_X], unsigned int row_count[SIZE_Y], unsigned int col_count[SIZE_X]) {
	unsigned long row_interactions = 0;
	unsigned long col_interactions = 0;
	unsigned long row_isa = 0;
	unsigned long col_isa = 0;

	const unsigned int END = (SIZE_X > SIZE_Y) ? SIZE_X : SIZE_Y;

	#pragma omp parallel for default(none) shared(M, END, row_count,col_count) reduction(+:row_interactions,row_isa,col_interactions,col_isa) schedule(static)
	for (int k = 0; k < END-1; k++) {
		for (int kk = k+1; kk < END; kk++) {
			if (kk < SIZE_Y) {
				for (int x = 0; x < SIZE_X; x++) {
					if (M[k][x] && M[kk][x]) {
						row_interactions++;
					}
				}
				row_isa += (row_count[k] < row_count[kk]) ? row_count[k] : row_count[kk];
			}
			if (kk < SIZE_X) {
				for (int y = 0; y < SIZE_Y; y++) {
					if (M[y][k] && M[y][kk]) {
						col_interactions++;
					}
				}
				col_isa += (col_count[k] < col_count[kk]) ? col_count[k] : col_count[kk];
			}
		}
	}
	double interactions = row_interactions + col_interactions;
	double isa = row_isa + col_isa;
	return interactions / isa;
}


void init(bool M[][SIZE_X]) {
	double d_x, d_y;
	const double s_x = SIZE_X;
	const double s_y = SIZE_Y;
	for (int y = 0; y < SIZE_Y; ++y) {
		for (int x = 0; x < SIZE_X; ++x) {
			M[y][x] = (rand() & 0b1);
			d_x = x / s_x;
			d_y = y / s_y;
			if ((d_x + d_y) <= 0.5) {
				M[y][x] |= (rand() & 0b1); // some extra nestedness on upper diagonals
				/*
				M[y][x] |= (rand() & 0b1);
				M[y][x] |= (rand() & 0b1);
				M[y][x] |= (rand() & 0b1);
				*/
			}
		}
	}
}


int main(int argc, char * argv[]) {
	bool (*M)[SIZE_Y][SIZE_X] = malloc(sizeof(bool[SIZE_Y][SIZE_X]));
	srand(0);

	init((*M));
	unsigned int row_count[SIZE_Y] = {0};
	unsigned int col_count[SIZE_X] = {0};
	calc_counts((*M), row_count, col_count);
	double start, stop, delta;
	printf("omp_get_max_threads() = %d\n", omp_get_max_threads());
	for (int k = 0; k < 10; ++k) {
		{
			start = omp_get_wtime();
			double n_aux = nestedness_aux((*M), row_count, col_count);
			stop = omp_get_wtime();
			delta = stop - start;
			printf("[%.2lf] nestedness_aux(M[%dx%d] @k = %d) = %lf\n", delta, SIZE_X, SIZE_Y, k, n_aux);
		}

		{
			// start = omp_get_wtime();
			// double n_single = single((*M), row_count, col_count);
			// stop = omp_get_wtime();
			// delta = stop - start;
			// printf("[%.2lf] single        (M[%dx%d] @k = %d) = %lf\n", delta, SIZE_X, SIZE_Y, k, n_single);
		}

		{
			// start = omp_get_wtime();
			// double n = nestedness((*M));
			// stop = omp_get_wtime();
			// delta = stop - start;
			// printf("[%.2lf] nestedness    (M[%dx%d] @k = %d) = %lf\n", delta, SIZE_X, SIZE_Y, k, n);
		}

		//assert(n == n_aux);
		//assert(n_aux == n_single);
	}
	return 0;
}
