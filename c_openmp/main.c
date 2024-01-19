#include <stdint.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>

#include <assert.h>
#include <omp.h>

#define SIZE_X 1024
#define SIZE_Y 512
// #define SIZE_X 1024*32
// #define SIZE_Y 1024*64


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
			if (M[y][x]) {
				col_count[x]++;
			}
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
	double interactions = row_interactions + col_interactions;
	unsigned long isa = row_isa + col_isa;
	return interactions / isa;
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

// using counts
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
	unsigned long isa = row_isa + col_isa;
	return interactions / isa;
}

// using counts and formula for interactions (triangular number [n])
double nestedness_aux_opti(unsigned int row_count[SIZE_Y], unsigned int col_count[SIZE_X]) {
	unsigned long row_interactions = 0;
	unsigned long col_interactions = 0;
	unsigned long row_isa = 0;
	unsigned long col_isa = 0;

	#pragma omp parallel for default(none) shared(row_count) reduction(+:col_interactions, row_isa) schedule(static)
	for (int i = 0; i < SIZE_Y; i++) {
		col_interactions += ((row_count[i])*(row_count[i]-1))/2;
		for (int ii = i+1; ii < SIZE_Y; ii++) {
			row_isa += (row_count[i] < row_count[ii]) ? row_count[i] : row_count[ii];
		}
	}
	#pragma omp parallel for default(none) shared(col_count) reduction(+:row_interactions, col_isa) schedule(static)
	for (int j = 0; j < SIZE_X; j++) {
		row_interactions += ((col_count[j])*(col_count[j]-1))/2;
		for (int jj = j+1; jj < SIZE_X; jj++) {
			col_isa += (col_count[j] < col_count[jj]) ? col_count[j] : col_count[jj];
		}
	}
	return (double) (row_interactions + col_interactions) / (row_isa + col_isa);
}

// function to find the partition position
// used for quicksort()
int partition(unsigned int array[], int low, int high) {
	unsigned int pivot = array[high];
	unsigned int _tmp;

	int i = (low - 1);
	for (int j = low; j < high; j++) {
		if (array[j] <= pivot) {
			i++;
			// swap [i] <-> [j]
			_tmp = array[i];
			array[i] = array[j];
			array[j] = _tmp;
		}
	}
	// swap [i+1] <-> [high]
	_tmp = array[i+1];
	array[i+1] = array[high];
	array[high] = _tmp;

	return (i + 1);
}

// https://www.programiz.com/dsa/quick-sort
void quicksort(unsigned int array[], int low, int high) {
	if (low < high) {
		int pivot = partition(array, low, high);
		quicksort(array, low, pivot - 1);
		quicksort(array, pivot + 1, high);
	}
}

// using formula for interactions & isa; needs count arrays to be sorted
double nestedness_sorted(unsigned int row_count_sorted[SIZE_Y], unsigned int col_count_sorted[SIZE_X]) {
	unsigned long interactions = 0;
	unsigned long isa = 0;

	//#pragma omp parallel for default(none) shared(row_count) reduction(+:col_interactions, row_isa) schedule(static)
	for (int i = 0; i < SIZE_Y; i++) {
		unsigned int v = row_count_sorted[i];
		interactions += (v * (v - 1)) / 2;

		unsigned int remaining = (SIZE_Y-1)-i;
		isa += v * remaining;
	}
	//#pragma omp parallel for default(none) shared(col_count) reduction(+:row_interactions, col_isa) schedule(static)
	for (int j = 0; j < SIZE_X; j++) {
		unsigned int v = col_count_sorted[j];
		interactions += ((v) * (v - 1)) / 2;
		unsigned int remaining = (SIZE_X-1)-j;
		isa += v * remaining;
	}
	return (double) interactions / isa;
}

// init M
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
				M[y][x] |= (rand() & 0b1);
				// M[y][x] |= (rand() & 0b1);
				// M[y][x] |= (rand() & 0b1);
			}
		}
	}
}

int bytes_to_human_readable(unsigned long n_bytes, char * buf, unsigned int len_buf) {
	const char * suffix[] = {"B", "KiB", "MiB", "GiB", "TiB"};
	const char len_suffix = sizeof(suffix) / sizeof(suffix[0]);
	int suffix_index = 0;
	double n = n_bytes;
	for (suffix_index = 0; ((n_bytes / 1024) > 0) && (suffix_index<(len_suffix-1)); suffix_index++) {
		n = n_bytes / 1024.0;
		n_bytes /= 1024;
	}
	return snprintf(buf, len_buf, "%.02lf %s", n, suffix[suffix_index]);
}

int main(int argc, char * argv[]) {
	unsigned long int n_bytes = sizeof(bool[SIZE_Y][SIZE_X]);
	assert(n_bytes < SIZE_MAX);

	bool (*M)[SIZE_Y][SIZE_X] = malloc(n_bytes);
	srand(0);

	char buf[256];
	bytes_to_human_readable(n_bytes, buf, sizeof(buf)/sizeof(buf[0]));
	printf("bytes(M) = %ld (%s)\n", n_bytes, buf);
	
	/// times (Seconds)
	double start, stop, delta;
	
	// init M:
	{
		start = omp_get_wtime();
		init((*M));
		stop = omp_get_wtime();
		delta = stop - start;
		printf("[%.3lfs] init(M[%dx%d])\n", delta, SIZE_X, SIZE_Y);
	}

	// get counts:
	unsigned int row_count[SIZE_Y] = {0};
	unsigned int col_count[SIZE_X] = {0};
	{
		start = omp_get_wtime();
		calc_counts((*M), row_count, col_count);
		stop = omp_get_wtime();
		delta = stop - start;
		printf("[%.3lfs] calc_counts(M[%dx%d])\n", delta, SIZE_X, SIZE_Y);
	}

	/*
	/// Real nestedness
	double n;
	{
		start = omp_get_wtime();
		n = nestedness((*M));
		stop = omp_get_wtime();
		delta = stop - start;
		printf("[%.3lfs] nestedness(M[%dx%d]) = %lf\n", delta, SIZE_X, SIZE_Y, n);
	}
	*/

	/// Real nestedness (using nestedness_aux_opti(...))
	double n;
	{
		start = omp_get_wtime();
		n = nestedness_aux_opti(row_count, col_count);
		stop = omp_get_wtime();
		delta = stop - start;
		printf("[%.3lfs] nestedness_aux_opti(M[%dx%d]) = %lf\n", delta, SIZE_X, SIZE_Y, n);
	}

	// sort counts:
	{
		//quicksort(row_count, 0, SIZE_Y-1);
		//quicksort(col_count, 0, SIZE_X-1);

		{
			start = omp_get_wtime();
			quicksort(row_count, 0, SIZE_Y-1);
			stop = omp_get_wtime();
			delta = stop - start;
			printf("[%.3lfs] quickSort(row_count[%d]) \n", delta, SIZE_Y);
		}
		{
			start = omp_get_wtime();
			quicksort(col_count, 0, SIZE_X-1);
			stop = omp_get_wtime();
			delta = stop - start;
			printf("[%.3lfs] quickSort(col_count[%d]) \n", delta, SIZE_X);
		}
	}

	// benchmarks:
	/// # of times to benchmark
	const unsigned int ITERATIONS = 8;
	/// total time (Seconds) spent doing benchmark work. Used for calculating delta_avg
	double total_delta = 0;
	printf("omp_get_max_threads() = %d\n", omp_get_max_threads());
	for (unsigned int iteration = 1; iteration <= ITERATIONS; ++iteration) {
		{
			start = omp_get_wtime();
			double n_aux = nestedness_sorted(row_count, col_count);
			stop = omp_get_wtime();
			delta = stop - start;
			assert(n == n_aux); // results must match real nestedness value
			printf("[%.3lfs] nestedness_sorted(M[%dx%d] @iteration = %d) = %lf\n", delta, SIZE_X, SIZE_Y, iteration, n_aux);
			total_delta += delta;
		}
	}
	double delta_avg = total_delta / ITERATIONS;
	printf("delta_avg = %.3lfs (%d iterations)\n", delta_avg, ITERATIONS);

	free(M);
	return 0;
}
