# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <stdbool.h>
// # include <time.h>

# define THRESHOLD 0.0001
# define NUM_ROWS 644
# define NUM_COLS 1056

typedef struct
{
    double nested_value;
    double p_value;
} Nested_elements;

// void abundances_matrix(double **matrix)
// {
    // double **abundances_matrix = (double **) malloc(rows * sizeof(double *));
    // for (i = 0; i < rows; i++) {
        // abundances_matrix[i] = (double *) malloc(cols * sizeof(double))
    // }
// }

void create_relative_abundances(double matrix_relative_abundances[][NUM_COLS], int matrix_absolute_abundances[][NUM_COLS], int num_bacterial_species_per_individual[])
{
    double relative_abundance;

    for (int row = 0; row < NUM_ROWS; row++) {
        for (int col = 0; col < NUM_COLS; col++) {
            relative_abundance = (double) matrix_absolute_abundances[row][col] / (double) num_bacterial_species_per_individual[row];
            matrix_relative_abundances[row][col] = relative_abundance;
        }
    }
}

void create_matrix_individuals(FILE *f_vertebrates, double matrix_individuals[][NUM_COLS])
{
    char line[10000], absolute_abundances_individual[10000], *absolute_abundance;
    int matrix_absolute_abundances[NUM_ROWS][NUM_COLS], num_bacterial_species_per_individual[NUM_ROWS] = {0};
    int row, col = 0;

    fgets(line, sizeof(line), f_vertebrates);       /* S'elimina la primera línia. */
    while (fgets(line, sizeof(line), f_vertebrates) != NULL) {
        sscanf(line, "%*s %[^\n]", absolute_abundances_individual);     /* S'elimina la primera columna. */

        absolute_abundance = strtok(absolute_abundances_individual, " ");
        row = 0;
        while (absolute_abundance != NULL) {
            matrix_absolute_abundances[row][col] = atoi(absolute_abundance);
            num_bacterial_species_per_individual[row] += atoi(absolute_abundance);
            absolute_abundance = strtok(NULL, " ");
            row++;
        }
        col++;
    }
    create_relative_abundances(matrix_individuals, matrix_absolute_abundances, num_bacterial_species_per_individual);
}

void discretize_matrix(double matrix[][NUM_COLS], int binary_matrix[][NUM_COLS], double threshold)
{
    int rows, columns;

    rows = 0;
    while (rows < NUM_ROWS) {
        columns = 0;
        while (columns < NUM_COLS) {
            if (matrix[rows][columns] > threshold) {
                binary_matrix[rows][columns] = 1;
            }
            else {
                binary_matrix[rows][columns] = 0;
            }
            columns++;
        }
        rows++;
    }
}

double calculate_nested_value(int matrix[][NUM_COLS])
{
    int first_isocline, second_isocline, third_isocline, fourth_isocline;
    int first_row, second_row, row, first_col, second_col, col, first_acum, second_acum;

    first_isocline = second_isocline = third_isocline = fourth_isocline = 0;

        /* Calculate the sum of the number of shared interactions between rows. */
    for (first_row = 0; first_row < NUM_ROWS; first_row++) {
        for (second_row = 0; second_row < NUM_ROWS; second_row++) {
            if (first_row < second_row) {
                for (col = 0; col < NUM_COLS; col++) {
                    if ((matrix[first_row][col] == 1) && (matrix[second_row][col] == 1)) {
                        first_isocline++;
                    }
                }
            }
        }
    }

        /* Calculate the sum of the number of shared interactions between columns. */
    for (first_col = 0; first_col < NUM_COLS; first_col++) {
        for (second_col = 0; second_col < NUM_COLS; second_col++) {
            if (first_col < second_col) {
                for (row = 0; row < NUM_ROWS; row++) {
                    if ((matrix[row][first_col] == 1) && (matrix[row][second_col] == 1)) {
                        second_isocline++;
                    }
                }
            }
        }
    }

        /* Calculate the sum of the number of interactions of rows. */
    for (first_row = 0; first_row < NUM_ROWS; first_row++) {
        for (second_row = 0; second_row < NUM_ROWS; second_row++) {
            if (first_row < second_row) {
                first_acum = second_acum = 0;
                for (col = 0; col < NUM_COLS; col++) {
                    first_acum += matrix[first_row][col];
                    second_acum += matrix[second_row][col];
                }
                if (first_acum < second_acum) {
                    third_isocline += first_acum;
                }
                else {
                    third_isocline += second_acum;
                }
            }
        }
    }

        /* Calculate the sum of the number of interactions of columns. */
    for (first_col = 0; first_col < NUM_COLS; first_col++) {
        for (second_col = 0; second_col < NUM_COLS; second_col++) {
            if (first_col < second_col) {
                first_acum = second_acum = 0;
                for (row = 0; row < NUM_ROWS; row++) {
                    first_acum += matrix[row][first_col];
                    second_acum += matrix[row][second_col];
                }
                if (first_acum < second_acum) {
                    fourth_isocline += first_acum;
                }
                else {
                    fourth_isocline += second_acum;
                }
            }
        }
    }

        /* Calculate and return the nested value of the matrix. */
    return ((double)(first_isocline + second_isocline) / (double)(third_isocline + fourth_isocline));
}

double calculate_nested_value_optimized(int matrix[][NUM_COLS])
{
    int sum_rows[NUM_ROWS] = {0};
    int sum_cols[NUM_COLS] = {0};
    int first_isocline, second_isocline, third_isocline, fourth_isocline, row, col;

        /* Calculate and save the number of interactions of every row. */
    for (row = 0; row < NUM_ROWS; row++) {
        for (col = 0; col < NUM_COLS; col++) {
            sum_rows[row] += matrix[row][col];
        }
    }

        /* Calculate and save the number of interactions of every column. */
    for (col = 0; col < NUM_COLS; col++) {
        for (row = 0; row < NUM_ROWS; row++) {
            sum_cols[col] += matrix[row][col];
        }
    }

    first_isocline = second_isocline = third_isocline = fourth_isocline = 0;

        /* Calculate the sum of the number of shared interactions between rows
           and the sum of the minimum of pairs of interactions of rows. */
    for (int first_row = 0; first_row < NUM_ROWS - 1; first_row++) {
        for (int second_row = first_row + 1; second_row < NUM_ROWS; second_row++) {
            for (col = 0; col < NUM_COLS; col++) {
                first_isocline += matrix[first_row][col] & matrix[second_row][col];
            }
            if (sum_rows[first_row] < sum_rows[second_row]) {
                third_isocline += sum_rows[first_row];
            }
            else {
                third_isocline += sum_rows[second_row];
            }
        }
    }

        /* Calculate the sum of the number of shared interactions between columns
           and the sum of the minimum of pairs of the number of interactions of columns. */
    for (int first_col = 0; first_col < NUM_COLS - 1; first_col++) {
        for (int second_col = first_col + 1; second_col < NUM_COLS; second_col++) {
            for (row = 0; row < NUM_ROWS; row++) {
                second_isocline += matrix[row][first_col] & matrix[row][second_col];
            }
            if (sum_cols[first_col] < sum_cols[second_col]) {
                fourth_isocline += sum_cols[first_col];
            }
            else {
                fourth_isocline += sum_cols[second_col];
            }
        }
    }

        /* Calculate and return the nested value of the matrix. */
    return ((double)(first_isocline + second_isocline) / (double)(third_isocline + fourth_isocline));
}

int count_ones_binary_matrix(int matrix[][NUM_COLS])
{
    int num_ones = 0;

    for (int row = 0; row < NUM_ROWS; row++) {
        for (int col = 0; col < NUM_COLS; col++) {
            num_ones += matrix[row][col];
        }
    }
    return num_ones;
}

void initialize_randomized_matrix(int randomized_matrix[][NUM_COLS])
{
    for (int row = 0; row < NUM_ROWS; row++) {
        for (int col = 0; col < NUM_COLS; col++) {
            randomized_matrix[row][col] = 0;
        }
    }
}

void generate_randomized_matrix(int randomized_matrix[][NUM_COLS], int num_ones)
{
    int cont_ones, pos, num_elements = NUM_ROWS * NUM_COLS;

    cont_ones = 0;
    while (cont_ones < num_ones) {
        pos = rand() % num_elements;
        if (randomized_matrix[pos / NUM_COLS][pos % NUM_COLS] != 1) {
            randomized_matrix[pos / NUM_COLS][pos % NUM_COLS] = 1;
            cont_ones++;
        }
    }
}

void generate_nested_values_randomized(int matrix[][NUM_COLS], double nested_values_randomized[], int num_randomized_matrices)
{
    int randomized_matrix[NUM_ROWS][NUM_COLS];
    int num_ones = count_ones_binary_matrix(matrix);

    for (int pos = 0; pos < num_randomized_matrices; pos++) {
        initialize_randomized_matrix(randomized_matrix);
        generate_randomized_matrix(randomized_matrix, num_ones);
        nested_values_randomized[pos] = calculate_nested_value(randomized_matrix);
        // nested_values_randomized[pos] = calculate_nested_value_optimized(randomized_matrix);
    }
}

int sort(double array[], int first, int last)
{
    int pivot = array[first], i = first + 1, j = last, aux;

    while (i < j) {
        if (array[i] > pivot) {
            if (array[j] <= pivot) {
                aux = array[i];
                array[i] = array[j];
                array[j] = aux;
                i++;
            }
            j--;
        }
        else {
            i++;
        }
    }
    if (array[i] < pivot) {
        aux = array[i];
        array[i] = pivot;
        array[first] = aux;
    }
    return i - 1;
}

void quicksort(double array[], int first, int last)
{
    int half;

        /* In direct case we do nothing. */
    if (first < last) {        /* Recursive case. */
        half = sort(array, first, last);
        quicksort(array, first, half);
        quicksort(array, half + 1, last);
    }
}

int get_index(double nested_values[], int num_elements, double nested_value)
{
    int pos = 0;
    bool found = false;

    while ((! found) && (pos < num_elements)) {
        if (nested_values[pos] == nested_value) {
            found = true;
        }
        else {
            pos++;
        }
    }
    return pos;
}

Nested_elements nested_test(int matrix[][NUM_COLS], int num_randomized_matrices)
{
    Nested_elements nested_elements;
    double nested_values[num_randomized_matrices + 1];

        /* Generate as many randomized matrices from the real matrix as it is specified and calculate their nested values. */
    generate_nested_values_randomized(matrix, nested_values, num_randomized_matrices);

        /* Calculate and store the nested value of the real matrix. */
    nested_elements.nested_value = calculate_nested_value(matrix);
    // nested_elements.nested_value = calculate_nested_value_optimized(matrix);
    nested_values[num_randomized_matrices] = nested_elements.nested_value;

        /* Sort the list of nested values. */
    quicksort(nested_values, 0, num_randomized_matrices);

        /* Calculate the fraction of randomized matrices that have a nested value greater than that of the real matrix. */
    nested_elements.p_value = ((double)(num_randomized_matrices - get_index(
            nested_values,num_randomized_matrices + 1, nested_elements.nested_value))
                    / (double) (num_randomized_matrices + 1));

    return nested_elements;
}

int main(int argc, char * argv[])
{
    FILE *f_vertebrates;
    Nested_elements nested_elements;
    double matrix_individuals[NUM_ROWS][NUM_COLS];
    // int binary_matrix[NUM_ROWS][NUM_COLS];
    // double nested_value;

    f_vertebrates = fopen(argv[1],"r");

    if (f_vertebrates == NULL) {
        printf("Error in opening the file.");
    }
    else {
        // srand(time(NULL));

        create_matrix_individuals(f_vertebrates, matrix_individuals);

        for (int i = 0; i < 1; i++) {
            for (int j = 0; j < NUM_COLS; j++) {
                printf(" %f ", matrix_individuals[i][j]);
            }
            printf("\n");
        }

        // discretize_matrix(matrix_individuals_relative_abundances, binary_matrix, THRESHOLD);

        // nested_value = calculate_nested_value(binary_matrix);
        // nested_value = calculate_nested_value_optimized(binary_matrix);
        // printf("\n%.2f\n", nested_value);

        // nested_elements = nested_test(binary_matrix, 1000);
        // printf("\nNested value: %f\nP-value: %f\n", nested_elements.nested_value, nested_elements.p_value);

        fclose(f_vertebrates);
    }
    return 0;
}
