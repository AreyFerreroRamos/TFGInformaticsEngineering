# include <stdio.h>

# define THRESHOLD 0.5
# define NUM_ROWS 4
# define NUM_COLS 4

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

double nestedness(int matrix[][NUM_COLS])
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

double nestedness_optimized(int matrix[][NUM_COLS])
{
    int sum_rows[NUM_ROWS] = {0, 0, 0, 0};
    int sum_cols[NUM_COLS] = {0, 0, 0, 0};
    int first_isocline, second_isocline, third_isocline, fourth_isocline;
    int row, col, first_row, second_row, first_col, second_col;

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
    for (first_row = 0; first_row < NUM_ROWS - 1; first_row++) {
        for (second_row = first_row + 1; second_row < NUM_ROWS; second_row++) {
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
    for (first_col = 0; first_col < NUM_COLS - 1; first_col++) {
        for (second_col = first_col + 1; second_col < NUM_COLS; second_col++) {
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
    int row, col, num_ones = 0;

    for (row = 0; row < NUM_ROWS; row++) {
        for (col = 0; col < NUM_COLS; col++) {
            num_ones += matrix[row][col];
        }
    }
    return num_ones;
}

void initialize_randomized_matrix(randomized_matrix[][NUM_COLS])
{
    int row, col;

    for (row = 0; row < NUM_ROWS; row++) {
        for (col = 0; col < NUM_COLS; col++) {
            randomized_matrix[row][col] = 0;
        }
    }
}

void randomize_matrix(int randomized_matrix[][NUM_COLS], int num_ones)
{
    
}

void generate_nested_values_randomized(int matrix[][NUM_COLS], double nested_values_randomized[], int num_randomized_matrices)
{
    int randomized_matrix[NUM_ROWS][NUM_COLS];
    int pos, num_ones = count_ones_binary_matrix(matrix);

    for (pos = 0; pos < num_randomized_matrices; pos++) {
        initialize_randomized_matrix(randomized_matrix);
        randomize_matrix(randomized_matrix, num_ones);
        nested_values_randomized[pos] = nestedness(randomized_matrix);
        // nested_values_randomized[pos] = nestedness_optimized(randomized_matrix);
    }
}

void sort(double array[])
{

}

int get_index(double nested_values[], double nested_value)
{

}

double nestedness_assesment(int matrix[][NUM_COLS], int num_randomized_matrices)
{
    Nested_elements nested_elements;
    double nested_values[num_randomized_matrices + 1];

    /* Generate as many randomized matrices from the real matrix as it is specified and calculate their nested values. */
    generate_nested_values_randomized(matrix, nested_values, num_randomized_matrices);

    /* Calculate and store the nested value of the real matrix. */
    nested_elements.nested_value = nestedness(matrix);
    // nested_elements.nested_value = nestedness_optimized(matrix);
    nested_values[num_randomized_matrices] = nested_elements.nested_value;

    /* Sort the list of nestedness values. */
    sort(nested_values);

    /* Calculate the fraction of randomized matrices that have a nested value greater than that of the real matrix. */
    nested_elements.p_value = ((double)(num_randomized_matrices - get_index(
            nested_values, nested_elements.nested_value)) / (double) (num_randomized_matrices + 1));

    return nested_elements;
}

int main(int argc, char * argv[])
{
    double abundances_matrix[NUM_ROWS][NUM_COLS] = {
            {0.9, 0.8, 0.7, 0.6},
            {0.8, 0.7, 0.6, 0.5},
            {0.7, 0.6, 0.5, 0.4},
            {0.6, 0.5, 0.4, 0.3}
    };
    int binary_matrix[NUM_ROWS][NUM_COLS];

    // Nested_elements nested_elements;
    double nested_value;
    int i, j;

    for (i = 0; i < NUM_ROWS; i++) {
        for (j = 0; j < NUM_COLS; j++) {
            printf("%.1f\t", abundances_matrix[i][j]);
        }
        printf("\n");
    }

    discretize_matrix(abundances_matrix, binary_matrix, THRESHOLD);

    printf("\n");
    for (i = 0; i < NUM_ROWS; i++) {
        for (j = 0; j < NUM_COLS; j++) {
            printf("%i\t", binary_matrix[i][j]);
        }
        printf("\n");
    }

    nested_value = nestedness(binary_matrix);
    // nested_value = nestedness_optimized(binary_matrix);
    printf("\n%.2f\n", nested_value);

    // nested_elements = nestedness_assesment(abundances_matrix, 1000);
    // printf("\nnested value: %.2f\np-value: %.2f", nested_elements.nested_value, nested_elements.p_value);
}
