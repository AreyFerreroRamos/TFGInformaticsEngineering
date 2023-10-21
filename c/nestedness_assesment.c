# include <stdio.h>

# define THRESHOLD 0.5
# define NUM_ROWS 4
# define NUM_COLS 4

// void abundances_matrix(double **matrix)
// {
    // double **abundances_matrix = (double **) malloc(rows * sizeof(double *));
    // for (i = 0; i < rows; i++) {
        // abundances_matrix[i] = (double *) malloc(cols * sizeof(double))
    // }
// }

void discretize_matrix(double matrix[][NUM_COLS], double threshold)
{
    int rows, columns;

    rows = 0;
    while (rows < NUM_ROWS) {
        columns = 0;
        while (columns < NUM_COLS) {
            if (matrix[rows][columns] > threshold) {
                matrix[rows][columns] = 1;
            }
            else {
                matrix[rows][columns] = 0;
            }
            columns++;
        }
        rows++;
    }
}

double nestedness(double matrix[][NUM_COLS])
{
    int first_isocline, second_isocline, third_isocline, fourth_isocline;
    int first_row, second_row, row, first_col, second_col, col, first_acum, second_acum;

    first_isocline = second_isocline = third_isocline = fourth_isocline = 0;

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

    return ((double)(first_isocline + second_isocline) / (double)(third_isocline + fourth_isocline));
}

int main(int argc, char * argv[])
{
    double nested_value;
    int i, j;

    double abundances_matrix[NUM_ROWS][NUM_COLS] = {
            {0.9, 0.8, 0.7, 0.6},
            {0.8, 0.7, 0.6, 0.5},
            {0.7, 0.6, 0.5, 0.4},
            {0.6, 0.5, 0.4, 0.3}
    };
    for (i = 0; i < NUM_ROWS; i++) {
        for (j = 0; j < NUM_COLS; j++) {
            printf("%.1f\t", abundances_matrix[i][j]);
        }
        printf("\n");
    }

    discretize_matrix(abundances_matrix, THRESHOLD);

    printf("\n");
    for (i = 0; i < NUM_ROWS; i++) {
        for (j = 0; j < NUM_COLS; j++) {
            printf("%.1f\t", abundances_matrix[i][j]);
        }
        printf("\n");
    }

    nested_value = nestedness(abundances_matrix);
    printf("\n%.2f\n", nested_value);

    // nestedness_assesment(abundances_matrix, atoi(argv[5]))
    // printf("\nnested value: %.2f\np-value: %.2f", nested_components.nested_value, nested_components.p_value)
}