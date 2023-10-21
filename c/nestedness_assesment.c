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

int main(int argc, char * argv[])
{
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

    // nestedness_assesment(abundances_matrix, atoi(argv[5]))
    // printf("\nnested value: %.2f\np-value: %.2f", nested_components.nested_value, nested_components.p_value)
}