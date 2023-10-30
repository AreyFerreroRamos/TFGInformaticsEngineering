# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <stdbool.h>
// # include <time.h>

# define NUM_INDIVIDUALS 644
# define NUM_VERTEBRATES 50
# define NUM_BACTERIAL_GENUS 1056
# define THRESHOLD 0.0001

typedef struct
{
   char *vertebrate;
   char *sample_type;
} Individual_type;

typedef struct
{
    double nested_value;
    double p_value;
} Nested_elements;

void create_relative_abundances(int num_rows, int num_cols, int matrix_absolute_abundances[][num_cols],
                                double matrix_relative_abundances[][num_cols], int num_bacterial_species_per_individual[])
{
    double relative_abundance;

    for (int row = 0; row < num_rows; row++) {
        for (int col = 0; col < num_cols; col++) {
            relative_abundance = (double) matrix_absolute_abundances[row][col]
                    /(double) num_bacterial_species_per_individual[row];
            matrix_relative_abundances[row][col] = relative_abundance;
        }
    }
}

void create_matrix_individuals(FILE *f_vertebrates, double matrix_individuals[][NUM_BACTERIAL_GENUS])
{
    char line[10000], absolute_abundances_individual[10000], *absolute_abundance;
    int matrix_absolute_abundances[NUM_INDIVIDUALS][NUM_BACTERIAL_GENUS];
    int num_bacterial_species_per_individual[NUM_INDIVIDUALS] = {0};
    int row, col = 0;

    fgets(line, sizeof(line), f_vertebrates);       /* S'elimina la primera línia. */
    while (fgets(line, sizeof(line), f_vertebrates) != NULL) {
        sscanf(line, "%*s %[^\n]", absolute_abundances_individual);     /* S'elimina la primera columna. */

        absolute_abundance = strtok(absolute_abundances_individual, " ");
        row = 0;
        while (absolute_abundance != NULL) {
            matrix_absolute_abundancesinformation[row][col] = atoi(absolute_abundance);
            num_bacterial_species_per_individual[row] += atoi(absolute_abundance);
            absolute_abundance = strtok(NULL, " ");
            row++;
        }
        col++;
    }
    create_relative_abundances(NUM_INDIVIDUALS, NUM_BACTERIAL_GENUS, matrix_absolute_abundances,
                               matrix_individuals, num_bacterial_species_per_individual);
}

individual_type get_individual_type(char *individual, FILE *f_metadata)
{
    Individual_type individual_type;
    char line[1000], *sample;
    bool found = false;

    fgets(line, sizeof(line), f_metadata);
    while ((fgets(line, sizeof(line), f_metadata) != NULL) && (! found) {
        strtok(line, ";");
        sample = strtok(line, ";");
        if (strcmp(individual, sample) == 0) {
            strtok(line, ";");
            individual_type.vertebrate = strtok(line, ";");
            strtok(line, ";");
            individual_type.sample_type = strtok(line, ";");
            found = true;
        }
    }
    return individual_type
}

void create_matrix_vertebrates(FILE *f_vertebrates, FILE *f_metadata, double matrix_vertebrates[][NUM_BACTERIAL_GENUS])
{
    char line[10000], absolute_abundances_individual[100000], *absolute_abundances, *individual;
    char *individuals[NUM_INDIVIDUALS];
    int matrix_absolute_abundances[NUM_VERTEBRATES][NUM_BACTERIAL_GENUS];
    int num_bacterial_species_per_individual[NUM_INDIVIDUALS] = {0};
    int row, col = 0;

    fgets(line, sizeof(line), f_vertebrates);
    individual = strtok(line, " ");
    while (individual != NULL) {
        individuals[col++] = individual;
        // strncpy(individual2, individual + 1, strlen(individual) - 2);
        // sscanf(individual, "%*c%[^\n]%*c", individuals[col]);
        individual = strtok(NULL, " ");
    }
    individuals[col++] = individual;
    // strncpy(individual2, individual + 1, strlen(individual) - 2);
    // sscanf(individual, "%*c %[^\n] %*c", individuals[col]);

    while(fgets(line, sizeof(line), f_vertebrates)) {

    }

    for (int i = 0; i < NUM_INDIVIDUALS; i++) {
       printf(" %s ", strtok(individuals[i], "\""));
    }
}

void discretize_matrix(int num_rows, int num_cols, double matrix[][num_cols],
                       int binary_matrix[][num_cols], double threshold)
{
    for (int row = 0; row < num_rows; row++) {
        for (int col = 0; col < num_cols; col++) {
            if (matrix[row][col] > threshold) {
                binary_matrix[row][col] = 1;
            }
            else {
                binary_matrix[row][col] = 0;
            }
        }
    }
}

double calculate_nested_value(int num_rows, int num_cols, int matrix[][num_cols])
{
    int first_isocline, second_isocline, third_isocline, fourth_isocline;
    int first_row, second_row, row, first_col, second_col, col, first_acum, second_acum;

    first_isocline = second_isocline = third_isocline = fourth_isocline = 0;

        /* Calculate the sum of the number of shared interactions between rows. */
    for (first_row = 0; first_row < num_rows; first_row++) {
        for (second_row = 0; second_row < num_rows; second_row++) {
            if (first_row < second_row) {
                for (col = 0; col < num_cols; col++) {
                    if ((matrix[first_row][col] == 1) && (matrix[second_row][col] == 1)) {
                        first_isocline++;
                    }
                }
            }
        }
    }

        /* Calculate the sum of the number of shared interactions between columns. */
    for (first_col = 0; first_col < num_cols; first_col++) {
        for (second_col = 0; second_col < num_cols; second_col++) {
            if (first_col < second_col) {
                for (row = 0; row < num_rows; row++) {
                    if ((matrix[row][first_col] == 1) && (matrix[row][second_col] == 1)) {
                        second_isocline++;
                    }
                }
            }
        }
    }

        /* Calculate the sum of the number of interactions of rows. */
    for (first_row = 0; first_row < num_rows; first_row++) {
        for (second_row = 0; second_row < num_rows; second_row++) {
            if (first_row < second_row) {
                first_acum = second_acum = 0;
                for (col = 0; col < num_cols; col++) {
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
    for (first_col = 0; first_col < num_cols; first_col++) {
        for (second_col = 0; second_col < num_cols; second_col++) {
            if (first_col < second_col) {
                first_acum = second_acum = 0;
                for (row = 0; row < num_rows; row++) {
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

double calculate_nested_value_optimized(int num_rows, int num_cols, int matrix[][num_cols])
{
    int sum_rows[NUM_INDIVIDUALS] = {0};
    // int sum_rows[NUM_VERTEBRATES] = {0};
    int sum_cols[NUM_BACTERIAL_GENUS] = {0};
    int first_isocline, second_isocline, third_isocline, fourth_isocline, row, col;

        /* Calculate and save the number of interactions of every row. */
    for (row = 0; row < num_rows; row++) {
        for (col = 0; col < num_cols; col++) {
            sum_rows[row] += matrix[row][col];
        }
    }

        /* Calculate and save the number of interactions of every column. */
    for (col = 0; col < num_cols; col++) {
        for (row = 0; row < num_rows; row++) {
            sum_cols[col] += matrix[row][col];
        }
    }

    first_isocline = second_isocline = third_isocline = fourth_isocline = 0;

        /* Calculate the sum of the number of shared interactions between rows
           and the sum of the minimum of pairs of interactions of rows. */
    for (int first_row = 0; first_row < num_rows - 1; first_row++) {
        for (int second_row = first_row + 1; second_row < num_rows; second_row++) {
            for (col = 0; col < num_cols; col++) {
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
    for (int first_col = 0; first_col < num_cols - 1; first_col++) {
        for (int second_col = first_col + 1; second_col < num_cols; second_col++) {
            for (row = 0; row < num_rows; row++) {
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

int count_ones_binary_matrix(int num_rows, int num_cols, int matrix[][num_cols])
{
    int num_ones = 0;

    for (int row = 0; row < num_rows; row++) {
        for (int col = 0; col < num_cols; col++) {
            num_ones += matrix[row][col];
        }
    }
    return num_ones;
}

void initialize_randomized_matrix(int num_rows, int num_cols, int randomized_matrix[][num_cols])
{
    for (int row = 0; row < num_rows; row++) {
        for (int col = 0; col < num_cols; col++) {
            randomized_matrix[row][col] = 0;
        }
    }
}

void generate_randomized_matrix(int num_rows, int num_cols, int num_ones, int randomized_matrix[][num_cols])
{
    int cont_ones, pos, num_elements = num_rows * num_cols;

    cont_ones = 0;
    while (cont_ones < num_ones) {
        pos = rand() % num_elements;
        if (randomized_matrix[pos / num_cols][pos % num_cols] != 1) {
            randomized_matrix[pos / num_cols][pos % num_cols] = 1;
            cont_ones++;
        }
    }
}

void generate_nested_values_randomized(int num_rows, int num_cols, int matrix[][num_cols],
                                       int num_randomized_matrices, double nested_values_randomized[])
{
    int randomized_matrix[num_rows][num_cols];
    int num_ones = count_ones_binary_matrix(num_rows, num_cols, matrix);

    for (int pos = 0; pos < num_randomized_matrices; pos++) {
        initialize_randomized_matrix(num_rows, num_cols, randomized_matrix);
        generate_randomized_matrix(num_rows, num_cols, num_ones, randomized_matrix);
        nested_values_randomized[pos] = calculate_nested_value(num_rows, num_cols, randomized_matrix);
        // nested_values_randomized[pos] = calculate_nested_value_optimized(num_rows, num_cols, randomized_matrix);
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

Nested_elements nested_test(int num_rows, int num_cols, int matrix[][num_cols], int num_randomized_matrices)
{
    Nested_elements nested_elements;
    double nested_values[num_randomized_matrices + 1];

        /* Generate as many randomized matrices from the real matrix as it is specified and calculate their nested values. */
    generate_nested_values_randomized(num_rows, num_cols, matrix, num_randomized_matrices,
                                      nested_values);

        /* Calculate and store the nested value of the real matrix. */
    nested_elements.nested_value = calculate_nested_value(num_rows, num_cols, matrix);
    // nested_elements.nested_value = calculate_nested_value_optimized(num_rows, num_cols, matrix);
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
    FILE *f_vertebrates, *f_metadata;
    Nested_elements nested_elements;
    double matrix_individuals[NUM_INDIVIDUALS][NUM_BACTERIAL_GENUS];
    // double matrix_vertebrates[NUM_VERTEBRATES][NUM_BACTERIAL_GENUS];
    // int binary_matrix_individuals[NUM_INDIVIDUALS][NUM_BACTERIAL_GENUS];
    // int binary_matrix_vertebrates[NUM_VERTEBRATES][NUM_BACTERIAL_GENUS];
    // double nested_value;

    f_vertebrates = fopen(argv[1],"r");
    f_metadata = fopen(argv[2], "r");

    if (f_vertebrates == NULL) {
        printf("Error in opening the file %s.", argv[1]);
    }
    else {
        if (f_metadata == NULL) {
            printf("Error in opening the file %s.", argv[2]);
        }
        else {
            // srand(time(NULL));

            create_matrix_individuals(f_vertebrates, matrix_individuals);
            // create_matrix_vertebrates(f_vertebrates, f_metadata, matrix_vertebrates);

            for (int i = 0; i < 1; i++) {
                for (int j = 0; j < NUM_BACTERIAL_GENUS; j++) {
                    // printf(" %f ", matrix_individuals[i][j]);
                    // printf(" %f ", matrix_vertebrates[i][j]);
                }
                printf("\n");
            }
            fflush(stdout);

            // discretize_matrix(NUM_INDIVIDUALS, NUM_BACTERIAL_GENUS, matrix_individuals, binary_matrix_individuals, THRESHOLD);
            // discretize_matrix(NUM_VERTEBRATES, NUM_BACTERIAL_GENUS, matrix_vertebrates, binary_matrix_vertebrates, THRESHOLD);


            /*for (int i = 0; i < 1; i++) {
                for (int j = 0; j < NUM_BACTERIAL_GENUS; j++) {
                    printf(" %i ", binary_matrix_individuals[i][j]);
                }
                printf("\n");
            }*/

            // nested_value = calculate_nested_value(NUM_INDIVIDUALS, NUM_BACTERIAL_GENUS, binary_matrix_individuals);
            // nested_value = calculate_nested_value(NUM_VERTEBRATES, NUM_BACTERIAL_GENUS, binary_matrix_vertebrates);
            // nested_value = calculate_nested_value_optimized(NUM_INDIVIDUALS, NUM_BACTERIAL_GENUS, binary_matrix_individuals);
            // nested_value = calculate_nested_value_optimized(NUM_VERTEBRATES, NUM_BACTERIAL_GENUS, binary_matrix_vertebrates);
            // printf("\n%.2f\n", nested_value);

            // nested_elements = nested_test(NUM_INDIVIDUALS, NUM_BACTERIAL_GENUS, binary_matrix_individuals, 1000);
            // nested_elements = nested_test(NUM_VERTEBRATES, NUM_BACTERIAL_GENUS, binary_matrix_vertebrates, 1000);
            // printf("\nNested value: %f\nP-value: %f\n", nested_elements.nested_value, nested_elements.p_value);

            fclose(f_metadata);
        }
        fclose(f_vertebrates);
    }
    return 0;
}
