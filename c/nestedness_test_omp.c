# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <stdbool.h>
# include <time.h>
# include <omp.h>

# define NUM_INDIVIDUALS 644
# define NUM_VERTEBRATES 50
# define NUM_BACTERIAL_GENUS 1056

typedef struct
{
    char code[7];
    char vertebrate[6];
    int sample_type;
} Individual;

typedef struct
{
    double nested_value;
    double p_value;
} Nested_elements;

void select_matrix(char *name_matrix, int *num_rows, int *num_cols)
{
    if (strcmp(name_matrix, "individuals") == 0) {
        *num_rows = NUM_INDIVIDUALS;
        *num_cols = NUM_BACTERIAL_GENUS;
    }
    else if (strcmp(name_matrix, "vertebrates") == 0) {
        *num_rows = NUM_VERTEBRATES;
        *num_cols = NUM_BACTERIAL_GENUS;
    }
}

void free_memory_doubles_matrix(double **matrix, int num_rows)
{
    for (int row = 0; row < num_rows; row++) {
        free(matrix[row]);
    }
    free(matrix);
}

void free_memory_integers_matrix(int **matrix, int num_rows)
{
    for (int row = 0; row < num_rows; row++) {
        free(matrix[row]);
    }
    free(matrix);
}

void free_memory_shorts_matrix(short **matrix, int num_rows)
{
    for (int row = 0; row < num_rows; row++) {
        free(matrix[row]);
    }
    free(matrix);
}

void free_memory_randomized_matrices(short **randomized_matrices[], int num_rows, int num_matrices)
{
    for (int pos = 0; pos < num_matrices; pos++) {
        free_memory_shorts_matrix(randomized_matrices[pos], num_rows);
    }
}

double** allocate_memory_doubles_matrix(int num_rows, int num_cols)
{
    double **matrix = (double **) malloc(num_rows * sizeof(double *));

    if (matrix != NULL) {
        for (int row = 0; row < num_rows; row++) {
            matrix[row] = (double *) malloc(num_cols * sizeof(double));

            if (matrix[row] == NULL) {
                free_memory_doubles_matrix(matrix, row);

                printf("Failed to allocate dynamic memory to the matrix row %i.\n", row);
                return NULL;
            }
        }
    }
    else {
        printf("Failed to allocate the dynamic memory to the matrix.\n");
        return NULL;
    }
    return matrix;
}

int** allocate_memory_integers_matrix(int num_rows, int num_cols) {
    int **matrix = (int **) malloc(num_rows * sizeof(int *));

    if (matrix != NULL) {
        for (int row = 0; row < num_rows; row++) {
            matrix[row] = (int *) malloc(num_cols * sizeof(int));

            if (matrix[row] == NULL) {
                free_memory_integers_matrix(matrix, row);

                printf("Failed to allocate the dynamic memory to the matrix row %i\n", row);
                return NULL;
            }
        }
    }
    else {
        printf("Failed to allocate the dynamic memory to the matrix.\n");
        return NULL;
    }
    return matrix;
}

short** allocate_memory_shorts_matrix(int num_rows, int num_cols)
{
    short **matrix = (short **) malloc(num_rows * sizeof(short *));

    if (matrix != NULL) {
        for (int row = 0; row < num_rows; row++) {
            matrix[row] = (short *) malloc(num_cols * sizeof(short));

            if (matrix[row] == NULL) {
                free_memory_shorts_matrix(matrix, row);

                printf("Failed to allocate the dynamic memory to the matrix row %i\n", row);
                return NULL;
            }
        }
    }
    else {
        printf("Failed to allocate the dynamic memory to the matrix.\n");
        return NULL;
    }
    return matrix;
}

void allocate_memory_randomized_matrices(short **randomized_matrices[], int num_rows, int num_cols, int num_matrices)
{
    for (int pos = 0; pos < num_matrices; pos++) {
        randomized_matrices[pos] = allocate_memory_shorts_matrix(num_rows, num_cols);
    }
}

void initialize_integers_matrix_zeros(int **matrix, int num_rows, int num_cols)
{
    for (int row = 0; row < num_rows; row++) {
        memset(matrix[row], 0, num_cols * sizeof(int));
    }
}

void initialize_matrix_shorts_zeros(short **matrix, int num_rows, int num_cols)
{
    for (int row = 0; row < num_rows; row++) {
        memset(matrix[row], 0, num_cols * sizeof(short));
    }
}

short** transpose_matrix(short **matrix, int num_rows, int num_cols)
{
    short **transposed_matrix = allocate_memory_shorts_matrix(num_cols, num_rows);

    for (int row = 0; row < num_rows; row++) {
        for (int col = 0; col < num_cols; col++) {
            transposed_matrix[col][row] = matrix[row][col];
        }
    }
    return transposed_matrix;
}

void create_relative_abundances(int **matrix_absolute_abundances, double **matrix_relative_abundances,
                                int num_rows, int num_cols, int num_bacterial_species[])
{
    for (int row = 0; row < num_rows; row++) {
        for (int col = 0; col < num_cols; col++) {
            if (num_bacterial_species[row] != 0) {
                matrix_relative_abundances[row][col] = (double) matrix_absolute_abundances[row][col]
                        / (double) num_bacterial_species[row];
            }
            else {
                matrix_relative_abundances[row][col] = 0;
            }
        }
    }
}

void create_matrix_individuals(char *vertebrates, double **matrix_individuals)
{
    FILE *f_vertebrates = fopen(vertebrates, "r");
    char line[10000], absolute_abundances_individual[10000], *absolute_abundance;
    int **matrix_absolute_abundances, num_bacterial_species_per_individual[NUM_INDIVIDUALS] = {0};
    int row, col = 0;

    if (f_vertebrates == NULL) {
        printf("Error in opening the file %s.", vertebrates);
    }
    else {
        matrix_absolute_abundances = allocate_memory_integers_matrix(NUM_INDIVIDUALS, NUM_BACTERIAL_GENUS);

        fgets(line, sizeof(line), f_vertebrates);       /* Removed from first row. */
        while (fgets(line, sizeof(line), f_vertebrates) != NULL) {
            sscanf(line, "%*s %[^\n]", absolute_abundances_individual);     /* Removed from first column. */

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
        create_relative_abundances(matrix_absolute_abundances,matrix_individuals,NUM_INDIVIDUALS,
                                   NUM_BACTERIAL_GENUS, num_bacterial_species_per_individual);

        free_memory_integers_matrix(matrix_absolute_abundances, NUM_INDIVIDUALS);

        fclose(f_vertebrates);
    }
}

void get_individuals(FILE *f_vertebrates, Individual *individuals, int num_individuals)
{
    char line[10000], *individual;
    int pos = 0;

    fgets(line, sizeof(line), f_vertebrates);
    individual = strtok(line, "\" \"");
    while (individual != NULL) {
        strcpy(individuals[pos++].code, individual);
        individual = strtok(NULL, "\" \"");
    }
}

int offset(char *sample_type)
{
    if (strcmp(sample_type, "Wild") == 0) {
        return 0;
    }
    else {
        return 1;
    }
}

void get_species_sample_types(char *metadata, Individual *individuals, int num_individuals)
{
    FILE *f_metadata = fopen(metadata, "r");
    char *sample, line[500];
    int pos;
    bool found;

    if (f_metadata == NULL) {
        printf("Error in opening the file %s.", metadata);
    }
    else {
        fgets(line, sizeof(line), f_metadata);          /* Removed from first row. */
        while (fgets(line, sizeof(line), f_metadata) != NULL) {
            strtok(line, ";");
            sample = strtok(NULL, ";");
            found = false;
            pos = 0;
            while ((!found) && (pos < num_individuals)) {
                if (strcmp(sample, individuals[pos].code) == 0) {
                    strtok(NULL, ";");
                    strcpy(individuals[pos].vertebrate, strtok(NULL, ";"));
                    strtok(NULL, ";");
                    individuals[pos].sample_type = offset(strtok(NULL, ";"));
                    found = true;
                }
                else {
                    pos++;
                }
            }
        }
        fclose(f_metadata);
    }
}

void create_matrix_vertebrates(char *vertebrates, char *metadata, double **matrix_vertebrates)
{
    FILE *f_vertebrates = fopen(vertebrates, "r");
    Individual *individuals = (Individual *) malloc(NUM_INDIVIDUALS * sizeof(Individual));
    char *absolute_abundance, absolute_abundances_genus[10000], line[10000];
    int **matrix_absolute_abundances, num_bacterial_species_per_vertebrate[NUM_VERTEBRATES] = {0};
    int pos_ant, pos, row, col = 0;

    if (f_vertebrates == NULL) {
        printf("Error in opening the file %s.", vertebrates);
    }
    else {
        get_individuals(f_vertebrates, individuals, NUM_INDIVIDUALS);
        get_species_sample_types(metadata, individuals, NUM_INDIVIDUALS);

        matrix_absolute_abundances = allocate_memory_integers_matrix(NUM_VERTEBRATES, NUM_BACTERIAL_GENUS);
        initialize_integers_matrix_zeros(matrix_absolute_abundances, NUM_VERTEBRATES, NUM_BACTERIAL_GENUS);

        while (fgets(line, sizeof(line), f_vertebrates) != NULL) {
            sscanf(line, "%*s %[^\n]", absolute_abundances_genus);     /* Removed from first column. */

            absolute_abundance = strtok(absolute_abundances_genus, " ");
            pos_ant = pos = row = 0;
            while (absolute_abundance != NULL) {
                if (strcmp(individuals[pos_ant].vertebrate, individuals[pos].vertebrate) != 0) {
                    row += 2;
                    pos_ant = pos;
                }

                matrix_absolute_abundances[row + individuals[pos].sample_type][col] += atoi(absolute_abundance);
                num_bacterial_species_per_vertebrate[row + individuals[pos++].sample_type] += atoi(absolute_abundance);

                absolute_abundance = strtok(NULL, " ");
            }
            col++;
        }
        free(individuals);
        create_relative_abundances(matrix_absolute_abundances,matrix_vertebrates, NUM_VERTEBRATES,
                                   NUM_BACTERIAL_GENUS, num_bacterial_species_per_vertebrate);

        free_memory_integers_matrix(matrix_absolute_abundances, NUM_VERTEBRATES);

        fclose(f_vertebrates);
    }
}

void create_matrix(char *name_matrix, char *vertebrates, char *metadata, double **abundances_matrix)
{
    if (strcmp(name_matrix, "individuals") == 0) {
        create_matrix_individuals(vertebrates, abundances_matrix);
    }
    else if (strcmp(name_matrix, "vertebrates") == 0) {
        create_matrix_vertebrates(vertebrates, metadata, abundances_matrix);
    }
}

void discretize_matrix(double **matrix, short **binary_matrix, int num_rows, int num_cols, double threshold)
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

double calculate_nested_value(int **matrix, int num_rows, int num_cols)
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

double calculate_nested_value_optimized(short **matrix, int num_rows, int num_cols)
{
    short **transposed_matrix = transpose_matrix(matrix, num_rows, num_cols);
    int *sum_rows = (int *) calloc(num_rows, sizeof(int));
    int *sum_cols = (int *) calloc(num_cols, sizeof(int));
    int first_isocline, second_isocline, third_isocline, fourth_isocline, row, col, first_col, second_col, total_cols;

        /* Calculate and save the number of interactions of every row and
           calculate and save the number of interactions of every column. */
    for (row = 0; row < num_rows; row++) {
        for (col = 0; col < num_cols; col++) {
            sum_rows[row] += matrix[row][col];
            sum_cols[col] += transposed_matrix[col][row];
        }
    }

    first_isocline = second_isocline = third_isocline = fourth_isocline = 0;
    total_cols = num_cols - 1;

        /* Calculate the sum of the number of shared interactions between rows
           and the sum of the minimum of pairs of interactions of rows. */
    for (int first_row = 0; first_row < num_rows - 1; first_row++) {
        for (int second_row = 0; second_row < num_rows; second_row++) {
            if (first_row < second_row) {
                for (col = 0; col < num_cols; col++) {
                    first_isocline += matrix[first_row][col] & matrix[second_row][col];
                }
                if (sum_rows[first_row] < sum_rows[second_row]) {
                    third_isocline += sum_rows[first_row];
                } else {
                    third_isocline += sum_rows[second_row];
                }
            }
        }
    }

        /* Calculate the sum of the number of shared interactions between columns
           and the sum of the minimum of pairs of the number of interactions of columns. */
    // #pragma omp parallel for private(first_col, second_col, row) shared(total_cols, num_cols, num_rows, transposed_matrix, sum_cols) reduction(+:second_isocline, fourth_isocline) default(none) schedule(dynamic)
    for (first_col = 0; first_col < total_cols; first_col++) {
        for (second_col = 0; second_col < num_cols; second_col++) {
            if (first_col < second_col) {
                for (row = 0; row < num_rows; row++) {
                    second_isocline += transposed_matrix[first_col][row] & transposed_matrix[second_col][row];
                }
                if (sum_cols[first_col] < sum_cols[second_col]) {
                    fourth_isocline += sum_cols[first_col];
                } else {
                    fourth_isocline += sum_cols[second_col];
                }
            }
        }
    }

    free_memory_shorts_matrix(transposed_matrix, num_cols);
    free(sum_rows);
    free(sum_cols);

        /* Calculate and return the nested value of the matrix. */
    return ((double)(first_isocline + second_isocline) / (double)(third_isocline + fourth_isocline));
}

int count_ones_binary_matrix(short **matrix, int num_rows, int num_cols)
{
    int num_ones = 0;

    for (int row = 0; row < num_rows; row++) {
        for (int col = 0; col < num_cols; col++) {
            num_ones += matrix[row][col];
        }
    }
    return num_ones;
}

void generate_randomized_matrix(short **randomized_matrix, int num_rows, int num_cols, int num_ones)
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

void generate_nested_values_randomized(short **matrix, int num_rows, int num_cols, int num_randomized_matrices,
                                       double nested_values_randomized[])
{
    short **randomized_matrices[omp_get_max_threads()];
    int pos, num_ones = count_ones_binary_matrix(matrix, num_rows, num_cols);

    allocate_memory_randomized_matrices(randomized_matrices, num_rows, num_cols, omp_get_max_threads());

    #pragma omp parallel for private(pos) shared(num_randomized_matrices, randomized_matrices, num_rows, num_cols, num_ones, nested_values_randomized) default(none) schedule(dynamic)
    for (pos = 0; pos < num_randomized_matrices; pos++) {
        initialize_matrix_shorts_zeros(randomized_matrices[omp_get_thread_num()], num_rows, num_cols);
        generate_randomized_matrix(randomized_matrices[omp_get_thread_num()], num_rows, num_cols, num_ones);
        nested_values_randomized[pos] = calculate_nested_value_optimized(randomized_matrices[omp_get_thread_num()], num_rows, num_cols);
    }
    free_memory_randomized_matrices(randomized_matrices, num_rows, omp_get_max_threads());
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
    int half, first = 0, last = num_elements - 1;
    bool found = false;

    while ((! found) && (first <= last)) {
        half = (first + last) / 2;

        if (nested_values[half] == nested_value) {
            found = true;       /* Found the value at the previous calculated index. */
        }
        else if (nested_values[half] < nested_value) {
            first = half + 1;       /* Search in the right half. */
        }
        else {
            last = half - 1;        /* Search in the left half. */
        }
    }
    if (found) {
        return half;
    }
    else {
        return -1;      /* Value not found. */
    }
}

Nested_elements nested_test(short **matrix, int num_rows, int num_cols, int num_randomized_matrices)
{
    Nested_elements nested_elements;
    double nested_values[num_randomized_matrices + 1];

        /* Generate as many randomized matrices from the real matrix as it is specified and calculate their nested values. */
    generate_nested_values_randomized(matrix, num_rows, num_cols, num_randomized_matrices,
                                      nested_values);

        /* Calculate and store the nested value of the real matrix. */
    nested_elements.nested_value = calculate_nested_value_optimized(matrix, num_rows, num_cols);
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
    double **abundances_matrix;
    short **binary_matrix;
    int num_rows, num_cols;
    Nested_elements nested_elements;
    // double nested_value;

    omp_set_num_threads(atoi(argv[5]));
    srand(time(NULL));

    select_matrix(argv[3], &num_rows, &num_cols);

    abundances_matrix = allocate_memory_doubles_matrix(num_rows, num_cols);

    create_matrix(argv[3], argv[1], argv[2], abundances_matrix);

    binary_matrix = allocate_memory_shorts_matrix(num_rows, num_cols);

    discretize_matrix(abundances_matrix, binary_matrix, num_rows, num_cols, atof(argv[4]));

    free_memory_doubles_matrix(abundances_matrix, num_rows);

    // nested_value = calculate_nested_value_optimized(binary_matrix, num_rows, num_cols);
    // printf("Nested value: %f\n", nested_value);

    nested_elements = nested_test(binary_matrix, num_rows, num_cols, 1000);
    printf("\nNested value: %f\nP-value: %f\n", nested_elements.nested_value, nested_elements.p_value);

    free_memory_shorts_matrix(binary_matrix, num_rows);

    return 0;
}
