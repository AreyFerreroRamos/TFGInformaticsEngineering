import pandas as pd
import numpy as np
import sys


def get_num_species_per_individual(individual, df_vertebrates):
    num_bacterial_species_per_individual = 0
    for num_bacterial_species_per_genus in df_vertebrates[individual]:
        num_bacterial_species_per_individual += num_bacterial_species_per_genus

    return num_bacterial_species_per_individual


def abundances_individuals_matrix(df_vertebrates):
    individuals_matrix = np.empty((df_vertebrates.shape[1], df_vertebrates.shape[0]))

    num_individuals = 0
    for individual in df_vertebrates:
        num_bacterial_species_per_individual = get_num_species_per_individual(individual, df_vertebrates)

        column_genus = 0
        for num_bacterial_species_per_genus in df_vertebrates[individual]:
            relative_abundance = num_bacterial_species_per_genus / num_bacterial_species_per_individual
            individuals_matrix[num_individuals][column_genus] = relative_abundance
            column_genus += 1

        num_individuals += 1

    return individuals_matrix


def get_specie_sample_type(individual, df_metadata):
    specie = sample_type = ""

    row = 1
    for sample in df_metadata[df_metadata.columns[0]]:
        if sample == individual:
            specie = df_metadata.loc[row, df_metadata.columns[2]]
            sample_type = df_metadata.loc[row, df_metadata.columns[4]]
        else:
            row += 1

    return specie, sample_type


def normalize_matrix(matrix, rows, columns, num_individuals):
    if num_individuals != 0:
        matrix[rows][columns] = matrix[rows][columns] / num_individuals
    else:
        matrix[rows][columns] = 0


def normalize_matrix_vertebrates(matrix_vertebrate_genus, num_specie, num_wild, num_captivity):
    num_genus = 0
    while num_genus < matrix_vertebrate_genus.shape[1]:
        normalize_matrix(matrix_vertebrate_genus, num_specie, num_genus, num_wild)
        normalize_matrix(matrix_vertebrate_genus, num_specie + 1, num_genus, num_captivity)
        num_genus += 1


def offset(sample_type):
    if sample_type == 'Wild':
        return 0
    else:
        return 1


def abundances_vertebrates_matrix(df_vertebrates, df_metadata):
    vertebrates_matrix = np.zeros((2, df_vertebrates.shape[0]))
    specie_ant, _ = get_specie_sample_type(df_vertebrates.columns[0], df_metadata)

    num_species = num_wild = num_captivity = 0
    for individual in df_vertebrates:
        specie, sample_type = get_specie_sample_type(individual, df_metadata)

        if specie_ant != specie:
            normalize_matrix_vertebrates(vertebrates_matrix, num_species, num_wild, num_captivity)
            num_species += 2
            num_wild = num_captivity = 0
            vertebrates_matrix.resize((vertebrates_matrix.shape[0] + 2, vertebrates_matrix.shape[1]), refcheck=False)
            specie_ant = specie

        num_bacterial_species_per_individual = get_num_species_per_individual(individual, df_vertebrates)

        column_genus = 0
        for num_bacterial_species_per_genus in df_vertebrates[individual]:
            relative_abundance = num_bacterial_species_per_genus / num_bacterial_species_per_individual
            vertebrates_matrix[num_species + offset(sample_type)][column_genus] += relative_abundance
            column_genus += 1

        if sample_type == "Wild":
            num_wild += 1
        else:
            num_captivity += 1

    normalize_matrix_vertebrates(vertebrates_matrix, num_species, num_wild, num_captivity)
    return vertebrates_matrix


def get_code_specie(name_specie, name_file_codes_vertebrates):
    f_codes_vertebrates = open(name_file_codes_vertebrates, 'r')
    code_specie = ""

    for vertebrate_specie in f_codes_vertebrates:
        if name_specie.replace(' ', '_', 1) == vertebrate_specie.split()[1]:
            code_specie = vertebrate_specie.split()[0]

    f_codes_vertebrates.close()
    return code_specie


def abundances_specie_matrix(vertebrate_specie, df_vertebrates, df_metadata):
    specie_matrix = np.empty((0, df_vertebrates.shape[0]))

    num_individuals = 0
    for individual in df_vertebrates:
        specie, _ = get_specie_sample_type(individual, df_metadata)

        if specie == vertebrate_specie:
            specie_matrix.resize((specie_matrix.shape[0] + 1, specie_matrix.shape[1]))

            num_bacterial_species_per_individual = get_num_species_per_individual(individual, df_vertebrates)

            column_genus = 0
            for num_bacterial_species_per_genus in df_vertebrates[individual]:
                relative_abundance = num_bacterial_species_per_genus / num_bacterial_species_per_individual
                specie_matrix[num_individuals][column_genus] = relative_abundance
                column_genus += 1

            num_individuals += 1

    return specie_matrix


def abundances_specie_sample_type_matrix(vertebrate_specie, vertebrate_sample_type, df_vertebrates, df_metadata):
    specie_sample_type_matrix = np.empty((0, df_vertebrates.shape[0]))

    num_individuals = 0
    for individual in df_vertebrates:
        specie, sample_type = get_specie_sample_type(individual, df_metadata)

        if specie == vertebrate_specie and sample_type == vertebrate_sample_type:
            specie_sample_type_matrix.resize((specie_sample_type_matrix.shape[0] + 1,
                                              specie_sample_type_matrix.shape[1]))

            num_bacterial_species_per_individual = get_num_species_per_individual(individual, df_vertebrates)

            column_genus = 0
            for num_bacterial_species_per_genus in df_vertebrates[individual]:
                relative_abundance = num_bacterial_species_per_genus / num_bacterial_species_per_individual
                specie_sample_type_matrix[num_individuals][column_genus] = relative_abundance
                column_genus += 1

            num_individuals += 1

    return specie_sample_type_matrix


def discretize_matrix(matrix, threshold):
    rows = 0
    while rows < matrix.shape[0]:
        columns = 0
        while columns < matrix.shape[1]:
            if matrix[rows][columns] > threshold:
                matrix[rows][columns] = 1
            else:
                matrix[rows][columns] = 0
            columns += 1
        rows += 1


def nestedness(matrix):
    first_isocline = 0
    second_isocline = 0
    third_isocline = 0
    fourth_isocline = 0

    # Calculate the sum of the number of shared interactions between rows.
    for first_row in range(matrix.shape[0]):
        for second_row in range(matrix.shape[0]):
            if first_row < second_row:
                for col in range(matrix.shape[1]):
                    if matrix[first_row][col] == 1 and matrix[second_row][col] == 1:
                        first_isocline += 1

    # Calculate the sum of the number of shared interactions between columns.
    for first_col in range(matrix.shape[1]):
        for second_col in range(matrix.shape[1]):
            if first_col < second_col:
                for row in range(matrix.shape[0]):
                    if matrix[row][first_col] == 1 and matrix[row][second_col] == 1:
                        second_isocline += 1

    # Calculate the sum of the number of interactions of rows.
    for first_row in range(matrix.shape[0]):
        for second_row in range(matrix.shape[0]):
            if first_row < second_row:
                first_acum = second_acum = 0
                for col in range(matrix.shape[1]):
                    first_acum += matrix[first_row][col]
                    second_acum += matrix[second_row][col]
                third_isocline += min(first_acum, second_acum)

    # Calculate the sum of the number of interactions of columns.
    for first_col in range(matrix.shape[1]):
        for second_col in range(matrix.shape[1]):
            if first_col < second_col:
                first_acum = second_acum = 0
                for row in range(matrix.shape[0]):
                    first_acum += matrix[row][first_col]
                    second_acum += matrix[row][second_col]
                fourth_isocline += min(first_acum, second_acum)

    # Calculate and return the nestedness of the matrix.
    return (first_isocline + second_isocline) / (third_isocline + fourth_isocline)


def nestedness_optimized(matrix):
    sum_rows = []
    sum_cols = []

    # Calculate and save the number of interactions of every row.
    for row in range(matrix.shape[0]):
        sum_rows.append(sum(matrix[row, :]))

    # Calculate and save the number of interactions of every column.
    for col in range(matrix.shape[1]):
        sum_cols.append(sum(matrix[:, col]))

    first_isocline = second_isocline = third_isocline = fourth_isocline = 0

    # Calculate the sum of the number of shared interactions between rows
    # and the sum of the minimum of pairs of interactions of rows.
    for first_row in range(matrix.shape[0] - 1):
        for second_row in range(first_row + 1, matrix.shape[0]):
            first_isocline += sum([first * second for first, second in zip(matrix[first_row, :], matrix[second_row, :])])
            third_isocline += min(sum_rows[first_row], sum_rows[second_row])

    # Calculate the sum of the number of shared interactions between columns
    # and the sum of the minimum of pairs of the number of interactions of columns.
    for first_col in range(matrix.shape[1] - 1):
        for second_col in range(first_col + 1, matrix.shape[1]):
            second_isocline += sum([first * second for first, second in zip(matrix[:, first_col], matrix[:, second_col])])
            fourth_isocline += min(sum_cols[first_col], sum_cols[second_col])

    # Calculate and return the nestedness of the matrix.
    return (first_isocline + second_isocline) / (third_isocline + fourth_isocline)


def count_ones_binary_matrix(matrix):
    num_ones = 0

    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            if matrix[row][column] == 1:
                num_ones += 1

    return num_ones


def generate_nested_values_randomized(matrix, num_randomized_matrices):
    nested_values_randomized = []
    num_ones = count_ones_binary_matrix(matrix)

    for i in range(num_randomized_matrices):
        randomized_matrix = np.zeros((matrix.shape[0], matrix.shape[1]), dtype=int)
        randomized_matrix.ravel()[np.random.choice(matrix.shape[0] * matrix.shape[1], num_ones, replace=False)] = 1
        nested_values_randomized.append(nestedness(randomized_matrix))

    return nested_values_randomized


def nestedness_assessment(matrix, num_randomized_matrices):
    if matrix.size == 0:
        return 0, 0

    # Generate as many randomized matrices from the real matrix as it is specified by parameter
    # and calculate their nestedness value.
    nested_values = generate_nested_values_randomized(matrix, num_randomized_matrices)

    # Calculate the nestedness value of the real matrix.
    nested_value = nestedness(matrix)
    nested_values.append(nested_value)

    # Sort the list of nestedness values.
    nested_values.sort()

    # Calculate the fraction of randomized matrices that have a nestedness value greater than that of the real matrix.
    p_value = (num_randomized_matrices - nested_values.index(nested_value)) / (num_randomized_matrices + 1)

    return nested_value, p_value


df_vertebrates = pd.read_table(sys.argv[1], delimiter=' ', header=0)
df_metadata = pd.read_table(sys.argv[2], delimiter=';', header=0)

if sys.argv[4] == "individuals":
    abundances_matrix = abundances_individuals_matrix(df_vertebrates)

elif sys.argv[4] == "vertebrates":
    abundances_matrix = abundances_vertebrates_matrix(df_vertebrates, df_metadata)

elif len(sys.argv[4].split()) >= 2:
    if len(sys.argv[4].split()) == 2:
        abundances_matrix = abundances_specie_matrix(get_code_specie(sys.argv[4], sys.argv[3]),
                                                       df_vertebrates, df_metadata)
    elif sys.argv[4].split()[2] == "Wild":
        abundances_matrix = abundances_specie_sample_type_matrix(get_code_specie(
            sys.argv[4].split()[0] + ' ' + sys.argv[4].split()[1], sys.argv[3]), 'Wild', df_vertebrates, df_metadata)
    elif sys.argv[4].split()[2] == "Captive":
        abundances_matrix = abundances_specie_sample_type_matrix(get_code_specie(
            sys.argv[4].split()[0] + ' ' + sys.argv[4].split()[1], sys.argv[3]), 'Captivity', df_vertebrates, df_metadata)

discretize_matrix(abundances_matrix, 0.0001)
nested_abundances_matrix, p_value = nestedness_assessment(abundances_matrix, int(sys.argv[5]))
print(nested_abundances_matrix, p_value)
