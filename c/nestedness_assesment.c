#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct
{
    double nested_value;
    double p_value;
} nested_components;

int abundances_individuals_matrix()
{
    return 0;
}

int abundances_vertebrates_matrix()
{
    return 0;
}

void discretize_matrix()
{

}

void nestedness_assesment()
{

}

int main(int argc, char * argv[])
{
    double **abundances_matrix;

    if (! strcmp(argv[4], "individuals"))
    {
        abundances_matrix = abundances_individuals_matrix();
    }
    else if (! strcmp(argv[4], "vertebrates"))
    {
        abundances_matrix = abundances_vertebrates_matrix()
    }

    discretize_matrix(abundances_matrix);
    nestedness_assesment(abundances_matrix, atoi(argv[5]))
    printf("\nnested value: %.2f\np-value: %.2f", nested_components.nested_value, nested_components.p_value)
}