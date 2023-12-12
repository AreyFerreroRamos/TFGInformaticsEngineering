nestedness_test: nestedness_test.c
	gcc -O3 nestedness_test.c -o nestedness_test

omp: nestedness_test_omp.c
	gcc -O3 -fopenmp nestedness_test_omp.c -o nestedness_test_parallel

mpi: nestedness_test_mpi.c
	gcc -O3 nestedness_test_mpi.c -o nestedness_test_parallel

mpi_omp: nestedness_test_mpi_omp.c
	gcc -O3 -fopenmp nestedness_test_mpi_omp.c -o nestedness_test_parallel

run: nestedness_test
	./nestedness_test ../input_files/count_Genus_all.tsv ../input_files/metadata.csv vertebrates 0.0001

parallel: nestedness_test_parallel
	./nestedness_test_parallel ../input_files/count_Genus_all.tsv ../input_files/metadata.csv vertebrates 0.0001 5

clean:
	rm -f nestedness_test
	rm -f nestedness_test_parallel
