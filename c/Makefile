nestedness_test: nestedness_test.c
	gcc -O3 nestedness_test.c -o nestedness_test

omp: nestedness_test_omp.c
	gcc -O3 -fopenmp nestedness_test_omp.c -o nestedness_test_omp

mpi: nestedness_test_mpi.c
	mpicc -O3 nestedness_test_mpi.c -o nestedness_test_mpi

mpi_omp: nestedness_test_mpi_omp.c
	mpicc -O3 -fopenmp nestedness_test_mpi_omp.c -o nestedness_test_mpi_omp

run: nestedness_test
	time ./nestedness_test ../input_files/count_Genus_all.tsv ../input_files/metadata.csv vertebrates 0.0001 1000

run_omp: nestedness_test_omp
	time ./nestedness_test_omp ../input_files/count_Genus_all.tsv ../input_files/metadata.csv vertebrates 0.0001 1000 5

run_mpi: nestedness_test_mpi
	time mpirun -np 5 ./nestedness_test_mpi ../input_files/count_Genus_all.tsv ../input_files/metadata.csv vertebrates 0.0001 1000

run_mpi_omp: nestedness_test_mpi_omp
	time mpirun -np 5 ./nestedness_test_mpi_omp ../input_files/count_Genus_all.tsv ../input_files/metadata.csv vertebrates 0.0001 1000 5

clean:
	rm -f nestedness_test
	rm -f nestedness_test_omp
	rm -f nestedness_test_mpi
	rm -f nestedness_test_mpi_omp
