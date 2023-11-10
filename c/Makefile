nestedness_test: nestedness_test.c
	gcc -pg nestedness_test.c -o nestedness_test

run: nestedness_test
	./nestedness_test ../input_files/count_Genus_all.tsv ../input_files/metadata.csv individuals 0.0001

clean:
	rm -f nestedness_test
