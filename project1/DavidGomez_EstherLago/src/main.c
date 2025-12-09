#include <stdio.h>
#include <trexio.h>
#include <stdlib.h>

trexio_t * open_fun(const char * filename);

void close_fun(trexio_t * trexio_file); 

double main() {
	trexio_t* trexio_file = open_fun("h2o.h5");
	trexio_exit_code rc;
	double Enn;
	int32_t n_up;
	int32_t mo_num;
	int64_t n_integrals;
	int64_t offset_file  = 0;              // To start in the first two-electron integral
	
	/* Read and write nuclear repulsion energy */
	rc = trexio_read_nucleus_repulsion(trexio_file, &Enn);
	if (rc != TREXIO_SUCCESS) {
		printf("TREXIO Error reading nuclear repulsion energy: \n%s\n", trexio_string_of_error(rc));
		exit(1);
	}
	printf("Nuclear Repulsion Energy (au)	= %lf\n", Enn);

	/*Read and write the number of occupied orbitals */
	rc = trexio_read_electron_up_num(trexio_file, &n_up);
	printf("Number of Occupied Orbitals	= %i\n", n_up);

	/*Read and write One-Electron Integrals*/
	rc = trexio_read_mo_num(trexio_file, &mo_num); //Number of MO
	printf("Number of Molecular Orbitals	= %i\n", mo_num);
	double * data = malloc(mo_num * mo_num * sizeof(double));
	rc = trexio_read_mo_1e_int_core_hamiltonian(trexio_file, data);
    	printf("\n\tOne-Electron Integrals (au)\n");
	for (int i=0 ; i<mo_num*mo_num ; i+=mo_num+1) { // 25 pq solo elementos diagonal ppal, y es matriz 24x24 
        printf("<i|h|i> (i=%d)       = %lf\n",i+1, data[i]);
    	}

	/*Read and write Two-Electron integrals*/
	rc = trexio_read_mo_2e_int_eri_size(trexio_file, &n_integrals);
	int64_t buffer_size  = n_integrals;    // Read all two-electron integrals
	int32_t* const index = malloc(4 * n_integrals * sizeof(int32_t));
	if (index == NULL) {
		fprintf(stderr, "Malloc failed for index");
		exit(1);
	}
	double* const value = malloc(n_integrals * sizeof(double));
	if (value == NULL) {
		fprintf(stderr, "Malloc failed for value");
		exit(1);
	}
	rc = trexio_read_mo_2e_int_eri(trexio_file, offset_file, &buffer_size, index, value);
	close_fun(trexio_file);
	/* Pruebas pero sobra
	printf("%ld\n", n_integrals);
	printf("%ld\n", offset_file);
	printf("%ld\n", buffer_size);
	for (int64_t i=0 ; i<n_integrals ; i++) {
		printf("integral[%ld]	=%lf\n", i, value[i]);
	}
        for (int64_t i=0 ; i<100 ; i++) {
                printf("integral[%ld]   =%d\n", i, index[i]);
        }
	*/

	return 0;
}
