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
	rc = trexio_read_nucleus_repulsion(trexio_file, &Enn);
	if (rc != TREXIO_SUCCESS) {
		printf("TREXIO Error reading nuclear repulsion energy: \n%s\n", trexio_string_of_error(rc));
		exit(1);
	}
	rc = trexio_read_electron_up_num(trexio_file, &n_up);
	rc = trexio_read_mo_num(trexio_file, &mo_num);
	double * data = malloc(mo_num * mo_num * sizeof(double));
	rc = trexio_read_mo_1e_int_core_hamiltonian(trexio_file, data);
	close_fun(trexio_file);
	printf("%lf\n", Enn);
	printf("%i\n", n_up);
	printf("%i\n", mo_num);
	printf("%lf\n", data[0]);
    	for (int i=0 ; i<mo_num*mo_num ; i+=25) { // 25 pq solo elementos diagonal ppal, y es matriz 24x24 
        printf("data[%d]       = %lf\n",i, data[i]);
    }
	return 0;
}
