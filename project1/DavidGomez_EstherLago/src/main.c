#include <stdio.h>
#include <trexio.h>
#include <stdlib.h>

trexio_t * open_fun(const char * filename);

trexio_t* close_fun(trexio_t * trexio_file); 

double main() {
	trexio_exit_code rc;
	double Enn;
	trexio_t* file = open_fun("h2o.h5");
	rc = trexio_read_nucleus_repulsion(file, &Enn);
	if (rc != TREXIO_SUCCESS) {
		printf("TREXIO Error reading nuclear repulsion energy: \n%s\n", trexio_string_of_error(rc));
		exit(1);
	}
	close_fun(file);
	printf("%lf\n", Enn);
}
