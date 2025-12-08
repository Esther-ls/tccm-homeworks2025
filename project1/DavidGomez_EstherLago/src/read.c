#include <trexio.h>
#include <stdio.h>
#include <stdlib.h>

trexio_t* open_fun(const char* filename){ 
 trexio_exit_code rc; // rc variable: save the exit code
 trexio_t* trexio_file = trexio_open(filename, 'r', TREXIO_AUTO, &rc); 
 /* trexio_file: open a trexio file
    filename: trexio filename
    r: read only
    TREXIO_AUTO: detects the extension
    &rc: address of the exit code
  */
 if (rc != TREXIO_SUCCESS) {
	 printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
	 exit(1);
 }
 return trexio_file;
}


trexio_t* close_fun(trexio_t * trexio_file){ 
 trexio_exit_code rc; // rc variable: save the exit code
 rc = trexio_close(trexio_file);
 if (rc != TREXIO_SUCCESS) {
	 printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
	 exit(1);
 }
 trexio_file = NULL;
}


/*trexio_t* NN_E(trexio_t* trexio_file){
	trexio_exit_code trexio_read_nucleus_repulsion(trexio_t* const trexio_file, double* const energy);
	return trexio_exit_code
}*/
