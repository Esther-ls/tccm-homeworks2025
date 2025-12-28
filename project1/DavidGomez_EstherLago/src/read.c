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

double Nuclear_repulsion_energy(trexio_t *  trexio_file){
	double Enn;
	trexio_exit_code rc;
	rc = trexio_read_nucleus_repulsion(trexio_file, &Enn);
        if (rc != TREXIO_SUCCESS) {
                printf("TREXIO Error reading nuclear repulsion energy: \n%s\n", trexio_string_of_error(rc));
                exit(1);
        }
	return Enn;
}

int Occupied_orbitals(trexio_t * trexio_file){
	int32_t n_up;
	trexio_exit_code rc;
	rc = trexio_read_electron_up_num(trexio_file, &n_up);
	if (rc != TREXIO_SUCCESS) {
                printf("TREXIO Error reading occupied orbitals: \n%s\n", trexio_string_of_error(rc));
                exit(1);
        }
	return n_up;
}



int Molecular_orbitals(trexio_t * trexio_file){
	int32_t mo_num;
	trexio_exit_code rc;
	rc = trexio_read_mo_num(trexio_file, &mo_num);
	if (rc != TREXIO_SUCCESS) {
                printf("TREXIO Error reading molecular orbitals: \n%s\n", trexio_string_of_error(rc));
                exit(1);
        }
	return mo_num;
}

double * one_e_integrals(trexio_t * trexio_file, int mo_num){
	double * data = malloc(mo_num * mo_num * sizeof(double));
	trexio_exit_code rc;
	rc = trexio_read_mo_1e_int_core_hamiltonian(trexio_file, data);
	if (rc != TREXIO_SUCCESS) {
                printf("TREXIO Error reading one electron integrals: \n%s\n", trexio_string_of_error(rc));
                exit(1);
        }
	return data;

}

double * mo_energies(trexio_t * trexio_file, int mo_num){
	double * mo_energy = malloc(mo_num*sizeof(double));
	trexio_exit_code rc;
	rc = trexio_read_mo_energy(trexio_file, mo_energy);
        if (rc != TREXIO_SUCCESS) {
                printf("TREXIO error reading molecular orbital energies: %s\n", trexio_string_of_error(rc));
                exit(1);
        }
	return mo_energy;
}


int number_2_e_integrals(trexio_t * trexio_file){
	int64_t n_integrals;
	trexio_exit_code rc;
	rc = trexio_read_mo_2e_int_eri_size(trexio_file, &n_integrals);
	if (rc != TREXIO_SUCCESS) {
                printf("TREXIO error reading the number of two electron integrals: %s\n", trexio_string_of_error(rc));
                exit(1);
        }
	return n_integrals;
}

void two_e_integrals(trexio_t * trexio_file, int n_integrals, int32_t** index_out, double** value_out){
	int64_t offset_file  = 0;
	int64_t buffer_size  = n_integrals;
	int32_t* index = malloc(4 * n_integrals * sizeof(int32_t));
        double* value = malloc(n_integrals * sizeof(double));
	trexio_exit_code rc;
	rc = trexio_read_mo_2e_int_eri(trexio_file, offset_file, &buffer_size, index, value);
	if (rc != TREXIO_SUCCESS) {
                printf("TREXIO error reading two electron integrals: %s\n", trexio_string_of_error(rc));
                exit(1);
        }
	*index_out=index;
	*value_out=value;
}


void close_fun(trexio_t * trexio_file){ 
 trexio_exit_code rc; // rc variable: save the exit code
 rc = trexio_close(trexio_file);
 if (rc != TREXIO_SUCCESS) {
	 printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
	 exit(1);
 }
 trexio_file = NULL;
}





