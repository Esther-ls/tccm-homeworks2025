/**
 * @file functions.c
 * @brief Functions to read quantum chemical data using the TREXIO library.
 *
 * @author David Gómez and Esther Lago
 */
#include <trexio.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * @brief Opens a TREXIO file.
 *
 * @param filename String containing the path to the TREXIO file (HDF5 format).
 * @return trexio_t* Pointer to the TREXIO file structure. Exits on error.
 */
trexio_t* open_fun(const char* filename){ 
 trexio_exit_code rc; // rc variable: save the exit code
 trexio_t* trexio_file = trexio_open(filename, 'r', TREXIO_AUTO, &rc); 
 if (rc != TREXIO_SUCCESS) {
	 printf("TREXIO Error opening the TREXIO file: %s\n", trexio_string_of_error(rc));
	 exit(1);
 }
 return trexio_file;
}

/**
 * @brief Reads the nuclear repulsion energy.
 *
 * @param trexio_file Pointer to the open TREXIO file.
 * @param Enn Pointer to a double where the energy value will be stored.
 */
void Nuclear_repulsion_energy(trexio_t *  trexio_file, double* Enn){
	trexio_exit_code rc;
	rc = trexio_read_nucleus_repulsion(trexio_file, Enn);
        if (rc != TREXIO_SUCCESS) {
                printf("TREXIO Error reading nuclear repulsion energy: \n%s\n", trexio_string_of_error(rc));
                exit(1);
        }
}

/**
 * @brief Reads the number of occupied orbitals.
 *
 * @param trexio_file Pointer to the open TREXIO file.
 * @param n_up Pointer to an integer where the number of up-spin electrons will be stored.
 */
void Occupied_orbitals(trexio_t * trexio_file, int32_t* n_up){
	trexio_exit_code rc;
	rc = trexio_read_electron_up_num(trexio_file, n_up);
	if (rc != TREXIO_SUCCESS) {
                printf("TREXIO Error reading occupied orbitals: \n%s\n", trexio_string_of_error(rc));
                exit(1);
        }
}

/**
 * @brief Reads the total number of molecular orbitals.
 *
 * @param trexio_file Pointer to the open TREXIO file.
 * @param mo_num Pointer to an integer where the number of orbitals will be stored.
 */
void Molecular_orbitals(trexio_t * trexio_file, int32_t* mo_num){
	trexio_exit_code rc;
	rc = trexio_read_mo_num(trexio_file, mo_num);
	if (rc != TREXIO_SUCCESS) {
                printf("TREXIO Error reading molecular orbitals: \n%s\n", trexio_string_of_error(rc));
                exit(1);
        }
}

/**
 * @brief Allocates memory and reads the one-electron integrals.
 *
 * @param trexio_file Pointer to the open TREXIO file.
 * @param mo_num Number of molecular orbitals (occupied + virtuals).
 * @param data Pointer to a double pointer where the one-electron integrals will be stored.
 */
void * one_e_integrals(trexio_t * trexio_file, int mo_num, double** data){
	* data = malloc(mo_num * mo_num * sizeof(double));
	trexio_exit_code rc;
	rc = trexio_read_mo_1e_int_core_hamiltonian(trexio_file,* data);
	if (rc != TREXIO_SUCCESS) {
                printf("TREXIO Error reading one electron integrals: \n%s\n", trexio_string_of_error(rc));
                exit(1);
        }
}

/**
 * @brief Allocates memory and reads the molecular orbital energies.
 *
 * @param trexio_file Pointer to the open TREXIO file.
 * @param mo_num Number of molecular orbitals (used for allocation size).
 * @param mo_energy Pointer to a double pointer where the molecular orbital energies will be stored.
 */
void * mo_energies(trexio_t * trexio_file, int mo_num, double** mo_energy){
	* mo_energy = malloc(mo_num*sizeof(double));
	trexio_exit_code rc;
	rc = trexio_read_mo_energy(trexio_file, *mo_energy);
        if (rc != TREXIO_SUCCESS) {
                printf("TREXIO error reading molecular orbital energies: %s\n", trexio_string_of_error(rc));
                exit(1);
        }
}

/**
 * @brief Reads the number of two-electron integrals.
 *
 * @param trexio_file Pointer to the open TREXIO file.
 * @param n_integrals Pointer to an int64 where the number of two-electron integrals will be stored.
 */
void number_2_e_integrals(trexio_t * trexio_file, int64_t* n_integrals){
	trexio_exit_code rc;
	rc = trexio_read_mo_2e_int_eri_size(trexio_file, n_integrals);
	if (rc != TREXIO_SUCCESS) {
                printf("TREXIO error reading the number of two electron integrals: %s\n", trexio_string_of_error(rc));
                exit(1);
        }
}

/**
 * @brief Reads the indices (i, j, a, b) and the value of the two-electron integrals.
 *
 * @param trexio_file Pointer to the open TREXIO file.
 * @param n_integrals Number of two-electron integrals.
 * @param index Pointer to an int32 pointer where the indices of the two-electron integrals will be stored. 
 * @param value Pointer to a double pointer where the value of the two-electron integrals will be stored. 
 */
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

/**
 * @brief Closes a TREXIO file.
 *
 * @param trexio_file Pointer to the TREXIO file structure to be closed.
 */
void close_fun(trexio_t * trexio_file){ 
 trexio_exit_code rc; 
 rc = trexio_close(trexio_file);
 if (rc != TREXIO_SUCCESS) {
	 printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
	 exit(1);
 }
 trexio_file = NULL;
}

