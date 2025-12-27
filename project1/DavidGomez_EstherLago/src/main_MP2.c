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
	double E_1e = 0;
	rc = trexio_read_mo_num(trexio_file, &mo_num); //Number of MO
	printf("Number of Molecular Orbitals	= %i\n", mo_num);
	double * data = malloc(mo_num * mo_num * sizeof(double));
	rc = trexio_read_mo_1e_int_core_hamiltonian(trexio_file, data);
    	printf("\n\tOne-Electron Integrals (au)\n");
	for (int x=0 ; x<n_up*mo_num ; x+=mo_num+1) { // 25 pq solo elementos diagonal ppal, y es matriz 24x24 
        E_1e = E_1e + data[x];
	printf("<i|h|i> (i=%d)       = %lf\n",(x%mo_num)+1, data[x]);
    	}
	
	//MO orbitals para MP2
	double * mo_energy = malloc(mo_num * sizeof(double));
	rc = trexio_read_mo_energy(trexio_file, mo_energy);
	if (rc != TREXIO_SUCCESS) {
    		printf("TREXIO error reading molecular orbital energies: %s\n", trexio_string_of_error(rc));
    		return 1;
	}
	/*Read and write Two-Electron integrals*/
	double E_2e_ijij = 0;
	double E_2e_ijji = 0;
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
	/* Pruebas pero sobra
	printf("%ld\n", n_integrals);
	printf("%ld\n", offset_file);
	printf("%ld\n", buffer_size);
	for (int64_t i=0 ; i<n_integrals ; i++) {
		printf("integral[%ld]	=%lf\n", i, value[i]);
	}
        for (int64_t i=0 ; i<100 ; i++) {
                printf("integral[%ld]   =%d\n", i, index[i]);
        } // aqui salen 4 valores seguidos, que representan los indices de la integral <ij||kl> */
	double division = 0;
	for (int64_t n=0; n<n_integrals; n++) {
		int i = index[4*n+0];
		int j = index[4*n+1];
		int a = index[4*n+2];
		int b = index[4*n+3];
		double integral = value[n];
		if (i<n_up && j<n_up && i==a && j==b) {
			printf("<ij|ij> = <%d %d|%d %d> =	%lf\n", i, j, a, b, integral); 
			E_2e_ijij = E_2e_ijij + integral;
			if (i!=j) {
				E_2e_ijij = E_2e_ijij + integral;
			}
			
		}	
		if (i<n_up && a<n_up && i==j && a==b) { //esto es lo mismo que que ij|ji por la conversion que se hace en la pagina 7 del pdf 
			printf("<ij|ji> = <ii|jj> = <%d %d|%d %d> =	%lf\n", i, j, a, b, integral);
			E_2e_ijji = E_2e_ijji + integral;
			if (i!=b) {
				E_2e_ijji = E_2e_ijji + integral;
			}
		}
		//Two-electron integrals para MP2
		
		
		
		if (a<n_up && b<n_up && i>=n_up && j>=n_up) {
			printf("uno<ij|ab> =	<%d %d|%d %d> = 	%lf\n", a, b, i, j, integral);//Estamos cogiendo <a b | i j>
			if (i==j) {
				double patatas = integral * integral ;
				double Filadelfia = mo_energy[a] + mo_energy[b] - 2*mo_energy[i];
				division = division + patatas/Filadelfia;
				printf("dos\n");
			} else {
				for (int64_t m=0; m<n_integrals; m++) {
		                	int i2 = index[4*m+0];
        	        		int j2 = index[4*m+1];
                			int a2 = index[4*m+2];
                			int b2 = index[4*m+3];
        	        		double exchange = value[m];
					if (a2==b && b2==a && i2==i && j2==j) {
						printf("dos<ij|ab> =	<%d %d|%d %d> =	%lf\n", a2, b2, i2, j2, exchange);
				
				double patatas = integral * (2*integral - exchange);
				double Filadelfia = mo_energy[a] + mo_energy[b] - mo_energy[i] - mo_energy[j];
				division = division + patatas/Filadelfia;
					}
				}
			}	
		}
	}
	//Final HF sum 
	double E;
	E = Enn + 2*E_1e + 2*E_2e_ijij - E_2e_ijji;
	printf("Enn=%lf\n E_1e=%lf\n 2*E_2e_ijij=%lf\n E_2e_ijji=%lf\n",Enn, E_1e, 2*E_2e_ijij, E_2e_ijji);
	printf("E = %lf\n", E);

	//Ahora vamos a hacer MP2
	/*
	double * mo_energy = malloc(mo_num * sizeof(double));
	rc = trexio_read_mo_energy(trexio_file, mo_energy);
	if (rc != TREXIO_SUCCESS) {
    		printf("TREXIO error reading molecular orbital energies: %s\n", trexio_string_of_error(rc));
    		return 1;
	}*/
	for (int64_t i=0; i<mo_num; i++) {
		printf("MP2[%ld]=%lf\n", i, mo_energy[i]);
	}
	close_fun(trexio_file);	
	printf("%lf\n", division);
		return 0;
}
