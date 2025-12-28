/**
 * @file main_final.c
 * @brief Program to read data from a TREXIO file using it to compute the 
 * Hartree-Fock and MP2 energies
 * @author David Gómez and Esther Lago
 */
#include <stdio.h>
#include <trexio.h>
#include <stdlib.h>
#define pointer(i,j,k,l) ((i)*mo_num*mo_num*mo_num + (j)*mo_num*mo_num + (k)*mo_num +(l))

trexio_t * open_fun(const char * filename);

double Nuclear_repulsion_energy(trexio_t * trexio_file);

int Occupied_orbitals(trexio_t * trexio_file);
int Molecular_orbitals(trexio_t * trexio_file);
double * one_e_integrals(trexio_t * trexio_file, int mo_num);
double * mo_energies(trexio_t * trexio_file, int mo_num);
int number_2_e_integrals(trexio_t * trexio_file);
double * two_e_integrals(trexio_t * trexio_file, int n_integrals, int32_t** index, double** value);
void close_fun(trexio_t * trexio_file); 

double main() {
	//Open input as TREXIO file
	char input[30];
	scanf("%s", input);
	trexio_t* trexio_file = open_fun(input);
	trexio_exit_code rc;

	double Enn;
	int32_t n_up;
	int32_t mo_num;
	int64_t n_integrals;
	
	/* Read and write nuclear repulsion energy */
	Enn = Nuclear_repulsion_energy(trexio_file);
	printf("Nuclear Repulsion Energy (au)	= %lf\n", Enn);

	/*Read and write the number of occupied orbitals */
	n_up = Occupied_orbitals(trexio_file);
	printf("Number of Occupied Orbitals	= %i\n", n_up);

	/*Read and write One-Electron Integrals*/
	double E_1e = 0;
	mo_num = Molecular_orbitals(trexio_file);
	printf("Number of Molecular Orbitals	= %i\n", mo_num);
	double * data;
	data = one_e_integrals(trexio_file, mo_num);
	
	for (int x=0 ; x<n_up*mo_num ; x+=mo_num+1) { // 25 pq solo elementos diagonal ppal, y es matriz 24x24 
        E_1e = E_1e + data[x];
    	}

	//MO orbitals para MP2
        double * mo_energy = malloc(mo_num * sizeof(double));
        mo_energy = mo_energies(trexio_file, mo_num);

	/*Read and write Two-Electron integrals*/
	double Two_e_HF = 0;
	n_integrals = number_2_e_integrals(trexio_file);
	int32_t* index = malloc(4 * n_integrals * sizeof(int32_t));
	
	double* value = malloc(n_integrals * sizeof(double));
	
	two_e_integrals(trexio_file, n_integrals, &index, &value);

	double* ERI = malloc(mo_num * mo_num * mo_num * mo_num * sizeof(double));//save n-fold integrals

        double MP2_corr = 0;
        for (int64_t n=0; n<n_integrals; n++) {
                int i = index[4*n+0];
                int j = index[4*n+1];
                int a = index[4*n+2];
                int b = index[4*n+3];
                double integral = value[n];
                ERI[pointer(i,j,a,b)] = integral;
                ERI[pointer(a,j,i,b)] = integral;
                ERI[pointer(i,b,a,j)] = integral;
                ERI[pointer(a,b,i,j)] = integral;
                ERI[pointer(j,i,b,a)] = integral;
                ERI[pointer(b,i,j,a)] = integral;
                ERI[pointer(j,a,b,i)] = integral;
                ERI[pointer(b,a,j,i)] = integral;
        }

	//Two-electron integrals para MP2
        for (int i=0; i<n_up; i++){
                for (int j=0; j<n_up; j++){
                        for (int a=0; a<mo_num; a++){
                                for (int b=0; b<mo_num; b++){
					if (i==a && j==b) {
						Two_e_HF = Two_e_HF + (2 * ERI[pointer(i,j,a,b)] - ERI[pointer(i,j,b,a)]);
					}


					if (a>=n_up && b>=n_up){
	                                        double num = ERI[pointer(i,j,a,b)] * (2*ERI[pointer(i,j,a,b)] - ERI[pointer(i,j,b,a)]);
        	                                double den = mo_energy[i] + mo_energy[j] - mo_energy[a] - mo_energy[b];
                	                        MP2_corr = MP2_corr + num/den;
                                	}
				}
                        }
                }
        }


	close_fun(trexio_file);	
	
	
	//Final HF sum 
	double E_HF;
	E_HF = Enn + 2*E_1e + Two_e_HF;
	printf("E_HF = %lf\n", E_HF);

	//Final MP2 energy
	printf("MP2 energy correction= %lf\n", MP2_corr);
	double E_MP2 = E_HF+MP2_corr;
	printf("MP2 energy=%lf\n", E_MP2);
		return 0;
}
