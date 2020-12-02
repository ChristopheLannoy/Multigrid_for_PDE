#include "multi_grid.h"
#include "main.h"
#include "umfpk.h"
#include <stdio.h>
#include <stdlib.h>

int multi_grid(int m_fine, int n_fine, int *ia, int *ja, double *a, double *b, double *u, double m_coarsest_level)

// Au = b ==> u= ?
{
	int smoothing_nb = 1; //number of pre-post smoothing
	int j;
/* générér le problème grossier */
	int m_coarse = ((m_fine-1) /2) +1;
	int n_coarse, *ia_coarse, *ja_coarse; 
	double *a_coarse, *b_osef;
	prob(m_coarse, &n_coarse, &ia_coarse, &ja_coarse, &a_coarse, &b_osef);
	//printf("\nPROBLEM: ");
	//printf("m_coarse = %5d   n_coarse = %8d  nnz_coarse = %9d\n", m_coarse, n_coarse, ia_coarse[n_coarse] );

/* génère les différents vecteurs pour les grilles fine et grosière   */	
	double *r_fine, *x_fine, *r_coarse, *x_coarse;	
    r_fine = malloc(n_fine * sizeof(double));
    x_fine = malloc(n_fine * sizeof(double));
    r_coarse = malloc(n_coarse * sizeof(double));
    x_coarse = malloc(n_coarse * sizeof(double));
    if (r_fine == NULL || x_fine == NULL || r_coarse == NULL || x_coarse == NULL) {
		printf("\n ERREUR : pas de mémoire dans le multi-grid (r/x, coarse/fine) \n\n");
		return 1;
    }
    
	smoothing(ia, ja, a, b, u, smoothing_nb, n_fine);

	// calcul du résidu apres 1 pré-smoothing
	mat_vec(u, r_fine, n_fine, ia, ja, a);   // r_fine = A*u
    for(j = 0; j < n_fine; j++){
		r_fine[j] = b[j] - r_fine[j];        // r = b - A*u_0
	}
	  
	fine_to_coarse(m_fine, r_fine, r_coarse);    // Restriction du résidu

/* solve A_c * x_c = r_c  */

	if(m_coarse > m_coarsest_level){    // > 8 pour le plus bas niveau
		for(j = 0; j < n_coarse; j++)      // aproximation initiale est un vecteur nulle
		    x_coarse[j] = 0;
	
		if(multi_grid(m_coarse, n_coarse, ia_coarse, ja_coarse, a_coarse, r_coarse, x_coarse, m_coarsest_level) )
			return 1;
	}
	
	else{
		if(solve_umfpack(n_coarse, ia_coarse, ja_coarse, a_coarse, r_coarse, x_coarse) ){ // x_c = A_c ^(-1)  *  r_c
			printf("\nPROBLEME avec umfpack dans multi_grid \n");	
			return 1;
		}
	}

    for(j = 0; j < n_fine; j++)
		x_fine[j] = 0;              //initialise x_fine a 0
	coarse_to_fine(m_coarse, x_fine, x_coarse);
	
	for(j = 0; j < n_fine; j++)
		u[j] = u[j] + x_fine[j];


	free(r_fine);free(r_coarse);free(x_coarse);free(x_fine);
	free(ia_coarse);free(ja_coarse);free(a_coarse);free(b_osef);
	
	smoothing(ia, ja, a, b, u, smoothing_nb, n_fine);
	
	return 0;
}
