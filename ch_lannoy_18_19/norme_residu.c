#include <stdlib.h>
#include <stdio.h>

#include "main.h"

double norm_residu( int n, double *b, double *x, int *ia,
                  int *ja, double *a)
/*
 * Calcule la norme résiduelle du systeme A*x = b 
 * Retourne ||r|| où r = b - A.x
 * 
 * INPUT:
 * 			x           : vecteur qui approche la solution du système
 * 			n           : taille des vecteur x et b  (A=(n x n))
 * 			ia, ja , a  : A au format CSR
 * OUTPUT:
 * 			nrmres      : norme résiduelle
 * 
 */

{
	double *r, nrmres;
	int i;
	r = malloc(n * sizeof(double) );
    if (r == NULL){ 
        printf("\n ERREUR : pas assez de mémoire pour générer la matrice (residu)\n\n");
        return 1;
    }
    
/* r =  A.x  */
	mat_vec(x, r, n, ia, ja, a);
/* r = A.x - b*/
	for (i=0; i < n; i++) {
		r[i] -= b[i];	
	}
/* norme 2 de r */
	nrmres = norme(n, r);
	free(r);
	
	return nrmres;
}
