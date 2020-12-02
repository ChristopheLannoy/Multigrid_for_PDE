#include <stdlib.h>

void mat_vec(double *vx, double *vy, int n, int *ia, int *ja, double *a)
/*
   But
   ===
   Calcule le produit matrice-vecteur
                              vy = A*vx 
   

   Arguments
   =========
   vx        (input) - vecteur(s) d'entrée
   n         (input) - taille des vecteur vx et vy
   ia, ja, a (input) - matrice A au format CSR
   vy       (output) - vecteur(s) de produit A*vx

*/
{
    int i, j;
	for(i = 0; i < n; i++){
		vy[i] = 0;
        for (j = ia[i]; j < ia[i + 1]; j++)
			vy[i] += a[j] * vx[ja[j]];
	}
}
