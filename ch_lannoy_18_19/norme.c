#include <math.h>

double norme(int n, double *x)

/*
 * Calcule la norme (2) du vecteur x
 * 
 * INPUT: 
 * 			n : taille du vecteur x
 * 			x : vecteur x
 * OUTPUT:
 * 			return ||x||
 */

{
	int i;
	double a = 0.0; 
	for (i=0; i < n; i++) {
		a += pow(x[i], 2);
	}
	
	return sqrt(a); 
}
