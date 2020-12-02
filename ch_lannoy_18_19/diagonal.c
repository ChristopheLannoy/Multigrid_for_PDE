#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int diagonal(int *ia, int *ja, double *a , double **d,int n)

/*
 * Return the diagonal of the matrix A in a vector d
 * 
 * INPUT :  
 * 			ia, ja, a : Matrix A in CSR format
 * 			n         : size of the matrix A (n*n) and vector d
 * 
 * OUPUT : 
 * 			d : vector d witch is the diagonal of A
 */

{		
	*d = malloc(n*sizeof(double));
	if (*d == NULL) {
        printf("\n ERREUR : pas assez de mémoire pour générer la diagonale d \n\n");
        return 1;
    }

	int i,j,nelem_a=0;
	for(i=0;i<n;i++){
		for(j= ia[i] ;j<ia[i+1];j++){
			if(ja[j]== i){              //if diagonal element
				(*d)[i]= a[nelem_a];
			}
		nelem_a++;
		}
	}

	return 0;
}
