#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include "smoothing.h"
#include "umfpk.h"

int smoothing(int *ia, int *ja,double *a,double *b, double *u_iter, int n_iter, int n)

{
	double *r, *temp;	
	r = malloc(n * sizeof(double));
	temp = malloc(n * sizeof(double));
	if (r == NULL || temp == NULL) {
		printf("\n ERREUR : pas de mémoire pour vecteur dans la fct smoothing n\n");
		return 1;
	}

	double *d;     //recupere la diagonale de A et la stocke dans le vecteur d
	diagonal(ia, ja, a , &d, n);
     
     // sim G_S
	
	int i, j;
	for(i = 0; i < n_iter; i++){
		
		mat_vec(u_iter, r, n, ia, ja, a);   // r = A*u_0
		for(j = 0; j < n; j++){
			r[j] = b[j] - r[j];  // r = b - A*u_0
		}

		solve_L(n, ia, ja, a, r, temp);     // temp = L^(-1) * r
		for(j = 0; j <n; j++){              // temp =   D *  L^(-1) * r
			temp[j] = d[j]*temp[j];
		}   
		solve_U(n, ia, ja, a, temp, r);      // r =(temp =)  U^(-1) *[ D *  L^(-1) * r]
	                                         // plus besoin de r ==> on peut stocker temp = B^(-1)*r dans r
		for(j = 0; j <n; j++){              // 
			u_iter[j] = u_iter[j] + r[j];
		}
		//r_k+1
		
	}
	free(d);
	free(r);free(temp);
	
      /* Jacobi   */
/*  
	for(i = 0; i < n_iter; i++){
		
		mat_vec(u_iter, r, n, ia, ja, a);   // r = A*u_0
		for(j = 0; j < n; j++){
			r[j] = b[j] - r[j];  // r = b - A*u_0
		}  
		
		for(j = 0; j <n; j++){              // temp =   D^(-1) * r
			temp[j] = (1/d[j])*r[j];
		}   
		for(j = 0; j <n; j++){              // 
		    u_iter[j] = u_iter[j] + temp[j];
		}	
	}
*/
	return 0;
}
