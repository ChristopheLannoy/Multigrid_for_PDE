#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "main.h"
#include "umfpk.h"
#include "time.h"



int main(int argc, char *argv[])
{   // Valeur max
	// x=7 coarsest = 40
	// x=9 coarsest = 8

	/* variables tunable by the user */
	double tau = 1;           // relaxation factor (=1 by default)	
	int print_level = 2;      // 0-2
	int x = 4;                // define the number of level in the multigrid 
	int coarsest_level = 8;	  // define the coarset level of the MG (must be a multiple of 8)
						   	  //(The mesh size is thus defined by x and coarsest_level)
	
	
	/* other variables (not tunable)*/
	int m = (pow(2,x)*coarsest_level) + 1;	// m correspond au nombre de point sur un coté du carré 
	int n, *ia, *ja; 	// n: nombre de points sur la membrane = nombre d'inconnue du problème
	double *a, *b;      //ia,ja,a matrice A au format CSR, (A*u=b)
	double t1, t2;
	
	/* générer le problème */
	if (prob(m, &n, &ia, &ja, &a, &b))
		return 1;

	printf("\nPROBLEM: ");
	printf("m = %5d   n = %8d  nnz = %9d\n \n", m, n, ia[n] );
 
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/////////////////      RESOLUTION MULTI-GRID       /////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
	
	int k;
	double *u, *u_p;      
	u = malloc(n *sizeof(double));    // allocate memory for solution vector 
	u_p = malloc(n * sizeof(double)); //previous solution (for the relaxation, tau)
	if (u == NULL || u_p == NULL) {
		printf("\n ERREUR : pas de mémoire pour vecteur solution u et u_p\n\n");
        return 1;
	}
	
	for(k = 0; k <n; k++)     // Approximation initiale nulle de la solution u_k
		u[k] = 0;
 
	int nb_it = 0;     //number of iteration
	int max_it = 50;   //maximum niumber of iteration
	double nrmres[max_it], conv_rate[max_it];   //residual norms
	nrmres[0] = norm_residu(n, b, u, ia, ja, a);
	if(print_level > 1)
		printf("0 iteration of MG: resiual norm = %1.3e \n", nrmres[0]);
	
	t1 = mytimer();
	
	while(nb_it < max_it){
		
		for(k = 0; k < n; k++)  //stock the previous solution
			u_p[k] = u[k];
	
		if(multi_grid(m, n, ia, ja, a, b, u, coarsest_level+1) )
			return 1;
		nb_it++;
		
		for(k = 0; k < n; k++)     //impementation of the relaxtion
			u[k] = u_p[k] + (u[k] - u_p[k])*tau;
		
		nrmres[nb_it] = norm_residu(n, b, u, ia, ja, a);
		conv_rate[nb_it - 1] = nrmres[nb_it]/nrmres[nb_it - 1];
		
		if(print_level > 1){
			printf("%d iterations of MG: residual norm = %1.3e || ",nb_it, nrmres[nb_it]);
			printf("Convergence rate  = %1.4f\n", conv_rate[nb_it - 1]);
		}
		if(nrmres[nb_it-1] < 1.1*nrmres[nb_it]){  //si la norme commence à stagner, on atteint la précision machine
			printf("\n\n------- END OF MULTIGRID : residual norm stagnates ----------------\n");
			break;					 //on arrete le multigrid 
		}
	}
	
	if(nb_it == max_it)
		printf("\n\n-------END OF MULTIGRID :max iteration limit (%d) reached ----- \n", max_it);

	t2 = mytimer();
	
	if(print_level > 0){
		printf("\n\nMULTIGRID ANALYSIS\n"
				"---------------------------------------------------\n");
		printf("  Solution time for multigrid (CPU): %5.4f sec\n",t2-t1);
		printf("  Number of iteration of multigrid: %d \n", nb_it);
		printf("  Residual norm = %1.3e \n",nrmres[nb_it]);
	}
	
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
///////////   Q4 Convergence analysis and relaxation   /////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


/* erreur relative  ||r||/||b|| */
	double rel_err = nrmres[nb_it]/(norme(n, b));
	double h2 = (L*L)/((m-1)*(m-1));
	double stability_nb = (h2*nrmres[nb_it])/(norme(n, u)*8);
	
	if(print_level > 0){
		printf("\n\nCONVERGENCE ANALYSIS\n"
				"---------------------------------------------------\n"
				"  Relative error (||r|| / ||b||)  = %e \n\n", rel_err);
		printf("  ||r|| / (||u||.||A||) = %3.3e ?<? u = 1.1e-16 \n", stability_nb);
		if(stability_nb < 1.1e-16)
			printf("      ==> Direct stability criterion verified \n\n");
	}
	
//Facteur de convergeance assymptotique (avant d'atteindre la limite de precision machine)
	double rho = 0;    
	for(k = nb_it - 1; k > nb_it - 10; k--){
		if(conv_rate[k] < 1.01 * conv_rate[k-1] ){
			rho = conv_rate[k];
			break;
		}
	}
	if(k == nb_it - 10){ //pas de stabilite du facteur de convergeance dans les 10 dernières itérations
		printf(" Asymptotic convergence rate not found \n");
	}else{
		printf("  Asymptotic convergence rate (with tau = %1.3f): rho = %1.3f \n\n",tau, rho); 
		if (tau == 1){
			printf("  Recommanded relaxtion factor: \n "
				   "     tau = 2/(lambda_max + lambda_min) = %5.3f \n \n",2.0/(2-rho));
			printf("  Theoritical estimate of new asymptotic convergence \n "
				   "  (with tau as recommanded)  rho = %5.3f \n \n", rho/(2-rho));
		}
	}
	
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//////////////////  Affichage graphique des solutions //////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////	
	
	printf("\n\nGRAPHIC ANALYSIS\n");
	printf("---------------------------------------------------\n");
	char c = ' ';
	while(c != 'y' && c != 'n'){
		printf("\n  Show the evolution of the residual norm (graphic) (y/n) ?  ");
		c = getchar();
		getchar();   //consomme le caractere <return>.
	}
	if(c == 'y')
		plot_vector(nb_it, nrmres);
	
	c = ' ';
	while(c != 'y' && c != 'n'){
		printf("\n  Show the graphic of the solution (y/n) ?  ");
		c = getchar();
		getchar();
	}
	if(c == 'y')
		plot(m, u);
	
 
////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////// 
///////////////////    Q5 VALEUR PROPRE MIN     ////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


	printf("\n\n\nEIGENVALUE PROBLEM RESOLVED WITH PRIMME\n"
		   "---------------------------------------------------\n");
  
    // primme - résolution 
	int nev = 1;           //nombre de pair vecteur/valeur propre désiré
	double *evals, *evecs;
  
    // allouer la memoire pour vecteurs & valeurs propres 
	evals = malloc(nev * sizeof(double));
	evecs = malloc(nev * n * sizeof(double));
	if (evals == NULL || evecs == NULL) {
		printf("\n ERREUR : pas assez de mémoire pour les vecteurs et valeurs propres\n\n");
		return 1;
	}  
    
	t1 = mytimer();
	if(primme(n, ia, ja, a, nev, evals, evecs, m))
		return 1;
	t2 = mytimer();
  
	// temps de solution 
	printf("\n  Temps de solution (CPU) (PRIMME): %5.1f sec\n",t2-t1);
	printf("  Valeur propre minimale calculée (PRIMME): %5.2f \n",evals[0]);
  
  
	c = ' ';
	while(c != 'y' && c != 'n'){
		printf("\n  Show the graphic of the solution of the eigenvalue problem (y/n) ?  ");
		c = getchar();
		getchar();
	}
	if(c == 'y')
		plot_eigen(m, evecs);
 
	/* libérér la mémoire */
	free(evals); free(evecs);
	free(ia); free(ja); free(a); free(b); free(u); free(u_p);
	
	
	
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////// Code servant a tester le bon fonctionnement du programme ////////
////////////////////////////////////////////////////////////////////////	
////////////////////////////////////////////////////////////////////////	
	
	
	                                     /* affichage de la matrice A */  
  
/*
	// on peut éventuellement afficher des matrices pour d'autre m (multiple de 8)
								    
	int jj;
  
	printf("\n MATRICE A \n");
	for (jj=0; jj<=n; jj++){
		printf("%d \n ", ia[jj]);
	}
	printf("\n \n");
	for (jj=0; jj<ia[n]; jj++){
		printf("%d \n ", ja[jj]);
	}
	printf("\n \n");
	for (jj=0; jj<ia[n]; jj++){
		printf("%f \n ", a[jj]);
	}
*/
	
		
	                                                /*  test mat_vec  */ 
/* 
	double a1[5] = {1,1,2,1,3}; //  1 0 1 
	int ja1[5] = {0,2,1,0,2};   //  0 2 0
	int ia1[4] = {0,2,3,5};     //  1 0 3
  
  
	printf("%p \n", a1);
	printf("%f \n", a1[2]);
	printf("%f \n", a1[0]);
	printf("%f \n \n", *a1);

	double vx[3] = {1,1,1};
	double vy[3];

	mat_vec(vx, vy, 3, ia1, ja1, a1);

	int jj;
	for (jj=0; jj<=2; jj++){
		printf("%f \n ", vy[jj]);
	}
*/
	
							  /*  résolution par smooting uniquement  */
/*
	smoothing(ia, ja, a, b, u, 1, n); 
	norm_residu( n, &nrmres, b, u, ia, ja, a);
	printf("norme résiduelle apres 117 it = %e \n", nrmres);
*/
 
 
                   /*  Comparason de la solution (via multigrid)
                    *  avec la solution fournie par le solveur exact  */
/*  
	double *x;
	x = malloc(n * sizeof(double)); // allouer la mémoire pour le vecteur de solution 
	if ( x == NULL ) {
		printf("\n ERREUR : pas de mémoire pour vecteur x\n\n");
		return 1;
	}
	getchar();

	plot(0,u,L);   // Aucune idée de pourquoi mais le soveur direct ne fonctionne pas sans cette ligne
 
	t1 = mytimer();
	//solve_umfpack(n, ia, ja, a, b, x);   //Pour une raison inconnue ce solveur est particulièrement loquace
	t2 = mytimer();
	printf("\nTemps de solution solveur direct(CPU): %5.1f sec\n",t2-t1);
  
	norm_residu( n, &nrmres, b, x, ia, ja, a);
  
	printf("norme résiduelle solveur direct = %e \n", nrmres);
 
	//plot(m,x,L,x_ratio,y_ratio);


//comparaison solveur direct et solveur itératif (u_k == x ?)
	int jj;
	for (jj=0; jj<n; jj+=100000){
		printf("%f \n ", x[jj]);
		printf("%f \n ", u[jj]);
        printf("\n ");
	}

	//plot(m,x,L,xratio,yratio);
	//plot(m,u,L,xratio,yratio);
	free(x);
*/ 
	

									  /* test solveur triangulaire L, U  */
/* 
	int n1= 3;                    //taille de la matrice 1
	double a1[5] = {1,1,2,1,3};   //  1 0 1 
	int ja1[5] = {0,2,1,0,2};     //  0 2 0
	int ia1[4] = {0,2,3,5};       //  1 0 3
 
	double temp[3] = {0,0,0};
	double r1[3] = {0,1,1};
    
	solve_U( n1, ia1, ja1, a1, r1, temp);
  
  
	int i;
	for(i=0; i<n1; i++){
	printf("temp[%d] = %e \n",i,temp[i]);
	}

	//Vecteur diagonal 
	double *d1;
	diagonal(ia1, ja1, a1 , &d1, n1);
	int r;
	for(r=0; r < n1; r++){
		printf("VECTEUR DIAG : d1[%d]= %e \n",r, d1[r]);
	}
*/
 

	return 0;
}

