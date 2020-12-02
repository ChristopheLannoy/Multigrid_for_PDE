#include <stdlib.h>
#include "primme.h"
#include "umfpk.h"
#include "main.h"

/* variables statiques -- accessibles  */
static double *a;
static int n, m, *ia, *ja;

void matvec_primme(void *vx, void *vy, int *blockSize, primme_params *primme)
/*
   But
   ===
   Calcule le produit matrice-vecteur
                              vy = A*vx 
   pour le solveur aux valeurs propres PRIMME. La matrice A doit être
   "stoquée" au préalable dans les variables statiques 'n', 'ia', 'ja' et 'a'
   en utilisant le format CSR (Compressed Sparse Rows). Par "stoquer" 
   on veut dire ici stoquer la valeur de 'n' et les pointeurs vers les 
   tableaux 'ia', 'ja' et 'a'.

   Arguments
   =========
   vx        (input) - vecteur(s) d'entrée
   vy       (output) - vecteur(s) de produit A*vx
   blockSize (input) - nombre de vecteurs d'entrée
   primme    (input) - paramètres fournis par primme pour optimiser le calcul
                       (pas utilisé) 
*/
{
    int i, j, b;
    double *x = vx, *y=vy;

    for(b = 0; b < (*blockSize)*n; b+=n)
        for(i = 0; i < n; i++){
            y[b+i] = 0;
            for (j = ia[i]; j < ia[i + 1]; j++)
                y[b+i] += a[j] * x[b+ja[j]];
        }
} 

void preconditioner_primme(void *vx, void *vy, int *blockSize, struct primme_params *primme)
/*
	Block preconditioner-multivector application, y = M^(−1) x where M is usually an approximation of A−σI
	or A − σB for finding eigenvalues close to σ. The function follows the convention of matrixMatvec.
*/
{
	double *x = vx, *y = vy;
	int k, l;
	
	for(k = 0; k < (*blockSize)*n; k+=n){
/*
		if(solve_umfpack(n, ia, ja, a, &x[k], &y[k]) ) 
			printf("\nPROBLEME avec umfpack dans PRIMME \n");
*/
		for(l = 0; l <n; l++)              // initial approx = nul vector
		    y[k + l] = 0;

		multi_grid(m, n, ia, ja, a, &x[k], &y[k], 9);  //Solve A*y = x ==> y = A^(-1) * x
	}
	//printf("Utilisation du Preconditioner dans PRIMME \n");
}

/**/

int primme(int primme_n, int *primme_ia, int *primme_ja, double *primme_a, 
           int nev, double *evals, double *evecs, int primme_m)

/*
   But
   ===
   Calcule les nve valeurs propres les plus basses de la matrice A stoquée 
   dans le format CSR à l'aide du scalaire primme_n et des vecteurs 
   primme_ia, primme_ja et primme_ia. 

  Arguments
  =========
  primme_n  (input) - le nombre d'inconus dans le système
  primme_ia (input) - le tableau 'ia' de la matrice A
  primme_ja (input) - le tableau 'ja' de la matrice A
  primme_a  (input) - le tableau 'a' de la matrice A
  nev       (input) - le nombre de valeurs propres recherchées
  evals    (output) - le tableau des valeurs propres
  evecs    (output) - le tableau des vecteurs propres
  max      (input)  - 0 si l'on cherche valeur propre min, 1 si max

  Retourne 0 si le caclul s'est bien déroulé, 1 si non.
*/
{
    int err;

    /* norme des résidus */
    double *resn = malloc(nev * sizeof(double));
    if (resn == NULL) {
        printf("\n ERREUR : pas assez de mémoire pour un vecteur auxilière dans la fonction primme\n\n");
        return 1;
    }

    /* sauvgarder les pointeurs dans des variables statiques */
    n = primme_n;
    m = primme_m;
    a = primme_a;
    ja = primme_ja;
    ia = primme_ia;

    /* encoder les paramètres de PRIMME */
    primme_params primme;
    primme_initialize (&primme);
    primme.matrixMatvec = matvec_primme; /* MV product */
    primme.n = primme_n;         /* Matrix dimensions */
    primme.numEvals = nev; /* Number of wanted eigenpairs */
    primme.printLevel = 2; /* 1-4 */
    
    primme.applyPreconditioner = preconditioner_primme;
    primme.correctionParams.precondition = 1;
    
    //Précise si l'on veut calculer les valeurs propres minimales ou maximales
	//primme.target = primme_largest;
	primme.target = primme_smallest;
    
    
    if((err = primme_set_method (DEFAULT_MIN_TIME, &primme))){
        printf("\nPRIMME: erreur N % dans le choix de la methode \n    (voir 'Error Codes' dans le guide d'utilisateur)\n",err);
        return 1;
    }
  
    /* afficher les papramètres de PRIMME */
    primme_display_params (primme);

    /* Caclul des valeurs et vecteurs propres */
    if( (err = dprimme (evals, evecs, resn, &primme)) ){
        printf("\nPRIMME: erreur N %d dans le calcul des valeurs propres \n    (voir 'Error Codes' dans le guide d'utilisateur)\n",err);
        return 1;
    }
    
    printf("\n\n  Norme residuelle du solveur primme %e", resn[0]);
    
    /* libérer la mémoire */
    primme_Free (&primme); free(resn);

    return 0;
}
