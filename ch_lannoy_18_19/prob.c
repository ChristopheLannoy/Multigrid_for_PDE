#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "main.h"

int prob(int m, int *n, int **ia, int **ja, double **a, double **b)
/*
   But
   ===
   Génère la matrice n x n qui correspond à la disrétisation sur une grille 
   cartesienne regulière m x m de l'operateur de Laplace à deux dimensions
              
            d    d        d    d
         - == ( == u ) - == ( == u )        sur [0,1] x [0,1]
           dx   dx       dy   dy

  avec la fonction u qui satisfait les conditions aux limites de Dirichlet
         
         u = 0  sur (0,y), (1,y), (x,0) et (x,1), avec 0 <= x,y <= 1 .
  
  La numérotation des inconnues est lexicographique, la direction x étént 
  parcourue avant celle de y. La matrice est retournée dans le format CRS
  qui est défini par le scalaire 'n' et les trois tableaux 'ia, 'ja' et 'a'.

  Arguments
  =========
  m (input)   - nombre de points par direction dans la grille 
  n  (output) - pointeur vers le nombre d'inconus dans le système
  ia (output) - pointeur vers le tableau 'ia' de la matrice A
  ja (output) - pointeur vers le tableau 'ja' de la matrice A
  a  (output) - pointeur vers le tableau 'a' de la matrice A
  b  (output) - pointeur vers le tableau 'b' 
  L (input)   - taille de la membrane        (via main.h)
  xratio, yratio (input) Définissent la forme de la membrane (via main.h)
  
   La matrice A à la dimension du nombre de points simulés potentiellement non nul
   (tant pour la hauteur du mode fondamental de vibration que pour la température)
  
   Par exemple:
   Pour m = 500, Simulation sur tout le domaine ==> dim(A) = 250 000 x 250 000
                 Simulation du problème 3  ==> dim(A) =environ  60 000 x  60 000 
*/
{
    int  nnz, ix, iy, nx, ind, ixlim, iylim, npts = 0;
    ind = 0;
    double invh2;

    nx = m - 2; /* noeuds de Dirichlet ne sont pas pris en compte */
    invh2 = (m-1)*(m-1)/(L*L); /* h^-2 pour L=1 *0.25 si L= 2m */
    *n  = nx * nx - ((m-1) * (m-1) * (1-x_ratio) * (1-y_ratio)) ;     /* nombre d'inconnues */
    nnz = 5 * (*n) - 4 * nx; /* nombre d'éléments non nuls */

    /* allocation des tableaux */

    *ia  = malloc((*n + 1) * sizeof(int)); 
    *ja  = malloc(nnz * sizeof(int));
    *a   = malloc(nnz * sizeof(double));
    *b   = malloc(*n * sizeof(double));

    /* allocation réussite? */

    if (*ia == NULL || *ja == NULL || *a == NULL || *b == NULL) {
        printf("\n ERREUR : pas assez de mémoire pour générer la matrice A\n\n");
        printf("n= %d , nnz = %d",*n,nnz);
        return 1;
    }

    /* partie principale : replissage de la matrice */
    nnz = 0;
    /* coordonnées du point caractérisant la membrane */
    ixlim = ((m-1) * x_ratio) - 1;
    iylim = ((m-1) * y_ratio) - 1;   
      
    for (iy = 0; iy < nx; iy++) {         // On parcourt tous les points
        for (ix = 0; ix < nx; ix++) {     // dans l'ordre lexicographique
            
            if(iy >= iylim && ix >= ixlim){  //Points hors de la menbrane
				npts++;      //On "saute" un point, on ne pose pas 
				continue;    //d'équation pour ce point et le compteur du
			}                //nombre de points ne donnant pas lieu à une équation augmente
            
            ind = ix + nx * iy -npts; /* numéro de l'équation */

            /* marquer le début de la ligne suivante dans le tableau 'ia' */
            (*ia)[ind] = nnz;
            
            /* calculer le membre de droite */
            (*b)[ind] = 0.0; 
   
            /* replissage de la ligne : voisin sud */
            if (iy > 0)  {
                (*a)[nnz] = -invh2; /* pour D=1 */
                
                if(iy > iylim){
					(*ja)[nnz] = ind - ixlim;
                }
                else{
					(*ja)[nnz] = ind - nx;
				}
                nnz++;
            }

            /* replissage de la ligne : voisin ouest */
            if (ix > 0)  {
                (*a)[nnz] = -invh2; /* pour D=1 */
                (*ja)[nnz] = ind - 1;
                nnz++;
            }
            /* replissage de la ligne : élém. diagonal */
            (*a)[nnz] = 4.0*invh2; /* pour D=1 */
            (*ja)[nnz] = ind;
            nnz++;

            /* replissage de la ligne : voisin est */
            if ( (iy < iylim && ix < nx - 1) || (iy >= iylim && ix < ixlim - 1) ){
                (*a)[nnz] = -invh2; /* pour D=1 */
                (*ja)[nnz] = ind + 1;
                nnz++;
            }

            /* replissage de la ligne : voisin nord */
            if ( (ix < ixlim && iy < nx - 1) ||(ix >= ixlim && iy < iylim - 1) ) {
                (*a)[nnz] = -invh2; /* pour D=1 */
                
                if(iy >= iylim){
					(*ja)[nnz] = ind + ixlim;
				}
				else{
					(*ja)[nnz] = ind + nx;
				}
                nnz++;
            }
            
            /* point a coté d'un bord du domaine non nul */
            if( (ix == ixlim-1  && iy >= iylim) || (iy == iylim-1 && ix >= ixlim) ){
				//printf("%d  %d \n ", ix, iy);
				(*b)[ind] += invh2; /* condition de Dirichlet */
			}    
        }
    }

    /* dernier élément du tableau 'ia' */
    (*ia)[ind + 1] = nnz;

    /* retour de fonction habituel */
    return 0;
}
