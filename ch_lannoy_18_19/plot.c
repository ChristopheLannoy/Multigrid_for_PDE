#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "main.h"

int plot(int m, double *evecs)

/*
 * 
 * Plot le vecteur evecs sur une grille de "m" points par coté
 * avec les condition au bord du problème de diffusion
 * 
 * 
 * 
 */ 

{
    
    double h = L/(m-1);
	int ix, iy, ind;
	
	int nx = m-2;
	
	int ixlim = ((m-1) * x_ratio) - 1;
    int iylim = ((m-1) * y_ratio) - 1;   /* coordonnées du point caractérisant la membrane */
    
    int npts = 0;
        
    FILE *data = fopen("data.txt", "w"); //ecriture de fichier

    for (iy = -1; iy <= nx; iy++) {
        for (ix = -1; ix <= nx; ix++) {
			
			
/* Points appartenant au bord du carré total  */
			if(ix == -1 || iy == -1 || (ix==nx && iy <= iylim)|| (iy == nx && ix <= ixlim)){
				fprintf(data, "%f %f %f\n",(ix+1)*h, (iy+1)*h , 0.0);
				continue;
			}
			
			if((ix==nx && iy > iylim)|| (iy == nx && ix > ixlim)){
				//fprintf(data, "%f %f %f\n",(ix+1)*h, (iy+1)*h , 1.0);
				continue;
			}

			
/* Reste les points appartenant à l'intérieur du carré */
            
         /* Points appartenant à l'intérieur du carré mais pas à l'intérieur de la membrane */
            if(iy >= iylim && ix >= ixlim){ 
				npts++;
				if( (iy == iylim && ix >= ixlim) || (iy >= iylim && ix == ixlim)){
					fprintf(data, "%f %f %f\n",(ix+1)*h, (iy+1)*h , 1.0);  //0.0 pour un vecteur correction !!!
				}	
				continue;
			}
		/*  Point intérieur à la membrane (donc seuls points potentiellement non nuls) */
			ind = ix + nx * iy -npts;   /* numéro de l'éauqtion */
			fprintf(data, "%f %f %f\n",(ix+1)*h, (iy+1)*h , evecs[ind]);
		}
		fprintf(data, "\n");
	}
	fclose(data);
	
	
	// fichier commande
	FILE *cmd = fopen("cmd.txt", "w");
	fprintf(cmd,
		"unset key\n"
		"set pm3d map \n"
		//"set cbrange [0:10] \n"
		"set title 'Graphique de la diffusion pour m = %d'\n" 
		"splot 'data.txt' using 1:2:3 \n",m
		 );
	fclose(cmd);

	//exécution
	system("gnuplot -persistent cmd.txt");
	
	return 0;
}
