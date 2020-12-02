#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "main.h"

int plot_eigen(int m, double *evecs)

/*
 * Plot le vecteur propre evecs carrespondant à une grille de "m" points par coté
 * avec les conditions nulles au bord
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
    
    
    FILE *data = fopen("data_eigen.txt", "w"); //ecriture de fichier


    for (iy = -1; iy <= nx; iy++) {
        for (ix = -1; ix <= nx; ix++) {
			
			
/* Points appartenant au bord du carré */
			if(ix == -1 || iy == -1 || ix == nx || iy == nx){
				fprintf(data, "%f %f %f\n",(ix+1)*h, (iy+1)*h , 0.0);
				continue;
			}
/* Reste les points appartenant à l'intérieur du carré */
            
         /* Points appartenant à l'intérieur du carré mais pas à l'intérieur de la membrane */
            if(iy >= iylim && ix >= ixlim){ 
				npts++;
				fprintf(data, "%f %f %f\n",(ix+1)*h, (iy+1)*h , 0.0);
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
		"set pm3d at b\n"  //remplacer "at b" par "map" pour avoir une présentation en 2D (utile pour les m élevé >256)
		"set title 'Graphique de la hauteur pour m = %d'\n" 
		"splot 'data_eigen.txt' using 1:2:3 with lines\n",m //enlever "with lines" pour representation 2D
		 );
	fclose(cmd);

	//exécution
	system("gnuplot -persistent cmd.txt");
	
	return 0;
}
