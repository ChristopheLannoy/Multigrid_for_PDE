#include "multi_grid.h"
#include "main.h"

void coarse_to_fine(int m_coarse, double *x_fine, double *x_coarse)

// x_fine must be initialised to 0 vector !!

/*
 * x_fine = P * x_coarse
 *  
 * Prolonge le vecteur x_coarse (solution du problème A_c * x_c = r_c) 
 * corespondant à une grille grossière de "m_coarse" points par coté de la grille
 * en un nouveau vecteur x_fine d'une grille fine
 * 
 * La grille fine est définie par h_fine (distance entre 2 points): 2.h_fine = h_coarse
 * 
 * 
 * INPUT :
 * 			m_coarse : paramètre de la grille coarse (nombre de point sur un coté)
 * 			x_coarse : vecteur initial de dimention correspondant à une grille m_coarse
 * 
 * OUTPUT : 
 * 			x_fine : vecteur prolongé de dimension correpondant à une grille m_fine
 *                    Attention il doit être initilisé comme vecteur nul! 
 * 
 */ 


{
	double ccpv;                          // curent coarse point value
	
	int m_fine = ((m_coarse-1) *2) +1;    // m_fine = (2 x multiple de 8 )+1
	int x, y, xlim, ylim;                 //coarses coordinates
	int x_f, y_f;                         //coordonnées de la grille fine coresspondant aux coordonnées (x,y) de la grille grossière
	xlim = ((m_coarse-1) * x_ratio) - 1;  
    ylim = ((m_coarse-1) * y_ratio) - 1;  //coordonnées grossières du point caractérisant la forme du problème 
	
	for (y = 0; y < m_coarse-2; y++){
		for (x = 0; x < m_coarse-2; x++){    // on parcourt tous les points de la grille grossière (x,y)
			if (y < ylim || x < xlim){       // si on est dans la membrane (et pas sur le bord)
				x_f = (2*x) +1;
				y_f = (2*y) +1;
				
				ccpv = x_coarse[elem_number( m_coarse, x, y) ]; // récupère la valeur du point(x,y) de la grille coarse
				
				x_fine[elem_number( m_fine, x_f   , y_f  )] = ccpv;
				
				x_fine[elem_number( m_fine, x_f+1 , y_f  )] += 0.5*ccpv;
				x_fine[elem_number( m_fine, x_f-1 , y_f  )] += 0.5*ccpv;
				x_fine[elem_number( m_fine, x_f   , y_f+1)] += 0.5*ccpv;
				x_fine[elem_number( m_fine, x_f   , y_f-1)] += 0.5*ccpv;
				
				x_fine[elem_number( m_fine, x_f+1 , y_f+1)] += 0.25*ccpv;
				x_fine[elem_number( m_fine, x_f+1 , y_f-1)] += 0.25*ccpv;
				x_fine[elem_number( m_fine, x_f-1 , y_f+1)] += 0.25*ccpv;
				x_fine[elem_number( m_fine, x_f-1 , y_f-1)] += 0.25*ccpv;
			} // pour chaque point de la grille grossière, on ajoute 
			  // leur contribution à la grille fine en fonction de leur emplacement
			  // ex : - un point commun aux deux grilles se verra attribuer
			  //      la même valeur que la grille grossière sur la grille fine
			  //
			  //      - un point fin situé au milieu d'un carré de 4 points 
			  //      grossiers recevra 4 fois (0.25 * la valeur de ces 4 points grossiers)
		}
	}
}
