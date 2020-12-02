#include <stdio.h>
#include "multi_grid.h"
#include "main.h"

void fine_to_coarse(int m_fine, double *r_fine, double *r_coarse)


/*
 * r_coarse = R * r_fine   (using full weighting)
 * 
 * FULL WEIGHTING: each coarse point is the weighted average of the 9
 * 					surrounding points
 *  
 * Restric the r_fine vector (r_f = b - A*u) corresponding to a fine grid
 * of "m_fine" points on each side to a new r_coarse vector corresponding to
 * a coarse grid (with "m_coarse" points on each side)
 * 
 * The coarse grid is defined by h_coarse (distance between 2 points): 2.h_fine = h_coarse
 * 
 * 
 * INPUT :
 * 			m_fine : parameter of the fine grid (number of point one one side)
 * 			r_fine : initial vector of dimension corresponding to a fine grid
 * 
 * OUTPUT : 
 * 			r_coarse : restricted vector of dimension corresponding to a coarse grid
 * 
 */ 

{
	int m_coarse = ((m_fine-1) /2) +1;    // m_fine = (2 x multiple of 8 )+1
	int count_coarse = 0;				  // Lexicographic number of the current coarse points
	int x, y, xlim, ylim;                 //coarse coordinates
	int x_f, y_f;                         //fine grid coordinates coresponding to the coarse coordinates (x,y)
	xlim = ((m_coarse-1) * x_ratio) - 1;  
    ylim = ((m_coarse-1) * y_ratio) - 1;  //coarse coordinates of the point representing the geometry of the problem
	
	for (y = 0; y < m_coarse-2; y++){
		for (x = 0; x < m_coarse-2; x++){    // we run through each corse point
			if (y < ylim || x < xlim){       // if we are inside the domain (and not on the edges)
				x_f = (2*x) +1;
				y_f = (2*y) +1;
				
				r_coarse[count_coarse] = 0.25 * r_fine[elem_number( m_fine, x_f, y_f) ];
				
				r_coarse[count_coarse] += 0.125 * r_fine[elem_number( m_fine, x_f-1, y_f) ];
				r_coarse[count_coarse] += 0.125 * r_fine[elem_number( m_fine, x_f+1, y_f) ];
				r_coarse[count_coarse] += 0.125 * r_fine[elem_number( m_fine, x_f, y_f-1) ];
				r_coarse[count_coarse] += 0.125 * r_fine[elem_number( m_fine, x_f, y_f+1) ];
				
				r_coarse[count_coarse] += 0.0625 * r_fine[elem_number( m_fine, x_f-1, y_f-1) ];
				r_coarse[count_coarse] += 0.0625 * r_fine[elem_number( m_fine, x_f+1, y_f-1) ];
				r_coarse[count_coarse] += 0.0625 * r_fine[elem_number( m_fine, x_f+1, y_f+1) ];
				r_coarse[count_coarse] += 0.0625 * r_fine[elem_number( m_fine, x_f-1, y_f+1) ];
				
				count_coarse +=1;
			}
		}
	}
	
	//verification
	//printf("count_coarse = %d \n ", count_coarse);
	
}
