#include "main.h"  

int elem_number(int m, int x, int y)

/* 
 * Return the number of a point in the lexicographic ordering 
 * giving its cartesian coordinates and the mesh size of the grid
 * 
 * The cartésian coordinates x and y are included in [O, m-2]  
 * 
 * INPUT: 
 * 			- x, y: cartesian coordinates of the point (x,y)
 * 			- m   : caracterise the mesh size (= number of point on one side)
 * 
 * 				(The geometric caracteristics of the grid (x_ratio, y_ratio)
 * 				 are given in main.h)
 * 
 * OUTPUT:
 * 			-res: position of the point (x,y) in lexicographic ordering
 */
 
 
{
	int xlim, ylim;   //point that caracterise the shape of the problem
	int res;
	int nx = m-2;     // exculdes the point on the boundary
	xlim = ((m-1) * x_ratio) - 1;
    ylim = ((m-1) * y_ratio) - 1;

	if(y < ylim){
		res = y*nx + x; 
	} 
	
	else{
		res = ylim * nx; //number of point in the inferior rectangle
		res += (y-ylim) * xlim;
		res += x;
	}
	
	return res;
}
