#include <stdio.h>
#include "init.h"

int readParameter(char *filename, char *simulationName, double *xlength, double *ylength, int *imax, int *jmax, 
	double *delx, double *dely, double *delt, double *del_vec, double *t_end, double *tau, int *itermax,
	double *eps, double *omg, double *alpha, double *Re, double *GX, double *GY, double *UI, double *VI, double *PI) {
	FILE *f = fopen(filename, "r");
	if (f == NULL) return -1;
	
	int error = 0;
	if (error == 0 && fscanf(f, "%*s = %s\n", simulationName) != 1)
		error = 1;
	if (error == 0 && fscanf(f, "%*s = %lf\n", xlength) != 1)
		error = 2;
	if (error == 0 && fscanf(f, "%*s = %lf\n", ylength) != 1)
		error =  3;
	if (error == 0 && fscanf(f, "%*s = %i\n", imax) != 1)
		error =  4;
	if (error == 0 && fscanf(f, "%*s = %i\n", jmax) != 1)
		error =  5;
	(*delx) = (*xlength) / (*imax);	
	(*dely) = (*ylength) / (*jmax);	
		
	if (error == 0 && fscanf(f, "%*s = %lf\n", delt) != 1)
		error =  6;
	if (error == 0 && fscanf(f, "%*s = %lf\n", t_end) != 1)
		error =  7;
	if (error == 0 && fscanf(f, "%*s = %lf\n", del_vec) != 1)
		error =  8;
	if (error == 0 && fscanf(f, "%*s = %lf\n", tau) != 1)
		error =  9;
	if (error == 0 && fscanf(f, "%*s = %i\n", itermax) != 1)
		error =  10;
	if (error == 0 && fscanf(f, "%*s = %lf\n", eps) != 1)
		error =  11;
	if (error == 0 && fscanf(f, "%*s = %lf\n", omg) != 1)
		error =  12;
	if (error == 0 && fscanf(f, "%*s = %lf\n", alpha) != 1)
		error =  13;
	if (error == 0 && fscanf(f, "%*s = %lf\n", Re) != 1)
		error =  14;
	if (error == 0 && fscanf(f, "%*s = %lf\n", GX) != 1)
		error =  15;
	if (error == 0 && fscanf(f, "%*s = %lf\n", GY) != 1)
		error =  16;
	if (error == 0 && fscanf(f, "%*s = %lf\n", UI) != 1)
		error =  17;
	if (error == 0 && fscanf(f, "%*s = %lf\n", VI) != 1)
		error =  18;
	if (error == 0 && fscanf(f, "%*s = %lf\n", PI) != 1)
		error =  19;
		
	fclose(f);
	
	return error;
}

void initField(double *field, int imax, int jmax, double value) {
	for (int j = 0; j <= jmax + 1; j++) 
		for (int i = 0; i <= imax + 1; i++)
			field[POS2D(i, j, imax + 2)] = value;
}
