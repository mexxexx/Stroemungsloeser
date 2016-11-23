#ifndef INIT_H
#define INIT_H

/* Rechnet eine 2D Position für den Zugriff auf ein 1D Array um */
#define POS2D(i, j, width) ((i) + ((j) * (width)))

void readParameter(char *filename, char *simulationName, double *xlength, double *ylength, int *imax, int *jmax, 
	double *delx, double *dely, double *delt, double *del_vec, double *t_end, double *tau, int *itermax,
	double *eps, double *omg, double *alpha, double *Re, double *GX, double *GY, double *UI, double *VI, double *PI);
	
void initField(double *field, int imax, int jmax, double value);

#endif
