#ifndef INIT_H
#define INIT_H

/* Rechnet eine 2D Position für den Zugriff auf ein 1D Array um */
#define POS2D(i, j, width) ((i) + ((j) * (width)))

void readParameter(char *filename, char *simulationName, char *heightMap, double *xlength, double *ylength, int *imax, int *jmax, 
	double *delx, double *dely, double *delt, double *del_vec, double *t_end, double *tau, int *itermax,
	double *eps, double *omg, double *alpha, double *Re, double *GX, double *GY, double *UI, double *VI, double *PI, int *wl, int *wr, int *wt, int *wb);
	
void initField(double *field, int imax, int jmax, double value);

int allocateVector(double **vector, int size);

void initFlag(char *heightMap, char *FLAG, int *imax, int *jmax);
#endif
