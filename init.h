#ifndef INIT_H
#define INIT_H

/* Rechnet eine 2D Position für den Zugriff auf ein 1D Array um */
#define POS2D(i, j, width) ((i) + ((j) * (width)))

void readParameter(char *filename, char *simulationName, char *heightMap, double *xlength, double *ylength, int *imax, int *jmax, 
	double *delx, double *dely, double *delt, double *del_vec, double *t_end, double *tau, int *itermax,
	double *eps, double *omg, double *alpha, double *Re, double*Pr, double *beta, double *GX, double *GY, double *UI, double *VI, double *PI,
	double *TI, int *wl, int *wr, int *wt, int *wb, double *posx1, double *posx2, double *posy1, double *posy2,
	int *tl, double *tl_value, int *tr, double *tr_value, int *tt, double *tt_value, int *tb, double *tb_value);
	
void initField(double *field, int imax, int jmax, double value);

int allocateVector(double **vector, int size);

void initFlag(char *heightMap, char *FLAG, int imax, int jmax, int *numFluidCells);
#endif
