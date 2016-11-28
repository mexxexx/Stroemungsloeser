#include <stdio.h>
#include <stdlib.h>
#include "init.h"

int allocateVector(double **vector, int size) {
	*vector = (double*)malloc(size * sizeof(double));
	if (*vector == NULL) {
		printf("Konnte keinen Speicher allokieren.\n");
		return 1;
	}

	return 0;
}

void parameterError(char *param) {
	printf("Error reading parameter %s\n", param);
	exit(EXIT_FAILURE);
}

void readParameter(char *filename, char *simulationName, double *xlength, double *ylength, int *imax, int *jmax, 
	double *delx, double *dely, double *delt, double *del_vec, double *t_end, double *tau, int *itermax,
	double *eps, double *omg, double *alpha, double *Re, double *GX, double *GY, double *UI, double *VI, double *PI) {
	FILE *f = fopen(filename, "r");
	if (f == NULL) {
		printf("Error opening parameter file\n");
		exit(EXIT_FAILURE);
	}
	
	if (fscanf(f, "%*s = %s\n", simulationName) != 1)
		parameterError("simulationName");
	if (fscanf(f, "%*s = %lf\n", xlength) != 1)
		parameterError("xlength");
	if (fscanf(f, "%*s = %lf\n", ylength) != 1)
		parameterError("ylength");
	if (fscanf(f, "%*s = %i\n", imax) != 1)
		parameterError("imax");
	if (fscanf(f, "%*s = %i\n", jmax) != 1)
		parameterError("jmax");
	(*delx) = (*xlength) / (*imax);	
	(*dely) = (*ylength) / (*jmax);	
		
	if (fscanf(f, "%*s = %lf\n", delt) != 1)
		parameterError("delt");
	if (fscanf(f, "%*s = %lf\n", t_end) != 1)
		parameterError("t_end");
	if (fscanf(f, "%*s = %lf\n", del_vec) != 1)
		parameterError("del_vec");
	if (fscanf(f, "%*s = %lf\n", tau) != 1)
		parameterError("tau");
	if (fscanf(f, "%*s = %i\n", itermax) != 1)
		parameterError("itermax");
	if (fscanf(f, "%*s = %lf\n", eps) != 1)
		parameterError("eps");
	if (fscanf(f, "%*s = %lf\n", omg) != 1)
		parameterError("omg");
	if (fscanf(f, "%*s = %lf\n", alpha) != 1)
		parameterError("alpha");
	if (fscanf(f, "%*s = %lf\n", Re) != 1)
		parameterError("Re");
	if (fscanf(f, "%*s = %lf\n", GX) != 1)
		parameterError("GX");
	if (fscanf(f, "%*s = %lf\n", GY) != 1)
		parameterError("GY");
	if (fscanf(f, "%*s = %lf\n", UI) != 1)
		parameterError("UI");
	if (fscanf(f, "%*s = %lf\n", VI) != 1)
		parameterError("VI");
	if (fscanf(f, "%*s = %lf\n", PI) != 1)
		parameterError("PI");
		
	fclose(f);
}

void initField(double *field, int imax, int jmax, double value) {
	for (int j = 0; j <= jmax + 1; j++) 
		for (int i = 0; i <= imax + 1; i++)
			field[POS2D(i, j, imax + 2)] = value;
}
