#include <stdio.h>
#include <stdlib.h>
#include "visual.h"
#include "init.h"

void printScalarField(double *field, int imax, int jmax, double xlength, double ylength, char *filename) {
	const double deltaX = xlength / imax;
	const double deltaY = ylength / jmax;
	
	FILE *f = fopen(filename, "w");
	if (f == NULL)
	{
		printf("Error opening file!\n");
		return;
	}
	
	fprintf(f, "# vtk DataFile Version 3.0\n");
	fprintf(f, "Scalar Field\n");
	fprintf(f, "ASCII\n");
	fprintf(f, "DATASET RECTILINEAR_GRID\n");
	fprintf(f, "DIMENSIONS %i %i 1\n", imax, jmax);
	fprintf(f, "X_COORDINATES %i double\n", imax);
	for (int i = 1; i <= imax; i++) 
		fprintf(f, "%f ", (i-0.5) * deltaX);
	fprintf(f, "\nY_COORDINATES %i double\n", jmax);
	for (int j = 1; j <= jmax; j++) 
		fprintf(f, "%f ", (j-0.5) * deltaY);
	fprintf(f, "\nZ_COORDINATES 1 double\n");
	fprintf(f, "0.0\n");
	fprintf(f, "POINT_DATA %i\n", imax * jmax);
	fprintf(f, "SCALARS Skalarfeld double 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int j = 1; j <= jmax; j++) 
		for (int i = 1; i <= imax; i++) 
			fprintf(f, "%f\n", field[POS2D(i, j, imax+2)]);
	
	fclose(f);
}

void printVectorField(double *U, double *V, int imax, int jmax, double xlength, double ylength, char *filename) {
	const double deltaX = xlength / imax;
	const double deltaY = ylength / jmax;
	
	FILE *f = fopen(filename, "w");
	if (f == NULL)
	{
		printf("Error opening %s!\n", filename);
		return;
	} 
	
	fprintf(f, "# vtk DataFile Version 3.0\n");
	fprintf(f, "Vector Field\n");
	fprintf(f, "ASCII\n");
	fprintf(f, "DATASET RECTILINEAR_GRID\n");
	fprintf(f, "DIMENSIONS %i %i 1\n", imax, jmax);
	fprintf(f, "X_COORDINATES %i double\n", imax);
	for (int i = 1; i <= imax; i++) 
		fprintf(f, "%f ", (i-0.5) * deltaX);
	fprintf(f, "\nY_COORDINATES %i double\n", jmax);
	for (int j = 1; j <= jmax; j++) 
		fprintf(f, "%f ", (j-0.5) * deltaY);
	fprintf(f, "\nZ_COORDINATES 1 double\n");
	fprintf(f, "0.0\n");
	fprintf(f, "POINT_DATA %i\n", imax * jmax);
	fprintf(f, "VECTORS Vektorfeld double\n");
	for (int j = 1; j <= jmax; j++) 
		for (int i = 1; i <= imax; i++) 
			fprintf(f, "%f %f 0.0\n", 
				0.5 * (U[POS2D(i, j, imax+2)]+U[POS2D(i-1, j, imax+2)]),
				0.5 * (V[POS2D(i, j, imax+2)]+V[POS2D(i, j-1, imax+2)]));
	
	fclose(f);
}

void printObstacles(char *field, int imax, int jmax, double xlength, double ylength, char *filename){
	const double deltaX = xlength / imax;
	const double deltaY = ylength / jmax;
	
	FILE *f = fopen(filename, "w");
	if (f == NULL)
	{
		printf("Error opening file!\n");
		return;
	}
	
	fprintf(f, "# vtk DataFile Version 3.0\n");
	fprintf(f, "Obstacles\n");
	fprintf(f, "ASCII\n");
	fprintf(f, "DATASET RECTILINEAR_GRID\n");
	fprintf(f, "DIMENSIONS %i %i 1\n", imax, jmax);
	fprintf(f, "X_COORDINATES %i double\n", imax);
	for (int i = 1; i <= imax; i++) 
		fprintf(f, "%f ", (i-0.5) * deltaX);
	fprintf(f, "\nY_COORDINATES %i double\n", jmax);
	for (int j = 1; j <= jmax; j++) 
		fprintf(f, "%f ", (j-0.5) * deltaY);
	fprintf(f, "\nZ_COORDINATES 1 double\n");
	fprintf(f, "0.0\n");
	fprintf(f, "POINT_DATA %i\n", imax * jmax);
	fprintf(f, "SCALARS Skalarfeld char 1\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	for (int j = 1; j <= jmax; j++) 
		for (int i = 1; i <= imax; i++) 
			fprintf(f, "%d\n", (field[POS2D(i, j, imax+2)]) ? 1 : 0);
	
	fclose(f);
}
