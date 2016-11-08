#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include "boundary.h"
#include "init.h"
#include "uvp.h"
#include "visual.h"

int allocateVector(double **vector, int size) {
	*vector = (double*)malloc(size * sizeof(double));
	if (*vector == NULL) {
		printf("Konnte keinen Speicher allokieren.\n");
		return 1;
	}

	return 0;
}

void printSnapshot(char *simulationName, int currentSnapshot, double *U, double *V, double *P, int imax, int jmax, 
					double xlength, double ylength, double t, double t_end) {
	char nameOfVelocity[256];
	char nameOfPressure[256];
	sprintf(nameOfVelocity, "%s/%s_velocity_%i.vtk", simulationName, simulationName, currentSnapshot);
	sprintf(nameOfPressure, "%s/%s_pressure_%i.vtk", simulationName, simulationName, currentSnapshot++);
	printVectorField(U, V, imax, jmax, xlength, ylength, nameOfVelocity); 
	printScalarField(P, imax, jmax, xlength, ylength, nameOfPressure);
	printf("Time: %f/%f (%3.1f%% completed)\n", t, t_end, t/t_end * 100);		
}

void initDrivenCavity(double *U, double *V, int imax, int jmax) {
	for (int i = 1; i <= imax; i++) {
		U[POS2D(i, jmax+1, imax+2)] = (10.0 - U[POS2D(i, jmax, imax+2)]);
		//U[POS2D(i, 0, imax+2)] = -(20.0 - U[POS2D(i, 1, imax+2)]);
	}
	/*for (int j = 1; j <= jmax; j++) {
		V[POS2D(0, j, imax+2)] = (20.0 - V[POS2D(1, j, imax+2)]);
		V[POS2D(imax+1, j, imax+2)] = -(20.0 - V[POS2D(imax, j, imax+2)]);
	}*/
}

void printMatrix(double *matrix, int imax, int jmax) {
	for (int j = jmax+1; j>=0; j--) {
		for (int i = 0; i < imax+2; i++) 
			printf("%f ", matrix[POS2D(i, j, imax+2)]);
		printf("\n");
	}
	printf("\n");
}

int calculateFluidDynamics(char* simulationName, double xlength, double ylength, int imax, int jmax, double delx, double dely, 
		double delt, double t_end, double del_vec, double tau, int itermax, double eps, double omg, 
		double alpha, double Re, double GX, double GY, double UI, double VI, double PI, double *U, double *V, double *P) {
	initField(U, imax, jmax, UI);
	initField(V, imax, jmax, VI);
	initField(P, imax, jmax, PI);
	
	double *F, *G, *rhs;
	if (allocateVector(&F, (imax+2) * (jmax+2)))
		return 1;
	if (allocateVector(&G, (imax+2) * (jmax+2))) {
		free(F);
		return 1;
	}
	if (allocateVector(&rhs, (imax+2) * (jmax+2))) {
		free(F);
		free(G);
		return 1;
	}		
			
	int currentSnapshot = 1;
	double t = 0;
	double frameDuration = 0;
	double umax = (fabs(UI) > eps) ? UI : 1;
	double vmax = (fabs(VI) > eps) ? VI : 1;
	//U[POS2D(imax/2, jmax/2, imax+2)]=1;
	while (t < t_end) {
		computeDelt(&delt, imax, jmax, delx, dely, umax, vmax, Re, tau);
		setBoundaryCont(U, V, imax, jmax);
		initDrivenCavity(U, V, imax, jmax);
		
		computeFG(U, V, F, G, imax, jmax, delt, delx, dely, GX, GY, alpha, Re);
		computeRHS(F, G, rhs, imax, jmax, delt, delx, dely);
		solvePoisson(P, rhs, omg, eps, itermax, delx, dely, imax, jmax);
		
		adapUV(U, V, F, G, P, imax, jmax, delt, delx, dely, &umax, &vmax); 
		
		if (frameDuration >= del_vec) {
			printSnapshot(simulationName, currentSnapshot++, U, V, P, imax, jmax, xlength, ylength, t, t_end);
			//printMatrix(P, imax, jmax);
			frameDuration -= del_vec;
		}
		t += delt;
		frameDuration += delt;
	}
	
	printSnapshot(simulationName, currentSnapshot, U, V, P, imax, jmax, xlength, ylength, t_end, t_end);
	free(F);
	free(G);
	free(rhs);
	return 0;
}

void createSimulationDirectory(char *simulationName) {
	struct stat st = {0};

	if (stat(simulationName, &st) == -1)
		mkdir(simulationName, 0700); 
}

int main() {
	char simulationName[256];
	int imax, jmax, itermax;
	double xlength, ylength, delx, dely, delt, t_end, del_vec, tau, eps, omg, alpha, Re, GX, GY, UI, VI, PI;
	int readingResult = readParameter("test.txt", simulationName, &xlength, &ylength, &imax, &jmax, &delx, &dely, &delt, 
							&del_vec, &t_end, &tau, &itermax, &eps, &omg, &alpha, &Re, &GX, &GY, &UI, &VI, &PI);
	
	createSimulationDirectory(simulationName);
	
	if (readingResult < 0) {
		printf("Fehler beim Öffnen der Datei.\n");
		return 1;
	}
	else if (readingResult > 0) {
		printf("Fehler beim Einlesen des %i-ten Parameters.\n", readingResult);
		return 1;
	}
	
	double *U, *V, *P;
	if (allocateVector(&U, (imax+2) * (jmax+2)))
		return 1;
	if (allocateVector(&V, (imax+2) * (jmax+2))) {
		free(U);
		return 1;
	}
	if (allocateVector(&P, (imax+2) * (jmax+2))) {
		free(U);
		free(V);
		return 1;
	}
	
	calculateFluidDynamics(simulationName, xlength, ylength, imax, jmax, delx, dely, delt, 
					t_end, del_vec, tau, itermax, eps, omg, alpha, Re, GX, GY, UI, VI, PI, U, V, P);
	
	free(U);
	free(V);
	free(P);
	return 0;
}
