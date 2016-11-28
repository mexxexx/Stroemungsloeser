#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "boundary.h"
#include "init.h"
#include "uvp.h"
#include "visual.h"

#define PROGRESS_DISPLAY 0.05
int currentProgressStep = 1;
time_t start;

char simulationName[256];

void createSimulationDirectory() {
	struct stat st = {0};
	if (stat(simulationName, &st) == -1) {
		int result;
		#ifdef __MSDOS__
		result = mkdir(simulationName);
		#else
		result = mkdir(simulationName, 0700);
		#endif
		
		if (result) {
			printf("Verzeichnis für die Simulation konnte nicht angelegt werden\n");
			exit(EXIT_FAILURE);
		}
	}
}

void printSnapshot(int currentSnapshot, double *U, double *V, double *P, int imax, int jmax, 
					double xlength, double ylength, double t, double t_end) {
	char nameOfVelocity[256];
	char nameOfPressure[256];
	sprintf(nameOfVelocity, "%s/%s_velocity_%05d.vtk", simulationName, simulationName, currentSnapshot);
	sprintf(nameOfPressure, "%s/%s_pressure_%05d.vtk", simulationName, simulationName, currentSnapshot++);
	printVectorField(U, V, imax, jmax, xlength, ylength, nameOfVelocity); 
	printScalarField(P, imax, jmax, xlength, ylength, nameOfPressure);
	
	double progress = t/t_end; 
	if (progress >= PROGRESS_DISPLAY * currentProgressStep) {
		double elapsedSeconds = (double)(time(NULL)-start);
		double secondsPerProgressDisplay = elapsedSeconds / currentProgressStep;
		double secondsLeft = (1.0/PROGRESS_DISPLAY - currentProgressStep) * secondsPerProgressDisplay;
		
		printf("%.1f%% completed ", progress * 100);		
		if (secondsLeft >= 120) 
			printf("(%.0f minutes remaining)\n", secondsLeft / 60);
		else 
			printf("(%.0f seconds remaining)\n", secondsLeft);
		
		currentProgressStep++;		
	}
}

void printMatrix(double *matrix, int imax, int jmax) {
	for (int j = jmax+1; j>=0; j--) {
		for (int i = 0; i < imax+2; i++) 
			printf("%f ", matrix[POS2D(i, j, imax+2)]);
		printf("\n");
	}
	printf("\n");
}

void initDrivenCavity(double *U, double *V, int imax, int jmax) {
	for (int i = 1; i <= imax; i++) {
		U[POS2D(i, jmax+1, imax+2)] = (2.0 - U[POS2D(i, jmax, imax+2)]);
		//U[POS2D(i, 0, imax+2)] = -(20.0 - U[POS2D(i, 1, imax+2)]);
	}
	/*for (int j = 1; j <= jmax; j++) {
		V[POS2D(0, j, imax+2)] = (20.0 - V[POS2D(1, j, imax+2)]);
		V[POS2D(imax+1, j, imax+2)] = -(20.0 - V[POS2D(imax, j, imax+2)]);
	}*/
}

int calculateFluidDynamics(double xlength, double ylength, int imax, int jmax, double delx, double dely, 
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
	double umax = (fabs(UI) > eps) ? UI : eps;
	double vmax = (fabs(VI) > eps) ? VI : eps;
	start = time(NULL);
	
	while (t < t_end) {
		computeDelt(&delt, imax, jmax, delx, dely, umax, vmax, Re, tau);
		setBoundaryCont(U, V, imax, jmax);
		initDrivenCavity(U, V, imax, jmax);
		
		computeFG(U, V, F, G, imax, jmax, delt, delx, dely, GX, GY, alpha, Re);
		computeRHS(F, G, rhs, imax, jmax, delt, delx, dely);
		solvePoisson(P, rhs, omg, eps, itermax, delx, dely, imax, jmax);
		
		adapUV(U, V, F, G, P, imax, jmax, delt, delx, dely, &umax, &vmax); 
		
		if (frameDuration >= del_vec) {
			printSnapshot(currentSnapshot++, U, V, P, imax, jmax, xlength, ylength, t, t_end);
			frameDuration -= del_vec;
		}
		t += delt;
		frameDuration += delt;
	}
	
	printSnapshot(currentSnapshot, U, V, P, imax, jmax, xlength, ylength, t_end, t_end);
	printf("%i frames taken in a time period of %.2f seconds (%.1f FPS)\n", currentSnapshot, t_end, currentSnapshot/t_end);
	printf("Duration: %.1f seconds\n", (double)(time(NULL)-start)); 
	free(F);
	free(G);
	free(rhs);
	return 0;
}

int main(int argc, char** argv) {
	char file[256];
	if (argc==1) {
		printf("Bitte eine Datei auswählen, you faggot: ");
		scanf("%s", file);
	}
	else
		sprintf(file, "%s", argv[1]);
	
	
	int imax, jmax, itermax;
	double xlength, ylength, delx, dely, delt, t_end, del_vec, tau, eps, omg, alpha, Re, GX, GY, UI, VI, PI;
	readParameter(file, simulationName, &xlength, &ylength, &imax, &jmax, &delx, &dely, &delt, 
							&del_vec, &t_end, &tau, &itermax, &eps, &omg, &alpha, &Re, &GX, &GY, &UI, &VI, &PI);
	
	createSimulationDirectory();
	
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
	//TESTCOMMENTAR
	calculateFluidDynamics(xlength, ylength, imax, jmax, delx, dely, delt, 
					t_end, del_vec, tau, itermax, eps, omg, alpha, Re, GX, GY, UI, VI, PI, U, V, P);
	
	free(U);
	free(V);
	free(P);
	return 0;
}
