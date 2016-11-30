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

int calculateFluidDynamics(double xlength, double ylength, int imax, int jmax, double delx, double dely, 
		double delt, double t_end, double del_vec, double tau, int itermax, double eps, double omg, 
		double alpha, double Re, double GX, double GY, double UI, double VI, double PI, double *U, double *V, double *P, char *FLAG, int numFluidCells, int wl, int wr, int wt, int wb, char* problem) {
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
		setBoundaryCond(U, V, FLAG, imax, jmax, wl, wr, wt, wb);
		setSpecialBoundaryCond(U, V, imax, jmax, problem);
		
		computeFG(U, V, F, G, FLAG, imax, jmax, delt, delx, dely, GX, GY, alpha, Re);
		computeRHS(F, G, rhs, imax, jmax, delt, delx, dely);
		solvePoisson(P, rhs, FLAG, omg, eps, itermax, delx, dely, imax, jmax, numFluidCells);
		adapUV(U, V, F, G, P, FLAG, imax, jmax, delt, delx, dely, &umax, &vmax); 
		
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
	char problem[256];
	char obstacelsMap[256];
	printf("test0\n");
	if (argc==1) {
		printf("Bitte eine Datei auswählen: ");
		scanf("%s", file);
	}
	else	{
		sprintf(file, "%s", argv[1]);
		if (argc == 3)
			sprintf(problem, "%s", argv[2]);	
	}
	
	printf("test1\n");
	int imax, jmax, itermax, wl, wr, wt, wb, numFluidCells;
	double xlength, ylength, delx, dely, delt, t_end, del_vec, tau, eps, omg, alpha, Re, GX, GY, UI, VI, PI;
	readParameter(file, simulationName, obstacelsMap, &xlength, &ylength, &imax, &jmax, &delx, &dely, &delt, 
							&del_vec, &t_end, &tau, &itermax, &eps, &omg, &alpha, &Re, &GX, &GY, &UI, &VI, &PI, &wl, &wr, &wt, &wb);
	
	createSimulationDirectory();
	
	printf("test2\n");
	char *FLAG = (char *)malloc((imax+2) * (jmax+2) * sizeof(char));
	initFlag(obstacelsMap, FLAG, imax, jmax, &numFluidCells);
	xlength = imax / 5.0;
	ylength = jmax / 5.0;
	char obstacleFile[256];
	sprintf(obstacleFile, "%s/obstacles.vtk", simulationName);
	printObstacles(FLAG, imax, jmax, xlength, ylength, obstacleFile);
	printf("%i\n", numFluidCells);
	double *U, *V, *P;
	if (allocateVector(&U, (imax+2) * (jmax+2))) {
		free(FLAG);
		return 1;
	}
	if (allocateVector(&V, (imax+2) * (jmax+2))) {
		free(FLAG);
		free(U);
		return 1;
	}
	if (allocateVector(&P, (imax+2) * (jmax+2))) {
		free(FLAG);
		free(U);
		free(V);
		return 1;
	}
	
	calculateFluidDynamics(xlength, ylength, imax, jmax, delx, dely, delt, 
					t_end, del_vec, tau, itermax, eps, omg, alpha, Re, GX, GY, UI, VI, PI, U, V, P, FLAG, numFluidCells, wl, wr, wt, wb, problem);
	
	free(FLAG);
	free(U);
	free(V);
	free(P);
	return 0;
}
