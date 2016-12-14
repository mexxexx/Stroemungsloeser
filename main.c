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
#include "particles.h"

#define PROGRESS_DISPLAY 0.05
int currentProgressStep = 1;
time_t start;

char simulationName[256], problem[256];
int imax, jmax, itermax, wl, wr, wt, wb, numFluidCells, partCount, anzahl;
double xlength, ylength, delx, dely, delt, t_end, del_vec, tau, eps, omg, alpha, Re, GX, GY, UI, VI, PI, posx1, posx2, posy1, posy2;

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

void printSnapshot(int currentSnapshot, double *U, double *V, double *P, Particle *particles, int partCount, double t) {
	char nameOfVelocity[256];
	char nameOfPressure[256];
	char nameOfParticles[256];
	sprintf(nameOfVelocity, "%s/%s_velocity_%05d.vtk", simulationName, simulationName, currentSnapshot);
	sprintf(nameOfParticles, "%s/%s_particles_%05d.vtk", simulationName, simulationName, currentSnapshot);
	sprintf(nameOfPressure, "%s/%s_pressure_%05d.vtk", simulationName, simulationName, currentSnapshot++);
	printVectorField(U, V, imax, jmax, xlength, ylength, nameOfVelocity); 
	printScalarField(P, imax, jmax, xlength, ylength, nameOfPressure);
	printParticles(particles, partCount, nameOfParticles);
	
	double progress = t/t_end; 
	if (progress >= PROGRESS_DISPLAY * currentProgressStep) {
		printf("%.1f%% completed\n", progress * 100);		
		
		currentProgressStep++;		
	}
}

void printField(double *field, int imax, int jmax) {
	for (int j = jmax+1; j >= 0; j--) {
		for (int i = 0; i <= imax+1; i++)
			printf("%f ", field[POS2D(i, j, imax+2)]);
		printf("\n");
	}
}

void printCharField(char *field, int imax, int jmax) {
	for (int j = jmax+1; j >= 0; j--) {
		for (int i = 0; i <= imax+1; i++)
			printf("%d ", field[POS2D(i, j, imax+2)]);
		printf("\n");
	}
}

int calculateFluidDynamics(double* U, double* V, double* P, char* FLAG, Particle *particles, int partCount) {
	initField(U, imax, jmax, UI);
	if (strcmp("Stufe", problem) == 0) {
		for (int i = 1; i <= imax; i++) 
			for (int j = 0; j <= jmax/2; j++) 
				U[POS2D(i,j,imax+2)]=0;
	}
	
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
	double seedTime = 0;
	double umax = (fabs(UI) > eps) ? UI : eps;
	double vmax = (fabs(VI) > eps) ? VI : eps;
	start = time(NULL);
	
	while (t < t_end) {
		computeDelt(&delt, imax, jmax, delx, dely, umax, vmax, Re, tau);
		setBoundaryCond(U, V, FLAG, imax, jmax, wl, wr, wt, wb);
		setSpecialBoundaryCond(U, V, imax, jmax, problem);
		computeFG(U, V, F, G, FLAG, imax, jmax, delt, delx, dely, GX, GY, alpha, Re);
		computeRHS(F, G, rhs, FLAG, imax, jmax, delt, delx, dely);
		solvePoisson(P, rhs, FLAG, omg, eps, itermax, delx, dely, imax, jmax, numFluidCells);
		adapUV(U, V, F, G, P, FLAG, imax, jmax, delt, delx, dely, &umax, &vmax); 

		if (seedTime >= del_vec/6) {
			particleSeed(particles, posx1, posx2, posy1, posy2, partCount, anzahl);
			seedTime -= del_vec/6;
		}
		particleVelocity(U, V, delx, dely, imax, jmax, particles, partCount);
		particleTransport(particles, delt, partCount, xlength, ylength);

		if (frameDuration >= del_vec) {
			printSnapshot(currentSnapshot++, U, V, P, particles, partCount, t);
			frameDuration -= del_vec;
		}
		t += delt;
		frameDuration += delt;
		seedTime += delt;
	}
	
	printSnapshot(currentSnapshot, U, V, P, particles, partCount, t_end);
	printf("%i frames taken in a time period of %.2f seconds (%.1f FPS)\n", currentSnapshot, t_end, currentSnapshot/t_end);
	printf("Duration: %.1f seconds\n", (double)(time(NULL)-start)); 
	free(F);
	free(G);
	free(rhs);
	return 0;
}

void readParams(int argc, char** argv, char* file) {
	if (argc==1) {
		printf("Bitte eine Datei auswählen: ");
		scanf("%s", file);
	}
	else
		sprintf(file, "%s", argv[1]);
	
	if (argc == 2){
		printf("Bitte ein Problem angeben: ");
		scanf("%s", problem);
	}	
	else 
		sprintf(problem, "%s", argv[2]);
		
	if (argc == 3)
		partCount = 10000;
	else 
		sscanf(argv[3], "%i", &partCount); 
		
	if (partCount < 200) 
		partCount = 200;
}

int main(int argc, char** argv) {
	char file[256];
	readParams(argc, argv, file);
	
	Particle *particles = malloc(sizeof(Particle)*partCount);
	if (particles == NULL) {
		printf("Konnte keinen Speicherplatz fuer particles allokieren");
		return 1;
	}
	
	particleInit(particles, partCount);
	char obstacelsMap[256];
	readParameter(file, simulationName, obstacelsMap, &xlength, &ylength, &imax, &jmax, &delx, &dely, &delt, 
							&del_vec, &t_end, &tau, &itermax, &eps, &omg, &alpha, &Re, &GX, &GY, &UI, &VI, &PI, &wl, &wr, &wt, &wb,
							&posx1, &posx2, &posy1, &posy2);
	
	createSimulationDirectory();
	anzahl = (int)(sqrt((posx1-posx2)*(posx1-posx2)+(posy1-posy2)*(posy1-posy2)))*30;
	
	char *FLAG = (char*)malloc((imax+2) * (jmax+2) * sizeof(char));
	if (FLAG == NULL) {
		free(particles);
		printf("Konnte keinen Speicherplatz fuer FLAG allokieren");
		return 1;
	}
	initFlag(obstacelsMap, FLAG, imax, jmax, &numFluidCells);
	
	char obstacleFile[256];
	sprintf(obstacleFile, "%s/obstacles.vtk", simulationName);
	printObstacles(FLAG, imax, jmax, xlength, ylength, obstacleFile);
	
	double *U, *V, *P;
	if (allocateVector(&U, (imax+2) * (jmax+2))) {
		free(particles);
		free(FLAG);
		return 1;
	}
	if (allocateVector(&V, (imax+2) * (jmax+2))) {
		free(particles);
		free(FLAG);
		free(U);
		return 1;
	}
	if (allocateVector(&P, (imax+2) * (jmax+2))) {
		free(particles);
		free(FLAG);
		free(U);
		free(V);
		return 1;
	}
	
	calculateFluidDynamics(U, V, P, FLAG, particles, partCount);
	
	free(particles);
	free(FLAG);
	free(U);
	free(V);
	free(P);
	return 0;
}
